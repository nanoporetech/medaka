from collections import defaultdict, Counter, OrderedDict
import concurrent.futures
from copy import deepcopy
import inspect
import itertools
import functools
from multiprocessing import Pool
import os
import sys
from timeit import default_timer as now

from Bio import SeqIO
import numpy as np
import pysam

import medaka.common
import medaka.datastore
import medaka.labels
import libmedaka


def pileup_counts(region, bam, dtype_prefixes=None, region_split=100000, workers=8, tag_name=None, tag_value=None, keep_missing=False):
    """Create pileup counts feature array for region.

    :param region: `medaka.common.Region` object
    :param bam: .bam file with alignments.
    :param dtype_prefixes: prefixes for query names which to separate counts.
        If `None` (or of length 1), counts are not split.
    :param tag_name: two letter tag name by which to filter reads.
    :param tag_value: integer value of tag for reads to keep.
    :param keep_missing: whether to keep reads when tag is missing.

    :returns: pileup counts array, reference positions, insertion postions
    """
    ffi, lib = libmedaka.ffi, libmedaka.lib
    logger = medaka.common.get_named_logger('PileUp')

    if dtype_prefixes is None or isinstance(dtype_prefixes, str) or len(dtype_prefixes) == 1:
        num_dtypes, dtypes = 1, ffi.NULL
    else:
        num_dtypes = len(dtype_prefixes)
        _dtypes = [ffi.new("char[]", d.encode()) for d in dtype_prefixes]
        dtypes = ffi.new("char *[]", _dtypes)

    if tag_name is None:
        tag_name = ffi.new("char[2]", "".encode())
        tag_value = 0
        keep_missing = False
    elif len(tag_name) != 2:
        raise ValueError("'tag_name' must be a length-2 string.")
    else:
        tag_name = ffi.new("char[2]", tag_name.encode())
    featlen = lib.featlen

    def _process_region(reg):
        # htslib start is 1-based, medaka.common.Region object is 0-based
        region_str = '{}:{}-{}'.format(reg.ref_name, reg.start + 1, reg.end)

        counts = lib.calculate_pileup(
            region_str.encode(), bam.encode(), num_dtypes, dtypes,
            tag_name, tag_value, keep_missing
        )

        # TODO: this should ALL probably not be hardcoded. Counts should return
        # all information needed about how to reconstruct the array
        size_sizet = np.dtype(np.uintp).itemsize
        np_counts = np.frombuffer(ffi.buffer(
            counts.counts, size_sizet * counts.n_cols * featlen * num_dtypes),
            dtype=np.uintp
        ).reshape(counts.n_cols, featlen * num_dtypes).copy()

        positions = np.empty(counts.n_cols, dtype=[('major', int), ('minor', int)])
        np.copyto(positions['major'],
            np.frombuffer(ffi.buffer(
            counts.major, size_sizet * counts.n_cols),
            dtype=np.uintp
        ))
        np.copyto(positions['minor'],
            np.frombuffer(ffi.buffer(
            counts.minor, size_sizet * counts.n_cols),
            dtype=np.uintp
        ))

        lib.destroy_plp_data(counts)
        return np_counts, positions

    # split large regions for performance
    regions = region.split(region_split)
    results = None
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        # First pass: need to check for discontinuities within chunks,
        # these show up as >2 changes in the major coordinate
        _results = list()
        for counts, positions in executor.map(_process_region, regions):
            move = np.ediff1d(positions['major'])
            gaps = np.where(move > 1)[0] + 1
            if len(gaps) == 0:
                _results.append((counts, positions))
            else:
                logger.info("Splitting discontiguous pileup region.")
                start = 0
                for i in gaps:
                    _results.append((counts[start:i], positions[start:i]))
                    start = i
                _results.append((counts[start:], positions[start:]))
        results = _results

    # Second pass: stitch abutting chunks together, anything not neighbouring
    # is kept separate whether it came from the same chunk originally or not
    def _finalize_chunk(c_buf, p_buf):
        chunk_counts = np.concatenate(c_buf)
        chunk_positions = np.concatenate(p_buf)
        return chunk_counts, chunk_positions
    counts_buffer, positions_buffer = list(), list()
    chunk_results = list()
    last = None
    for counts, positions in results:
        if len(positions) == 0:
            continue
        first = positions['major'][0]
        if len(counts_buffer) == 0 or first - last == 1:
            # new or contiguous
            counts_buffer.append(counts)
            positions_buffer.append(positions)
            last = positions['major'][-1]
        else:
            # discontinuity
            chunk_results.append(_finalize_chunk(counts_buffer, positions_buffer))
            counts_buffer = [counts]
            positions_buffer = [positions]
            last = positions['major'][-1]
    if len(counts_buffer) != 0:
        chunk_results.append(_finalize_chunk(counts_buffer, positions_buffer))

    return chunk_results


class FeatureEncoder(object):
    _norm_modes_ = ['total', 'fwd_rev', None]


    def __init__(self,
            normalise:str='total', dtypes=('',),
            tag_name=None, tag_value=None, tag_keep_missing=False):
        """Class to generate neural network input features.

        :param normalise: str, how to normalise the data.
        :param dtypes: iterable of str, read id prefixes of distinct data types that should be counted separately.
        :param tag_name: two letter tag name by which to filter reads.
        :param tag_value: integer value of tag for reads to keep.
        :param tag_keep_missing: whether to keep reads when tag is missing.

        """
        self.logger = medaka.common.get_named_logger('Feature')
        self.normalise = normalise
        self.feature_dtype = np.float32
        self.dtypes = dtypes
        self.tag_name = tag_name
        self.tag_value = tag_value
        self.tag_keep_missing = tag_keep_missing

        if self.normalise not in self._norm_modes_:
            raise ValueError('normalise={} is not one of {}'.format(self.normalise, self._norm_modes_))

        opts = inspect.signature(FeatureEncoder.__init__).parameters.keys()
        opts = {k:getattr(self, k) for k in opts if k != 'self'}
        self.logger.debug("Creating features with: {}".format(opts))

        #TODO: this is defined by the C-code, it should be obtained from there
        read_decoding = []
        for dtype in self.dtypes:
            # dtype, rev/fwd, base/gap
            read_decoding += [
                (dtype,) + (direction, base)
                for (direction, base) in itertools.product((True, False), medaka.common._alphabet_)
            ]
            # forward and reverse gaps
            read_decoding += [(dtype, True, None), (dtype, False, None)]
        self.encoding = OrderedDict(((a, i) for i, a in enumerate(read_decoding)))
        self.logger.debug("Feature decoding is:\n{}".format('\n'.join(
            '{}: {}'.format(i, x) for i, x in enumerate(read_decoding)
        )))


    @property
    def feature_indices(self):
        """Location of feature vector components.

        :returns: dictionary mapping read data type and strand to feature vector
            indices (a list over all bases)

        """
        #TODO: should be defined along with self.encoding (is it just read_decoding from __init__)?
        return {
            (dt, strand):
                [v for k, v in self.encoding.items()
                    if k[0] == dt and strand == k[1]]
                for dt, strand in itertools.product(self.dtypes, (True, False))
        }


    def bam_to_sample(self, reads_bam, region):
        """Converts a section of an alignment pileup (as shown
        by e.g. samtools tview) to a base frequency feature array

        :param reads_bam: (sorted indexed) bam with read alignment to reference
        :param region: `medaka.common.Region` object with ref_name, start and end attributes.
        :param start: starting position within reference
        :param end: ending position within reference

        :returns: `medaka.common.Sample` object

        """
        pileup = pileup_counts(
            region, reads_bam, dtype_prefixes=self.dtypes,
            tag_name=self.tag_name, tag_value=self.tag_value,
            keep_missing=self.tag_keep_missing
        )
        samples = list()
        for counts, positions in pileup:
            if len(counts) == 0:
                msg = 'Pileup-feature is zero-length for {} indicating no reads in this region.'.format(region)
                self.logger.warning(msg)
                samples.append(
                    medaka.common.Sample(
                        ref_name=region.ref_name, features=None,
                        labels=None, ref_seq=None,
                        postions=positions, label_probs=None
                    )
                )
                continue

            start, end = positions['major'][0], positions['major'][-1]
            if start != region.start or end + 1 != region.end: # TODO investigate off-by-one
                self.logger.warning(
                    'Pileup counts do not span requested region, requested {}, '
                    'received {}-{}.'.format(region, start, end)
                )

            # find the position index for parent major position of all minor positions
            minor_inds = np.where(positions['minor'] > 0)
            major_pos_at_minor_inds = positions['major'][minor_inds]
            major_ind_at_minor_inds = np.searchsorted(positions['major'], major_pos_at_minor_inds, side='left')

            depth = np.sum(counts, axis=1)
            depth[minor_inds] = depth[major_ind_at_minor_inds]

            if self.normalise == 'total':
                # normalize counts by total depth at major position, since the
                # counts include deletions this is a count of spanning reads
                feature_array = counts / np.maximum(1, depth).reshape((-1, 1)) # max just to avoid div error
            elif self.normalise == 'fwd_rev':
                # normalize forward and reverse and by dtype
                feature_array = np.empty_like(counts, dtype=float)
                for (dt, is_rev), inds in self.feature_indices.items():
                    dt_depth = np.sum(counts[:, inds], axis=1)
                    dt_depth[minor_inds] = dt_depth[major_ind_at_minor_inds]
                    feature_array[: , inds] = counts[:, inds] / np.maximum(1, dt_depth).reshape((-1, 1)) # max just to avoid div error
            else:
                feature_array = counts

            feature_array = feature_array.astype(self.feature_dtype)

            sample = medaka.common.Sample(
                ref_name=region.ref_name, features=feature_array,
                labels=None, ref_seq=None,
                positions=positions, label_probs=None
            )
            samples.append(sample)
            self.logger.info('Processed {} (median depth {})'.format(sample.name, np.median(depth)))
        return samples


    def bams_to_training_samples(self, truth_bam, bam, region, reference=None, truth_haplotag=None):
        """Prepare training data chunks.

        :param truth_bam: .bam file of truth aligned to ref to generate labels.
        :param bam: input .bam file.
        :param region: `medaka.common.Region` obj.
            the reference will be parsed.
        :param reference: reference `.fasta`, should correspond to `bam`.
        :param truth_haplotag: two letter tag name used for grouping truth labels by haplotype.

        :returns: tuple of `medaka.common.Sample` objects.

        .. note:: Chunks might be missing if `truth_bam` is provided and
            regions with multiple mappings were encountered.

        """
        # pick function to get pairs for labels, should be modified for RLE
        aln_to_pairs = medaka.common.get_pairs

        # filter truth alignments to restrict ourselves to regions of the ref where the truth
        # in unambiguous
        alns = medaka.labels.TruthAlignment.bam_to_alignments(truth_bam, region, haplotag=truth_haplotag)
        if len(alns) == 0:
            self.logger.info("Filtering and grouping removed all alignments of truth to ref from {}.".format(region))

        samples = []
        pad = (medaka.common.encoding[medaka.common._gap_], 1)
        for aln in alns:
            # truth_labels should be shape (pos, ploidy) and dtype (base, run_length)
            truth_pos, truth_labels = medaka.labels.TruthAlignment.get_positions_and_labels(aln, aln_to_pairs)
            aln_samples = self.bam_to_sample(bam, medaka.common.Region(region.ref_name, aln[0].start, aln[0].end))
            ploidy = truth_labels.shape[-1]
            for sample in aln_samples:
                # Create labels according to positions in pileup
                padded_labels = np.empty((len(sample.positions), ploidy), dtype=truth_labels.dtype)
                # fill with pad so that insertions not present in labels have correct gap-label
                padded_labels.fill(pad)
                truth_inds = np.where(np.in1d(truth_pos, sample.positions))
                sample_inds = np.where(np.in1d(sample.positions, truth_pos))
                assert len(truth_inds[0]) == len(sample_inds[0])
                assert np.alltrue(truth_pos[truth_inds] == sample.positions[sample_inds])

                padded_labels[sample_inds] = truth_labels[truth_inds]

                sample = sample._asdict()
                sample['labels'] = padded_labels
                samples.append(medaka.common.Sample(**sample))
        return tuple(samples)


def alphabet_filter(sample_gen, alphabet=None, filter_labels=True, filter_ref_seq=True):
    """Skip chunks in which labels and/or ref_seq contain bases not in `alphabet`.

    :param sample_gen: generator of `medaka.common.Sample` named tuples.
    :param alphabet: set of str of allowed bases. If None, automatically generated from decoding.
    :param filter_labels: bool, whether to filter on labels.
    :param filter_ref_seq: bool, whether to filter on ref_seq.

    :yields: `medaka.common.Sample` named tuples.
    """
    if alphabet is None:
        alphabet = set([c for c in medaka.common._alphabet_ + medaka.common._gap_])
    logger = medaka.common.get_named_logger('AlphaFilter')
    logger.debug("alphabet: {}".format(alphabet))

    alphabet = set([medaka.common.encoding[c] for c in alphabet])

    def _find_bad_bases(s, field, alphabet):
        seq_rle = getattr(s, field)
        bases = set(np.unique(seq_rle['base']))
        if not bases.issubset(alphabet):
            diff = [medaka.common.decoding[i] for i in bases - alphabet]
            msg = "Skipping {}:{}-{} ({} bases) due to {} {}"
            pos = s.positions
            logger.info(msg.format(s.ref_name, pos['major'][0], pos['major'][-1],
                                   len(pos), field, diff))
            return True

    for s in sample_gen:
        if filter_labels and s.labels is not None and _find_bad_bases(s, 'labels', alphabet):
            continue
        if filter_ref_seq and s.ref_seq is not None and _find_bad_bases(s, 'ref_seq', alphabet):
            continue
        yield s


class SampleGenerator(object):

    def __init__(self, bam, region, model, truth_bam=None,
                 truth_haplotag=None, chunk_len=1000,
                 chunk_overlap=200, tag_name=None, tag_value=None,
                 tag_keep_missing=False, enable_chunking=True):
        """Generate chunked inference (or training) samples.

        :param bam: `.bam` containing alignments from which to generate samples.
        :param region: a `medaka.common.Region` for which to generate samples.
        :param model: a medaka model.
        :param truth_bam: a `.bam` containing alignment of truth sequence to
            `reference` sequence. Required only for creating training chunks.
        :param truth_haplotag: two letter tag name used for grouping truth labels by haplotype.
        :param reference: reference `.fasta`, should correspond to `bam`.
        :param tag_name: two letter tag name by which to filter reads.
        :param tag_value: integer value of tag for reads to keep.
        :param tag_keep_missing: whether to keep reads when tag is missing.
        :param enable_chunking: when yielding samples, do so in chunks.

        """
        self.logger = medaka.common.get_named_logger("Sampler")
        self.sample_type = "training" if truth_bam is not None else "consensus"
        self.logger.info("Initializing sampler for {} of region {}.".format(self.sample_type, region))

        #TODO: this will need changing when we switch to simply saving a class
        with medaka.datastore.DataStore(model) as ds:
            self.fencoder_args = ds.meta['medaka_features_kwargs']
        dtypes = ('',) if 'dtypes' not in self.fencoder_args else self.fencoder_args['dtypes']
        self.fencoder = FeatureEncoder(
            normalise=self.fencoder_args['normalise'], dtypes=dtypes,
            tag_name=tag_name, tag_value=tag_value, tag_keep_missing=tag_keep_missing,
            )

        self.bam = bam
        self.region = region
        self.model = model
        self.truth_bam = truth_bam
        self.truth_haplotag = truth_haplotag
        self.chunk_len = chunk_len
        self.chunk_overlap = chunk_overlap
        self.enable_chunking = enable_chunking
        self._source = None # the base data to be chunked
        self._quarantined = list() # samples which are shorter than chunk size

        #TODO: check reference has been given if model/feature encoder requires it


    def _fill_features(self):
        if self._source is None:
            self._quarantined = None
            t0 = now()
            if self.truth_bam is not None:
                self._source = self.fencoder.bams_to_training_samples(
                    self.truth_bam, self.bam, self.region,
                    truth_haplotag=self.truth_haplotag)
            else:
                self._source = self.fencoder.bam_to_sample(
                    self.bam, self.region)
            t1 = now()
            self.logger.info("Took {:.2f}s to make features.".format(t1-t0))


    def _quarantine_sample(self, sample):
        """Add sample name and pileup width to a list."""
        # Note: the below assumes we haven't split a pileup on minor positions.
        # This should be the case: chunking on minor positions only occurs for
        # larger regions.
        start, _ = sample.first_pos
        end, _ = sample.last_pos
        end += 1 # end exclusive
        self._quarantined.append((
            medaka.common.Region(sample.ref_name, start, end), sample.size
        ))


    @property
    def samples(self):
        """List of (possibly) chunked samples."""
        self._fill_features()
        self._quarantined = list()
        all_data = []
        for source in self._source:
            if source.is_empty:
                continue
            if not self.enable_chunking:
                chunks = [source]
            else:
                if source.size < self.chunk_len:
                    msg = "Region {} ({} positions) is smaller than inference chunk length {}, quarantining.".format(
                        source.name, source.size, self.chunk_len)
                    self.logger.warning(msg)
                    self._quarantine_sample(source)
                    continue
                self.logger.debug(
                    "Chunking pileup data into {} columns with overlap of {}.".format(
                    self.chunk_len, self.chunk_overlap))
                chunks = source.chunks(chunk_len=self.chunk_len, overlap=self.chunk_overlap)
            self.logger.info("Pileup for {} is of width {}".format(source.name, source.size))
            all_data.extend(alphabet_filter(chunks))
        return all_data


def create_samples(args):
    raise NotImplementedError('Creation of unlabelled samples is currently disabled')


def _labelled_samples_worker(args, region):
    logger = medaka.common.get_named_logger('PrepWork')
    logger.info("Processing region {}.".format(region))
    data_gen = SampleGenerator(
        args.bam, region, args.model, truth_bam=args.truth, truth_haplotag=args.truth_haplotag,
        chunk_len=args.chunk_len, chunk_overlap=args.chunk_ovlp)
    return list(data_gen.samples), region, deepcopy(data_gen.fencoder_args), deepcopy(data_gen.fencoder.decoding)


def create_labelled_samples(args):
    logger = medaka.common.get_named_logger('Prepare')
    if args.chunk_ovlp >= args.chunk_len:
        raise ValueError('chunk_ovlp {} is not smaller than chunk_len {}'.format(args.chunk_ovlp, args.chunk_len))
    regions = medaka.common.get_regions(args.bam, args.regions)
    reg_str = '\n'.join(['\t\t\t{}'.format(r) for r in regions])
    logger.info('Got regions:\n{}'.format(reg_str))

    labels_counter = Counter()

    no_data = False
    with medaka.datastore.DataStore(args.output, 'w') as ds:
        # write feature options to file
        logger.info("Writing meta data to file.")
        with medaka.datastore.DataStore(args.model) as model:
            meta = { k: model.meta[k] for k in ('medaka_features_kwargs', 'medaka_feature_decoding')}
        ds.update_meta(meta)
        # TODO: this parallelism would be better in `SampleGenerator.bams_to_training_samples`
        #       since training alignments are usually chunked.
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
            # break up overly long chunks
            MAX_SIZE= int(1e6)
            regions = itertools.chain(*(r.split(MAX_SIZE) for r in regions))
            futures = [executor.submit(_labelled_samples_worker, args, reg) for reg in regions]
            for fut in concurrent.futures.as_completed(futures):
                if fut.exception() is None:
                    samples, region, fencoder_args, fencoder_decoder = fut.result()
                    logger.info("Writing {} samples for region {}".format(len(samples), region))
                    for sample in samples:
                        ds.write_sample(sample)
                else:
                    logger.info(fut.exception())
                    logger.info(fut.result())
                fut._result = None  # python issue 27144
        no_data = ds.n_samples == 0

    if no_data:
        logger.critical("Warning: No training data was written to file, deleting output.")
        os.remove(args.output)


# read / reference fasta / fastq compression
_scores_ = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'


def get_runs_from_fastq(fastq, ref_name):
    """Get run length encoding of a read.
    :param fastq: str, path to fastq file or `SeqIO.index` object.
    :param ref_name: str, name of read with fastq.
    :returns: structured `numpy.ndarray` with dtype=[('start', int), ('length', int)].
    """
    if isinstance(fastq, str):
        fastq = SeqIO.index(fastq, 'fastq')
    read = fastq[ref_name]
    lengths = read.letter_annotations['phred_quality']
    return medaka.common.lengths_to_rle(lengths)


def compress_seq(read):
    """Compress homopolymers within a basecall, encoding their lengths in qscores.

    :param read: `Bio.SeqRecord.SeqRecord` object
    :returns: (str read.description,
               str compressed sequence,
               str qscores,
               structured array with fields `start`, `length`, and `value`).
    """
    seq_array = np.fromiter(read.seq, dtype='U1', count=len(read.seq))
    runs = medaka.common.rle(seq_array, low_mem=True)
    compressed_seq = ''.join(runs['value'])
    # we can only encode up to a homopolymer length 93
    # if we want to the score decoding in pysam to correspond to counts
    # (we could get to 94 if we make score zero correspond to count 1).
    is_too_long = runs['length'] >= len(_scores_)
    if np.any(is_too_long):
        logger.warning('Some homopolymers in {} are longer than the longest supported length\n'.format(read.name))
        inds = np.where(is_too_long)
        runs['length'][inds] = len(_scores_) - 1
    compressed_scores = ''.join([ _scores_[run['length']] for run in runs])
    return read.description, compressed_seq, compressed_scores, runs


def compress(args):
    if args.output is None:
        fh = sys.stdout
    else:
        fh = open(args.output, 'w')

    formats = {'a': 'fasta', 'q': 'fastq'}
    if args.input[-1] not in formats:
        msg='Could not guess file format of {}, rename to .f(ast)a/.f(ast)q'
        raise KeyError(msg.format(args.input))

    reads = SeqIO.parse(args.input, formats[args.input[-1]])
    if args.threads > 1:
        pool = Pool(args.threads)
        compressed = pool.imap(compress_seq, reads)
    else:
        compressed = (compress_seq(r) for r in reads)

    t0 = now()
    for description, compressed_seq, compressed_scores, runs in compressed:
        fh.write('@{}\n{}\n'.format(description, compressed_seq))
        fh.write('{}\n{}\n'.format('+', compressed_scores))
    t1 = now()
    logger.info('Compressing {} took {:.3f}s.'.format(args.input, t1 - t0))

    if args.output is not None:
        fh.close()
