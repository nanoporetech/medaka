from collections import defaultdict, Counter, OrderedDict
import concurrent.futures
from copy import deepcopy
import inspect
import itertools
from functools import partial
import logging
from multiprocessing import Pool
import sys
from timeit import default_timer as now

from Bio import SeqIO
import h5py
import numpy as np
import pysam

#TODO: from medaka import common
from medaka.common import (yield_compressed_pairs, Sample, lengths_to_rle, rle,
                           Region,
                           decoding, encoding, get_regions,
                           _gap_, _alphabet_, _feature_opt_path_, _feature_decoding_path_,
                           get_pairs, get_pairs_with_hp_len, seq_to_hp_lens,
                           write_yaml_data, load_yaml_data, gen_train_batch, serial_gen_train_batch,
                           _feature_batches_path_, _label_batches_path_, _label_counts_path_,
                           _label_decod_path_,)

from medaka.labels import TruthAlignment
import libmedaka


def pileup_counts(region, bam, dtype_prefixes=None):
    """Create pileup counts feature array for region.

    :param region: `Region` object
    :param bam: .bam file with alignments.
    :param dtype_prefixes: prefixes for query names which to separate counts.
        If `None` (or of length 1), counts are not split.

    :returns: pileup counts array, reference positions, insertion postions
    """
    ffi, lib = libmedaka.ffi, libmedaka.lib
    
    # htslib start is 1-based, Region object is 0-based
    region_str = '{}:{}-{}'.format(region.ref_name, region.start + 1, region.end)

    num_dtypes, dtypes = 1, ffi.NULL
    if isinstance(dtype_prefixes, str):
        dtype_prefixes = [dtype_prefixes]
    if dtype_prefixes is not None and len(dtype_prefixes) > 1:
        num_dtypes = len(dtype_prefixes)
        _dtypes = [ffi.new("char[]", d.encode()) for d in dtype_prefixes]
        dtypes = ffi.new("char *[]", _dtypes)

    counts = lib.calculate_pileup(
        region_str.encode(), bam.encode(), num_dtypes, dtypes
    )

    featlen = lib.featlen

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

    # get rid of 'first' counts row for each datatype (counts of alternative bases)
    mask = np.ones(np_counts.shape[1], dtype=bool)
    mask[[x * featlen for x in range(0, num_dtypes)]] = False
    np_counts = np_counts[:, mask]

    return np_counts, positions


class FeatureEncoder(object):
    _ref_modes_ = ['onehot', 'base_length', 'index', None]
    _norm_modes_ = ['total', 'fwd_rev', None]


    def __init__(self, ref_mode:str=None, max_hp_len:int=10,
                 log_min:int=None, normalise:str='total', with_depth:bool=False,
                 consensus_as_ref:bool=False, is_compressed:bool=True, dtypes=('',)):
        """Class to support multiple feature encodings

        :param ref_mode: str, how to represent the reference.
        :param max_hp_len: int, longest homopolymer run which can be represented, longer runs will be truncated.
        :param log_min: int, take log10 of counts/fractions and set zeros to 10**log_min.
        :param normalise: str, how to normalise the data.
        :param with_depth: bool, whether to include a feature describing the total depth.
        :param consensus_as_ref: bool, whether to use a naive max-count consensus instead of the reference.
        :param is_compressed: bool, whether to use HP compression. If false, treat as uncompressed.
        :param dtypes: iterable of str, read id prefixes of distinct data types that should be counted separately.

        """
        self.ref_mode = ref_mode
        self.consensus_as_ref = consensus_as_ref
        self.max_hp_len = max_hp_len
        self.log_min = log_min
        self.normalise = normalise
        self.feature_dtype = np.float32 if (self.normalise is not None or self.log_min is not None) else np.uint64
        self.with_depth = with_depth
        self.is_compressed = is_compressed
        self.logger = logging.getLogger(__package__)
        self.logger.name = 'Feature'
        self.dtypes = dtypes

        if self.ref_mode not in self._ref_modes_:
            raise ValueError('ref_mode={} is not one of {}'.format(self.ref_mode, self._ref_modes_))
        if self.normalise not in self._norm_modes_:
            raise ValueError('normalise={} is not one of {}'.format(self.normalise, self._norm_modes_))

        opts = inspect.signature(FeatureEncoder.__init__).parameters.keys()
        opts = {k:getattr(self, k) for k in opts if k != 'self'}

        read_decoding = []
        for dtype in self.dtypes:
            # set up one-hot encoding of read run lengths for each dtype
            read_decoding += [
                (dtype,) + k
                for k in itertools.product(
                    (True, False), _alphabet_, range(1, max_hp_len + 1)
                )
            ]

            # forward and reverse gaps
            read_decoding += [(dtype, True, None, 1), (dtype, False, None, 1)]

        if self.ref_mode == 'onehot':
            ref_decoding = [('ref', b, l) for b, l in itertools.product(
                alphabet, range(1, max_hp_len + 1))]
            ref_decoding.append(('ref', _gap_, 1))  # gaps
        elif self.ref_mode == 'base_length':
            ref_decoding = ['ref_base', 'ref_length']
            self.ref_base_encoding = {b:i for i, b in enumerate(_alphabet_ + _gap_)}
        elif self.ref_mode == 'index':
            ref_decoding = ['ref_index']
        else:
            ref_decoding = []

        self.decoding = tuple(read_decoding + ref_decoding)
        if self.with_depth:
            self.decoding = self.decoding + ('depth',)
        self.encoding = OrderedDict(((a, i) for i, a in enumerate(self.decoding)))
        self.logger.info("Creating features with: {}".format(opts))

        self.logger.debug("Label decoding is:\n{}".format('\n'.join(
            '{}: {}'.format(i, x) for i, x in enumerate(self.decoding)
        )))


    def process_ref_seq(self, ref_name, ref_fq):
        """Translate a sequence to RLE or truncated homopolymer form if required.

        :param ref_name: str, name of contig within `ref_fq`.
        :param ref_fq: path to reference fastq, or `SeqIO.index` obj.
        """
        if self.is_compressed:  # get (compressed) rle_encoded HP len from fq
            if ref_fq is None:
                raise ValueError('If homopolymers have been compressed, ref_fq must be provided')
            return get_runs_from_fastq(ref_fq, ref_name)
        elif self.max_hp_len > 1:  # get (uncompressed) hp_lens from fq
            if ref_fq is None:
                raise ValueError('If max_hp_len > 1, ref_fq must be provided')
            if isinstance(ref_fq, str):
                ref_fq = SeqIO.index(ref_fq, 'fastq')
            read = ref_fq[ref_name]
            return seq_to_hp_lens(read.seq)
        else:
            return None


    @property
    def feature_indices(self):
        """Location of feature vector components.
        
        :returns: dictionary mapping read data type and strand to feature vector
        indices (a list over all bases)

        """
        return {
            (dt, strand):
                [self.encoding[k] for k in self.encoding.keys()
                    if k[0] == dt and strand == k[1]]
                for dt, strand in itertools.product(self.dtypes, (True, False))
        }


    def bam_to_sample_c(self, reads_bam, region):
        """Converts a section of an alignment pileup (as shown
        by e.g. samtools tview) to a base frequency feature array

        :param reads_bam: (sorted indexed) bam with read alignment to reference
        :param region: `Region` object with ref_name, start and end attributes.
        :param start: starting position within reference
        :param end: ending position within reference
        :returns: `Sample` object
        """
        assert self.ref_mode is None
        assert not self.consensus_as_ref
        assert self.max_hp_len == 1
        assert self.log_min is None
        assert self.normalise == 'total' or self.normalise == 'fwd_rev' or self.normalise is None
        assert not self.with_depth
        assert not self.is_compressed

        counts, positions = pileup_counts(region, reads_bam, dtype_prefixes=self.dtypes)
        if len(counts) == 0:
            msg = 'Pileup-feature is zero-length for {} indicating no reads in this region.'.format(region)
            self.logger.warning(msg)
            return Sample(ref_name=region.ref_name, features=None,
                          labels=None, ref_seq=None,
                          positions=positions, label_probs=None)

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

        sample = Sample(ref_name=region.ref_name, features=feature_array,
                        labels=None, ref_seq=None,
                        positions=positions, label_probs=None)

        self.logger.info('Processed {} (median depth {})'.format(sample.name, np.median(depth)))
        return sample


    def bam_to_sample(self, reads_bam, region, reference=None, read_fraction=None, force_py=False):
        """Converts a section of an alignment pileup (as shown
        by e.g. samtools tview) to a base frequency feature array

        :param reads_bam: (sorted indexed) bam with read alignment to reference
        :param region: `Region` object with ref_name, start and end attributes.
        :param reference: reference `.fasta`, should correspond to `bam`.
            Required only for run length encoded references and reads.
        :param read_fraction: fraction of reads to use, if `None` use all.
        :param force_py: bool, if True, force use of python code (rather than c library).
        :returns: `Sample` object
        """

        ref_rle = self.process_ref_seq(region.ref_name, reference)

        # Try to use fast c function if we can, else fall back on this function
        if not force_py and (ref_rle is None and read_fraction is None):
            try:
                return self.bam_to_sample_c(reads_bam, region)
            except Exception as e:
                self.logger.info('Could not process sample with bam_to_sample_c, using python code instead.\n({}).'.format(e))
                pass

        if self.is_compressed:
            aln_to_pairs = partial(yield_compressed_pairs, ref_rle=ref_rle)
        elif self.max_hp_len == 1:
            aln_to_pairs = get_pairs
        else:
            aln_to_pairs = partial(get_pairs_with_hp_len, ref_seq=ref_rle)

        # accumulate data in dicts
        aln_counters = defaultdict(Counter)
        ref_bases = dict()
        with pysam.AlignmentFile(reads_bam, 'rb') as bamfile:
            aln_reads = bamfile.fetch(region.ref_name, region.start, region.end)
            if read_fraction is not None:
                low, high = read_fraction
                np.random.seed((int(now()) * region.start) % 2**32)
                fraction = ((high - low) * np.random.random_sample(1) + low)[0]
                aln_reads = [a for a in aln_reads]
                n_reads = len(aln_reads)
                n_reads_to_keep = max(int(fraction * n_reads), 1)
                replace = n_reads_to_keep > n_reads
                msg = "Resampling (replace {}) from {} to {} ({:.3f}) for {}"
                self.logger.debug(msg.format(replace, n_reads, n_reads_to_keep, fraction, region))
                aln_reads = np.random.choice(aln_reads, n_reads_to_keep, replace=replace)

            start = region.start
            end = region.end
            if start is None:
                start = 0
            if end is None:
                end = float('Inf')

            for aln in aln_reads:
                # get the dtype from the prefix of the query name
                try:
                    dtype = self.dtypes[np.where([aln.query_name.startswith(dt) for dt in self.dtypes])[0][0]]
                except:
                    msg = "Skipping read {} as dtype not in {}"
                    logging.info(msg.format(aln.query_name, self.dtypes))
                    continue

                reverse = aln.is_reverse
                pairs = aln_to_pairs(aln)
                ins_count = 0
                for pair in itertools.dropwhile(lambda x: (x.rpos is None)
                                                or (x.rpos < start), pairs):
                    if ((pair.rpos == aln.reference_end - 1) or
                        (pair.rpos is not None and pair.rpos >= end)):
                        break
                    if pair.rpos is None:
                        ins_count += 1
                    else:
                        ins_count = 0
                        current_pos = pair.rpos

                    (aln_counters[(current_pos, ins_count)]
                    [self.encoding[dtype, reverse, pair.qbase, min(pair.qlen, self.max_hp_len)]]) += 1

                    ref_base = pair.rbase.upper() if pair.rbase is not None else '*'
                    ref_bases[(current_pos, ins_count)] = (ref_base, pair.rlen)

            # create feature array
            aln_cols = len(aln_counters)
            if aln_cols == 0:
                msg = 'Pileup-feature is zero-length for {} indicating no reads in this region.'.format(region)
                self.logger.warning(msg)
                return Sample(ref_name=region.ref_name, features=None,
                              labels=None, ref_seq=None,
                              positions=positions, label_probs=None)

            feature_len = len(self.encoding)
            feature_array = np.zeros(shape=(aln_cols, feature_len), dtype=self.feature_dtype)
            if self.log_min is not None:
                feature_array.fill(np.nan)
            ref_array = np.empty(shape=(aln_cols), dtype=[('base', int), ('run_length', int)])
            positions = np.empty(aln_cols, dtype=[('major', int),
                                                ('minor', int)])
            depth_array = np.empty(shape=(aln_cols), dtype=int)

            # keep track of which features are for fwd/rev reads of each dtype
            inds_by_type = self.feature_indices

            #TODO: refactor so common combinations of options can be handled as in C-function
            for i, ((pos, counts), (_, (ref_base, ref_len))) in \
                    enumerate(zip(sorted(aln_counters.items()),
                                sorted(ref_bases.items()))):
                positions[i] = pos
                ref_array[i] = (encoding[ref_base], ref_len)
                for j in counts.keys():
                    feature_array[i, j] = counts[j]

                if self.consensus_as_ref:
                    cons_i = np.argmax(feature_array[i])
                    cons_is_reverse, cons_base, cons_length = self.decoding[cons_i]
                    ref_base = cons_base if cons_base is not None else _gap_
                    ref_len = cons_length

                if positions[i]['minor'] == 0:
                    major_depth = sum(counts.values())
                    # get the depth of each fwd and rev dtype
                    major_depths_by_type = {t: sum((counts[i] for i in inds_by_type[t])) for t in inds_by_type}
                    assert sum(major_depths_by_type.values()) == major_depth

                if self.normalise is not None:
                    if self.normalise == 'total':
                        feature_array[i, :] /= max(major_depth, 1)
                    elif self.normalise == 'fwd_rev':
                        # normalize fwd and reverse seperately for each dtype
                        for dt, inds in inds_by_type.items():
                            feature_array[i, inds] /= max(major_depths_by_type[dt], 1)

                depth_array[i] = major_depth

                if self.with_depth:
                    feature_array[i, self.encoding['depth']] = depth_array[i]

                if self.log_min is not None:  # counts/proportions and depth will be normalised
                    # when we take log of probs, make it easier for network by keeping all log of
                    # probs positive. add self.log_min to any log probs so they are positive,
                    # if self.log_min is 10, we can cope with depth up to 10**9
                    feature_array[i, :] = np.log10(feature_array[i, :],
                                                   out=feature_array[i, :])
                    feature_array[i, :] += self.log_min
                    feature_array[i, :] = np.nan_to_num(feature_array[i, :], copy=False)

                if self.ref_mode == 'onehot':
                    feature_array[i, self.encoding[('ref', str(ref_base), int(ref_len))]] = 1
                elif self.ref_mode == 'base_length':
                    feature_array[i, self.encoding['ref_base']] = self.ref_base_encoding[ref_base]
                    feature_array[i, self.encoding['ref_length']] = ref_len
                elif self.ref_mode == 'index':
                    # index of count which ref would contribute to were it a read
                    feature_array[i, self.encoding['ref_index']] = self.encoding[(False, min(ref_len, self.max_hp_len), ref_base)]

            sample = Sample(ref_name=region.ref_name, features=feature_array,
                            labels=None, ref_seq=ref_array,
                            positions=positions, label_probs=None)
            self.logger.info('Processed {} (median depth {})'.format(sample.name, np.median(depth_array)))

            return sample


    def bams_to_training_samples(self, truth_bam, bam, region, reference=None, read_fraction=None):
        """Prepare training data chunks.

        :param truth_bam: .bam file of truth aligned to ref to generate labels.
        :param bam: input .bam file.
        :param region: `Region` obj.
            the reference will be parsed.
        :param reference: reference `.fasta`, should correspond to `bam`.

        :returns: tuple of `Sample` objects (one item for each input bam) for
            each chunk.

        .. note:: Chunks might be missing if `truth_bam` is provided and
            regions with multiple mappings were encountered.

        """
        ref_rle = self.process_ref_seq(region.ref_name, reference)

        # filter truth alignments to restrict ourselves to regions of the ref where the truth
        # in unambiguous
        alignments = TruthAlignment.bam_to_alignments(truth_bam, region.ref_name, start=region.start, end=region.end)
        filtered_alignments = TruthAlignment.filter_alignments(alignments, start=region.start, end=region.end)
        if len(filtered_alignments) == 0:
            self.logger.info("Filtering removed all alignments of truth to ref from {}.".format(region))

        samples = []
        for aln in filtered_alignments:
            mock_compr = self.max_hp_len > 1 and not self.is_compressed
            truth_pos, truth_labels = aln.get_positions_and_labels(ref_compr_rle=ref_rle, mock_compr=mock_compr,
                                                                   is_compressed=self.is_compressed, rle_dtype=True)
            sample = self.bam_to_sample(bam, Region(region.ref_name, aln.start, aln.end), ref_rle, read_fraction=read_fraction)
            # Create labels according to positions in pileup
            pad = (encoding[_gap_], 1) if len(truth_labels.dtype) > 0 else encoding[_gap_]
            padder = itertools.repeat(pad)
            position_to_label = defaultdict(padder.__next__,
                                            zip([tuple(p) for p in truth_pos],
                                                [a for a in truth_labels]))
            padded_labels = np.fromiter((position_to_label[tuple(p)] for p in sample.positions),
                                            dtype=truth_labels.dtype, count=len(sample.positions))

            sample = sample._asdict()
            sample['labels'] = padded_labels
            samples.append(Sample(**sample))
        return tuple(samples)


def alphabet_filter(sample_gen, alphabet=None, filter_labels=True, filter_ref_seq=True):
    """Skip chunks in which labels and/or ref_seq contain bases not in `alphabet`.

    :param sample_gen: generator of `Sample` named tuples.
    :param alphabet: set of str of allowed bases. If None, automatically generated from decoding.
    :param filter_labels: bool, whether to filter on labels.
    :param filter_ref_seq: bool, whether to filter on ref_seq.

    :yields: `Sample` named tuples.
    """
    if alphabet is None:
        alphabet = set([c for c in _alphabet_ + _gap_])
    logger = logging.getLogger('AlphaFilter')
    logger.debug("alphabet: {}".format(alphabet))

    alphabet = set([encoding[c] for c in alphabet])

    def _find_bad_bases(s, field, alphabet):
        seq_rle = getattr(s, field)
        bases = set(np.unique(seq_rle['base']))
        if not bases.issubset(alphabet):
            diff = [decoding[i] for i in bases - alphabet]
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

    def __init__(self, bam, region, model, rle_ref=None, truth_bam=None, read_fraction=None, chunk_len=1000, chunk_overlap=200):
        """Generate chunked inference (or training) samples.

        :param bam: `.bam` containing alignments from which to generate samples.
        :param region: a `Region` for which to generate samples.
        :param model: a medaka model.
        :param truth_bam: a `.bam` containing alignment of truth sequence to
            `reference` sequence. Required only for creating training chunks.
        :param reference: reference `.fasta`, should correspond to `bam`.

        """
        self.logger = logging.getLogger("Sampler")
        self.sample_type = "training" if truth_bam is not None else "consensus"
        self.logger.info("Initializing sampler for {} or region {}.".format(self.sample_type, region))
        self.fencoder_args = load_yaml_data(model, _feature_opt_path_)
        self.fencoder = FeatureEncoder(**self.fencoder_args)

        self.bam = bam
        self.region = region
        self.model = model
        self.rle_ref = rle_ref
        self.truth_bam = truth_bam
        self.read_fraction = read_fraction
        self.chunk_len = chunk_len
        self.chunk_overlap = chunk_overlap
        self._source = None # the base data to be chunked

        #TODO: check reference has been given if model/feature encoder requires it


    def _fill_features(self):
        if self._source is None:
            t0 = now()
            if self.truth_bam is not None:
                self._source = self.fencoder.bams_to_training_samples(
                    self.truth_bam, self.bam, self.region, self.rle_ref,
                    self.read_fraction)
            else:
                self._source = self.fencoder.bam_to_sample(
                    self.bam, self.region, self.rle_ref, self.read_fraction)
                self._source = (self._source,) # wrap to be the same as above
            t1 = now()
            self.logger.info("Took {:.2f}s to make features.".format(t1-t0))

    @property
    def n_samples(self):
        """The approximate number of samples that will be yielded by `.samples`."""
        self._fill_features()
        return 1 + sum(s.size // (self.chunk_len - self.chunk_overlap) for s in self._source)


    @property
    def samples(self):
        """Iterator over chunked samples."""
        self._fill_features()
        for source in self._source:
            if source.is_empty:
                continue
            if source.size < self.chunk_len:
                msg = "Region {} ({} positions) is smaller than inference chunk length {}.".format(
                    source.name, source.size, self.chunk_len)
                self.logger.warning(msg)
                continue

            self.logger.debug(
                "Chunking pileup data into {} columns with overlap of {}.".format(
                self.chunk_len, self.chunk_overlap))
            chunks = source.chunks(chunk_len=self.chunk_len, overlap=self.chunk_overlap)
            yield from alphabet_filter(chunks)


    def training_samples(self, max_label_len):
        """Iterator of (feature, label) pairs for training."""
        self.logger.info("Maxlabellen: {}".format(max_label_len))
        if self.truth_bam is None:
            raise ValueError("Cannot iterate over training pairs when truth bam has not been given.""")
        label_encoding, label_decoding = get_label_encoding(max_label_len)
        for s in self.samples:
            if s.labels is None: # this shouldn't happen
                raise ValueError("Cannot train without labels.")
            x = s.features
            # labels can either be unicode strings or (base, length) integer tuples
            if isinstance(s.labels[0], np.unicode):
                y = np.fromiter(
                    (label_encoding[l[:min(max_label_len, len(l))]] for l in s.labels),
                    dtype=int, count=len(s.labels))
            else:
                y = np.fromiter(
                    (label_encoding[tuple((l['base'], min(max_label_len, l['run_length'])))]
                        for l in s.labels),
                    dtype=int, count=len(s.labels))
            y = y.reshape(y.shape + (1,))
            yield x, y


def get_label_encoding(max_label_len):
    """Get label encodings for a given maximum label length.

    :param max_label_len: int, maximum label length.

    :returns: (label_encoding, label_decoding_strs)
        label_encoding: {(int encoded base, int run length): int label encoding}.
        label_decoding_strs: list of str of label decodings.

    >>> get_label_encoding(2)
    ({(5, 1): 0,
      (5, 2): 1,
      (6, 1): 2,
      (6, 2): 3,
      (7, 1): 4,
      (7, 2): 5,
      (8, 1): 6,
      (8, 2): 7,
      (0, 1): 8,
      },
      ['A', 'AA', 'C', 'CC', 'G', 'GG', 'T', 'TT', '*']
    )

    """
    encoded_bases = [encoding[b] for b in _alphabet_]
    label_decoding = [(b, l) for b, l in itertools.product(
        encoded_bases, range(1, max_label_len + 1))]
    label_decoding.append((encoding[_gap_], 1))  # gaps
    label_encoding = {t: i for i, t in enumerate(label_decoding)}
    label_decoding_strs = [l * decoding[b] for (b,l) in label_decoding]
    return label_encoding, label_decoding_strs


def create_samples(args):
    raise NotImplementedError('Creation of unlabelled samples is currently disabled')


def _labelled_samples_worker(args, region):
    logger = logging.getLogger('PrepWork')
    logger.info("Processing region {}.".format(region))
    fencoder_args, fencoder_decoder = None, None
    data_gen = SampleGenerator(
        args.bam, region, args.model, args.rle_ref, truth_bam = args.truth,
        read_fraction=args.read_fraction, chunk_len=args.chunk_len, chunk_overlap=args.chunk_ovlp)
    if fencoder_args is None:
        fencoder_args = deepcopy(data_gen.fencoder_args)
        fencoder_decoder = deepcopy(data_gen.fencoder.decoding)
    batches = serial_gen_train_batch(
        data_gen.training_samples(args.max_label_len),
        args.batch_size)
    #TODO: shuffle
    return list(batches)


def create_labelled_samples(args):
    logger = logging.getLogger('Prepare')
    regions = get_regions(args.bam, args.regions)
    reg_str = '\n'.join(['\t\t\t{}'.format(r) for r in regions])
    logger.info('Got regions:\n{}'.format(reg_str))


    labels_counter = Counter()
    with h5py.File(args.output, 'w') as hdf:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            worker = partial(_labelled_samples_worker, args)
            futures = executor.map(worker, regions)
            for fut in concurrent.futures.as_completed(futures):
                if fut.exception is None:
                    batches = fut.result()
                    for i, (x_batch, y_batch) in enumerate(batches):
                        unique, counts = np.lib.arraysetops.unique(y_batch, return_counts=True)
                        labels_counter.update(dict(zip(unique, counts)))
                        hdf['{}/{}_{}'.format(_feature_batches_path_, region, i)] = x_batch
                        hdf['{}/{}_{}'.format(_label_batches_path_, region, i)] = y_batch

    # write feature options to file, so we can later check model and features
    # are compatible and label counts, so we can weight labels by inverse counts.
    logger.info("Writing meta data to file.")
    label_encoding, label_decoding = get_label_encoding(args.max_label_len)
    to_save = {_feature_opt_path_: fencoder_args,
               _feature_decoding_path_: fencoder_decoder,
               _label_counts_path_: labels_counter,
               _label_decod_path_: label_decoding}
    for fname in (args.output, args.output + '.yml'):
        write_yaml_data(fname, to_save)


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
    return lengths_to_rle(lengths)


def compress_seq(read):
    """Compress homopolymers within a basecall, encoding their lengths in qscores.

    :param read: `Bio.SeqRecord.SeqRecord` object
    :returns: (str read.description,
               str compressed sequence,
               str qscores,
               structured array with fields `start`, `length`, and `value`).
    """
    seq_array = np.fromiter(read.seq, dtype='U1', count=len(read.seq))
    runs = rle(seq_array, low_mem=True)
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

