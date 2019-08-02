from collections import defaultdict, Counter, OrderedDict
import concurrent.futures
from copy import deepcopy
import inspect
import itertools
from functools import partial
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


def pileup_counts(region, bam, dtype_prefixes=None, region_split=100000, workers=12, tag_name=None, tag_value=None, keep_missing=False):
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

    num_dtypes, dtypes = 1, ffi.NULL
    if isinstance(dtype_prefixes, str):
        dtype_prefixes = [dtype_prefixes]
    if dtype_prefixes is not None and len(dtype_prefixes) > 1:
        num_dtypes = len(dtype_prefixes)
        _dtypes = [ffi.new("char[]", d.encode()) for d in dtype_prefixes]
        dtypes = ffi.new("char *[]", _dtypes)
    if tag_name is None:
        tag_name = ffi.new("char[2]", "".encode())
        tag_value = 0
        keep_missing = False
    else:
        if len(tag_name) > 2:
            raise ValueError("'tag_value' must be a length-2 string.")
        tag_name = ffi.new("char[2]", tag_name.encode())

    featlen = lib.featlen

    def _process_region(reg):
        # htslib start is 1-based, medaka.common.Region object is 0-based
        region_str = '{}:{}-{}'.format(reg.ref_name, reg.start + 1, reg.end)

        counts = lib.calculate_pileup(
            region_str.encode(), bam.encode(), num_dtypes, dtypes,
            tag_name, tag_value, keep_missing
        )

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
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(_process_region, regions))

    # First pass: need to check for discontinuities within chunks,
    # these show up as >2 changes in the major coordinate
    _results = list()
    for counts, positions in results:
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
        # get rid of 'first' counts row for each datatype (counts of
        # alternative bases)
        mask = np.ones(chunk_counts.shape[1], dtype=bool)
        mask[[x * featlen for x in range(0, num_dtypes)]] = False
        chunk_counts = chunk_counts[:, mask]
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
    _ref_modes_ = ['onehot', 'base_length', 'index', None]
    _norm_modes_ = ['total', 'fwd_rev', None]


    def __init__(self, ref_mode:str=None, max_hp_len:int=10,
                 log_min:int=None, normalise:str='total', with_depth:bool=False,
                 consensus_as_ref:bool=False, is_compressed:bool=True, dtypes=('',),
                 tag_name=None, tag_value=None, tag_keep_missing=False, sym_indels=False,
                 ):
        """Class to support multiple feature encodings

        :param ref_mode: str, how to represent the reference.
        :param max_hp_len: int, longest homopolymer run which can be represented, longer runs will be truncated.
        :param log_min: int, take log10 of counts/fractions and set zeros to 10**log_min.
        :param normalise: str, how to normalise the data.
        :param with_depth: bool, whether to include a feature describing the total depth.
        :param consensus_as_ref: bool, whether to use a naive max-count consensus instead of the reference.
        :param is_compressed: bool, whether to use HP compression. If false, treat as uncompressed.
        :param dtypes: iterable of str, read id prefixes of distinct data types that should be counted separately.
        :param tag_name: two letter tag name by which to filter reads.
        :param tag_value: integer value of tag for reads to keep.
        :param tag_keep_missing: whether to keep reads when tag is missing.
        :param sym_indels: bool, whether to count a lack of an insertion as a deletion.

        """
        self.ref_mode = ref_mode
        self.consensus_as_ref = consensus_as_ref
        self.max_hp_len = max_hp_len
        self.log_min = log_min
        self.normalise = normalise
        self.feature_dtype = np.float32 if (self.normalise is not None or self.log_min is not None) else np.uint64
        self.with_depth = with_depth
        self.is_compressed = is_compressed
        self.logger = medaka.common.get_named_logger('Feature')
        self.dtypes = dtypes
        self.tag_name = tag_name
        self.tag_value = tag_value
        self.tag_keep_missing = tag_keep_missing
        self.sym_indels = sym_indels

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
                    (True, False), medaka.common._alphabet_, range(1, max_hp_len + 1)
                )
            ]

            # forward and reverse gaps
            read_decoding += [(dtype, True, None, 1), (dtype, False, None, 1)]

        if self.ref_mode == 'onehot':
            ref_decoding = [('ref', b, l) for b, l in itertools.product(
                alphabet, range(1, max_hp_len + 1))]
            ref_decoding.append(('ref', medaka.common._gap_, 1))  # gaps
        elif self.ref_mode == 'base_length':
            ref_decoding = ['ref_base', 'ref_length']
            self.ref_base_encoding = {b:i for i, b in enumerate(medaka.common._alphabet_ + medaka.common._gap_)}
        elif self.ref_mode == 'index':
            ref_decoding = ['ref_index']
        else:
            ref_decoding = []

        self.decoding = tuple(read_decoding + ref_decoding)
        if self.with_depth:
            self.decoding = self.decoding + ('depth',)
        self.encoding = OrderedDict(((a, i) for i, a in enumerate(self.decoding)))
        self.logger.debug("Creating features with: {}".format(opts))

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
            return medaka.common.seq_to_hp_lens(read.seq)
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
        :param region: `medaka.common.Region` object with ref_name, start and end attributes.
        :param start: starting position within reference
        :param end: ending position within reference
        :returns: `medaka.common.Sample` object
        """
        assert self.ref_mode is None
        assert not self.consensus_as_ref
        assert self.max_hp_len == 1
        assert self.log_min is None
        assert self.normalise == 'total' or self.normalise == 'fwd_rev' or self.normalise is None
        assert not self.with_depth
        assert not self.is_compressed

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
                    medaka.common.Sample(ref_name=region.ref_name, features=None,
                           labels=None, ref_seq=None,
                           postions=positions, label_probs=None
                    ))
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

            if self.sym_indels:
                # make indels at ref and non-ref positions symmetric.
                # major columns otherwise have counts of reads with and without a
                # deletion, whilst minor (inserted) columns only have counts of
                # the reads with an isertion.
                # To make ref and non-ref positions symmetric, fill in counts of reads which don't have insertions
                # i.e. depth_del = depth_major - depth_ins
                for (dt, is_rev), inds in self.feature_indices.items():
                    dt_depth = np.sum(counts[:, inds], axis=1)
                    del_ind = self.encoding[(dt, is_rev, None, 1)]
                    counts[minor_inds, del_ind] = dt_depth[major_ind_at_minor_inds] - dt_depth[minor_inds]
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


    def bam_to_sample(self, reads_bam, region, reference=None, read_fraction=None, force_py=False):
        """Converts a section of an alignment pileup (as shown
        by e.g. samtools tview) to a base frequency feature array

        :param reads_bam: (sorted indexed) bam with read alignment to reference
        :param region: `medaka.common.Region` object with ref_name, start and end attributes.
        :param reference: reference `.fasta`, should correspond to `bam`.
            Required only for run length encoded references and reads.
        :param read_fraction: fraction of reads to use, if `None` use all.
        :param force_py: bool, if True, force use of python code (rather than c library).
        :returns: iterable of `medaka.common.Sample` objects
        """

        ref_rle = self.process_ref_seq(region.ref_name, reference)

        # Try to use fast c function if we can, else fall back on this function
        if not force_py and (ref_rle is None and read_fraction is None):
            try:
                return self.bam_to_sample_c(reads_bam, region)
            except Exception as e:
                self.logger.info('Could not process sample with bam_to_sample_c, using python code instead.\n({}).'.format(e))
                pass
        if self.tag_name is not None:
            raise NotImplementedError("Filtering alignments by tag is not supported in python code.")

        #TODO: The code below will abut discontiguous regions in a pileup i.e.
        #      where no reads span a reference position the major position
        #      is dropped from the pileup. The correct behaviour would be to
        #      split apart the sub-regions and return them separately.
        #      The C implementation does this splitting.

        if self.is_compressed:
            aln_to_pairs = partial(medaka.common.yield_compressed_pairs, ref_rle=ref_rle)
        elif self.max_hp_len == 1:
            aln_to_pairs = medaka.common.get_pairs
        else:
            aln_to_pairs = partial(medaka.common.get_pairs_with_hp_len, ref_seq=ref_rle)

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
                    self.logger.info(msg.format(aln.query_name, self.dtypes))
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
            feature_len = len(self.encoding)
            feature_array = np.zeros(shape=(aln_cols, feature_len), dtype=self.feature_dtype)
            if self.log_min is not None:
                feature_array.fill(np.nan)
            ref_array = np.empty(shape=(aln_cols), dtype=[('base', int), ('run_length', int)])
            positions = np.empty(aln_cols, dtype=[('major', int),
                                                ('minor', int)])

            if aln_cols == 0:
                msg = 'Pileup-feature is zero-length for {} indicating no reads in this region.'.format(region)
                self.logger.warning(msg)
                return [medaka.common.Sample(ref_name=region.ref_name, features=None,
                              labels=None, ref_seq=None,
                              positions=positions, label_probs=None)]

            depth_array = np.empty(shape=(aln_cols), dtype=int)

            # keep track of which features are for fwd/rev reads of each dtype
            inds_by_type = self.feature_indices

            #TODO: refactor so common combinations of options can be handled as in C-function
            for i, ((pos, counts), (_, (ref_base, ref_len))) in \
                    enumerate(zip(sorted(aln_counters.items()),
                                sorted(ref_bases.items()))):
                positions[i] = pos
                ref_array[i] = (medaka.common.encoding[ref_base], ref_len)
                for j in counts.keys():
                    feature_array[i, j] = counts[j]

                if self.consensus_as_ref:
                    cons_i = np.argmax(feature_array[i])
                    cons_is_reverse, cons_base, cons_length = self.decoding[cons_i]
                    ref_base = cons_base if cons_base is not None else medaka.common._gap_
                    ref_len = cons_length

                if positions[i]['minor'] == 0:
                    major_depth = sum(counts.values())
                    # get the depth of each fwd and rev dtype
                    major_depths_by_type = {t: sum((counts[i] for i in inds_by_type[t])) for t in inds_by_type}
                    assert sum(major_depths_by_type.values()) == major_depth

                if self.sym_indels and positions[i]['minor'] > 0:
                    # make indels at ref and non-ref positions symmetric (see comment in bam_to_sample_c).
                    for (dtype, is_rev), inds in inds_by_type.items():
                        del_ind = self.encoding[(dtype, is_rev, None, 1)]
                        assert feature_array[i, del_ind] == 0
                        feature_array[i, del_ind] = major_depths_by_type[(dtype, is_rev)] - feature_array[i, inds].sum()

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

            sample = medaka.common.Sample(ref_name=region.ref_name, features=feature_array,
                            labels=None, ref_seq=ref_array,
                            positions=positions, label_probs=None)
            self.logger.info('Processed {} (median depth {})'.format(sample.name, np.median(depth_array)))

            return [sample]


    def bams_to_training_samples(self, truth_bam, bam, region, reference=None, read_fraction=None, truth_haplotag=None):
        """Prepare training data chunks.

        :param truth_bam: .bam file of truth aligned to ref to generate labels.
        :param bam: input .bam file.
        :param region: `medaka.common.Region` obj.
            the reference will be parsed.
        :param reference: reference `.fasta`, should correspond to `bam`.
        :param read_fraction: fraction of reads to use, if `None` use all.
        :param truth_haplotag: two letter tag name used for grouping truth labels by haplotype.

        :returns: tuple of `medaka.common.Sample` objects.

        .. note:: Chunks might be missing if `truth_bam` is provided and
            regions with multiple mappings were encountered.

        """
        ref_rle = self.process_ref_seq(region.ref_name, reference)

        # pick function to get pairs for labels
        mock_compr = self.max_hp_len > 1 and not self.is_compressed

        if self.is_compressed:
            aln_to_pairs = medaka.common.yield_compressed_pairs
        elif mock_compr:
            aln_to_pairs = medaka.common.get_pairs_with_hp_len
        else:
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
            aln_samples = self.bam_to_sample(bam, medaka.common.Region(region.ref_name, aln[0].start, aln[0].end),
                                             ref_rle, read_fraction=read_fraction)
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

    def __init__(self, bam, region, model, rle_ref=None, truth_bam=None,
                 truth_haplotag=None, read_fraction=None, chunk_len=1000,
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
        with medaka.datastore.DataStore(model) as ds:
            self.fencoder_args = ds.meta['medaka_features_kwargs']
        self.fencoder = FeatureEncoder(
            tag_name=tag_name, tag_value=tag_value, tag_keep_missing=tag_keep_missing,
            **self.fencoder_args)

        self.bam = bam
        self.region = region
        self.model = model
        self.rle_ref = rle_ref
        self.truth_bam = truth_bam
        self.truth_haplotag = truth_haplotag
        self.read_fraction = read_fraction
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
                    self.truth_bam, self.bam, self.region, self.rle_ref,
                    self.read_fraction, truth_haplotag=self.truth_haplotag)
            else:
                self._source = self.fencoder.bam_to_sample(
                    self.bam, self.region, self.rle_ref, self.read_fraction)
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
        """Iterator over chunked samples."""
        self._fill_features()
        self._quarantined = list()
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
            yield from alphabet_filter(chunks)


def create_samples(args):
    raise NotImplementedError('Creation of unlabelled samples is currently disabled')


def _labelled_samples_worker(args, region):
    logger = medaka.common.get_named_logger('PrepWork')
    logger.info("Processing region {}.".format(region))
    data_gen = SampleGenerator(
        args.bam, region, args.model, args.rle_ref, truth_bam=args.truth, truth_haplotag=args.truth_haplotag,
        read_fraction=args.read_fraction, chunk_len=args.chunk_len, chunk_overlap=args.chunk_ovlp)
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

