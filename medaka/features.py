"""Creation of neural network input features."""
import abc
from collections import defaultdict
import concurrent.futures
import functools
import inspect
import itertools
import os
import time
from timeit import default_timer as now

import numpy as np

import libmedaka
import medaka.common
import medaka.datastore
import medaka.labels


def _plp_data_to_numpy(plp_data, n_rows):
    """Create numpy representation of feature data.

    Copy the feature matrix and alignment column names from a
    `plp_data` structure returned from C library function calls.

    :param plp_data: a cffi proxy to a `plp_data*` pointer
    :param nrows: the number of rows in the plp_data.matrix (the number
        of elements in the feature per pileup column).

    :returns: pileup counts numpy array, reference positions

    """
    ffi = libmedaka.ffi
    size_sizet = np.dtype(np.uintp).itemsize
    np_counts = np.frombuffer(ffi.buffer(
        plp_data.matrix, size_sizet * plp_data.n_cols * n_rows),
        dtype=np.uintp
    ).reshape(plp_data.n_cols, n_rows).copy()

    positions = np.empty(plp_data.n_cols, dtype=[
        ('major', int), ('minor', int)])
    np.copyto(
        positions['major'], np.frombuffer(
            ffi.buffer(plp_data.major, size_sizet * plp_data.n_cols),
            dtype=np.uintp))
    np.copyto(
        positions['minor'],
        np.frombuffer(ffi.buffer(
            plp_data.minor, size_sizet * plp_data.n_cols), dtype=np.uintp))
    return np_counts, positions


def __enforce_pileup_chunk_contiguity(pileups):
    """Split and join ordered pileup chunks to ensure contiguity.

    :param pileups: iterable of (counts, pileups) as constructed by
        `_plp_data_to_numpy`.

    :returns: a list of reconstituted (counts, pileups) where discontinuities
        in the inputs cause breaks and abutting inputs are joined.

    """
    split_results = list()
    # First pass: need to check for discontinuities within chunks,
    # these show up as >1 changes in the major coordinate
    for counts, positions in pileups:
        move = np.ediff1d(positions['major'])
        gaps = np.where(move > 1)[0] + 1
        if len(gaps) == 0:
            split_results.append((counts, positions))
        else:
            start = 0
            for i in gaps:
                split_results.append((counts[start:i], positions[start:i]))
                start = i
            split_results.append((counts[start:], positions[start:]))

    # Second pass: stitch abutting chunks together, anything not neighbouring
    # is kept separate whether it came from the same chunk originally or not
    def _finalize_chunk(c_buf, p_buf):
        chunk_counts = np.concatenate(c_buf)
        chunk_positions = np.concatenate(p_buf)
        return chunk_counts, chunk_positions

    counts_buffer, positions_buffer = list(), list()
    chunk_results = list()
    last = None
    for counts, positions in split_results:
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
            chunk_results.append(_finalize_chunk(
                counts_buffer, positions_buffer))
            counts_buffer = [counts]
            positions_buffer = [positions]
            last = positions['major'][-1]
    if len(counts_buffer) != 0:
        chunk_results.append(_finalize_chunk(counts_buffer, positions_buffer))
    return chunk_results


def _tidy_libfunc_args(
        dtype_prefixes, tag_name, tag_value, keep_missing, read_group):
    ffi = libmedaka.ffi

    if dtype_prefixes is None or isinstance(dtype_prefixes, str) \
            or len(dtype_prefixes) == 1:
        num_dtypes, dtypes, _dtypes = 1, ffi.NULL, [ffi.NULL]
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
    if read_group is None:
        read_group = ffi.NULL
    else:
        read_group = ffi.new("char[]", read_group.encode())

    return (num_dtypes, dtypes, _dtypes, tag_name, tag_value,
            keep_missing, read_group)


def pileup_counts(
        region, bam, dtype_prefixes=None, region_split=100000, workers=8,
        tag_name=None, tag_value=None, keep_missing=False, num_qstrat=1,
        weibull_summation=False, read_group=None):
    """Create pileup counts feature array for region.

    :param region: `medaka.common.Region` object
    :param bam: .bam file with alignments.
    :param dtype_prefixes: prefixes for query names which to separate counts.
        If `None` (or of length 1), counts are not split.
    :param region_split: largest region to process in single thread.
    :param workers: worker threads for calculating pileup.
    :param tag_name: two letter tag name by which to filter reads.
    :param tag_value: integer value of tag for reads to keep.
    :param keep_missing: whether to keep reads when tag is missing.
    :param num_qstrat: number of layers for qscore stratification.
    :param weibull_summation: use a Weibull partial-counts approach,
        requires 'WL' and 'WK' float-array tags.

    :returns: pileup counts array, reference positions, insertion positions
    """
    lib = libmedaka.lib
    featlen = lib.featlen
    (num_dtypes, dtypes, _dtypes, tag_name, tag_value,
        keep_missing, read_group) = _tidy_libfunc_args(
            dtype_prefixes, tag_name, tag_value, keep_missing, read_group)

    def _process_region(reg):
        # htslib start is 1-based, medaka.common.Region object is 0-based
        region_str = '{}:{}-{}'.format(reg.ref_name, reg.start + 1, reg.end)

        counts = lib.calculate_pileup(
            region_str.encode(), bam.encode(), num_dtypes, dtypes, num_qstrat,
            tag_name, tag_value, keep_missing, weibull_summation, read_group
        )
        np_counts, positions = _plp_data_to_numpy(
            counts, featlen * num_dtypes * num_qstrat)
        lib.destroy_plp_data(counts)
        return np_counts, positions

    # split large regions for performance
    regions = region.split(region_split, fixed_size=False)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) \
            as executor:
        results = executor.map(_process_region, regions)
        chunk_results = __enforce_pileup_chunk_contiguity(results)

    return chunk_results


def get_trimmed_reads(
        region, bam, dtype_prefixes=None, region_split=750, chunk_overlap=150,
        workers=8, tag_name=None, tag_value=None, keep_missing=False,
        partial=True, num_qstrat=1, read_group=None):
    """Fetch reads trimmed to a region.

    Overlapping chunks of the imput region will be produced, with each chunk
    having its reads trimmed to the reference sequence coordinates of th
    chunk.

    :param region: `medaka.common.Region` object
    :param bam: .bam file with alignments.
    :param dtype_prefixes: prefixes for query names which to separate counts.
        If `None` (or of length 1), counts are not split.
    :param region_split: largest region to process in single thread.
    :param chunk_overlap: overlap between chunks.
    :param workers: worker threads for calculating pileup.
    :param tag_name: two letter tag name by which to filter reads.
    :param tag_value: integer value of tag for reads to keep.
    :param keep_missing: whether to keep reads when tag is missing.
    :param partial: whether to keep reads which don't fully span the region.
    :param num_qstrat: number of layers for qscore stratification.

    :returns: lists of trimmed reads.
    """
    ffi, lib = libmedaka.ffi, libmedaka.lib
    (num_dtypes, dtypes, _dtypes, tag_name, tag_value,
        keep_missing, read_group) = _tidy_libfunc_args(
            dtype_prefixes, tag_name, tag_value, keep_missing, read_group)

    def _process_region(reg):
        # htslib start is 1-based, medaka.common.Region object is 0-based
        region_str = '{}:{}-{}'.format(reg.ref_name, reg.start + 1, reg.end)
        stuff = lib.PY_retrieve_trimmed_reads(
            region_str.encode(), bam.encode(), num_dtypes, dtypes,
            tag_name, tag_value, keep_missing, partial, read_group,
        )
        # last string is reference
        seqs = [(False, ffi.string(stuff.seqs[stuff.n_seqs - 1]).decode())]
        for i in range(stuff.n_seqs - 1):
            seqs.append((stuff.is_rev[i], ffi.string(stuff.seqs[i]).decode()))
        lib.PY_destroy_reads(stuff)
        return reg, seqs

    # split large regions for performance
    regions = region.split(region_split, chunk_overlap)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) \
            as executor:
        results = executor.map(_process_region, regions)

    return results


def pileup_counts_norm_indices(dtypes, num_qstrat=1):
    """Calculate feature vector normalization groups.

    Calculates (per-datatype, per-read-orientation) locations of bases in
    `pileup_counts` output for the purpose of various normalization
    strategies.

    e.g. For two datatype the feature vector is counts of bases:
    (acgtACGTdDacgtACGTdD), i.e. the first datatype followed by
    the second, and where lowercase denotes reverse reads and the
    base ordering is defined in the library consts `plp_bases`. The
    resultant dictionary contains entries such as:
    (datatype1, False):[4,5,6,7,9], for the datatype1 forward bases.

    :param dtypes: list of datatype names.
    :param num_qstrat: number of layers for qscore stratification.

    :returns: a dictionary of the form
        `{(datatype, is_rev): [base1_index, base2_index, ...]}`.
    """
    ffi, lib = libmedaka.ffi, libmedaka.lib
    plp_bases = lib.plp_bases
    featlen = lib.featlen

    indices = defaultdict(list)
    codes = ffi.string(plp_bases).decode()
    assert len(codes) == featlen
    # from the C code:

    #     pileup->matrix[major_col + featlen * dtype * num_qstrat + \
    #     featlen * qstrat + base_i] += 1;
    # where j is the insert index, so feature is ordered: datatype > base
    for dti, dt in enumerate(dtypes):
        for qindex in range(num_qstrat):
            for base_i, code in enumerate(codes):
                is_rev = code.islower()
                indices[dt, is_rev].append(
                    base_i + dti * num_qstrat * len(codes) +
                    qindex * len(codes))

    return dict(indices)


feature_encoders = dict()


class FeatureEncoderRegistrar(type):
    """Class for registering feature encoders."""

    def __new__(cls, clsname, bases, attrs):
        """Register class to `feature_encoders` dict upon instantiation."""
        newclass = super(FeatureEncoderRegistrar, cls).__new__(
            cls, clsname, bases, attrs)
        cls.register_feature_encoder(clsname, newclass)
        return newclass

    def register_feature_encoder(clsname, cls):
        """Add `FeatureEncoder` to `feature_encoders` dict."""
        # do not display base class as command line option
        if clsname != 'BaseFeatureEncoder':
            feature_encoders[clsname] = cls


class FeatureEncoderMeta(abc.ABC, FeatureEncoderRegistrar):
    """Metaclass facilitating registration of `FeatureEncoder` s."""

    pass


class BaseFeatureEncoder(metaclass=FeatureEncoderMeta):
    """Base class for creation of feature arrays from a `.bam` file."""

    @abc.abstractmethod
    def __init__(*args, **kwargs):
        """Initialize feature encoder."""
        raise NotImplementedError

    @abc.abstractmethod
    def _pileup_function(self, region, reads_bam):
        # Called to create pileup matrix and position arrays
        raise NotImplementedError

    @abc.abstractmethod
    def _post_process_pileup(self, matrix, positions, region):
        # A chance to mutate the output of pileup_function
        raise NotImplementedError

    @abc.abstractmethod
    def bams_to_training_samples(
            self, truth_bam, bam, region, label_scheme, truth_haplotag=None,
            min_length=1000):
        """Create labelled samples, should internally call bam_to_sample."""
        # TODO: should be moved outside this class?
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def feature_vector_length(self):
        """Return size of a single feature vector.

        The length of a neural network input at a single time point.
        """
        raise NotImplementedError

    # this shouldn't be overridden
    def bam_to_sample(self, reads_bam, region):
        """Convert a section of an alignment pileup to a sample.

        :param reads_bam: (sorted indexed) bam with read alignment to reference
        :param region: `medaka.common.Region` object with ref_name, start and
            end attributes.

        :returns: `medaka.common.Sample` object

        """
        pileups = self._pileup_function(region, reads_bam)
        samples = list()
        for counts, positions in pileups:
            if len(counts) == 0:
                msg = (
                    'Pileup-feature is zero-length for {} indicating no '
                    'reads in this region.').format(region)
                self.logger.warning(msg)
                samples.append(
                    medaka.common.Sample(
                        ref_name=region.ref_name, features=None,
                        labels=None, ref_seq=None,
                        positions=positions, label_probs=None
                    )
                )
                continue
            samples.append(self._post_process_pileup(
                counts, positions, region))
        return samples

    def __getstate__(self):
        """Modify object so it is pickleable."""
        # Logs are not picklable in python < 3.7
        state = self.__dict__.copy()
        del state['logger']
        return state

    def __setstate__(self, state):
        """Modify object after unpickling."""
        self.__dict__.update(state)
        self.logger = medaka.common.get_named_logger('Feature')


class CountsFeatureEncoder(BaseFeatureEncoder):
    """Create a pileup array of counts of observed bases."""

    _norm_modes_ = ['total', 'fwd_rev', None]
    feature_dtype = np.float32

    def __init__(
            self, normalise='total', dtypes=('',),
            tag_name=None, tag_value=None, tag_keep_missing=False,
            read_group=None):
        """Initialize creation of neural network input features.

        :param normalise: str, how to normalise the data.
        :param dtypes: iterable of str, read id prefixes of distinct data types
            that should be counted separately.
        :param tag_name: two letter tag name by which to filter reads.
        :param tag_value: integer value of tag for reads to keep.
        :param tag_keep_missing: whether to keep reads when tag is missing.
        :param read_group: value of RG tag to which to filter reads.

        """
        self.logger = medaka.common.get_named_logger('Feature')
        self.normalise = normalise
        self.dtypes = dtypes
        self.feature_indices = pileup_counts_norm_indices(self.dtypes)
        self.tag_name = tag_name
        self.tag_value = tag_value
        self.tag_keep_missing = tag_keep_missing
        self.read_group = read_group

        if self.normalise not in self._norm_modes_:
            raise ValueError('normalise={} is not one of {}'.format(
                self.normalise, self._norm_modes_))

        # TODO: refactor/remove/move to base class.
        opts = inspect.signature(
            CountsFeatureEncoder.__init__).parameters.keys()
        opts = {k: getattr(self, k) for k in opts if k != 'self'}
        self.logger.debug("Creating features with: {}".format(opts))

    @property
    def feature_vector_length(self):
        """Return size of a single feature vector.

        The length of a neural network input at a single time point.
        """
        featlen = libmedaka.lib.featlen
        return len(self.dtypes) * featlen

    def _pileup_function(self, region, bam):
        return pileup_counts(
            region, bam,
            dtype_prefixes=self.dtypes,
            tag_name=self.tag_name, tag_value=self.tag_value,
            keep_missing=self.tag_keep_missing, read_group=self.read_group)

    def _post_process_pileup(self, counts, positions, region):
        # Normalise produced counts using chosen method.

        start, end = positions['major'][0], positions['major'][-1]
        # TODO investigate off-by-one
        if start != region.start or end + 1 != region.end:
            self.logger.warning(
                'Pileup counts do not span requested region, requested {}, '
                'received {}-{}.'.format(region, start, end)
            )

        # find the position index for parent major position of all minor
        # positions
        minor_inds = np.where(positions['minor'] > 0)
        major_pos_at_minor_inds = positions['major'][minor_inds]
        major_ind_at_minor_inds = np.searchsorted(
            positions['major'], major_pos_at_minor_inds, side='left')

        depth = np.sum(counts, axis=1)
        depth[minor_inds] = depth[major_ind_at_minor_inds]

        if self.normalise == 'total':
            # normalize counts by total depth at major position, since the
            # counts include deletions this is a count of spanning reads
            # max just to avoid div error
            feature_array = counts / np.maximum(1, depth).reshape((-1, 1))
        elif self.normalise == 'fwd_rev':
            # normalize forward and reverse and by dtype
            feature_array = np.empty_like(counts, dtype=self.feature_dtype)
            for (dt, is_rev), inds in self.feature_indices.items():
                dt_depth = np.sum(counts[:, inds], axis=1)
                dt_depth[minor_inds] = dt_depth[major_ind_at_minor_inds]
                # max just to avoid div err
                feature_array[:, inds] = \
                    counts[:, inds] / np.maximum(1, dt_depth).reshape((-1, 1))
        else:
            feature_array = counts
        feature_array = feature_array.astype(self.feature_dtype)

        sample = medaka.common.Sample(
            ref_name=region.ref_name, features=feature_array,
            labels=None, ref_seq=None,
            positions=positions, label_probs=None
        )
        self.logger.info('Processed {} (median depth {})'.format(
            sample.name, np.median(depth)))
        return sample

    def bams_to_training_samples(
            self, truth_bam, bam, region, label_scheme, truth_haplotag=None,
            min_length=1000):
        """Prepare training data chunks.

        :param truth_bam: .bam file of truth aligned to ref to generate labels.
        :param bam: input .bam file.
        :param region: `medaka.common.Region` instance for region to process.
        :param label_scheme: a `LabelScheme` describing network outputs.
        :param truth_haplotag: two letter tag name used for grouping truth
            labels by haplotype.
        :param min_length: minimum length for valid alignments.

        :returns: tuple of `medaka.common.Sample` objects.

        .. note:: Chunks might be missing if `truth_bam` is provided and
            regions with multiple mappings were encountered.

        """
        # Find truth alignments (with some filtering).
        alns = medaka.labels.TruthAlignment.bam_to_alignments(
            truth_bam, region, haplotag=truth_haplotag,
            min_length=min_length)
        if len(alns) == 0:
            self.logger.info(
                "Filtering and grouping removed all alignments "
                "of truth to ref from {}.".format(region))

        samples = []
        for aln in alns:
            # get labels from truth alignments.
            truth_pos, truth_labels = label_scheme.encode(aln)

            # get features from read alignment data
            aln_samples = self.bam_to_sample(bam, medaka.common.Region(
                region.ref_name, aln[0].start, aln[0].end))

            for sample in aln_samples:
                # create a label array that respects positions present
                # in the feature's position array (where single reads
                # may have inserted bases, creating minor postions absent
                # from the labels position array)
                shape = list(truth_labels.shape)
                shape[0] = len(sample.positions)
                padded_labels = np.full(shape, label_scheme.padding_vector,
                                        dtype=truth_labels.dtype)
                truth_inds = np.where(np.in1d(truth_pos, sample.positions))
                sample_inds = np.where(np.in1d(sample.positions, truth_pos))
                assert len(truth_inds[0]) == len(sample_inds[0])
                assert np.alltrue(
                    truth_pos[truth_inds] == sample.positions[sample_inds])

                padded_labels[sample_inds] = truth_labels[truth_inds]

                sample = sample.amend(labels=padded_labels)
                samples.append(sample)
        return tuple(samples)


class HardRLEFeatureEncoder(CountsFeatureEncoder):
    """Create a pileup array of counts of observed bases.

    Counts are segregated according to run lengths stored in the quality
    information of reads.
    """

    def __init__(
            self, normalise='total', dtypes=('', ), tag_name=None,
            tag_value=None, tag_keep_missing=False, num_qstrat=15,
            read_group=None):
        """Class to generate neural network input features.

        :param normalise: str, how to normalise the data.
        :param dtypes: iterable of str, read id prefixes of distinct data
            types that should be counted separately.
        :param tag_name: two letter tag name by which to filter reads.
        :param tag_value: integer value of tag for reads to keep.
        :param tag_keep_missing: whether to keep reads when tag is missing.
        :param num_qstrat: number of layers for qscore stratification.
        :param read_group: value of RG tag to which to filter reads.

        """
        self.num_qstrat = num_qstrat
        super().__init__(
            normalise, dtypes=dtypes, tag_name=tag_name, tag_value=tag_value,
            tag_keep_missing=tag_keep_missing, read_group=read_group)
        self.feature_indices = pileup_counts_norm_indices(
            self.dtypes, num_qstrat=self.num_qstrat)

    def _pileup_function(self, region, bam):
        return pileup_counts(
            region, bam,
            dtype_prefixes=self.dtypes,
            tag_name=self.tag_name, tag_value=self.tag_value,
            keep_missing=self.tag_keep_missing, num_qstrat=self.num_qstrat,
            read_group=self.read_group)

    @property
    def feature_vector_length(self):
        """Return size of a single feature vector.

        The length of a neural network input at a single time point.
        """
        featlen = libmedaka.lib.featlen
        return len(self.dtypes) * featlen * self.num_qstrat


class SymHardRLEFeatureEncoder(HardRLEFeatureEncoder):
    """HRLE encoder where lack of insertion == deletion.

    In minor positions, when a read spans an insertion but
    does not contain the inserted base, this will be counted
    as a deletion.
    """

    def _pileup_function(self, region, bam):
        [(counts, positions)] = super()._pileup_function(region, bam)

        minor_inds = np.where(positions['minor'] > 0)
        major_pos_at_minor_inds = positions['major'][minor_inds]
        major_ind_at_minor_inds = np.searchsorted(
            positions['major'], major_pos_at_minor_inds, side='left')

        # Correct count of indels in minor positions, where reads that
        # do not contain the insertion are not counted as deletion in
        # minor positions, unlike in major positions
        for (dt, is_rev), inds in self.feature_indices.items():
            dt_depth = np.sum(counts[:, inds], axis=1)

            # Every dtype needs a vector with size featlen x num_qstrat,
            # the elements divided in two (split fwd/reverse). Also, by design,
            # indels in RLE pileups are accumulated in  the first layer of
            # stratification available to that dtype. E.g., for num_qstrat=2
            # and 2 dtypes, 'R1' and 'R2':
            # {('R1', False): [4, 5, 6, 7, 9, 14, 15, 16, 17, 19],
            #  ('R1', True): [0, 1, 2, 3, 8, 10, 11, 12, 13, 18],
            #  ('R2', False): [24, 25, 26, 27, 29, 34, 35, 36, 37, 39],
            #  ('R2', True): [20, 21, 22, 23, 28, 30, 31, 32, 33, 38]}
            featlen_index = libmedaka.lib.rev_del if is_rev else \
                libmedaka.lib.fwd_del
            dtype_size = libmedaka.lib.featlen * self.num_qstrat
            del_ind = [x for x in inds if x % dtype_size == featlen_index][0]
            counts[minor_inds, del_ind] = dt_depth[major_ind_at_minor_inds] -\
                dt_depth[minor_inds]

        return [(counts, positions)]


class SoftRLEFeatureEncoder(HardRLEFeatureEncoder):
    """Create pileups using soft RLE calls."""

    def _pileup_function(self, region, bam):
        return pileup_counts(
            region, bam, dtype_prefixes=self.dtypes, tag_name=self.tag_name,
            tag_value=self.tag_value, keep_missing=self.tag_keep_missing,
            num_qstrat=self.num_qstrat, weibull_summation=True,
            read_group=self.read_group)


class SampleGenerator(object):
    """Orchestration of neural network inputs with chunking."""

    def __init__(self, bam, region, feature_encoder, truth_bam=None,
                 label_scheme=None, truth_haplotag=None, chunk_len=1000,
                 chunk_overlap=200, enable_chunking=True,
                 min_truth_length=1000):
        """Generate chunked inference (or training) samples.

        :param bam: `.bam` containing alignments from which to generate
            samples.
        :param region: a `medaka.common.Region` for which to generate samples.
        :param feature_encoder: a `FeatureEncoder` object used for
            constructing training features.
        :param truth_bam: a `.bam` containing alignment of truth sequence to
            `reference` sequence. Required only for creating training chunks.
        :param label_scheme: a `LabelScheme` used for deriving training truth
            labels.
        :param truth_haplotag: two letter tag name used for grouping truth
            labels by haplotype.
        :param reference: reference `.fasta`, should correspond to `bam`.
        :param enable_chunking: when yielding samples, do so in chunks.
        :param min_truth_length: minimum length of truth alignments to accepts
            for generation of truth labels.

        """
        self.logger = medaka.common.get_named_logger("Sampler")
        self.sample_type = "training" if truth_bam is not None else "consensus"
        self.logger.info("Initializing sampler for {} of region {}.".format(
            self.sample_type, region))
        self.fencoder = feature_encoder
        self.bam = bam
        self.region = region
        self.truth_bam = truth_bam
        self.label_scheme = label_scheme
        self.truth_haplotag = truth_haplotag
        self.chunk_len = chunk_len
        self.chunk_overlap = chunk_overlap
        self.enable_chunking = enable_chunking
        self.min_truth_length = min_truth_length
        self._source = None  # the base data to be chunked
        self._quarantined = list()  # samples which are shorter than chunk size

        if self.truth_bam is not None and self.label_scheme is None:
            raise ValueError(
                "A `LabelScheme` must be given to create training data.")

    def _fill_features(self):
        if self._source is None:
            self._quarantined = None
            t0 = now()
            if self.truth_bam is not None:
                self._source = self.fencoder.bams_to_training_samples(
                    self.truth_bam, self.bam, self.region, self.label_scheme,
                    truth_haplotag=self.truth_haplotag,
                    min_length=self.min_truth_length)
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
        end += 1  # end exclusive
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
                    msg = (
                        "Region {} ({} positions) is smaller than "
                        "inference chunk length {}, quarantining.").format(
                            source.name, source.size, self.chunk_len)
                    self.logger.warning(msg)
                    self._quarantine_sample(source)
                    continue

                self.logger.debug(
                    "Chunking pileup data into {} columns with overlap "
                    "of {}.".format(self.chunk_len, self.chunk_overlap))
                chunks = source.chunks(
                    chunk_len=self.chunk_len, overlap=self.chunk_overlap)
            self.logger.debug("Pileup for {} is of width {}".format(
                source.name, source.size))
            all_data.extend(chunks)

        return all_data


def _samples_worker(args, region, feature_encoder, label_scheme):
    logger = medaka.common.get_named_logger('PrepWork')
    logger.info("Processing region {}.".format(region))
    data_gen = SampleGenerator(
        args.bam, region, feature_encoder, truth_bam=args.truth,
        label_scheme=label_scheme, truth_haplotag=args.truth_haplotag,
        chunk_len=args.chunk_len, chunk_overlap=args.chunk_ovlp)

    return list(data_gen.samples), region


def create_samples(args):
    """Entry point for creation of feature .hdfs (labelled or unlabelled)."""
    logger = medaka.common.get_named_logger('Prepare')
    if args.chunk_ovlp >= args.chunk_len:
        raise ValueError(
            'chunk_ovlp {} is not smaller than chunk_len {}'.format(
                args.chunk_ovlp, args.chunk_len))
    regions = medaka.common.get_regions(args.bam, args.regions)
    reg_str = '\n'.join(['\t\t\t{}'.format(r) for r in regions])
    logger.info('Got regions:\n{}'.format(reg_str))
    if args.truth is None:
        logger.warn(
            'Running medaka features without a truth bam, '
            'unlabelled data will be produced. Is this intended?')
        time.sleep(3)

    no_data = False
    with medaka.datastore.DataStore(args.output, 'w') as ds:
        # write feature options to file
        logger.info("Writing meta data to file.")

        num_qstrat = args.feature_encoder_args.get('num_qstrat')
        max_run = args.label_scheme_args.get('max_run')
        # If one of them is set, set the other to agree.
        # If both are set, force them to agree.
        # If none is set or they are the same continue merrily.
        if max_run is None and num_qstrat is not None:
            args.label_scheme_args['max_run'] = num_qstrat
        elif max_run is not None and num_qstrat is None:
            args.feature_encoder_args['num_qstrat'] = max_run
        elif max_run != num_qstrat:
            raise ValueError(
                 'num_qstrat in feature_encoder_args must agree '
                 'with max_run in feature_encoder_args')

        # Create and serialise to file model ancilliaries
        feature_encoder = feature_encoders[args.feature_encoder](
            **args.feature_encoder_args)
        ds.set_meta(feature_encoder, 'feature_encoder')

        label_scheme = medaka.labels.label_schemes[args.label_scheme](
            **args.label_scheme_args)
        ds.set_meta(label_scheme, 'label_scheme')

        model_function = functools.partial(
            medaka.models.build_model,
            feature_encoder.feature_vector_length,
            len(label_scheme._decoding))
        ds.set_meta(model_function, 'model_function')

        # TODO: this parallelism would be better in
        # `SampleGenerator.bams_to_training_samples` since training
        # alignments are usually chunked.
        ExecutorClass = concurrent.futures.ProcessPoolExecutor
        with ExecutorClass(max_workers=args.threads) as executor:
            # break up overly long chunks
            MAX_SIZE = int(1e6)
            regions = itertools.chain(*(r.split(MAX_SIZE) for r in regions))
            futures = [executor.submit(
                _samples_worker, args, reg,
                feature_encoder, label_scheme) for reg in regions]

            for fut in concurrent.futures.as_completed(futures):
                if fut.exception() is None:
                    samples, region = fut.result()
                    logger.info(
                        "Writing {} samples for region {}".format(
                            len(samples), region))
                    for sample in samples:
                        ds.write_sample(sample)
                else:
                    logger.info(fut.exception())
                    logger.info(fut.result())
                fut._result = None  # python issue 27144
        no_data = ds.n_samples == 0

    if no_data:
        logger.critical(
            "Warning: No training data was written to file, "
            "deleting output.")
        os.remove(args.output)
