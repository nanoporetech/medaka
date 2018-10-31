from collections import defaultdict, Counter, OrderedDict
import inspect
import itertools
from functools import partial
from multiprocessing import Pool
import os
import sys
from timeit import default_timer as now

from Bio import SeqIO
import h5py
import numpy as np
import pysam

from medaka.common import (yield_compressed_pairs, Sample, lengths_to_rle, rle,
                           Region, parse_regions, chunk_samples, write_samples_to_hdf,
                           segment_limits, encode_sample_name, decoding, encoding,
                           _gap_, _alphabet_, _feature_opt_path_, _feature_decoding_path_,
                           get_pairs, get_pairs_with_hp_len, seq_to_hp_lens,
                           write_yaml_data, load_yaml_data, sample_to_x_y, gen_train_batch,
                           _feature_batches_path_, _label_batches_path_, _label_counts_path_,
                           _label_decod_path_,)

from medaka.labels import TruthAlignment
import libmedaka


import logging
logger = logging.getLogger(__name__)


def pileup_counts(region, bam):
    """Create pileup counts feature array for region.

    :param region: `Region` object
    :param bam: .bam file with alignments.

    :returns: pileup counts array, reference positions, insertion postions
    """
    # htslib start is 1-based, Region object is 0-based
    region_str = '{}:{}-{}'.format(region.ref_name, region.start + 1, region.end)

    ffi, lib = libmedaka.ffi, libmedaka.lib
    counts = lib.calculate_pileup(
        region_str.encode(), bam.encode()
    )

    size_sizet = np.dtype(np.uintp).itemsize
    np_counts = np.frombuffer(ffi.buffer(
        counts.counts, size_sizet * counts.n_cols * 11),
        dtype=np.uintp
    ).reshape(counts.n_cols, 11).copy()

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


class FeatureEncoder(object):
    """Class to support multiple feature encodings.

    """
    _ref_modes_ = ['onehot', 'base_length', 'index', None]
    _norm_modes_ = ['total', 'fwd_rev', None]


    def __init__(self, ref_mode:str=None, max_hp_len:int=10,
                 log_min:int=None, normalise:str='total', with_depth:bool=False,
                 consensus_as_ref:bool=False, is_compressed:bool=True):
        """Class to support multiple feature encodings

        :param ref_mode: str, how to represent the reference.
        :param max_hp_len: int, longest homopolymer run which can be represented, longer runs will be truncated.
        :param log_min: int, take log10 of counts/fractions and set zeros to 10**log_min.
        :param normalise: str, how to normalise the data.
        :param with_depth: bool, whether to include a feature describing the total depth.
        :param consensus_as_ref: bool, whether to use a naive max-count consensus instead of the reference.
        :param is_compressed: bool, whether to use HP compression. If false, treat as uncompressed.

        :returns: `Sample` object
        """
        self.ref_mode = ref_mode
        self.consensus_as_ref = consensus_as_ref
        self.max_hp_len = max_hp_len
        self.log_min = log_min
        self.normalise = normalise
        # No point in using default np.float64 when most GPUs want max np.float32
        # TODO: any benefits going to np.float16 or smaller?
        self.feature_dtype = np.float32 if (self.normalise is not None or self.log_min is not None) else np.uint64
        self.with_depth = with_depth
        self.is_compressed = is_compressed

        if self.ref_mode not in self._ref_modes_:
            raise ValueError('ref_mode={} is not one of {}'.format(self.ref_mode, self._ref_modes_))
        if self.normalise not in self._norm_modes_:
            raise ValueError('normalise={} is not one of {}'.format(self.normalise, self._norm_modes_))

        # set up one-hot encoding of read run lengths
        read_decoding = [ k for k in itertools.product((True, False),
                                                            _alphabet_,
                                                            range(1, max_hp_len + 1))]
        # forward and reverse gaps
        read_decoding += [(True, None, 1), (False, None, 1)]

        if self.ref_mode == 'onehot':
            ref_decoding = [ ('ref', b, l) for b, l in itertools.product(alphabet, range(1, max_hp_len + 1))]
            ref_decoding.append(('ref', _gap_, 1))  # gaps
        elif self.ref_mode == 'base_length':
            ref_decoding = ['ref_base', 'ref_length']
            self.ref_base_encoding = { b:i for i, b in enumerate(_alphabet_ + _gap_)}
        elif self.ref_mode == 'index':
            ref_decoding = ['ref_index']
        else:
            ref_decoding = []

        self.decoding = tuple(read_decoding + ref_decoding)
        if self.with_depth:
            self.decoding = self.decoding + ('depth',)
        self.encoding = OrderedDict(((a, i) for i, a in enumerate(self.decoding)))


    def process_ref_seq(self, ref_name, ref_fq):
        """Get encoded reference

        :param ref_name: str, name of contig within `ref_fq`.
        :param ref_fq: path to reference fastq, or `SeqIO.index` obj.
        """
        if self.is_compressed:  # get (compressed) rle_encoded HP len from fq
            if ref_fq is None:
                raise ValueError('If homopolymers have been compressed, ref_rle_fq must be provided')
            return get_runs_from_fastq(ref_fq, ref_name)
        elif self.max_hp_len > 1:  # get (uncompressed) hp_lens from fq
            if ref_fq is None:
                raise ValueError('If max_hp_len > 1, ref_rle_fq must be provided')
            if isinstance(ref_fq, str):
                ref_fq = SeqIO.index(ref_fq, 'fastq')
            read = ref_fq[ref_name]
            return seq_to_hp_lens(read.seq)
        else:
            return None


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
        assert self.normalise == 'total' or self.normalise is None
        assert not self.with_depth
        assert not self.is_compressed

        counts, positions = pileup_counts(region, reads_bam)
        # get rid of first counts row (counts of alternative bases)
        counts = counts[:, 1:]

        depth = np.sum(counts, axis=1)

        # we don't have counts for reads which don't have an insertion
        # so we have to take the depth at major positions
        minor_inds = np.where(positions['minor'] > 0)
        major_pos_at_minor_inds = positions['major'][minor_inds]
        major_ind_at_minor_inds = np.searchsorted(positions['major'], major_pos_at_minor_inds, side='left')
        depth[minor_inds] = depth[major_ind_at_minor_inds]

        if self.normalise == 'total':
            # normalise counts by total depth
            feature_array = counts / depth.reshape((-1, 1))
        else:
            feature_array = counts

        feature_array = feature_array.astype(self.feature_dtype)

        sample = Sample(ref_name=region.ref_name, features=feature_array,
                        labels=None, ref_seq=None,
                        positions=positions, label_probs=None)

        logging.info('Processed {} (median depth {})'.format(encode_sample_name(sample), np.median(depth)))
        return sample


    def bam_to_sample(self, reads_bam, region, ref_rle_fq, read_fraction=None, force_py=False):
        """Converts a section of an alignment pileup (as shown
        by e.g. samtools tview) to a base frequency feature array

        :param reads_bam: (sorted indexed) bam with read alignment to reference
        :param region: `Region` object with ref_name, start and end attributes.
        :param ref_rle_fq: result of `process_ref_seq`
        :param force_py: bool, if True, force use of python code (rather than c library).
        :returns: `Sample` object
        """
        # Try to use fast c function if we can, else fall back on this function
        if not force_py and (ref_rle_fq is None and read_fraction is None):
            try:
                return self.bam_to_sample_c(reads_bam, region)
            except Exception as e:
                logging.debug('Could not process sample with bam_to_sample_c, using python code instead.\n({}).'.format(e))
                pass

        if self.is_compressed:
            aln_to_pairs = partial(yield_compressed_pairs, ref_rle=ref_rle_fq)
        elif self.max_hp_len == 1:
            aln_to_pairs = get_pairs
        else:
            aln_to_pairs = partial(get_pairs_with_hp_len, ref_seq=ref_rle_fq)

        # accumulate data in dicts
        aln_counters = defaultdict(Counter)
        ref_bases = dict()
        with pysam.AlignmentFile(reads_bam, 'rb') as bamfile:
            n_aln = bamfile.count(region.ref_name, region.start, region.end)
            #print('{} reads aligned to ref segment.'.format(n_aln))
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
                logging.debug(msg.format(replace, n_reads, n_reads_to_keep, fraction, region))
                aln_reads = np.random.choice(aln_reads, n_reads_to_keep, replace=replace)

            start = region.start
            end = region.end
            if start is None:
                start = 0
            if end is None:
                end = float('Inf')

            for aln in aln_reads:
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
                    [self.encoding[reverse, pair.qbase, min(pair.qlen, self.max_hp_len)]]) += 1

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
            depth_array = np.empty(shape=(aln_cols), dtype=int)

            # keep track of which features are for fwd/rev reads
            fwd_feat_inds = [self.encoding[k] for k in self.encoding.keys() if not k[0]]
            rev_feat_inds = [self.encoding[k] for k in self.encoding.keys() if k[0]]

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
                    major_depth_fwd = sum((counts[i] for i in fwd_feat_inds))
                    major_depth_rev = sum((counts[i] for i in rev_feat_inds))
                    assert major_depth_fwd + major_depth_rev == major_depth

                if self.normalise is not None:
                    if self.normalise == 'total':
                        feature_array[i, :] /= max(major_depth, 1)
                    elif self.normalise == 'fwd_rev':
                        # normalize fwd and reverse seperately
                        feature_array[i, fwd_feat_inds] /= max(major_depth_fwd, 1)
                        feature_array[i, rev_feat_inds] /= max(major_depth_rev, 1)

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
            logging.info('Processed {} (median depth {}) (in python)'.format(encode_sample_name(sample), np.median(depth_array)))

            return sample


    def bams_to_training_samples(self, truth_bam, bam, region, ref_rle_fq,
                                chunk_len=1000, overlap=0, read_fraction=None):
        """Prepare training data chunks.

        :param truth_bam: .bam file of truth aligned to ref to generate labels.
        :param bam: input .bam file.
        :param region: `Region` obj.
            the reference will be parsed.
        :param ref_rle_fq: path to fastq rle-compressed reference, or result
            of `get_runs_from_fastq`
        :param chunk_len: int, chunk length in reference coordinates.
        :param overlap: int, overlap of adjacent chunks in reference
            coordinates.

        :returns: tuple of `Sample` objects (one item for each input bam) for
            each chunk.

        .. note:: Chunks might be missing if `truth_bam` is provided and
            regions with multiple mappings were encountered.

        """

        # filter truth alignments to restrict ourselves to regions of the ref where the truth
        # in unambiguous
        alignments = TruthAlignment.bam_to_alignments(truth_bam, region.ref_name, start=region.start, end=region.end)
        filtered_alignments = TruthAlignment.filter_alignments(alignments, start=region.start, end=region.end)
        if len(filtered_alignments) == 0:
            logging.info("Filtering removed all alignments of truth to ref from {}.".format(region))

        samples = []
        for aln in filtered_alignments:
            mock_compr = self.max_hp_len > 1 and not self.is_compressed
            truth_pos, truth_labels = aln.get_positions_and_labels(ref_compr_rle=ref_rle_fq, mock_compr=mock_compr,
                                                                   is_compressed=self.is_compressed, rle_dtype=True)
            sample = self.bam_to_sample(bam, Region(region.ref_name, aln.start, aln.end), ref_rle_fq, read_fraction=read_fraction)
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
    logging.debug("Alphabet filter: alphabet: {}".format(alphabet))

    alphabet = set([encoding[c] for c in alphabet])

    def _find_bad_bases(s, field, alphabet):
        seq_rle = getattr(s, field)
        bases = set(np.unique(seq_rle['base']))
        #logging.info("{}: {} bases {} alpha {}".format(encode_sample_name(s), field, bases, alphabet))
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


def min_positions_filter(sample_gen, min_positions):
    """Filter samples with fewer positions than a minimum

    :param sample_gen: generator of `Sample` named tuples.
    :param min_positions: int, samples with fewer positions will be skipped.
    :yields: `Sample` named tuples
    """
    for s in sample_gen:
        if len(s.positions) < min_positions:
            msg = 'Skipping sample {} which has {} columns < min {}.'
            logger.info(msg.format(encode_sample_name(s), len(s.positions), min_positions))
        else:
            yield s


def _get_fe_kwargs(model=None):
    if model is not None:
        fe_kwargs = load_yaml_data(model, _feature_opt_path_)
    else:
        fe_kwargs = { k:v.default for (k,v) in inspect.signature(FeatureEncoder.__init__).parameters.items()
                      if v.default is not inspect.Parameter.empty}
    return fe_kwargs


def get_feature_gen(args):
    mapper = Pool(args.threads).imap if args.threads > 1 else map

    regions = _get_regions(args.bam, args.regions)

    fe_kwargs = _get_fe_kwargs(args.model)
    opt_str = '\n'.join(['{}: {}'.format(k,v) for k, v in fe_kwargs.items()])
    logging.info('FeatureEncoder options: \n{}'.format(opt_str))

    fe = FeatureEncoder(**fe_kwargs)

    truth = hasattr(args, 'truth') and args.truth is not None

    if truth:
        logging.info('Creating training features.')
    else:
        logging.info('Creating consensus features.')

    reg_str = '\n'.join(['\t\t\t{ref_name}:{start}-{end}'.format(**r._asdict()) for r in regions])
    logging.info('Got regions:\n{}'.format(reg_str))
    logging.info('Got read fraction {}'.format(args.read_fraction))

    sample_gens = []
    for region in regions:
        ref_rle = fe.process_ref_seq(region.ref_name, args.ref_fastq)
        if truth:
            f = partial(fe.bams_to_training_samples, args.truth, args.bam, ref_rle_fq=ref_rle, read_fraction=args.read_fraction)
        else:
            f = partial(fe.bam_to_sample, args.bam, ref_rle_fq=ref_rle, read_fraction=args.read_fraction)
        seg_regions = [Region(region.ref_name, start, end)
                        for start, end in segment_limits(region.start, region.end,
                                            segment_len=10*args.chunk_len, overlap_len=5*args.chunk_ovlp)]
        sample_gens.append(mapper(f, seg_regions))
        logging.debug("Created sample gen for region {}".format(region))

    logging.debug("Created list of sample generators")
    if truth:
        # generators yield tuples of samples
        samples = (s for g in sample_gens for t in g for s in t)
    else:
        # generatores yield samples
        samples = (s for g in sample_gens for s in g)

    logging.debug("Created chained sample generators")

    return alphabet_filter(chunk_samples(min_positions_filter(samples, args.chunk_len),
                           chunk_len=args.chunk_len, overlap=args.chunk_ovlp))


def features(args):

    write_samples_to_hdf(args.output, get_feature_gen(args))
    # write feature options to file, so we can later check model and features
    # are compatible.
    fe_kwargs = _get_fe_kwargs(args.model)
    fe = FeatureEncoder(**fe_kwargs)
    to_save = {_feature_opt_path_: fe_kwargs,
               _feature_decoding_path_: fe.decoding}
    for fname in (args.output, args.output + '.yml'):
        write_yaml_data(fname, to_save)


def _process_region(region_getter, filter, chunker, s2xy, region_getter_kwargs):
    region = region_getter(**region_getter_kwargs)
    gen = filter(chunker(region))
    return [s2xy(s) for s in gen]


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
    label_decoding = [ (b, l) for b, l in itertools.product(encoded_bases,
                                                      range(1, max_label_len + 1))]
    label_decoding.append((encoding[_gap_], 1))  # gaps
    label_encoding = {t: i for i, t in enumerate(label_decoding)}
    label_decoding_strs = [l * decoding[b] for (b,l) in label_decoding]
    return label_encoding, label_decoding_strs


def training_batches(args):
    if args.truth is None:
        raise ValueError('If creating training batches, a truth bam must be provided.')
    mapper = Pool(args.threads).imap_unordered if args.threads > 1 else map

    regions = _get_regions(args.bam, args.regions)


    if args.model is not None:
        fe_kwargs = load_yaml_data(args.model, _feature_opt_path_)
    else:
        fe_kwargs = { k:v.default for (k,v) in inspect.signature(FeatureEncoder.__init__).parameters.items()
                      if v.default is not inspect.Parameter.empty}
    opt_str = '\n'.join(['{}: {}'.format(k,v) for k, v in fe_kwargs.items()])
    logging.info('FeatureEncoder options: \n{}'.format(opt_str))

    fe = FeatureEncoder(**fe_kwargs)

    reg_str = '\n'.join(['\t\t\t{ref_name}:{start}-{end}'.format(**r._asdict()) for r in regions])
    logging.info('Got regions:\n{}'.format(reg_str))

    chunked_regions = [Region(region.ref_name, start, end)
                       for region in regions
                       for start, end in segment_limits(region.start, region.end,
                                                        segment_len=10*args.chunk_len,
                                                        overlap_len=5*args.chunk_ovlp)]

    unique_refs = set((region.ref_name for region in regions))
    ref_rles = {r: fe.process_ref_seq(r, args.ref_fastq) for r in unique_refs}

    rg = partial(fe.bams_to_training_samples, args.truth, args.bam, read_fraction=args.read_fraction)
    chunker = partial(chunk_samples, chunk_len=args.chunk_len, overlap=args.chunk_ovlp)

    label_encoding, label_decoding = get_label_encoding(args.max_label_len)
    s2xy = partial(sample_to_x_y, encoding=label_encoding, max_label_len=args.max_label_len)

    f = partial(_process_region, rg, alphabet_filter, chunker, s2xy)

    np.random.shuffle(chunked_regions)  # shuffle in place
    kwargs_gen = ({'region': r, 'ref_rle_fq': ref_rles[r.ref_name]} for r in chunked_regions)
    region_gen = mapper(f, kwargs_gen)
    xy_tups = (xy for xy_chunk in region_gen for xy in xy_chunk)

    labels_counter = Counter()

    with h5py.File(args.output, 'w') as hdf:
        for i, (x_batch, y_batch) in enumerate(gen_train_batch(xy_tups, batch_size=args.batch_size)):
            unique, counts = np.lib.arraysetops.unique(y_batch, return_counts=True)
            labels_counter.update(dict(zip(unique, counts)))
            hdf['{}/{}'.format(_feature_batches_path_, i)] = x_batch
            hdf['{}/{}'.format(_label_batches_path_, i)] = y_batch

    # write feature options to file, so we can later check model and features
    # are compatible and label counts, so we can weight labels by inverse counts.
    to_save = {_feature_opt_path_: fe_kwargs,
               _feature_decoding_path_: fe.decoding,
               _label_counts_path_: labels_counter,
               _label_decod_path_: label_decoding}
    for fname in (args.output, args.output + '.yml'):
        write_yaml_data(fname, to_save)


def _get_regions(bam, region_strs=None):
    """Create Region objects from a bam and region strings.
    :param bam: str, path to bam file.
    :param region_strs: iterable of str in samtools region format e.g. ref:start-end
        or filepath containing a region str per line.
    :returns: list of `Region` objects.
    """
    with pysam.AlignmentFile(bam) as bam_fh:
        ref_lengths = dict(zip(bam_fh.references, bam_fh.lengths))
    if region_strs is not None:
        if os.path.isfile(region_strs[0]):
            with open(region_strs[0]) as fh:
                region_strs = [ l.strip() for l in fh.readlines()]

        regions = []
        for r in parse_regions(region_strs):
            start = r.start if r.start is not None else 0
            end = r.end if r.end is not None else ref_lengths[r.ref_name]
            regions.append(Region(r.ref_name, start, end))
    else:
        regions = [Region(ref_name, 0, end) for ref_name, end in ref_lengths.items()]

    return regions


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
        logging.warning('Some homopolymers in {} are longer than the longest supported length\n'.format(read.name))
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
    logging.info('Compressing {} took {:.3f}s.'.format(args.input, t1 - t0))

    if args.output is not None:
        fh.close()


def choose_feature_func(args):
    if args.batch_size is not None:
        training_batches(args)
    else:
        features(args)
