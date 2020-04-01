"""Creation of variant call files from network outputs."""
import collections
import itertools

import intervaltree
import numpy as np
import pysam

import medaka.common
import medaka.datastore
import medaka.vcf


def trim_samples(sample_gen):
    """Generate trimmed samples.

    Samples are trimmed to remove overlap between adjacent samples.

    :param sample_gen: generator yielding `medaka.common.Sample` s

    :yields: (`medaka.common.Sample` view, bool is_last_in_contig)

    """
    logger = medaka.common.get_named_logger('TrimOverlap')
    s1 = next(sample_gen)
    # do not trim beginning of s1
    start_1 = None
    # initialise in case we have one sample
    start_2 = None

    for s2 in itertools.chain(sample_gen, (None,)):
        s1_name = 'Unknown' if s1 is None else s1.name
        s2_name = 'Unknown' if s2 is None else s2.name

        is_last_in_contig = False
        # s1 is last chunk
        if s2 is None:
            # go to end of s1
            end_1 = None
            is_last_in_contig = True
        else:
            rel = medaka.common.Sample.relative_position(s1, s2)
            # skip s2 if it is contained within s1
            if rel is medaka.common.Relationship.s2_within_s1:
                logger.info('{} is contained within {}, skipping.'.format(
                    s2_name, s1_name))
                continue
            elif rel is medaka.common.Relationship.forward_overlap:
                end_1, start_2, _ = medaka.common.Sample.overlap_indices(
                    s1, s2)
            elif rel is medaka.common.Relationship.forward_gapped:
                is_last_in_contig = True
                end_1, start_2 = (None, None)
                msg = '{} and {} cannot be concatenated as there is ' + \
                      'no overlap and they do not abut.'
                logger.info(msg.format(s1_name, s2_name))
            else:
                raise RuntimeError(
                    'Unexpected sample relationship {} '
                    'between {} and {}'.format(repr(rel), s1.name, s2.name))

        yield s1.slice(slice(start_1, end_1)), is_last_in_contig

        s1 = s2
        start_1 = start_2


def join_samples(sample_gen, ref_seq, label_scheme):
    """Process stream of trimmed Samples.

    Care is taken to ensure a variant is not split across
    multiple Samples.

    :param sample_gen: stream of (`medaka.common.Sample`,
        bool is_last_in_contig)
    :param ref_seq: str, reference sequence

    :yields: `medaka.common.Sample` s

    """
    queue = []
    for s, is_last_in_contig in sample_gen:
        if is_last_in_contig:
            # there are no further samples in this contig
            # all remaining variants must be in this contig
            # or in the queue, so process and empty queue
            queue.append(s)
            yield medaka.common.Sample.from_samples(queue)
            queue = []
            continue

        # find last non-variant call, i.e. a major position which
        # has no minor positions and is the same as the ref

        # call retaining gaps cast to array
        call_with_gaps = np.array(list(
            label_scheme.decode_consensus(s, with_gaps=True)), dtype='|U1')

        def get_symbol(p):
            return ref_seq[p['major']] if p['minor'] == 0 else '*'

        ref_seq_with_gaps = np.fromiter(
            (get_symbol(p) for p in s.positions),
            dtype='|U1',
            count=len(s.positions))

        assert len(call_with_gaps) == len(ref_seq_with_gaps)

        is_diff = call_with_gaps != ref_seq_with_gaps
        # don't count a minor position called as gap as a match
        both_gap = np.logical_and(call_with_gaps == '*',
                                  ref_seq_with_gaps == '*')
        is_var = np.logical_or(is_diff, both_gap)

        # if the entire slice of Sample is a variant, add it to queue
        if np.all(is_var):
            queue.append(s)
            continue

        major_inds = np.where(s.positions['minor'] == 0)
        major_pos = s.positions['major'][major_inds]
        major_call = call_with_gaps[major_inds]
        ref_seq_as_array = ref_seq_with_gaps[major_inds]

        assert len(major_call) == len(ref_seq_as_array)

        is_diff = major_call != ref_seq_as_array

        for offset, is_diff_i in enumerate(is_diff[::-1]):
            # Even if the last position does not look like a variant,
            # we can't know if the first pos in next chunk is not an
            # insertion => we need last non-variant position before
            # last position.
            if not is_diff_i:
                break
        last_non_var_pos = major_pos[-1 - offset]
        last_non_var_start = np.searchsorted(
            s.positions['major'],
            last_non_var_pos, side='left')
        to_yield = s.slice(slice(None, last_non_var_start))
        to_queue = s.slice(slice(last_non_var_start, None))

        yield medaka.common.Sample.from_samples(queue + [to_yield])
        queue = [to_queue]
    if len(queue) > 0:
        raise ValueError(
            'Reached end of generator at {} without '
            'is_last_in_contig being True'.format(s.name))


def snps_from_hdf(args):
    """Entry point for SNP calling from HDF5 files.

    A `LabelScheme` read from file is used to decode SNPs. All `LabelScheme` s
    define a `decode_snps` public method. We do not need to use `join_samples`
    to look for variants overlapping sample slice boundaries because we only
    analyse a single locus at a time. This means that `LabelScheme` s that do
    not define the `decode_consensus` method (called within `join_samples`) can
    be used.

    """
    logger = medaka.common.get_named_logger('SNPs')

    index = medaka.datastore.DataIndex(args.inputs)

    if args.regions is None:
        args.regions = sorted(index.index)
    regions = [
        medaka.common.Region.from_string(r)
        for r in args.regions]

    # lookup LabelScheme stored in HDF5
    try:
        label_scheme = index.metadata['label_scheme']
    except KeyError:
        logger.debug(
            "Could not find `label_scheme` metadata in input file, "
            "assuming HaploidLabelScheme.")
        label_scheme = medaka.labels.HaploidLabelScheme()

    # tell label_scheme whether we want verbose info fields
    label_scheme.verbose = args.verbose

    logger.debug("Label decoding is:\n{}".format(
        '\n'.join('{}: {}'.format(k, v)
                  for k, v in label_scheme._decoding.items())))

    meta_info = label_scheme.snp_metainfo

    ref_names = [r.ref_name for r in regions]
    with medaka.vcf.VCFWriter(
            args.output, 'w', version='4.1',
            contigs=ref_names, meta_info=meta_info) as vcf_writer:
        for reg in regions:
            logger.info("Processing {}.".format(reg))
            ref_seq = pysam.FastaFile(args.ref_fasta).fetch(
                reference=reg.ref_name).upper()

            samples = index.yield_from_feature_files(regions=[reg])
            trimmed_samples = trim_samples(samples)

            for sample, is_last in trimmed_samples:
                snps = label_scheme.decode_snps(
                    sample, ref_seq, threshold=args.threshold)
                vcf_writer.write_variants(snps, sort=True)


def variants_from_hdf(args):
    """Entry point for variant calling from HDF5 files.

    A `LabelScheme` read from HDF must define both a `decode_variants`
    and `decode_consnesus` method. The latter is used with `join_samples`
    to detect multi-locus variants spanning `Sample` slice boundaries.

    """
    logger = medaka.common.get_named_logger('Variants')

    index = medaka.datastore.DataIndex(args.inputs)

    if args.regions is None:
        args.regions = sorted(index.index)
    regions = [
        medaka.common.Region.from_string(r)
        for r in args.regions]

    # lookup LabelScheme stored in HDF5
    try:
        label_scheme = index.metadata['label_scheme']
    except KeyError:
        logger.debug(
            "Could not find `label_scheme` metadata in input file, "
            "assuming HaploidLabelScheme.")
        label_scheme = medaka.labels.HaploidLabelScheme()

    logger.debug("Label decoding is:\n{}".format(
        '\n'.join('{}: {}'.format(k, v)
                  for k, v in label_scheme._decoding.items())))

    if not hasattr(label_scheme, 'decode_variants'):
        raise AttributeError(
            '{} does not support decoding of variants'.format(label_scheme))

    if not hasattr(label_scheme, 'decode_consensus'):
        raise AttributeError(
            '{} does not support consensus decoding required '
            'for variant calling.'.format(label_scheme))

    # tell label_scheme whether we want verbose info fields
    label_scheme.verbose = args.verbose

    meta_info = label_scheme.variant_metainfo

    ref_names = [r.ref_name for r in regions]
    with medaka.vcf.VCFWriter(
            args.output, 'w', version='4.1',
            contigs=ref_names, meta_info=meta_info) as vcf_writer:
        for reg in regions:
            logger.info("Processing {}.".format(reg))
            ref_seq = pysam.FastaFile(args.ref_fasta).fetch(
                reference=reg.ref_name).upper()

            samples = index.yield_from_feature_files([reg])
            trimmed_samples = trim_samples(samples)
            joined_samples = join_samples(
                trimmed_samples, ref_seq, label_scheme)

            for sample in joined_samples:
                variants = label_scheme.decode_variants(sample, ref_seq)
                vcf_writer.write_variants(variants, sort=True)


def samples_to_bed(args):
    """Write a bed file from samples in a datastore file."""
    logger = medaka.common.get_named_logger('Variants')

    index = medaka.datastore.DataIndex(args.inputs)

    trees = collections.defaultdict(intervaltree.IntervalTree)
    logger.info("Building interval tree")
    for s, f in index.samples:
        d = medaka.common.Sample.decode_sample_name(s)
        # start and end are string repr of floats (major.minor coordinates)
        start, end = int(float(d['start'])), int(float(d['end']))
        # add one to end of interval, as intervaltree intervals and bed file
        # intervals are end-exclusive (i.e. they don't contain the last
        # coordinate), whilst the last position in a sample is included in that
        # sample.
        trees[d['ref_name']].add(intervaltree.Interval(start, end + 1))

    with open(args.output, 'w') as fh:
        for contig, tree in trees.items():
            # strict=False as consecutive samples can start and end on the same
            # major (overlap is in minor) hence if samples are abutting but not
            # overlapping in major coords, merge them
            tree.merge_overlaps(strict=False)
            logger.info("Writing intervals for {}".format(contig))
            for i in sorted(tree.all_intervals):
                fh.write("{}\t{}\t{}\n".format(contig, i.begin, i.end))

    logger.info("All done, bed file written to {}".format(args.output))
