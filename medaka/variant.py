"""Creation of vcf files from network outputs or consensus sequences."""
import collections
import itertools
from timeit import default_timer as now

import intervaltree
import numpy as np
import pysam

import medaka.common
import medaka.datastore
import medaka.vcf


def apply_variants(variants, ref_seq):
    """Apply variants (like bcftools consensus does).

    :param variants: iterator of `medaka.vcf.Variant` objs.
    :param ref_seq: str, reference sequence.

    :returns: str, sequence with variants applied.
    """
    query_from_vcf = list(ref_seq)
    for v in variants:
        query_from_vcf[v.pos: v.pos + len(v.ref)] = len(v.ref) * ['']
        query_from_vcf[v.pos] = v.alt[0]
    return ''.join(query_from_vcf)


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
    for s, is_last_in_contig, _ in sample_gen:
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
        both_gap = np.logical_and(
            call_with_gaps == '*', ref_seq_with_gaps == '*')
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

        left_slice = s.slice(slice(None, last_non_var_start))
        right_slice = s.slice(slice(last_non_var_start, None))

        to_yield = queue
        if last_non_var_start > 0:
            # left_slice only appended if it contains at least one position
            to_yield += [left_slice]

        if to_yield:
            # Cowardly refuse to call from_samples on empty list
            yield medaka.common.Sample.from_samples(to_yield)

        queue = [right_slice]

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
        args.regions = index.regions

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

    with pysam.FastaFile(args.ref_fasta) as fa:
        lengths = dict(zip(fa.references, fa.lengths))

    with medaka.vcf.VCFWriter(
            args.output, 'w', version='4.1',
            contigs=['{},length={}'.format(r.ref_name, lengths[r.ref_name])
                     for r in args.regions],
            meta_info=meta_info) as vcf_writer:
        for reg in args.regions:
            logger.info("Processing {}.".format(reg))
            ref_seq = pysam.FastaFile(args.ref_fasta).fetch(
                reference=reg.ref_name).upper()

            samples = index.yield_from_feature_files(regions=[reg])
            trimmed_samples = medaka.common.Sample.trim_samples(samples)

            for sample, is_last, _ in trimmed_samples:
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
        args.regions = index.regions

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

    with pysam.FastaFile(args.ref_fasta) as fa:
        lengths = dict(zip(fa.references, fa.lengths))

    with medaka.vcf.VCFWriter(
            args.output, 'w', version='4.1',
            contigs=['{},length={}'.format(r.ref_name, lengths[r.ref_name])
                     for r in args.regions],
            meta_info=meta_info) as vcf_writer:
        for reg in args.regions:
            logger.info("Processing {}.".format(reg))
            ref_seq = pysam.FastaFile(args.ref_fasta).fetch(
                reference=reg.ref_name).upper()

            samples = index.yield_from_feature_files([reg])
            trimmed_samples = medaka.common.Sample.trim_samples(samples)
            joined_samples = join_samples(
                trimmed_samples, ref_seq, label_scheme)

            for sample in joined_samples:
                variants = label_scheme.decode_variants(
                    sample, ref_seq, ambig_ref=args.ambig_ref,
                    return_all=args.gvcf)
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


AlignPos = collections.namedtuple('AlignPos',  ('rpos', 'rbase', 'qbase'))


def yield_variants_from_aln(aln, ref_seq):
    """Yield variants in an alignment.

    :param aln: `pysam.AlignedSegment` obj.
    :param ref_seq: str, refernece sequence.
    :yields: `medaka.vcf.Variant` objs.
    """
    if aln.is_unmapped or aln.is_secondary or dict(aln.tags)['NM'] == 0:
        return ()
    last_match = None
    seq = aln.query_sequence
    rstart, rend = aln.reference_start, aln.reference_end
    qstart, qend = aln.query_alignment_start, aln.query_alignment_end
    gt = {'GT': '1'}
    chrm = aln.reference_name
    queue = []

    def is_not_end(x):
        return x[1] != rend and x[0] != qend

    def is_not_start(x):
        return x[1] != rstart and x[0] != qstart

    def decode_variant(queue):
        def pos_is_none(x):
            return x.rpos is None
        pos = next(itertools.dropwhile(pos_is_none, queue)).rpos
        ref, alt = (
            ''.join([getattr(p, a) for p in queue]).replace('-', '').upper()
            for a in ('rbase', 'qbase')
        )
        return medaka.vcf.Variant(chrm, pos, ref, alt, genotype_data=gt).trim()

    pairs = itertools.takewhile(
        is_not_end, itertools.dropwhile(is_not_start,
                                        aln.get_aligned_pairs()))
    for qp, rp in pairs:
        qb = seq[qp] if qp is not None else '-'
        rb = ref_seq[rp] if rp is not None else '-'
        p = AlignPos(rp, rb, qb)

        if qb == rb:  # match
            if queue:  # decode queue to variant
                # append match to queue so variants are padded by a match on
                # both side so we can easily handle indels at start/end of ref
                queue.append(p)
                yield decode_variant(queue)
                queue = []
            last_match = p
        else:  # variant
            if not queue and last_match is not None:
                queue.append(last_match)
            queue.append(p)
    if queue:
        # we have a variant at the end of the alignment, there is no next match
        yield decode_variant(queue)


def edlib_chunked_align_fastas(query_fasta, ref_fasta, contigs=None, **kwargs):
    """Perfom chunked alignment of query to reference using edlib.

    :param query_fasta: str, path to query fasta file.
    :param ref_fasta: str, path to reference fasta file.
    :param contigs: iterable of str, if provided, limit alignments to these
        contigs. Contig names should be the same in query and reference fastas.
    :param kwargs: `medaka.align.chunked_edlib_align` keyword arguments.
    :yields: `pysam.AlignedSegment` objs.
    """
    contigs = medaka.common.common_fasta_contigs(
        [ref_fasta, query_fasta], contigs)
    if len(contigs) == 0:
        raise KeyError('Reference and query contig names should match')

    ref_fasta = pysam.FastaFile(ref_fasta)
    query_fasta = pysam.FastaFile(query_fasta)
    for contig in contigs:
        ref = ref_fasta.fetch(contig)
        query = query_fasta.fetch(contig)
        yield from medaka.align.chunked_edlib_align(query, ref, contig,
                                                    **kwargs)


def vcf_from_fasta(args):
    """Entry point for calling variants by consensus sequence alignment."""
    logger = medaka.common.get_named_logger('CONS2VCF')

    with pysam.FastaFile(args.ref_fasta) as fasta:
        ref_seqs = {name: fasta.fetch(name) for name in fasta.references}
        contig_lengths = dict(zip(fasta.references, fasta.lengths))
        total_bp = sum(fasta.lengths)
        ref_contigs = fasta.references
        h = pysam.AlignmentHeader().from_references(
            fasta.references, fasta.lengths)

    if args.bam is not None:
        alns = pysam.AlignmentFile(args.bam)
        out_bam = None
    else:
        out_bam = pysam.AlignmentFile(args.out_prefix + '.bam', 'wb', header=h)
        if args.regions is not None:
            contigs = [r.ref_name for r in args.regions]
        else:
            contigs = None
        alns = edlib_chunked_align_fastas(
            args.consensus, args.ref_fasta, contigs,
            chunk_size=args.chunk_size, pad=args.pad, mode=args.mode,
            header=h)
    vcf_fp = args.out_prefix + '.vcf'
    trees = collections.defaultdict(intervaltree.IntervalTree)
    t_log = now()
    log_interval = 5
    msg = 'Processed {:.2%} of reference.'
    bp_done = collections.Counter()

    header_contigs = [
        '{},length={}'.format(c, contig_lengths[c])
        for c in ref_contigs]
    meta_info = [medaka.vcf.MetaInfo(
        'FORMAT', 'GT', 1, 'String', 'Medaka genotype.')]
    with medaka.vcf.VCFWriter(
            vcf_fp, contigs=header_contigs, meta_info=meta_info) as writer:
        for aln in alns:
            # reference_start is 0 based, reference_end points to one past
            # the last aligned residue, i.e. same as bed file
            ref = aln.reference_name
            rstart, rend = aln.reference_start, aln.reference_end
            if trees[ref].overlaps(rstart, rend) and args.bam is not None:
                # We expect edlib alignments to overlap by 1 match so only
                # apply this check for a user-provided bam.
                logger.warning(
                    ('WARNING: alignment {}:{}-{} overlaps another ' +
                        'alignment, which could cause overlapping variants.' +
                        '\nCheck output bam and vcf for details.').format(
                        ref, rstart, rend))
            trees[ref].add(intervaltree.Interval(rstart, rend))
            for v in yield_variants_from_aln(aln, ref_seqs[ref]):
                if 'N' in v.ref:
                    continue
                writer.write_variant(v)
                if now() - t_log > log_interval:
                    done = bp_done[ref] + v.pos - rstart
                    logger.info(msg.format(done / total_bp))
                    t_log = now()
            bp_done[ref] += rend - rstart
            if out_bam is not None:
                out_bam.write(aln)

    if out_bam is not None:
        out_bam.close()
        pysam.index(out_bam.filename.decode("utf-8"))

    bed_fp = args.out_prefix + '_coverage.bed'
    gap_bed_fp = args.out_prefix + '_coverage_gaps.bed'
    for tree in trees.values():
        # strict=False to merge abutting alignments.
        tree.merge_overlaps(strict=False)
    medaka.common.write_intervaltrees_to_bed(trees, bed_fp)
    gap_trees = medaka.common.complement_intervaltrees(trees, contig_lengths)
    medaka.common.write_intervaltrees_to_bed(gap_trees, gap_bed_fp)
    # loop over contigs for which we have alignments checking for gaps
    for contig in trees:
        if len(gap_trees[contig]):
            logger.info((
                'WARNING: There are alignment gaps for ref contig'
                ' {}, see bed files for details.').format(contig))
    if len(ref_contigs) != len(trees):
        logger.info(
            'WARNING: Some contigs have no alignments, see bed files'
            ' for details.')
    # bp_done calculated above does not take account of overlapping alignments
    # hence recalculate here based on merged alignment intervals.
    aligned_bp = sum((i.length() for tree in trees.values() for i in tree))
    msg = 'Alignments spanned {:%} of the reference.'
    logger.info(msg.format(aligned_bp / total_bp))
    msg = 'Check bed files {} and {} for alignment coverage and gaps.'
    logger.info(msg.format(bed_fp, gap_bed_fp))
    logger.info('All done. VCF written to {}.'.format(vcf_fp))
