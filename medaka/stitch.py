"""Creation of contiguous consensus sequences from chunked network outputs."""
import collections
import concurrent.futures
import functools
import operator

import intervaltree
import pysam

import medaka.common
import medaka.datastore
import medaka.labels


def write_fasta(filename, contigs):
    """Write a fasta file from tuples of (name, sequence).

    :param filename: output filename.
    :param contigs: tuples of the form (sequence name, base sequence).

    """
    with open(filename, 'w') as fasta:
        for name, seq in contigs:
            fasta.write('>{}\n{}\n'.format(name, seq))


def stitch_from_probs(h5_fp, regions):
    """Join overlapping label probabilities from HDF5 files.

     Network outputs from multiple samples stored within a file are spliced
     together into a logically contiguous array and decoded to generate
     contiguous sequence(s).

    :param h5_fp: iterable of HDF5 filepaths.
    :param regions: iterable of `medaka.common.Region` objs to process.

    :returns: list of (region string, sequence).
    """
    logger = medaka.common.get_named_logger('Stitch')
    if isinstance(regions, medaka.common.Region):
        regions = [regions]
    logger.info("Stitching regions: {}".format([str(r) for r in regions]))
    index = medaka.datastore.DataIndex(h5_fp)
    label_scheme = index.metadata['label_scheme']
    logger.debug("Label decoding is:\n{}".format(
        '\n'.join('{}: {}'.format(k, v)
                  for k, v in label_scheme._decoding.items())))

    def get_pos(sample, i):
        return '{}.{}'.format(
            sample.positions[i]['major'],
            sample.positions[i]['minor'])

    ref_assemblies = []
    for reg in regions:
        logger.info("Processing {}.".format(reg))
        data_gen = medaka.common.Sample.trim_samples(
            index.yield_from_feature_files(regions=[reg]))
        cur_ref_name = None
        heuristic_use = 0
        seq_parts = list()
        cur_segment = 0
        start = None
        for s, is_last_in_contig, heuristic in data_gen:
            cur_ref_name = s.ref_name if cur_ref_name is None else cur_ref_name
            start = get_pos(s, 0) if start is None else start
            seq_parts.append(label_scheme.decode_consensus(s))
            if is_last_in_contig:
                ref_assemblies.append((
                    '{}_segment{}'.format(s.ref_name, cur_segment),
                    '{}:{}-{}'.format(s.ref_name, start, get_pos(s, -1)),
                    ''.join(seq_parts)))
                seq_parts = list()
                start = None
                cur_segment += 1
            heuristic_use += heuristic
        logger.info("Used heuristic {} times for {}.".format(
            heuristic_use, reg))
    return ref_assemblies


def fill_gaps(contigs, draft):
    """Fill gaps between polished contigs with draft sequence.

    :param contigs: iterable of (name, info, seq)
    :param draft: `pysam.FastaFile` or filepath of draft sequence.

    :returns: iterable of (name, info, seq)
    """
    if isinstance(draft, str):
        draft = pysam.FastaFile(draft)

    def parse_info(info):
        d = medaka.common.Sample.decode_sample_name(info)
        # start and end are string repr of floats (major.minor coordinates)
        start, end = int(float(d['start'])), int(float(d['end']))
        return d['ref_name'], start, end

    contig_trees = collections.defaultdict(intervaltree.IntervalTree)
    ordered_contigs = collections.OrderedDict()
    for name, info, seq in contigs:
        ref_name, start, end = parse_info(info)
        # add one to end of interval, as intervaltree intervals and bed file
        # intervals are end-exclusive (i.e. they don't contain the last
        # coordinate), whilst the last position in a sample is included.
        contig_trees[ref_name].addi(start, end + 1, data=seq)
        ordered_contigs[ref_name] = None

    contig_lengths = dict(zip(draft.references, draft.lengths))
    contig_lengths = {k: v for k, v in contig_lengths.items() if k in
                      contig_trees}
    gap_trees = medaka.common.complement_intervaltrees(contig_trees,
                                                       contig_lengths)
    stitched_contigs = []
    for ref_name in ordered_contigs:
        contig_trees[ref_name].update(gap_trees[ref_name])
        draft_seq = draft.fetch(ref_name)
        pieces = []
        for i in sorted(contig_trees[ref_name],
                        key=operator.attrgetter('begin')):
            if i.data is None:  # this is a gap
                seq = draft_seq[i.begin: i.end]
            else:
                seq = i.data
            pieces.append(seq)
        stitched_contigs.append((ref_name, ref_name, ''.join(pieces)))
    return stitched_contigs, gap_trees


def _stitcher(inputs, draft, ref_names):
    contigs = stitch_from_probs(inputs, ref_names)
    if draft is not None:
        return fill_gaps(contigs, draft)
    else:
        return contigs, None


def stitch(args):
    """Entry point for stitching program."""
    index = medaka.datastore.DataIndex(args.inputs)
    if args.regions is None:
        args.regions = index.regions

    # batch size is a simple empirical heuristic
    rgrps = list(medaka.common.grouper(args.regions,
                 batch_size=max(1, len(args.regions) // (2 * args.threads))))
    gap_trees = {}
    with open(args.output, 'w') as fasta:
        Executor = concurrent.futures.ProcessPoolExecutor
        with Executor(max_workers=args.threads) as executor:
            worker = functools.partial(_stitcher, args.inputs, args.draft)
            for contigs, gap_tree in executor.map(worker, rgrps):
                for name, info, seq in contigs:
                    print('Writing {} {}'.format(name, info))
                    fasta.write('>{} {}\n{}\n'.format(name, info, seq))
                if gap_tree is not None:
                    gap_trees.update(gap_tree)

    if args.draft is not None:
        bed_out = args.output + '.gaps_in_draft_coords.bed'
        medaka.common.write_intervaltrees_to_bed(gap_trees, bed_out)
