"""Creation of contiguous consensus sequences from chunked network outputs."""
import collections
import concurrent.futures
import functools
import itertools
import logging
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


def stitch_from_probs(h5_fp, region):
    """Join overlapping label probabilities from HDF5 files.

     Network outputs from multiple samples stored within a file are spliced
     together into a logically contiguous array and decoded to generate
     contiguous sequence(s).

    :param h5_fp: iterable of HDF5 filepaths.
    :param region: `medaka.common.Region` instance

    :returns: list of (region string, sequence).
    """
    logger = medaka.common.get_named_logger('Stitch')
    logger.debug("Stitching region: {}".format(str(region)))
    index = medaka.datastore.DataIndex(h5_fp)
    label_scheme = index.metadata['label_scheme']
    logger.debug("Label decoding is:\n{}".format(
        '\n'.join('{}: {}'.format(k, v)
                  for k, v in label_scheme._decoding.items())))

    samples = index.yield_from_feature_files(regions=[region])
    data_gen = medaka.common.Sample.trim_samples_to_region(
        samples, start=region.start, end=region.end)
    cur_ref_name = None
    heuristic_use = 0
    seq_parts = list()
    start = None
    contigs = []
    for s, is_last_in_contig, heuristic in data_gen:
        cur_ref_name = s.ref_name if cur_ref_name is None else cur_ref_name
        start = s.positions[0]['major'] if start is None else start
        seq_parts.append(label_scheme.decode_consensus(s))
        if is_last_in_contig:
            contigs.append((
                (s.ref_name, start, s.positions[-1]['major']),
                seq_parts))
            seq_parts = list()
            start = None
        heuristic_use += heuristic
    if len(seq_parts) > 0:
        contigs.append((
            (s.ref_name, start, s.positions[-1]['major']),
            seq_parts))

    logger.debug("Used heuristic {} times for {}.".format(
        heuristic_use, region))
    return contigs


def fill_gaps(contigs, draft):
    """Fill gaps between polished contigs with draft sequence.

    :param contigs: iterable of ((ref_name, start, stop), sequence parts)
    :param draft: `pysam.FastaFile` or filepath of draft sequence.

    :returns: iterable of (name, info, seq)
    """
    if isinstance(draft, str):
        draft = pysam.FastaFile(draft)

    contig_trees = collections.defaultdict(intervaltree.IntervalTree)
    ordered_contigs = collections.OrderedDict()
    for info, sequence_parts in contigs:
        ref_name, start, stop = info
        # add one to end of interval, as intervaltree intervals and bed file
        # intervals are end-exclusive (i.e. they don't contain the last
        # coordinate), whilst the last position in a sample is included.
        contig_trees[ref_name].addi(start, stop + 1, data=sequence_parts)
        ordered_contigs[ref_name] = None

    contig_lengths = dict(zip(draft.references, draft.lengths))
    contig_lengths = {
        k: v for k, v in contig_lengths.items()
        if k in contig_trees}
    gap_trees = medaka.common.complement_intervaltrees(
        contig_trees, contig_lengths)

    stitched_contigs = []
    for ref_name in ordered_contigs:
        contig_trees[ref_name].update(gap_trees[ref_name])
        draft_seq = draft.fetch(ref_name)
        pieces = []
        for i in sorted(
                contig_trees[ref_name], key=operator.attrgetter('begin')):
            if i.data is None:  # this is a gap
                seq = [draft_seq[i.begin: i.end]]
            else:
                seq = i.data
            pieces.extend(seq)
        # return as input
        stitched_contigs.append((
            (ref_name, 0, contig_lengths[ref_name]), pieces))
    return stitched_contigs, gap_trees


def collapse_neighbours(contigs):
    """Build larger contigs by joining neighbours.

    :param contigs: a stream of ordered (partial)-contigs.
    """
    try:
        contig = next(contigs)
    except StopIteration:
        return
    ref_name, start, stop = contig[0]
    buffer = contig[1]
    for contig in contigs:
        c_rn, c_start, c_stop, = contig[0]
        if c_rn == ref_name and c_start == stop + 1:
            stop = c_stop
            buffer.extend(contig[1])
        else:
            # clear buffer, start anew
            yield (
                (ref_name, start, stop), buffer)
            ref_name, start, stop = contig[0]
            buffer = contig[1]
    yield ((ref_name, start, stop), buffer)


def stitch(args):
    """Entry point for stitching program."""
    index_log = medaka.common.get_named_logger('DataIndex')
    index_log.setLevel(logging.WARNING)
    index = medaka.datastore.DataIndex(args.inputs)
    if args.regions is None:
        args.regions = index.regions

    # split up draft contigs into chunks for parallelism
    MAX_REGION_SIZE = int(1e6)
    draft = pysam.FastaFile(args.draft)
    draft_lengths = dict(zip(draft.references, draft.lengths))
    regions = list()
    for ref_name, start, end in args.regions:
        if start is None:
            start = 0
        if end is None:
            end = draft_lengths[ref_name]
        regions.append(medaka.common.Region(ref_name, start, end))
    args.regions = regions
    regions = itertools.chain.from_iterable((
        r.split(MAX_REGION_SIZE, overlap=0, fixed_size=False)
        for r in args.regions))

    gap_trees = {}
    with open(args.output, 'w') as fasta:
        Executor = concurrent.futures.ProcessPoolExecutor
        with Executor(max_workers=args.threads) as executor:
            worker = functools.partial(stitch_from_probs, args.inputs)
            # we rely on map being ordered
            pieces = itertools.chain.from_iterable(
                executor.map(worker, regions))
            contigs = collapse_neighbours(pieces)
            if args.fillgaps:
                # TODO: ideally fill_gaps would be a generator
                contigs, gt = fill_gaps(contigs, args.draft)
                gap_trees.update(gt)
                for (ref_name, start, stop), seq_parts in contigs:
                    fasta.write(">{}\n".format(ref_name))
                    for s in seq_parts:
                        fasta.write(s)
                    fasta.write("\n")
            else:
                ref_name = None
                counter = 0
                for (rname, start, stop), seq_parts in contigs:
                    if ref_name == rname:
                        counter += 1
                    else:
                        counter = 0
                    fasta.write(">{}_{} {}-{}\n".format(
                        rname, counter, start, stop + 1))
                    for s in seq_parts:
                        fasta.write(s)
                    fasta.write("\n")
                    ref_name = rname

    if args.draft is not None:
        bed_out = args.output + '.gaps_in_draft_coords.bed'
        medaka.common.write_intervaltrees_to_bed(gap_trees, bed_out)
