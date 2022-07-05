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


def write_fastx_segment(fh, contig, qualities=True):
    """Write a fastx file from a tuple of (name, sequence, qualities).

    :param filename: output filehandle.
    :param contigs: tuple of the form (sequence name, base sequence,
        base qualities).

    """
    if qualities:
        fastx_prefix = '@'
    else:
        fastx_prefix = '>'
    fh.write('{}{}\n{}\n'.format(fastx_prefix, contig[0], ''.join(contig[1])))
    if qualities:
        fh.write('+\n{}\n'.format(''.join(contig[2])))


def stitch_from_probs(h5_fp, region, min_depth):
    """Join overlapping label probabilities from HDF5 files.

     Network outputs from multiple samples stored within a file are spliced
     together into a logically contiguous array and decoded to generate
     contiguous sequence(s).

    :param h5_fp: iterable of HDF5 filepaths.
    :param region: `medaka.common.Region` instance

    :returns: list of (region string, sequence, qualities).
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
    if min_depth:
        data_gen = medaka.common.Sample.filter_samples(
            data_gen, min_depth=min_depth)
    cur_ref_name = None
    heuristic_use = 0
    seq_parts = list()
    qualities = list()
    start = None
    contigs = []
    for s, is_last_in_contig, heuristic in data_gen:
        cur_ref_name = s.ref_name if cur_ref_name is None else cur_ref_name
        start = s.positions[0]['major'] if start is None else start
        seq, qual = label_scheme.decode_consensus(s, with_qualities=True)
        seq_parts.append(seq)
        qualities.append(qual)
        if is_last_in_contig:
            contigs.append((
                (s.ref_name, start, s.positions[-1]['major']),
                seq_parts,
                qualities))
            seq_parts = list()
            qualities = list()
            start = None
        heuristic_use += heuristic
    if len(seq_parts) > 0:
        contigs.append((
            (s.ref_name, start, s.positions[-1]['major']),
            seq_parts,
            qualities))

    logger.debug("Used heuristic {} times for {}.".format(
        heuristic_use, region))
    return contigs


def fill_gaps(contigs, draft, fill_char=None):
    """Fill gaps between contigs with draft sequence or designated char.

    :param contigs: iterable of ((ref_name, start, stop), sequence parts,
        qualities)
    :param draft: `pysam.FastaFile` or filepath of draft sequence.
    :param fill_char: character to fill gaps (default: None, fills gaps
        with content from the draft sequence)

    :returns: iterable of (name, info, seq)
    """
    if isinstance(draft, str):
        draft = pysam.FastaFile(draft)
    # interpret an empty string for fill_char as "fill-with-draft"
    fill_char = None if fill_char in (None, "") else str(fill_char)[0]

    contig_trees = collections.defaultdict(intervaltree.IntervalTree)
    ordered_contigs = collections.OrderedDict()
    for info, sequence_parts, qualities in contigs:
        ref_name, start, stop = info
        # add one to end of interval, as intervaltree intervals and bed file
        # intervals are end-exclusive (i.e. they don't contain the last
        # coordinate), whilst the last position in a sample is included.
        contig_trees[ref_name].addi(start, stop + 1, data=(
            sequence_parts, qualities))
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
        seq_pieces = []
        qual_pieces = []
        for i in sorted(
                contig_trees[ref_name], key=operator.attrgetter('begin')):
            if i.data is None:  # this is a gap
                if fill_char is None:
                    seq = [draft_seq[i.begin: i.end]]
                else:
                    seq = [fill_char * (i.end - i.begin)]
                qual = ['!' * (i.end - i.begin)]  # set Q=0 for padding bases
            else:
                seq = i.data[0]
                qual = i.data[1]
            seq_pieces.extend(seq)
            qual_pieces.extend(qual)
        # return as input
        stitched_contigs.append((
            (ref_name, 0, contig_lengths[ref_name]), seq_pieces,
            qual_pieces))
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
    qual = contig[2]
    for contig in contigs:
        c_rn, c_start, c_stop, = contig[0]
        if c_rn == ref_name and c_start == stop + 1:
            stop = c_stop
            buffer.extend(contig[1])
            qual.extend(contig[2])
        else:
            # clear buffer, start anew
            yield (
                (ref_name, start, stop), buffer, qual)
            ref_name, start, stop = contig[0]
            buffer = contig[1]
            qual = contig[2]
    yield ((ref_name, start, stop), buffer, qual)


def stitch(args):
    """Entry point for stitching program."""
    logger = medaka.common.get_named_logger("Stitcher")
    index_log = medaka.common.get_named_logger("DataIndex")
    index_log.setLevel(logging.WARNING)
    index = medaka.datastore.DataIndex(args.inputs)
    draft = pysam.FastaFile(args.draft)
    draft_lengths = dict(zip(draft.references, draft.lengths))

    # regions requested for processing
    req_regions = args.regions
    if req_regions is None:
        req_regions = {
            medaka.common.Region.from_string(r) for r in draft.references}

    # split up draft contigs into chunks for parallelism
    MAX_REGION_SIZE = int(1e6)
    regions_to_process = list()
    index_region_names = {r.ref_name for r in index.regions}
    for ref_name, start, end in req_regions:
        if ref_name not in index_region_names:
            # for clarity: don't try to process data which doesn't exist.
            # This shouldn't be strictly necessary, just for sanity.
            continue
        if start is None:
            start = 0
        if end is None:
            end = draft_lengths[ref_name]
        regions_to_process.append(medaka.common.Region(ref_name, start, end))
    # regions will be sections for which we have at least some data for
    # the corresponding contig
    regions = itertools.chain.from_iterable((
        r.split(MAX_REGION_SIZE, overlap=0, fixed_size=False)
        for r in regions_to_process))

    gap_trees = {}
    with open(args.output, 'w') as fastx:
        Executor = concurrent.futures.ProcessPoolExecutor
        with Executor(max_workers=args.threads) as executor:
            worker = functools.partial(
                stitch_from_probs, args.inputs,
                min_depth=args.min_depth)
            # we rely on map being ordered
            pieces = itertools.chain.from_iterable(
                executor.map(worker, regions))
            contigs = collapse_neighbours(pieces)
            if args.fillgaps:
                # TODO: ideally fill_gaps would be a generator
                contigs, gt = fill_gaps(contigs, args.draft, args.fill_char)
                gap_trees.update(gt)
                for (ref_name, start, stop), seq_parts, qualities in contigs:
                    write_fastx_segment(
                        fastx,
                        (ref_name, seq_parts, qualities),
                        qualities=args.qualities)
            else:
                ref_name = None
                counter = 0
                for (rname, start, stop), seq_parts, qualities in contigs:
                    if ref_name == rname:
                        counter += 1
                    else:
                        counter = 0
                    write_fastx_segment(
                        fastx,
                        ("{}_{} {}-{}".format(
                            rname, counter, start, stop + 1),
                         seq_parts, qualities),
                        qualities=args.qualities)
                    ref_name = rname
        # the data index doesn't contain entries for contigs that
        # were entirely skipped (and we explicitely removed these
        # from consideration above). If desired, go back and copy
        # theses across verbatim, and fill in an entry in the gap_trees
        if args.fillgaps:
            processed_regions = {r.ref_name for r in regions_to_process}
            required_regions = {r.ref_name for r in req_regions}
            missing = required_regions - processed_regions
            for reg in missing:
                logger.info(
                    "Copying contig '{}' verbatim from input.".format(reg))
                seq = draft.fetch(reg)
                write_fastx_segment(
                    fastx,
                    (reg, seq, '!' * len(seq)),
                    qualities=args.qualities)
                tree = intervaltree.IntervalTree()
                tree.addi(0, draft_lengths[reg])
                gap_trees[reg] = tree

    if args.fillgaps:
        bed_out = args.output + ".gaps_in_draft_coords.bed"
        medaka.common.write_intervaltrees_to_bed(gap_trees, bed_out)
