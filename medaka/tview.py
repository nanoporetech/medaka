import functools
import h5py
import itertools
import numpy as np
import os
import pandas as pd
import pysam
import re
import subprocess

from collections import OrderedDict, namedtuple, defaultdict
from medaka.common import segment_limits, rle, sliding_window, get_common_index
from timeit import default_timer as now

import logging
logger = logging.getLogger(__name__)

samtools = 'samtools'  # we place into the virtual environment

# Codec for converting tview output to ints. Note that '#' here is used
#   for read positions which a both a gap ('*') in a read and reference.
_ref_gap_ = '#'
_gap_ = '*'
_read_sep_ = ' '
decoding = _ref_gap_ + _gap_ + 'acgtACGT' + _read_sep_
# store encoding in ordered dict as the order will always be the same
# (we need order to be the same for training and inference)
encoding = OrderedDict(((a, i) for i, a in enumerate(decoding)))

# String output of running samtools tview
TView = namedtuple('TView', ['coords', 'ref_seq', 'consensus', 'pileup', 'start', 'end', 'end_col'])

# Encoded (and possibly reindexed) pileup
Pileup = namedtuple('Pileup', ['bam', 'ref_name', 'reads', 'positions'])


def run_tview(bam, ref_fasta, ref_name, start, columns):
    """Run samtools tview for a given region and columns.

    :param bam: str, path to bam file of reads aligned to a common reference.
    :param ref_fasta: str, path to fasta file of reference to which reads are aligned.
    :param ref_name: str, name of reference to which reads are aligned.
    :param start: int, start position in reference coordinates.
    :param columns: int, number of columns to output.

    :returns: `TView` object.
    """

    def _find_coords(line):
        # Find the columns of the first and last coordinate markers from tview line.
        coords_iter = re.finditer(r'\S+', line)
        first = next(coords_iter)
        last = None
        for last in coords_iter:
            pass
        return int(first.group()) - 1, int(last.group()) - 1, last.start()

    os.environ['COLUMNS'] = str(columns)
    cmd = [samtools, 'tview', bam, ref_fasta, '-d', 'T', '-p', '{}:{}'.format(ref_name, start + 1)]
    t0 = now()
    output = subprocess.check_output(cmd).decode()
    t1 = now()
    logger.info("Ran tview for {}:{}, COLUMNS={} \t({:.3f}s)".format(bam, start, columns, t1 - t0))

    # first line contains coords, then reference, then consensus
    coords, ref, consensus, remainder = output.split("\n", 3)
    start, end, end_col = _find_coords(coords)
    pileup = remainder.split("\n")
    return TView(coords, ref, consensus, pileup, start, end, end_col)


def tview_to_numpy(tview, end):
    """From `TView` text output create a numpy array encoding and
    position index.

    :param tview: a `TView` object.
    :param end: int, last reference coordinate to parse.

    :returns: pileup, positions;
        pileup: numpy array shape (depth, len(positions)),
        positions: structured numpy array with `major` and `minor` reference position columns.

    """

    first_coord, last_coord, last_pos = tview.start, tview.end, tview.end_col

    coords = tview.coords[:last_pos]
    ref_seq = tview.ref_seq[:last_pos]

    # Record positions in reference in a (major, minor) format. Minor
    #   positions do not occur in the reference. Record positions of gaps
    #   in ref for use below - any matching gaps in reads will
    #   be encoded distinctly from deletions wrt to ref.
    pos = []
    ref_gaps = []
    major, minor = first_coord - 1, 0  # -1, so first position will be 0
    for i, c in enumerate(ref_seq):
        if c == _gap_:
            minor += 1
            ref_gaps.append(i)
        else:
            major += 1
            minor = 0
        pos.append((major, minor))
    ref_gaps = np.array(ref_gaps)
    ref_pos = np.array(pos, dtype=[('major', int), ('minor', int)])

    # at the end of the reference, if columns exceeds end-start, tview
    # will create fake positions and pad the reference with 'N's until
    # the number of requested columns are returned. Hence the need to
    # trim back to end.
    end_i = np.searchsorted(ref_pos['major'], end, side='right')
    ref_pos = ref_pos[:end_i]
    ref_seq = ref_seq[:end_i]
    ref_gaps = ref_gaps[np.where(ref_gaps < end_i)].astype(int)

    reads = []
    for r in (x[:end_i] for x in tview.pileup):
        if set(r) == set(_read_sep_) or len(r) == 0:  # after trimming we could have a blank line
            continue
        read = np.fromiter(
            # (encoding[b] for b in r.upper()),  # if we want both forward and reverse strands in caps.
            (encoding[b] for b in r),
            dtype='uint32', count=len(r)
        )
        # swap out gaps which also occur in the ref
        swap = np.where(read[ref_gaps] == encoding[_gap_])[0]
        read[ref_gaps[swap]] = encoding[_ref_gap_]
        reads.append(read)

    reads = np.stack(reads)
    return reads, ref_pos


def bam_to_pileup(bam, ref_fasta, ref_name, start, end, bases_per_ref_base=3.0, margin=1.2):
    """Obtain the read pileup spanning a reference region as a numpy array,
    along with some meta information.

    :param bam: str, path to bam file of reads aligned to a common reference.
    :param ref_fasta: str, path to fasta file of reference to which reads are aligned.
    :param ref_name: str, name of reference to which reads are aligned.
    :param start: int, start position in reference coordinates.
    :param end: int, end position in reference coordinates.
    :param bases_per_ref_base: float, used to estimate number of tview columns will be needed to span end-start.
    :param margin: float, fraction by which to increase estimated `bases_per_ref_base`.
        If `bases_per_ref_base` is underestimated, it is set to a value calculated from the shortfall in ref positions.
        `margin` is a multiplicative factor by which this calculated `bases_per_ref_base` is increased.

    :returns: (`Pileup` object, `bases_per_ref_base`)
    """
    actual_end = end - 1
    while actual_end < end:
        columns = int(bases_per_ref_base * (end - start))
        tview = run_tview(bam, ref_fasta, ref_name, start, columns)
        actual_end = tview.end
        bases_per_ref_base = margin * columns / float(actual_end - start)
        columns = int(bases_per_ref_base * (end - start))

    t0 = now()
    reads, positions = tview_to_numpy(tview, end)
    t1 = now()
    logger.info("Parsed tview output for {}:{}-{} \t({:.3f}s)".format(bam, start, end, t1 - t0))

    return Pileup(bam, ref_name, reads, positions), bases_per_ref_base


def reindex_pileup(pileup, new_positions):
    """Reindex pileup and fix any gaps.

    :param pileup: Pileup obj with `reads` and `positions` attributes.
    :param new_positions: iterable of tuples (e.g. structured numpy array) of new minor and major positions.

    :returns: Pileup obj

    """

    reindexed = np.empty((pileup.reads.shape[0], len(new_positions)))
    reindexed.fill(np.nan)
    inds = np.searchsorted(new_positions, pileup.positions)
    reindexed[:, inds] = pileup.reads

    # replace all np.nan in minor indices with _ref_gap_ and all other nans with _gap_
    # col is a 'row' of the tview pileup (1 read or more with spaces between).
    is_minor_pos = new_positions['minor'] > 0

    for i, row in enumerate(reindexed):
        is_minor_pos_and_nan = np.logical_and(is_minor_pos, np.isnan(row))
        reindexed[i, np.where(is_minor_pos_and_nan)] = encoding[_ref_gap_]
        # fill all remaining np.nan with _gap_
        reindexed[i, np.where(np.isnan(reindexed[i]))] = encoding[_gap_]
        # join up / expand read separators
        reindexed[i] = join_read_seps(reindexed[i], encoding[_read_sep_], encoding[_gap_])

    # after the merge, all encodings seem to be in float, convert back to pileup.dtype
    reindexed = reindexed.astype(pileup.reads.dtype)
    return Pileup(pileup.bam, pileup.ref_name, reindexed, new_positions)


def join_read_seps(pileup_row, read_sep, gap_val):
    """Join up blanks converting any surrounding `gap_val`s with `read_sep`s.
    i.e. in substrings containing only `read_sep` and `gap_val`, replace all
    `gap_val` with `read_sep`.

    :param pileup_row: np.array of input data.
    :param read_sep: permitted gap symbol.
    :param gap_val: symbol to be replaced by `read_sep`.
    :returns: np.array of modified pileup_row

    >>> ''.join(join_read_seps(np.array([ s for s in '*C*AA***  *  *** GGA']), ' ', '*'))
    '*C*AA            GGA'
    >>> ''.join(join_read_seps(np.array([ s for s in '*C*AA*** ']), ' ', '*'))
    '*C*AA    '
    >>> ''.join(join_read_seps(np.array([ s for s in '*C*AA***']), ' ', '*'))
    '*C*AA***'

    """
    # TODO: move to common, generalise variable names in line with description
    #      of algorithm in docstring.

    # can't pad with unicode str arrays, so use np.concatenate instead
    pad = np.array([gap_val], dtype=pileup_row.dtype)

    if read_sep in pileup_row and set(pileup_row) != set((read_sep,)):
        pileup_row = np.concatenate((pad, pileup_row, pad))
        blocks = rle(pileup_row)
        # create rolling view of blocks before, central, after (width 3)
        shape = (len(blocks) - 3 + 1, 3)
        strides = blocks.strides + (blocks.strides[-1],)
        sliding = np.lib.stride_tricks.as_strided(blocks, shape=shape, strides=strides)
        central_is_blank = sliding[:, 1]['value'] == read_sep
        for before, central, after in sliding[np.where(central_is_blank)]:
            for contx in before, after:
                if contx['value'] == gap_val:
                    pileup_row[contx['start']: contx['start'] + contx['length']] = read_sep

        pileup_row = pileup_row[len(pad):-len(pad)]
    return pileup_row


def _process_labels(truth, pileups):
    """Process labels
    :param truth: `Pileup` object
    :param pileups: list of `Pileup` objects
    :returns: np.array of labels
    """
    # get common index of truth and reads
    all_pos = get_common_index((pileups[0].positions, truth.positions)) 

    # reindex reads to common index, taking care of _gap_ and _read_sep_
    truth = reindex_pileup(truth, all_pos)
    labels = [decoding[i].upper() for i in truth.reads[0]]

    # Where we have no evidence in the pileup (and hence ref draft assembly)
    # for a base which is present in the truth, we want the label for the
    # proceeding base to include the missing bases, so e.g. [G] [*] would be e.g. [GC] [*].
    # When this happens, we will have a truth minor position not present in the reads.

    # create index of pileup (without truth, which has already been removed from pileups)
    pileup_index = pd.Index(pileups[0].positions)
    labels = pd.DataFrame(labels, index=all_pos)
    # loop over minor positions in truth
    for (major, minor) in truth.positions[np.where(truth.positions['minor'] > 0)]:
        if (major, minor) not in pileup_index:
            labels.iloc[labels.index.get_loc((major, 0))] += labels.iloc[labels.index.get_loc((major, minor))]
            labels.drop((major, minor), inplace=True)  # we don't want this position to appear in the labels
    # convert any _ref_gap_ labels to _gap_
    # TODO: do we need to do this?
    labels[labels[labels.columns[0]] == _ref_gap_] = _gap_

    # convert labels to numpy array
    labels = np.char.decode(labels.to_records(index=False).astype('S'))
    assert len(labels) == len(pileups[0].positions)
    return labels


def bams_to_pileup(bams, ref_fasta, ref_name, start, end, truth_bam=None, base_ratios=None):
    """Merge and optionally label multiple bams yielding chunks of pileup.

    :param bams: iterable of input .bam files.
    :param ref_fasta: input .fasta file.
    :param ref_name: name of reference within .bam to parse.
    :param start: start reference coordinate.
    :param end: end reference coordinate.
    :param truth_bam: .bam file of truth aligned to ref to generate labels
    :param base_ratios: {bam: bases_per_ref_pos} used to guess how many columns to run tview with.
    :returns: (pileups, labels, base_ratios);
        pileups: a list of `Pileup` objects
        labels: array of truth labels or `None`.
        base_ratios: dict of updated base ratios.
    """

    labels = None

    ratios = defaultdict(lambda: 1.0)
    if base_ratios is not None:
        ratios.update(base_ratios)

    pileups = []  # list of tview pileups objects for each bam

    # TODO: run this in parallel?
    for bam in bams:
        p, ratios[bam] = bam_to_pileup(bam, ref_fasta, ref_name, start, end,
                                       bases_per_ref_base=ratios[bam])
        pileups.append(p)

    if len(set([p.positions['major'][0] for p in pileups])) > 1:
        raise ValueError("Something went wrong with the tview call")

    # reindex reads to common index, taking care of _gap_ and _read_sep_
    if len(pileups) > 1:  # otherwise we don't need to reindex
        positions = get_common_index([p.positions for p in pileups])
        t0 = now()
        pileups = [reindex_pileup(p, positions) for p in pileups]
        t1 = now()
        logger.info("Re-indexed pileups for {}-{} \t({:.3f}s)".format(start, end, t1 - t0))

    if truth_bam is not None:
        truth, ratios[truth_bam] = bam_to_pileup(truth_bam, ref_fasta, ref_name, start, end,
                                                 bases_per_ref_base=ratios[bam])
        # TODO make less wasteful
        if truth.reads.shape[0] > 1:
            logging.warn("Unwilling to process truthset region with multiple mappings.")
            return None, ratios

        # reindex labels to common index, taking care of _gap_ and _read_sep_
        t0 = now()
        labels = _process_labels(truth, pileups)
        t1 = now()
        logger.info("Processed labels for {}-{} \t({:.3f}s)".format(start, end, t1 - t0))

    return pileups, labels, ratios


def generate_pileup_chunks(bams, ref_fasta, ref_name, start=0, end=None,
                           truth_bam=None, base_ratios=None, chunk_len=20000, overlap=500):
    """Yield of overlapping chunks of pileup generated by `bams_to_pileup`.

    :param bams: iterable of input .bam files.
    :param ref_fasta: input .fasta file.
    :param ref_name: name of reference within .bam to parse.
    :param start: start reference coordinate.
    :param end: end reference coordinate. If `None` the full extent of
        the reference will be parsed.
    :param truth_bam: .bam file of truth aligned to ref to generate labels.
    :param base_ratios: {bam: bases_per_ref_pos} used to guess how many columns to run tview with.
    :param chunk_len: int, chunk length in reference coordinates.
    :param overlap: int, overlap of adjacent chunks in reference coordinates.

    :yields: `(pileups, labels)` for each chunk;
        `pileups`: a list of `Pileup` objects
        `labels`: array of truth labels or `None`.

        Chunks might be missing if `truth_bam` is provided and
        regions with multiple mappings were encountered.

    """
    for bam in bams:
        if ref_name not in pysam.AlignmentFile(bam).references:
            raise ValueError("{} is not a reference within {}".format(ref_name, bam))

    if end is None:
        align_fhs = (pysam.AlignmentFile(bam) for bam in bams)
        end = min([dict(zip(p.references, p.lengths))[ref_name] for p in align_fhs])

    logger.info("Chunking {}: {}-{} in chunks of {} overlapping by {}".format(ref_name, start, end, chunk_len, overlap))

    for chunk_start, chunk_end in segment_limits(start, end, segment_len=chunk_len, overlap_len=overlap):
        *chunk, base_ratios = bams_to_pileup(bams, ref_fasta, ref_name, chunk_start, chunk_end,
                                            truth_bam=truth_bam, base_ratios=base_ratios)
        if chunk[0] is not None:  # None if we encountered multiple truth mappings
            yield chunk


def bam_pileup_to_hdf(hdf, bams, ref_fasta, ref_name, start=0, end=None, truth_bam=None,
                      chunk_len=20000, overlap=500, prefix='pileup'):
    """Create an .hdf file containing chunks of a pileup.

    :param hdf: output .hdf file (will be overwritten).
    :param bams: iterable of input .bam files.
    :param ref_fasta: input .fasta file.
    :param ref_name: name of reference within .bam to parse.
    :param start: start reference coordinate.
    :param end: end reference coordinate. If `None` the full extent of
        the reference will be parsed.
    :param truth_bam: .bam file of truth aligned to ref to generate labels.
    :param chunk_len: int, chunk length in reference coordinates.
    :param overlap: int, overlap of adjacent chunks in reference coordinates.
    :param prefix: prefix for hdf datasets. Datasets will be name as
        `prefix_{start}_{end}`, indicating the reference coordinates spanned.
    """
    gen = generate_pileup_chunks(bams, ref_fasta, ref_name=ref_name, start=start, end=end,
                                 truth_bam=truth_bam, chunk_len=chunk_len, overlap=overlap)

    with h5py.File(hdf, 'w') as hdf:
        for pileups, labels in gen:
            pos = pileups[0].positions
            chunk_str = '{}_{}_{}'.format(prefix, pos['major'][0], pos['major'][-1] + 1)
            if labels is not None:
                # save labels
                hdf['{}/labels'.format(chunk_str)] = np.char.encode(labels)

            for bam_count, p in enumerate(pileups):
                grp = '{}/datatype{}'.format(chunk_str, bam_count)
                hdf['{}/positions'.format(grp)] = p.positions
                hdf[grp].attrs['rname'] = p.ref_name
                hdf[grp].attrs['start'] = p.positions['major'][0]
                hdf[grp].attrs['end'] = p.positions['major'][-1]
                hdf[grp].attrs['bam'] = p.bam
                hdf['{}/data'.format(grp)] = p.reads


def load_pileup(pileup_hdf):
    """Load hdf pileup yielding chunks of pileup.

    :param hdf: input .hdf file.

    :yields: `(pileups, labels)` for each chunk;
        `pileups`: a list of `Pileup` objects
        `labels`: array of truth labels or `None`.
    """
    with h5py.File(pileup_hdf, 'r') as hdf:
        # hdf file should have pileups for each data type in /pileup_<start>_<end</datatype<N>/data
        # and labels in /pileup_<start>_<end</labels
        for chunk_str in hdf:  # pileup_<start>_<end>

            if 'labels' not in hdf[chunk_str]:
                labels = None
            else:
                labels = np.char.decode(hdf[chunk_str]['labels'][()])

            data_types = [d for d in hdf[chunk_str] if d.startswith('datatype')]

            pileups = []
            for bam_count, datatype in enumerate(data_types):
                grp = '{}/datatype{}'.format(chunk_str, bam_count)
                ref_name = hdf[grp].attrs['rname']
                bam = hdf[grp].attrs['bam']
                reads = hdf['{}/data'.format(grp)][()]
                positions = hdf['{}/positions'.format(grp)][()]
                pileups.append(Pileup(bam, ref_name, reads, positions))

            yield pileups, labels


def rechunk(pileup_gen, chunk_size=1000, overlap=200):
    """Rechunk chunks of pileup into smaller overlapping chunks.

    :param pileup_gen: generator of (pileups, labels);
        pileups: a list of `Pileup` objects
        labels: array of truth labels or `None`.
    :yields: (pileups, labels)
        pileups: a list of `Pileup` objects
        labels: array of truth labels or `None`.
    """
    chunker = functools.partial(sliding_window, window=chunk_size, step=chunk_size-overlap, axis=-1)

    logger.info("Rechunking into chunks of {} columns each overlapping by {} columns.".format(chunk_size, overlap))

    for pileups, labels in pileup_gen:
        read_chunkers = [chunker(p.reads) for p in pileups]
        # all pileups should have same positions so we just use first
        pos_chunker = chunker(pileups[0].positions)

        if labels is None:
            labels_gen = itertools.repeat(None)
        else:
            labels_gen = chunker(labels)

        for pos in pos_chunker:
            # logger.info("Chunking positions {}:{} ({} bases)".format(pos['major'][0], pos['major'][-1], len(pos)))
            yield [Pileup(p.bam, p.ref_name, next(rc), pos) for p, rc in zip(pileups, read_chunkers)], next(labels_gen)


def prepare(args):
    bam_pileup_to_hdf(
        args.output, (args.bam,), args.ref_fasta, args.ref_name, chunk_len=args.chunk_len, start=args.start, end=args.end,
        truth_bam=args.truth
    )
