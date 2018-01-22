import functools
import h5py
import itertools
import numpy as np
import os
import pandas as pd
import pysam
import re
import subprocess
from typing import Tuple

from collections import namedtuple
from medaka.common import segment_limits, rle, sliding_window, get_common_index
from medaka.common import _gap_, _ref_gap_, _read_sep_, decoding, encoding, Pileup, LabelledPileup
from medaka.labels import TruthAlignment

from timeit import default_timer as now

import logging
logger = logging.getLogger(__name__)

samtools = 'samtools'  # we place into the virtual environment

# String output of running samtools tview
RawTView = namedtuple('TView', ['coords', 'ref_seq', 'consensus', 'pileup', 'start', 'end'])

class TViewException(Exception):
    pass


class TView(object):

    def __init__(self, bam:str, reference:str, columns:int=20000, cache:int=None, minor_gaps:bool=False):
        """Interface to samtools tview functionality.

        :param bam: path to bam file of reads aligned to a common
            reference.
        :param reference: path to fasta file of reference to which reads
            are aligned.
        :param columns: default number of tview columns to read in one go.
        :param cache: Currently not implemented.
        :param minor_gaps: encoded gaps in non-reference positions distinctly.
        """
        self.logger = logging.getLogger('TView')
        self._bam = bam
        self._reference = reference
        self._columns = columns
        self._cache = cache
        self._minor_gaps = minor_gaps
        if self._minor_gaps:
            raise NotImplementedError
        self._columns_per_base = 3 # reasonable estimate
        self._column_pad = 1.1 # multiplier when _columns_per_base is updated

        # crosscheck the inputs
        with pysam.AlignmentFile(bam, 'rb') as b:
            brefs = b.references
            blens = b.lengths
        rrefs = set(pysam.FastaFile(reference).references)
        sbrefs = set(brefs)
        if not sbrefs.issuperset(rrefs):
            raise ValueError("input bam and reference do not appear consistent: {} vs {}.".format(sbrefs, rrefs))
        self._refs = set(brefs)
        self._reflens = dict(zip(brefs, blens))

        self._pileup = None # buffered `Pileup` object


    def _run_tview(self, ref_name, start, columns):
        """Run samtools tview."""
        def _find_coords(coord_line, ref_line):
            # Find the start and end ref coordinate of a tview chunk
            # from the coordinates line and reference line.
            coords_iter = re.finditer(r'\S+', coord_line)
            first = next(coords_iter)
            last = None
            for last in coords_iter:
                pass
            # return 0-based numbers coords
            first_coord = int(first.group()) - len(ref_line[:first.start()].replace('*', '')) - 1
            last_coord = int(last.group()) + len(ref_line[last.start() + 1:].replace('*', '')) - 1
            return first_coord, last_coord


        # tview is 1 based
        os.environ['COLUMNS'] = str(columns)
        cmd = [samtools, 'tview', self._bam, self._reference, '-d', 'T', '-p', '{}:{}'.format(ref_name, start + 1)]
        t0 = now()
        output = subprocess.check_output(cmd).decode()
        t1 = now()
        self.logger.info("Ran tview for {}, COLUMNS={} \t({:.3f}s)".format(start, columns, t1 - t0))

        # first line contains coords, then reference, then consensus
        try:
            coords, ref, consensus, remainder = output.split("\n", 3)
        except:
            raise TViewException("Tview did not output any reads.")
        start, end = _find_coords(coords, ref)
        pileup = remainder.split("\n")
        return RawTView(coords, ref, consensus, pileup, start, end)


    def _enqueue_region(self, ref_name:str, start:int=0, end:int=None):
        """Update internal buffer of tview output.

        :param ref_name: name of reference to which reads are aligned.
        :param start: zero-based start position in reference coordinates.
        :param end: zero-based end position in reference coordinates.

        :returns: self, updates `self._pileup` with fresh `Pileup` instance
        """
        if end is None:
            end = self._rlens[ref_name] - 1

        actual_end = end - 1
        it = 0
        tview = None
        while actual_end < end:
            it += 1
            columns = int(self._columns_per_base * (end - start))
            tview = self._run_tview(ref_name, start, columns)
            actual_end = tview.end # actually this is the last number in the header line
            self._columns_per_base = self._column_pad * columns / float(actual_end - start)
        self.logger.debug("Grabbing tview region required {} iterations.".format(it))

        t0 = now()
        reads, positions, ref_seq = self._tview_to_numpy(tview, end)
        t1 = now()
        self.logger.debug(
            "Parsed tview output for {}:{}-{}. Shape: {}.\t({:.3f}s)".format(
            self._bam, start, end, reads.shape, t1 - t0)
        )
        assert len(ref_seq) == len(positions)
        self._ref_seq = ref_seq
        self._pileup = Pileup(self._bam, ref_name, reads, positions)
        return self


    def _tview_to_numpy(self, tview, end:int):
        """From `TView` text output create a numpy array encoding and
        position index.

        :param tview: a `TView` object.
        :param end: last reference coordinate to parse.

        :returns: pileup, positions;
            pileup: numpy array shape (depth, len(positions)),
            positions: structured numpy array with `major` and `minor` reference position columns.

        """

        first_coord, last_coord = tview.start, tview.end
        ref_seq = tview.ref_seq

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
        end_i = np.searchsorted(ref_pos['major'], end, side='left')
        ref_seq_ar = np.array([i for i in ref_seq], dtype='U1')
        ref_pos = ref_pos[:end_i]
        ref_seq = tview.ref_seq[:end_i]
        ref_seq_ar = ref_seq_ar[:end_i]
        ref_gaps = ref_gaps[np.where(ref_gaps < end_i)].astype(int)

        reads = []
        for r in (x[:end_i] for x in tview.pileup):
            r_set = set(r)
            if r_set == set(_read_sep_) or len(r) == 0:  # after trimming we could have a blank line
                continue
            if not r_set.issubset(decoding):
                msg = ( "Got unexpected tview read characters: {}. " +
                        "\nYou might want to try removing non-primary reads (samtools view -F 2038)." )
                raise ValueError(msg.format(r_set.difference(decoding)))
            read = np.fromiter(
                # (encoding[b] for b in r.upper()),  # if we want both forward and reverse strands in caps.
                (encoding[b] for b in r),
                dtype='uint32', count=len(r)
            )
            # swap out gaps which also occur in the ref
            if self._minor_gaps:
                swap = np.where(read[ref_gaps] == encoding[_gap_])[0]
                read[ref_gaps[swap]] = encoding[_ref_gap_]
            reads.append(read)

        try:
            reads = np.stack(reads)
        except:
            raise TViewException("No reads (lines) after conversion to numpy.")
        return reads, ref_pos, ref_seq_ar


    def pileup(self, ref_name:str, start:int=0, end:int=None):
        """Obtain the read pileup spanning a reference region as a numpy array,
        along with some meta information.

        :param ref_name: name of reference to which reads are aligned.
        :param start: start position in reference coordinates.
        :param end: end position in reference coordinates.

        :returns: `Pileup` tuples
        """
        #TODO: enable caching
        self._enqueue_region(ref_name, start=start, end=end)
        return self._pileup, self._ref_seq


class MultiTView(object):

    def __init__(self, bams:Tuple[str], reference:str, columns:int=20000, cache:int=None):
        """Produce/merge pileup data from multiple bams.

        :param bam: path to bam file(s) of reads aligned to a common
            reference.
        :param reference: path to fasta file(s) of reference to which reads
            are aligned.
        :param columns: default number of tview columns to read in one go.
        """
        self.tviews = tuple(TView(b, reference, columns) for b in bams)
        self.n_inputs = len(bams)


    def pileup(self, ref_name:str, start:int=0, end:int=None):
        """Produce aligned pileups spanning a reference region as a numpy
        array, along with some meta information.

        :param ref_name: name of reference to which reads are aligned.
        :param start: start position in reference coordinates.
        :param end: end position in reference coordinates.

        :returns: a tuple of `Pileup` tuples

        """
        # TODO: caching the reindexed regions
        pileups, ref_seqs = zip(*(tv.pileup(ref_name, start=start, end=end) for tv in self.tviews))

        if self.n_inputs > 1:
            positions = get_common_index([p.positions for p in pileups])
            # all pileups should be aligned to same ref, so just use the first
            ref_p = Pileup(reads=ref_seqs[0], positions=pileups[0].positions)
            ref_seq = reindex_labels(ref_p, positions)
            pileups = (reindex_pileup(p, positions, self._minor_gaps) for p in pileups)
        else:
            ref_seq = ref_seqs[0]

        pileups = tuple(pileups)
        return pileups, ref_seq


def reindex_pileup(pileup, new_positions, minor_gaps=False):
    """Reindex pileup and fix any gaps.

    :param pileup: Pileup obj with `reads` and `positions` attributes.
    :param new_positions: iterable of tuples (e.g. structured numpy array) of
        new minor and major positions.
    :param minor_gaps: used distinct encoding for gaps in minor positions.

    :returns: Pileup obj

    """

    reindexed = np.empty((pileup.reads.shape[0], len(new_positions)))
    reindexed.fill(np.nan)
    inds = np.searchsorted(new_positions, pileup.positions)
    reindexed[:, inds] = pileup.reads

    # replace all np.nan with _gap_
    # col is a 'row' of the tview pileup (1 read or more with spaces between).
    is_minor_pos = new_positions['minor'] > 0

    for i, row in enumerate(reindexed):
        if minor_gaps:
            is_minor_pos_and_nan = np.logical_and(is_minor_pos, np.isnan(row))
            reindexed[i, np.where(is_minor_pos_and_nan)] = encoding[_ref_gap_]
        # fill all remaining np.nan with _gap_
        reindexed[i, np.where(np.isnan(reindexed[i]))] = encoding[_gap_]
    # convert back to pileup.dtype (we used nan above)
    reindexed = reindexed.astype(pileup.reads.dtype)

    # join up / expand read separators (relies on equality, hence done after casting back)
    for i, row in enumerate(reindexed):
        reindexed[i] = join_read_seps(reindexed[i], encoding[_read_sep_], encoding[_gap_])

    return Pileup(pileup.bam, pileup.ref_name, reindexed, new_positions)


def join_read_seps(pileup_row, read_sep, gap_val):
    """Join up `read_sep`s converting any surrounding `gap_val`s with `read_sep`s.
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

    # TODO: support multiple gap_vals in case we want to clean up _gap_ and _ref_gap_

    # can't pad with unicode str arrays, so use np.concatenate instead
    pad = np.array([gap_val], dtype=pileup_row.dtype)

    if read_sep in pileup_row and set(pileup_row) != set((read_sep,)):
        pileup_row = np.concatenate((pad, pileup_row, pad))
        blocks = rle(pileup_row)
        # create rolling view of blocks before, central, after (width 3)
        shape = (len(blocks) - 3 + 1, 3)
        strides = blocks.strides + (blocks.strides[-1],)
        sliding = np.lib.stride_tricks.as_strided(blocks, shape=shape, strides=strides)
        # find windows with central blank...
        central_is_blank = sliding[:, 1]['value'] == read_sep
        for contx in (0, 2):
            # ...and where before (/after) is gap...
            contx_is_gap = sliding[:, contx]['value'] == gap_val
            # ...replace the intersection
            replace = np.logical_and(central_is_blank, contx_is_gap)
            for window in sliding[replace]:
                con = window[contx]
                pileup_row[con['start']: con['start'] + con['length']] = read_sep

        pileup_row = pileup_row[len(pad):-len(pad)]
    return pileup_row


def reindex_labels(truth, new_positions, max_label_len=None):
    """Generate labels from a truth pileup and positions (of some other
    pileup).

    :param truth: `Pileup` object
    :param pileups: list of `Pileup` objects
    :param new_positions: structured array of new positions with 'major' and
        'minor' fields.
    :param max_label_len: int, where we have no evidence in the pileup (and
        hence ref) for a base which is present in the truth, we want the
        label for the previous base to include the missing bases up to
        `max_label_len`.

    :returns: np.array of labels
    """
    # get common index of truth and reads
    all_pos = get_common_index((new_positions, truth.positions))

    # reindex truth to common index
    truth = reindex_pileup(truth, all_pos)
    labels = truth.reads[0]
    if not isinstance(labels[0], str):
        labels = [decoding[i].upper() for i in labels]

    # Where we have no evidence in the pileup (and hence ref draft assembly)
    # for a base which is present in the truth, we want the label for the
    # previous base to include the missing bases. When this happens,
    # we will have a truth minor position not present in the reads.

    # Create index of pileup
    pileup_index = pd.Index(new_positions)
    labels = pd.DataFrame(labels, index=all_pos)
    # loop over minor positions in truth
    for (major, minor) in truth.positions[np.where(truth.positions['minor'] > 0)]:
        if (major, minor) not in pileup_index:
            # Append this (additional) label to previous and remove
            labels.iloc[labels.index.get_loc((major, 0))] += labels.iloc[labels.index.get_loc((major, minor))]
            labels.drop((major, minor), inplace=True)

    # convert labels to numpy array, truncating labels to max_label_len.
    if max_label_len is None:
        max_label_len = ''
    labels = labels.to_records(index=False).astype('U{}'.format(max_label_len))
    assert len(labels) == len(new_positions)
    return labels


def generate_pileup_chunks(bams, ref_fasta, ref_name, start=0, end=None,
                           chunk_len=50000, overlap=1000):
    """Yield chunks of pileup(s).

    :param bams: iterable of input .bam files.
    :param ref_fasta: input .fasta file.
    :param ref_name: name of reference within .bam to parse.
    :param start: start reference coordinate.
    :param end: end reference coordinate. If `None` the full extent of
        the reference will be parsed.
    :param chunk_len: int, chunk length in reference coordinates.
    :param overlap: int, overlap of adjacent chunks in reference coordinates.

    :yields: (tuple `Pileup` objects, None)

    ..note:: the None in the yielded values is for compatibility with an
        old API and will be removed.
    """
    tview = MultiTView(bams, ref_fasta, columns=chunk_len)

    if end is None:
        align_fhs = (pysam.AlignmentFile(bam) for bam in bams)
        end = min([dict(zip(p.references, p.lengths))[ref_name] for p in align_fhs])

    logger.info("Chunking {}: {}-{} in chunks of {} overlapping by {}".format(ref_name, start, end, chunk_len, overlap))

    #TODO: we could change how this is done since the class could cache
    #      everything in the entire region of interest.
    for chunk_start, chunk_end in segment_limits(start, end, segment_len=chunk_len, overlap_len=overlap):
        try:
            pileups, ref_seq = tview.pileup(ref_name, chunk_start, chunk_end)
        except TViewException:
            logger.info("Skipping region {}-{} as TView did not find any reads".format(chunk_start, chunk_end))
        else:
            yield LabelledPileup(pileups=pileups, labels=None, ref_seq=ref_seq)


def generate_training_chunks(truth_bam, bams, ref_fasta, ref_name, start=0, end=None,
                             chunk_len=50000, overlap=0):
    """Prepare training data chunks.

    :param truth_bam: .bam file of truth aligned to ref to generate labels.
    :param bams: iterable of input .bam files.
    :param ref_fasta: input .fasta file.
    :param ref_name: name of reference within .bam to parse.
    :param start: start reference coordinate.
    :param end: end reference coordinate. If `None` the full extent of
        the reference will be parsed.
    :param chunk_len: int, chunk length in reference coordinates.
    :param overlap: int, overlap of adjacent chunks in reference coordinates.

    :yields: tuple of `Pileup` objects (one item for each input bam) for each
        chunk.

    .. note:: Chunks might be missing if `truth_bam` is provided and
        regions with multiple mappings were encountered.

    """
    tview = MultiTView(bams, ref_fasta, columns=chunk_len)

    # filter truth alignments to restrict ourselves to regions of the ref where the truth
    # in unambiguous
    alignments = TruthAlignment.bam_to_alignments(truth_bam, ref_name, start=start, end=end)
    filtered_alignments = TruthAlignment.filter_alignments(alignments, start=start, end=end)
    if len(filtered_alignments) == 0:
        raise RuntimeError("Filtering removed all alignments of truth to ref, cannot continue.")

    for aln in filtered_alignments:
        msg = "Chunking {}: {}-{} in chunks of {} overlapping by {}"
        logger.info(msg.format(ref_name, aln.start, aln.end, chunk_len, overlap))

        for chunk_start, chunk_end in segment_limits(aln.start, aln.end, segment_len=chunk_len, overlap_len=overlap):
            truth_chunk = aln.get_positions_and_labels(start=chunk_start, end=chunk_end)
            try:
                pileups, ref_seq = tview.pileup(ref_name, chunk_start, chunk_end)
            except TViewException:
                logger.info("Skipping region {}-{} as TView did not find any reads".format(chunk_start, chunk_end))
                continue

            assert truth_chunk.positions[0]['major'] == pileups[0].positions[0]['major']
            assert truth_chunk.positions[-1]['major'] == pileups[0].positions[-1]['major']

            # Create labels according to positions in pileup
            t0 = now()
            labels = reindex_labels(truth_chunk, pileups[0].positions)
            t1 = now()
            logger.info("Processed labels for {}-{} \t({:.3f}s)".format(chunk_start, chunk_end, t1 - t0))

            yield LabelledPileup(pileups=pileups, labels=labels, ref_seq=ref_seq)


def load_pileup(pileup_hdf):
    """Load hdf pileup yielding chunks of pileup.

    :param hdf: input .hdf file.

    :yields: `(pileups, labels)` for each chunk;
        `pileups`: a list of `Pileup` objects
        `labels`: array of truth labels or `None`.
    """
    #TODO: decide whether this should be for training only (i.e. always with truth labels)

    data_types = None
    with h5py.File(pileup_hdf, 'r') as hdf:
        # hdf file should have pileups for each data type in /<ref>_<start>_<end</datatype<N>/data
        # and labels in /<ref>_<start>_<end</labels
        for chunk_str in hdf:  # pileup_<start>_<end>

            if 'labels' not in hdf[chunk_str]:
                labels = None
            else:
                labels = np.char.decode(hdf[chunk_str]['labels'][()])

            if 'ref_seq' not in hdf[chunk_str]:
                ref_seq = None
            else:
                ref_seq = np.char.decode(hdf[chunk_str]['ref_seq'][()])

            if data_types is None:
                data_types = [d for d in hdf[chunk_str] if d.startswith('datatype')]

            pileups = []
            for dtype in data_types:
                grp = hdf[chunk_str][dtype]
                ref_name = grp.attrs['rname']
                bam = grp.attrs['bam']
                reads = grp['data'][()]
                positions = grp['positions'][()]
                pileups.append(Pileup(bam, ref_name, reads, positions))

            yield LabelledPileup(pileups=pileups, labels=labels, ref_seq=ref_seq)


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

    for c in pileup_gen:

        read_chunkers = [chunker(p.reads) for p in c.pileups]
        # all pileups should have same positions so we just use first
        pos_chunker = chunker(c.pileups[0].positions)

        if c.labels is None:
            labels_gen = itertools.repeat(None)
        else:
            labels_gen = chunker(c.labels)

        if c.ref_seq is None:
            ref_seq_gen = itertools.repeat(None)
        else:
            ref_seq_gen = chunker(c.ref_seq)

        for pos in pos_chunker:
            chunk_labels = next(labels_gen)
            chunk_ref_seq = next(ref_seq_gen)
            chunk_pileups = tuple([Pileup(p.bam, p.ref_name, next(rc), pos) for p, rc in zip(c.pileups, read_chunkers)])
            msg = "Chunking positions {}:{} ({} bases)"
            logger.debug(msg.format(pos['major'][0], pos['major'][-1], len(pos)))
            yield LabelledPileup(pileups=chunk_pileups, labels=chunk_labels, ref_seq=chunk_ref_seq)
