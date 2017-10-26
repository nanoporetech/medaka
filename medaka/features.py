from collections import defaultdict
import numpy as np

from medaka.common import _read_sep_, decoding, encoding, Pileup, LabelledPileup
import logging
logger = logging.getLogger(__name__)


def pcnt_gap_filter(pileup_gen, depth=15, max_pcnt_gap=20.0):
    """Sort rows by number of alignment separators (fewest at top), remove any
    with more than `max_pcnt_gap` percent gap values, and return
    `depth` rows (or all rows if `depth` is None), skipping chunks with insufficient depth.

    :param pileup_gen: generator of `LabelledPileup` named tuples.
    :param depth: int, number of rows to return (if None, return all rows).
    :param min_depth: int, chunks with insufficient depth at any position will not be yielded.
    :param max_pcnt_gap: float, upper limit of percentage gaps per row tolerated.

    :yields: `LabelledPileup` named tuples.
    """
    logging.info("Coverage filter: depth: {} max_pcnt_gap: {}".format(depth, max_pcnt_gap))
    sep = encoding[_read_sep_]
    for c in pileup_gen:
        deep_enough = True
        filtered = []  # list of filtered pileups
        for p in c.pileups:
            row_cnts = [defaultdict(int, zip(*np.unique(row, return_counts=True))) for row in p.reads]
            gaps_per_row = np.array([c[sep] for c in row_cnts])
            sort_inds = np.argsort(gaps_per_row)
            reads = p.reads[sort_inds]
            if max_pcnt_gap < 100.0:
                max_n_gaps = int(p.reads.shape[1] * max_pcnt_gap/100.0)
                gaps_per_row = gaps_per_row[sort_inds]
                cut_i = np.searchsorted(gaps_per_row, max_n_gaps, side='right')
                reads = reads[:cut_i]
            if depth is not None and reads.shape[0] < depth:
                logging.info("Skipping {}-{} due to lack of depth.".format(p.positions[0], p.positions[-1]))
                deep_enough = False
                break
            else:
                reads = reads[:depth]
            filtered.append(Pileup(p.bam, p.ref_name, reads, p.positions))
        if deep_enough:
            yield LabelledPileup(pileups=filtered, labels=c.labels, ref_seq=c.ref_seq)


def depth_filter(pileup_gen, min_depth=15, max_depth=np.inf):
    """Calculate depth (excluding gaps) at every position and skip chunks with depth < `min_depth`.

    :param pileup_gen: generator of `LabelledPileup` named tuples.
    :param min_depth: int, chunks with insufficient depth at any position will not be yielded.
    :param max_depth: int, chunks with excessive depth at any position will not be yielded.

    :yields: `LabelledPileup` named tuples.
    """
    logging.info("Depth filter: min_depth: {} max_depth: {}".format(min_depth, max_depth))
    for c in pileup_gen:
        depths = np.sum(counts(c.pileups), axis=1)
        is_too_deep = np.any(depths > max_depth)
        is_too_shallow = np.any(depths < min_depth)
        if not (is_too_deep or is_too_shallow):
            yield c
        else:
            pos = c.pileups[0].positions
            if is_too_deep:
                msg = "Skipping {}:{}-{} due to excessive depth ({})."
                d = np.max(depths)
            else:
                msg = "Skipping {}:{}-{} due to low depth ({})."
                d = np.min(depths)
            logging.info(msg.format(c.pileups[0].ref_name, pos[0], pos[-1], d))


def alphabet_filter(pileup_gen, alphabet=None, filter_labels=True, filter_ref_seq=True):
    """Skip chunks in which labels and/or ref_seq contain bases not in `alphabet`.

    :param pileup_gen: generator of `LabelledPileup` named tuples.
    :param alphabet: set of str of allowed bases. If None, automatically generated from decoding.
    :param filter_labels: bool, whether to filter on labels.
    :param filter_ref_seq: bool, whether to filter on ref_seq.

    :yields: `LabelledPileup` named tuples.
    """
    if alphabet is None:
        alphabet = set([c for c in decoding.upper().replace(_read_sep_,'')])
    logging.info("Alphabet filter: alphabet: {}".format(alphabet))

    def _find_bad_bases(c, field, alphabet):
        # Labels could contain multi-labels e.g. ['G', 'TG', 'CCA']
        field_bases = set((i for i in ''.join(i for i in getattr(c, field))))
        if not field_bases.issubset(alphabet):
            msg = "Skipping {}:{} ({} bases) due to {} {}"
            pos = c.pileups[0].positions
            logger.info(msg.format(pos['major'][0], pos['major'][-1],
                                    len(pos), field, field_bases.difference(alphabet)))
            return True

    for c in pileup_gen:
        if filter_labels and c.labels is not None and _find_bad_bases(c, 'labels', alphabet):
            continue
        if filter_ref_seq and c.ref_seq is not None and _find_bad_bases(c, 'ref_seq', alphabet):
            continue
        yield c


def counts(pileups):
    """Get counts features for a pileup.

    :param pileups: iterable of `Pileup` objects

    :returns: np.array of counts features
    """
    # we don't want to a count _read_sep_ in our counts feature
    symbols = [v for k, v in encoding.items() if k != _read_sep_]
    n_features = len(pileups) * len(symbols)
    chunk_len = pileups[0].reads.shape[1]
    counts = np.zeros((chunk_len, n_features), dtype=pileups[0].reads.dtype)
    for i, p in enumerate(pileups):
        col_cnts = [defaultdict(int, zip(*np.unique(col, return_counts=True))) for col in p.reads.T]
        for j, s in enumerate(symbols):
            counts[:, i * len(symbols) + j] = [c[s] for c in col_cnts]
    return counts


def pileup_and_counts(pileups):
    """Get features comprising pileup and counts.

    :param pileups: iterable of `Pileup` objects

    :returns: np.array of counts and pileup features.
    """
    pileup_feature = np.hstack([p.reads.T for p in pileups])
    counts_feature = counts(pileups)
    return np.hstack((pileup_feature, counts_feature))
