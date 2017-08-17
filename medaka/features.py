from collections import defaultdict
import numpy as np

from medaka.tview import _read_sep_, encoding, Pileup
import logging
logger = logging.getLogger(__name__)


def coverage_filter(pileup_gen, depth=50, max_pcnt_gap=20.0):
    """Sort rows by number of alignment separators (fewest at top), remove any
    with more than `max_pcnt_gap` percent gap values, and return
    `depth` rows (or all rows if `depth` is -1), skipping chunks with insufficient depth.

    :param pileup_gen: generator of (pileups, labels) tuples;
        pileups: list of `Pileup` objects
        labels: np.array of labels or None.
    :param depth: int, number of rows to return (if None, return all rows).
    :param max_pcnt_gap: float, upper limit of percentage gaps per row tolerated.

    :yields:
    """
    logging.info("Coverage filter: depth: {} max_pcnt_gap: {}".format(depth, max_pcnt_gap))
    sep = encoding[_read_sep_]
    for pileups, labels in pileup_gen:
        deep_enough = True
        filtered = []  # list of filtered pileups
        for p in pileups:
            row_cnts = [defaultdict(int, zip(*np.unique(row, return_counts=True))) for row in p.reads]
            gaps_per_row = np.array([c[sep] for c in row_cnts])
            sort_inds = np.argsort(gaps_per_row)
            reads = p.reads[sort_inds]
            if max_pcnt_gap < 100.0:
                max_n_gaps = int(p.reads.shape[1] * max_pcnt_gap/100.0)
                gaps_per_row = gaps_per_row[sort_inds]
                cut_i = np.searchsorted(gaps_per_row, max_n_gaps, side='right')
                reads = reads[:cut_i]
            if reads.shape[0] < depth:
                logging.info("Skipping {}-{} due to lack of depth.".format(p.positions[0], p.positions[-1]))
                deep_enough = False
                break
            else:
                reads = reads[:depth]
            filtered.append(Pileup(p.bam, p.ref_name, reads, p.positions))
        if deep_enough:
            yield filtered, labels


def counts(pileups):
    """Get counts features for a pileup.

    :param pileups: iterable of `Pileup` objects

    :returns: np.array of counts features
    """
    symbols = encoding.values()
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
