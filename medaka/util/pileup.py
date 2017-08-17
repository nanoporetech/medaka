from collections import defaultdict, Counter, namedtuple
import itertools
import re
import sys
import numpy as np
import pysam


def get_reference_names(bam):
    """Obtain reference names

    :param bam: bam file
    :returns: list of reference names
    """
    idx_stats = pysam.idxstats(bam)
    ref_stats = [r.split('\t') for r in idx_stats.split('\n')[:-2]]
    ref_names = [r[0] for r in ref_stats]
    return ref_names


def get_reference_length(bam, ref):
    """Obtain length of a reference

    :param bam: bam file
    :param ref: name of reference
    :returns: int length of reference
    """
    idx_stats = pysam.idxstats(bam)
    ref_stats = [r.split('\t') for r in idx_stats.split('\n')[:-2]]
    ref_length = {r[0]: r[1] for r in ref_stats}
    return int(ref_length[ref])


def get_coverage(bam, ref, start=None, end=None):
    """Generate coverage stats across a region.

    :param bam: bam file
    :param ref: name of reference to process
    :param start: starting position within reference
    :param end: ending position within reference
    :yields: tuples of (position, coverage)
    """
    if start is None and end is None:
        depth = pysam.depth('-aa', bam)
    else:
        assert start is not None and end is not None, \
            "must provide neither or both start and end coords"
        pos = '{}:{}-{}'.format(ref, start, end)
        depth = pysam.depth('-aa', bam, '-r', pos)
    for line in (line.split('\t') for line in depth.split('\n')[:-1]):
        yield (int(p[1]) - 1, int(p[2]))


def no_coverage(bam, ref, start=None, end=None):
    """Obtain reference positions with 0 aligned query sequences

    :param bam: bam file
    :param ref: name of reference to process
    :param start: starting position within reference
    :param end: ending position within reference
    :returns: list of int reference positions without query sequence coverage
    """
    coverage = get_coverage(bam, ref, start=None, end=None)
    return [p[0] for p in coverage if p[1] == 0]


def multiple_coverage(bam, ref, start=None, end=None):
    """Returns list of reference positions with 2 or more
    aligned query sequences

    :param bam: bam file
    :param ref: name of reference to process
    :param start: starting position within reference
    :param end: ending position within reference
    :returns: list of int reference positions with multiple query
        sequence coverage
    """
    coverage = get_coverage(bam, ref, start=None, end=None)
    return [p[0] for p in coverage if p[1] > 1]


AlignPos = namedtuple('AlignPos', ('qpos', 'qbase', 'rpos', 'rbase'))


def feature_index(rev, base):
    """Defines integer encoding of base and del symbols

    :param rev: boolean read maps to reverse orientation
    :param base: {'A', 'C', 'G', 'T', None}. None indicates deletion.
    :returns: int encoded base
    """
    fwd_base_mapping = {'A': 1, 'C': 2, 'G': 3, 'T': 4, None: 5}
    rev_base_mapping = {'A': 6, 'C': 7, 'G': 8, 'T': 9, None: 10}
    if rev:
        return rev_base_mapping[base]
    return fwd_base_mapping[base]

# convert oriented base symbol to integer code
orient_base_to_int = {x: feature_index(*x)
                      for x in itertools.product((True, False),
                                                 list('ACGT') + [None])}

# convert integer code back to base symbols
int_to_orient_base = {code: ('d' if is_rev and base is None
                             else 'D' if base is None else
                             base.lower() if is_rev else
                             base)
                      for (is_rev, base), code in orient_base_to_int.items()}


def get_pairs(aln):
    """Return generator yielding AlignPos objects for
    aligned pairs of an AlignedSegment.

    :param aln: pysam AlignedSegment
    """
    seq = aln.query_sequence
    pairs = (AlignPos(qp, seq[qp], rp, rb) if qp is not None
             else AlignPos(qp, None, rp, rb)
             for qp, rp, rb in aln.get_aligned_pairs(with_seq=True))
    return pairs


def iterate_bam(bam, ref, start, end):
    """Iterate though alignment positions in the region of a bam.

    :param reads_bam: (sorted indexed) bam with read alignment to reference
    :param ref: name of reference to process
    :param start: starting position within reference
    :param end: ending position within reference
    :yields: tuple(`pysam.aligned_pair`, position, insertion index, is_reverse)
    """

    with pysam.AlignmentFile(bam, 'rb') as bamfile:
        aln_reads = bamfile.fetch(ref, start, end)

        if start is None:
            start = 0
        if end is None:
            end = float('Inf')

        for aln in aln_reads:
            reverse = aln.is_reverse
            pairs = get_pairs(aln)
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
                yield pair, current_pos, ins_count, reverse


def bam_to_feature_array(reads_bam, ref, start=None, end=None):
    """Converts a section of an alignment pileup (as shown
    by e.g. samtools tview) to a base frequency feature array

    :param reads_bam: (sorted indexed) bam with read alignment to reference
    :param ref: name of reference to process
    :param start: starting position within reference
    :param end: ending position within reference
    :returns: tuple(positions, feature_array)

        - positions: numpy structured array with 'ref_major'
          (reference position index) and 'ref_minor'
          (trailing insertion index) fields.

        - feature_array: 2D numpy array where each row is an
          11-element feature array for a pileup position:
          array([(encoded) reference base, nA+, nC+, nG+, 
          nT+, nDel+, nA-, nC-, nG-, nT-, nDel-])

          (e.g. nA+ = num. forward reads aligning with A
          to ref pos) (reference bases are encoded by
          :func:`feature_index`)

    :Example:

    ::

         Outputs for a minimal pileup with 3 aligned reads:

         reference sequence: AC**T
                      read1: ACG*T  (forward alignment)
                      read2: ACGGT  (forward alignment)
                      read3: ac**t  (reverse alignment)

         positions:

             ref_major  ref_minor
             1          0
             2          0
             2          1
             2          2
             3          0

         feature_array:

             array([[1, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0]   (2xA+, 1xA-)
                    [2  0, 2, 0, 0, 0, 0, 1, 0, 0, 0]   (2xC+, 1xC-)
                    [5  0, 0, 2, 0, 0, 0, 0, 0, 0, 1]   (2xG+, 1xDel-)
                    [5  0, 0, 1, 0, 1, 0, 0, 0, 0, 1]   (1xG+, 1xDel+, 1xDel-)
                    [4  0, 0, 0, 2, 0, 0, 0, 0, 1, 0]]) (1xT+, 1xT-)

    """
    aln_counters = defaultdict(Counter)
    ref_bases = dict()

    aln_info = iterate_bam(bam, ref, start, stop)
    for (pair, current_pos, ins_count, reverse) in aln_info:
        (aln_counters[(current_pos, ins_count)]
         [orient_base_to_int[reverse, pair.qbase]]) += 1
        ref_base = (orient_base_to_int[False, pair.rbase.upper()]
                     if pair.rbase else orient_base_to_int[False, None])
        ref_bases[(current_pos, ins_count)] = ref_base

    aln_cols = len(aln_counters)
    feature_len = len(orient_base_to_int) + 1
    feature_array = np.zeros(shape=(aln_cols, feature_len))
    positions = np.empty(aln_cols, dtype=[('ref_major', int),
                                          ('ref_minor', int)])

    for i, ((pos, counts), (_, ref_base)) in \
            enumerate(zip(sorted(aln_counters.items()),
                          sorted(ref_bases.items()))):
        positions[i] = pos
        feature_array[i, 0] = ref_base
        for j in counts.keys():
            feature_array[i, j] = counts[j]

    return positions, feature_array


def bam_to_label(truth_bam, ref, start=None, end=None):
    """Create encoded truth array.

    :param truth_bam: (sorted indexed) bam with true sequence aligned to reference
    :param ref: name of reference to process
    :param start: starting position within reference
    :param end: ending position within reference
    :returns: tuple(positions, encoded_label_array)

        - positions: numpy structured array with 'ref_major'
          (reference position index) and 'ref_minor'
          (trailing insertion index) fields.

        - feature_array: 1D numpy array of encoded labels
    """
    position_to_label = {}

    aln_info = iterate_bam(bam, ref, start, stop)
    for (pair, current_pos, ins_count, reverse) in aln_info:
        # Note: fwd/rev is collapsed here
        position_to_label[(current_pos, ins_count)] = (
            orient_base_to_int[False, pair.qbase.upper()]
            if pair.qbase else orient_base_to_int[False, None]
        )

    aln_cols = len(position_to_label)
    label_array = np.zeros(shape=(aln_cols, 1))
    positions = np.empty(aln_cols, dtype=[('ref_major', int),
                                          ('ref_minor', int)])

    for i, (pos, label) in enumerate(sorted(position_to_label.items())):
        positions[i] = pos
        label_array[i] = label

    return positions, label_array


def prepare_training_data(reads_bam, truth_bam, ref, limits=(None, None)):
    """Prepare data for model training

    :param reads_bam: (sorted indexed) bam with read alignment to reference
    :param truth_bam: (sorted indexed) bam with true sequence aligned to reference
    :param ref: name of reference to process
    :param limits: start and end points of reference
    :returns: tuple (filtered_features, labels)

        - filtered_features: 2D numpy array with base frequency features
          as described in :func:`bam_to_feature_array` for valid reference
          positions (where there is an unambiguous 'truth' label provided by
          the known sequence.

        - labels: 1D numpy array of integer encoded labels for each valid
          reference position.
    """
    start, end = limits
    read_pos, read_features = bam_to_feature_array(reads_bam, ref,
                                                   start, end)
    truth_pos, label_array = bam_to_label(truth_bam, ref, start, end)

    position_to_label = defaultdict(itertools.repeat(5).__next__,
                                    zip([tuple(p) for p in truth_pos],
                                        [int(a) for a in label_array]))

    no_label = no_coverage(truth_bam, ref, start, end)
    multi_label = multiple_coverage(truth_bam, ref, start, end)
    blacklist = np.array(no_label + multi_label)

    good_pos_mask = np.in1d(read_pos['ref_major'], blacklist, invert=True)
    filtered_features = read_features[good_pos_mask]
    filtered_read_pos = read_pos[good_pos_mask]

    valid_positions = ((pos[0], pos[1]) for pos in filtered_read_pos)
    labels = np.array([position_to_label[vp] for vp in valid_positions]).reshape(-1,1)

    return filtered_features, labels, filtered_read_pos


