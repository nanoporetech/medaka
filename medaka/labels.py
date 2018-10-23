import itertools
import logging
import numpy as np
import pysam
from copy import copy
from operator import attrgetter
from medaka.common import _gap_, encoding
from medaka.common import get_pairs, yield_compressed_pairs, get_pairs_with_hp_len


class TruthAlignment(object):
    def __init__(self, alignment):
        """Create a TruthAlignment oblist from a `pysam.libcalignedsegment.AlignedSegment` object.

        :param alignment: `pysam.libcalignedsegment.AlignedSegment` object.
        """
        self.aln = alignment  # so we can get positions and labels later
        # initialise start and end (which might be moved)
        self.start = self.aln.reference_start  # zero-based
        self.end = self.aln.reference_end
        self.is_kept = True

    def get_overlap_with(self, other):
        first, second = sorted((self, other), key=attrgetter('aln.reference_start'))
        overlap = first.aln.reference_end > second.aln.reference_start
        if overlap:
            # return positions of start and end of overlapping region
            return second.aln.reference_start, first.aln.reference_end
        else:
            return None

    @staticmethod
    def filter_alignments(alignments, min_length=1000, length_ratio=2.0, overlap_fraction=0.5,
                         start=0, end=None):
        """Filter alignments to yield only segments suitable for training.

        :param alignments: iterable of `TruthAlignment` objects.
        :param min_length: int, minimum length of alignment segment.
        :param length_ratio: float, in cases of overlap, ratio of longer segment / shorter segment
            above which the longer segment is assumed to be more reliable.
        :param overlap_fraction: float, length_ratio: float, in cases of overlap, fraction of
            shorter segment overlapping with the longer segment above which the segments are
            considered highly overlapping.
        :param start: zero-based start position in reference coordinates. Any alignments which
            start before this value will take this new starting coordinate.
        :param end: zero-based end position in reference coordinates. and alignments
            which end after this value will take this new ending coordinate.
        :returns: list of `TruthAlignment` objects

        """
        filtered_alignments = [ copy(a) for a in alignments]  # don't want to modify original alignments
        for al_i, al_j in itertools.combinations(filtered_alignments, 2):
            first, second = sorted((al_i, al_j), key=attrgetter('aln.reference_start'))
            overlap = first.get_overlap_with(second)
            if overlap is None:
                continue
            ovlp_start, ovlp_end = overlap
            shorter, longer = sorted((al_i, al_j), key=attrgetter('aln.reference_length'))
            length_ratio_ij = longer.aln.reference_length / shorter.aln.reference_length
            overlap_fraction_ij = (ovlp_end - ovlp_start) / shorter.aln.reference_length
            # 4 cases
            if length_ratio_ij < length_ratio:  # we don't trust one more than the other
                if overlap_fraction_ij >= overlap_fraction:
                    # 1) they overlap a lot; we have significant ambiguity, discard both
                    shorter.is_kept = False
                    longer.is_kept = False
                else:
                    # 2) they overlap a little; just trim overlapping portions of both alignments
                    first.end = ovlp_start
                    second.start = ovlp_end
            else:  # we trust the longer one more than the shorter one
                if overlap_fraction_ij >= overlap_fraction:
                    # 3) they overlap a lot; discard shorter alignment
                    shorter.is_kept = False
                else:
                    # 4) they overlap a little; trim overlapping portion of shorter alignment
                    second.start = ovlp_end

        # trim starts and ends if needed
        if start > 0 or end is not None:
            for al in filtered_alignments:
                if start > 0:
                    al.start = max(start, al.start)
                if end is not None:
                    al.end = min(end, al.end)
        # do filtering
        filtered_alignments = [al for al in filtered_alignments
                               if (al.is_kept and al.end - al.start >= min_length)]
        filtered_alignments.sort(key=attrgetter('start'))
        return filtered_alignments

    @staticmethod
    def bam_to_alignments(truth_bam, ref_name, start=None, end=None):
        """Create list of TruthAlignment objects from a bam of Truth aligned to ref.

        :param truth_bam: (sorted indexed) bam with true sequence aligned to reference
        :param ref: name of reference to process
        :param start: starting position within reference
        :param end: ending position within reference
            (all alignments with any overlap with the interval start:end will be retrieved)
        :returns: tuple(positions, encoded_label_array)

            - positions: numpy structured array with 'ref_major'
              (reference position index) and 'ref_minor'
              (trailing insertion index) fields.

            - feature_array: 1D numpy array of encoded labels
        """
        with pysam.AlignmentFile(truth_bam, 'rb') as bamfile:
            aln_reads = bamfile.fetch(reference=ref_name, start=start, end=end)
            alignments = [TruthAlignment(r) for r in aln_reads if not (r.is_unmapped or r.is_secondary)]
            alignments.sort(key=attrgetter('start'))
        return alignments

    def get_positions_and_labels(self, start=None, end=None, ref_compr_rle=None, mock_compr=False, is_compressed=False, rle_dtype=False):
        """Create labels and positions array.

        :param start: starting position within reference
        :param end: ending position within reference
        :param ref_compr_rle: np.ndarray, dtype [('start', int), ('end', int)]
            containing uncompressed start position and length for each compressed position.
        :param mock_compr: bool, labels will be encoded as for a compressed sequence, even if ref_compr_rle is None.
        :param is_compressed: bool, whether the sequence has been run length encoded.
        :param rle_dtype: bool, whether to force dtype of labels to be (int base , int length).
        :returns: tuple(positions, encoded_label_array)

            - positions: numpy structured array with 'ref_major'
              (reference position index) and 'ref_minor'
              (trailing insertion index) fields.

            - label_array: 1D numpy array of labels
        """
        if start is None:
            start = 0
        if end is None:
            end = float('Inf')
        # ensure start and end fall within filtered region
        start = max(start, self.start)
        end = min(end, self.end)

        positions_labels = []
        if is_compressed:
            pairs = yield_compressed_pairs(self.aln, ref_compr_rle)
        elif mock_compr:
            pairs = get_pairs_with_hp_len(self.aln, ref_compr_rle)
        else:
            pairs = get_pairs(self.aln)

        ins_count = 0

        for pair in itertools.dropwhile(lambda x: (x.rpos is None)
                                        or (x.rpos < start), pairs):
            if ((pair.rpos == self.aln.reference_end ) or
                (pair.rpos is not None and pair.rpos >= end)):
                break
            if pair.rpos is None:
                ins_count += 1
            else:
                ins_count = 0
                current_pos = pair.rpos
            pos = (current_pos, ins_count)
            label = pair.qbase.upper() if pair.qbase else _gap_
            if label == 'N':
                logging.info('Found {} at pos {}'.format(label, pos))
            label = encoding[label]
            if is_compressed or mock_compr:
                label = (label, pair.qlen)
            elif rle_dtype:
                label = (label, 1)
            positions_labels.append((pos, label))

        if is_compressed or mock_compr or rle_dtype:
            dtype = [('base', int), ('run_length', int)]
        else:
            dtype = int
        label_array = np.empty(shape=(1, len(positions_labels)), dtype=dtype)
        positions = np.empty(len(positions_labels),
                             dtype=[('major', int), ('minor', int)])

        for i, (pos, label) in enumerate(positions_labels):
            positions[i] = pos
            label_array[0, i] = label

        return positions, label_array
