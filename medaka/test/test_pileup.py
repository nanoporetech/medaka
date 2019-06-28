import os
import unittest
from collections import namedtuple
from medaka.labels import TruthAlignment

__truth_bam__ = os.path.join(os.path.dirname(__file__), 'data', 'truth_to_ref.bam')
__ref_fasta__ = os.path.join(os.path.dirname(__file__), 'data', 'draft_ref.fasta')
__ref_name__ = 'Consensus_Consensus_Consensus_Consensus_utg000001l'

MockAlignment = namedtuple('MockAlignment', ['reference_start', 'reference_end', 'reference_length'])

class TruthAlignmentTest(unittest.TestCase):

    def test_case1(self):
        # case 1: longer < 2 x len shorter and >= 50% of shorter overlaps longer both should be removed
        starts_ends = [(2000, 2999), (2500, 3000)]
        expected = []

        alignments = [TruthAlignment(MockAlignment(start, end, end-start)) for start, end in starts_ends]
        filtered = [(f.start, f.end) for f in TruthAlignment.filter_alignments(alignments)]
        if filtered != expected:
            raise AssertionError('got {}, expected {}'.format(filtered, expected))


    def test_case2(self):
        # case 2: longer < 2 x len shorter and < 50% shorter overlaps longer, trim both
        starts_ends = [(5000, 5999), (5900, 6401)]

        alignments = [TruthAlignment(MockAlignment(start, end, end-start)) for start, end in starts_ends]

        filtered = [(f.start, f.end) for f in TruthAlignment.filter_alignments(alignments, min_length=401)]
        expected = [(5000, 5900), (5999, 6401)]
        if filtered != expected:
            raise AssertionError('got {}, expected {}'.format(filtered, expected))

        filtered = [(f.start, f.end) for f in TruthAlignment.filter_alignments(alignments, min_length=899)]
        expected = [(5000, 5900)]
        if filtered != expected:
            raise AssertionError('got {}, expected {}'.format(filtered, expected))

        filtered = [(f.start, f.end) for f in TruthAlignment.filter_alignments(alignments, min_length=1000)]
        expected = []
        if filtered != expected:
            raise AssertionError('got {}, expected {}'.format(filtered, expected))


    def test_case3(self):
        # case 3: longer >= 2 x len shorter and < 50% shorter overlaps longer, remove shorter
        starts_ends = [(7000, 8000), (7501, 8000)]
        expected = [(7000, 8000)]

        alignments = [TruthAlignment(MockAlignment(start, end, end-start)) for start, end in starts_ends]
        filtered = [(f.start, f.end) for f in TruthAlignment.filter_alignments(alignments)]
        if filtered != expected:
            raise AssertionError('got {}, expected {}'.format(filtered, expected))

        # case 3 (contained): shorter contained within longer, shorter should be removed
        starts_ends = [(0, 1000), (100, 200)]
        expected = [(0, 1000)]

        alignments = [TruthAlignment(MockAlignment(start, end, end-start)) for start, end in starts_ends]
        filtered = [(f.start, f.end) for f in TruthAlignment.filter_alignments(alignments)]
        if filtered != expected:
            raise AssertionError('got {}, expected {}'.format(filtered, expected))


    def test_case4(self):
        # case 4: longer >= 2 x len shorter and < 50% shorter overlaps longer, trim shorter
        starts_ends = [(3000, 4000), (3800, 4299)]
        expected = [(3000, 4000), (4000, 4299)]

        alignments = [TruthAlignment(MockAlignment(start, end, end-start)) for start, end in starts_ends]
        filtered = [(f.start, f.end) for f in TruthAlignment.filter_alignments(alignments, min_length=298)]
        if filtered != expected:
            raise AssertionError('got {}, expected {}'.format(filtered, expected))


    def test_many(self):
        starts_ends = [(0, 1000), (100, 200),
                       (3000, 4000), (3800, 4299),
                       (7000, 8000), (7501, 8000),
                       (2000, 2999), (2500, 3000),
                       (5000, 5999), (5900, 6401),]

        alignments = [TruthAlignment(MockAlignment(start, end, end-start)) for start, end in starts_ends]
        filtered = [(f.start, f.end) for f in TruthAlignment.filter_alignments(alignments, min_length=1)]
        expected = [(0, 1000), (3000, 4000), (4000, 4299), (5000, 5900), (5999, 6401), (7000, 8000)]
        if filtered != expected:
            raise AssertionError('got {}, expected {}'.format(filtered, expected))

        starts_ends_m = starts_ends + [(0, 5000)]
        alignments = [TruthAlignment(MockAlignment(start, end, end-start)) for start, end in starts_ends_m]
        filtered = [(f.start, f.end) for f in TruthAlignment.filter_alignments(alignments, min_length=1)]
        expected = [(0, 5000), (5000, 5900), (5999, 6401), (7000, 8000)]
        if filtered != expected:
            raise AssertionError('got {}, expected {}'.format(filtered, expected))

        starts_ends_m = starts_ends + [(0, 5200)]
        alignments = [TruthAlignment(MockAlignment(start, end, end-start)) for start, end in starts_ends_m]
        filtered = [(f.start, f.end) for f in TruthAlignment.filter_alignments(alignments, min_length=1)]
        expected = [(0, 5200), (5200, 5900), (5999, 6401), (7000, 8000)]
        if filtered != expected:
            raise AssertionError('got {}, expected {}'.format(filtered, expected))

        starts_ends_m = starts_ends + [(0, 10000)]
        alignments = [TruthAlignment(MockAlignment(start, end, end-start)) for start, end in starts_ends_m]
        filtered = [(f.start, f.end) for f in TruthAlignment.filter_alignments(alignments, min_length=1)]
        expected = [(0, 10000)]
        if filtered != expected:
            raise AssertionError('got {}, expected {}'.format(filtered, expected))


    def test_labels_trimmed_back(self):
        # we should have two alignments which partially overlap
        # (318288, 417741)
        # (417732, 422799)
        # in this case, the first is >2 x longer than the second, so we trim back the second
        # check resulting positions and labels are non-overlapping
        alignments = TruthAlignment.bam_to_alignments(__truth_bam__, __ref_name__, start=318288, end=422799)
        assert alignments[0].start == 318288
        assert alignments[0].end == 417741
        assert alignments[1].start == 417732
        assert alignments[1].end == 422799
        filtered_alignments = TruthAlignment.filter_alignments(alignments)
        assert filtered_alignments[0].start == 318288
        assert filtered_alignments[0].end == 417741
        assert filtered_alignments[1].start == 417741
        assert filtered_alignments[1].end == 422799

        p1_positions, p1_labels = filtered_alignments[0].get_positions_and_labels()
        p2_positions, p2_labels = filtered_alignments[1].get_positions_and_labels()
        assert p1_positions[0]['major'] == filtered_alignments[0].start
        assert p1_positions[-1]['major'] == filtered_alignments[0].end - 1
        assert p2_positions[0]['major'] == filtered_alignments[1].start
        assert p2_positions[-1]['major'] == filtered_alignments[1].end - 1


if __name__ == '__main__':
    unittest.main()
