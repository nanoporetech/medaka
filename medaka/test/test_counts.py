import array
import numpy as np
import os
import tempfile
import unittest

import pysam

import libmedaka
import medaka.features
from medaka.common import Region

__reads_bam__ = os.path.join(os.path.dirname(__file__), 'data', 'test_reads.bam')
__two_type_bam__ = os.path.join(os.path.dirname(__file__), 'data', 'test_two_type.bam')
__gapped_bam__ = os.path.join(os.path.dirname(__file__), 'data', 'reads_gapped.bam')
__region__ = Region('Consensus_Consensus_Consensus_Consensus_utg000001l', start=50000, end=100000)
__region_start__ = Region('Consensus_Consensus_Consensus_Consensus_utg000001l', start=0, end=200)


def create_rle_bam(fname):
    """ Create a small bam file with RLE encoding coded in the qscores.

     Ref          A    C    A    T    *    G    A    T    G

    Basecall1:   2A   1C   4A   5T        1G   1A   2T   1G
    Basecall2:   3A   1C   4A    *        1G   1A   1T   2G
    Basecall3:   2A   1C   4A   5T   1A   1G   1A   2T   1G
    Basecall4:   2A   1C   4A   1C        1G   1A   2T   1G

    """

    ref = 'ACATGATG'
    basecall1 = {
        'seq': 'ACATGATG',  # exactly ref
        'quality': array.array('B', [2, 1, 4, 5, 1, 1, 2, 1]),
        'cigarstring': '8='}

    basecall2 = {
        'seq': 'ACAGATG',  # deletion of T in the middle
        'quality': array.array('B', [3, 1, 4, 1, 1, 1, 2]),
        'cigarstring': '3=1D4='}

    basecall3 = {
        'seq': 'ACATAGATG',  # insertion of A in the middle
        'quality':  array.array('B', [2, 1, 4, 5, 1, 1, 1, 2, 1]),
        'cigarstring': '4=1I4='}

    basecall4 = {
        'seq': 'ACACGATG',  # substitution T->C
        'quality': array.array('B', [2, 1, 4, 1, 1, 1, 2, 1]),
        'cigarstring': '3=1X4='}

    header = {'HD': {'VN': '1.0'},
              'SQ': [{'LN': 8, 'SN': 'ref'}]}

    tmp_file = '{}.tmp'.format(fname)
    with pysam.AlignmentFile(
            tmp_file, 'wb', reference_names=['ref', ],
            reference_lengths=[8, ], header=header) as output:
        for index, basecall in enumerate(
                (basecall1, basecall2, basecall3, basecall4)):
            a = pysam.AlignedSegment()
            a.query_name = "basecall_{}".format(index)
            a.reference_id = 0
            a.reference_start = 0
            a.query_sequence = basecall['seq']
            a.cigarstring = basecall['cigarstring']

            a.flag = 0
            a.mapping_quality = 50
            a.query_qualities = basecall['quality']
            output.write(a)

    pysam.sort("-o", fname, tmp_file)
    os.remove(tmp_file)
    pysam.index(fname)


class CountsTest(unittest.TestCase):

    def test_001_basic_counting(self):
        kwargs = {'normalise': None}
        encoder = medaka.features.CountsFeatureEncoder(**kwargs)
        sample = encoder.bam_to_sample(__reads_bam__, __region__)
        sample = sample[0]
        assert tuple(sample.positions.shape) == (81730,)
        assert tuple(sample.positions[0]) == (50000, 0)
        assert tuple(sample.positions[-1]) == (99999, 1)
        assert sample.features.shape == (81730, 10)
        # test counts
        np.testing.assert_array_equal(sample.features[0], np.array([ 0, 21, 0, 1, 0, 14, 0, 0, 0, 0]))
        # test mean depth
        np.testing.assert_almost_equal(np.mean(np.sum(sample.features, axis=1)), 19.83996)


class CountsSplittingTest(unittest.TestCase):

    def test_000_split_gap(self):
        # The gapped bam has:
        # @SQ    SN:ref    LN:30
        # seq1    0    ref    1    7    10M
        # seq2    0    ref    15    13    16M
        # so an alignment from [0:10] and one from [14:30] without insertions
        chunk_lengths = [10, 16]

        region = Region.from_string('ref:0-30')
        results = medaka.features.pileup_counts(region, __gapped_bam__)
        self.assertEqual(len(results), 2, 'Number of chunks from gapped alignment')
        for exp_len, chunk in zip(chunk_lengths, results):
            for i in (0, 1):
                # check both pileup and positions
                self.assertEqual(exp_len, len(chunk[i]))


class CountsQscoreStratification(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        (counts_strat1, positions_strat1) = medaka.features.pileup_counts(
            __region__, __reads_bam__, num_qstrat=1)[0]
        (counts_strat2, positions_strat2) = medaka.features.pileup_counts(
            __region__, __reads_bam__, num_qstrat=2)[0]

        cls.counts_strat1 = counts_strat1
        cls.positions_strat1 = positions_strat1
        cls.counts_strat2 = counts_strat2
        cls.positions_strat2 = positions_strat2

    def test_000_positions_indifferent_to_qstrat(self):
        """The positions in a bam is indifferent to the stratification by qscores."""

        self.assertSequenceEqual(self.positions_strat1.tolist(), self.positions_strat2.tolist())

    def test_001_num_columns(self):
        """The number of columns does not depend on the stratification
        of qscores."""

        rows_strat1 = self.counts_strat1.shape[0]
        rows_strat2 = self.counts_strat2.shape[0]

        self.assertEqual(rows_strat1, rows_strat2)

    def test_002_length_features(self):
        """The length in the feature doubles for a stratification with 2 layers."""

        columns_strat1 = self.counts_strat1.shape[1]
        columns_strat2 = self.counts_strat2.shape[1]

        self.assertEqual(columns_strat1 * 2, columns_strat2)

    def test_003_marginalise_qstrat(self):
        """When you add all sections of the stratified counts, you need
            to get the non-stratified one."""
        strat2_q0 = self.counts_strat2[:,0:10]
        strat2_q1 = self.counts_strat2[:,10:]

        added = strat2_q0 + strat2_q1
        self.assertTrue(np.all(added == self.counts_strat1))


class PileupCountsNormIndices(unittest.TestCase):

    def test_000_single_dtype_no_qstrat(self):
        """ If there is a single dtype, no qscore stratification the codes are: acgtACGTdD. """
        expected = {('', True): [0, 1, 2, 3, 8],
                    ('', False): [4, 5, 6, 7, 9]}
        got = medaka.features.pileup_counts_norm_indices(['',], num_qstrat=1)

        self.assertDictEqual(expected, got)

    def test_001_two_dtypes_no_qstrat(self):
        """With 2 qtypes, the codes are: acgtACGTdDacgtACGTdD."""
        expected = {('1', True): [0, 1, 2, 3, 8],
                    ('2', True): [10, 11, 12, 13, 18],
                    ('1', False): [4, 5, 6, 7, 9],
                    ('2', False): [14, 15, 16, 17, 19]}
        got = medaka.features.pileup_counts_norm_indices(['1', '2'], num_qstrat=1)

        self.assertDictEqual(expected, got)

    def test_002_one_dtype_two_qstrat(self):
        """Single qtype, but two layers of qscore stratification:  acgtACGTdDacgtACGTdD
        """
        expected = {('', True): [0, 1, 2, 3, 8, 10, 11, 12, 13, 18],
                    ('', False): [4, 5, 6, 7, 9, 14, 15, 16, 17, 19]}
        got = medaka.features.pileup_counts_norm_indices([''], num_qstrat=2)
        self.assertDictEqual(expected, got)


class HardRLEFeatureEncoder(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        kwargs = {'normalise': None}
        encoder = medaka.features.HardRLEFeatureEncoder(**kwargs)

        # Create a bam file where we know the alignments
        RLE_bam = tempfile.NamedTemporaryFile(suffix='.bam').name
        create_rle_bam(RLE_bam)
        sample = encoder.bam_to_sample(RLE_bam, Region('ref', 0, 8))
        cls.sample = sample[0]
        cls.num_qstrat = encoder.num_qstrat

    def test_001_check_number_positions(self):
        """Check number of positions returned."""
        self.assertEqual(len(self.sample.positions), 9)

    def test_002_check_start_end_region(self):
        """Check the positions start and end where they should."""
        self.assertSequenceEqual(tuple(self.sample.positions[0]), (0, 0))
        self.assertSequenceEqual(tuple(self.sample.positions[-1]), (7, 0))

    def test_003_check_counts_shape(self):
        """Check shape of counts is correct:

            feature size will be 10 numbers per qscore stratification layer
        """
        featlen = libmedaka.lib.featlen
        self.assertEqual(self.sample.features.shape, (9, featlen * self.num_qstrat))

    def test_004_check_counts(self):
        """Check the counts themselves. See above create_rle_bam to understand
         the values in the counts. For example, 3 basecalls believe the first
         column of the alignment is `2A`, while one has `3A`. Thus, index
         [0, 14] should contain 3 and [0, 24] should show a 1.
        """

        # dictionary with value: index encoding
        values_for_positions = {
            1: ([0, 24], [3, 9], [3, 5], [4, 4], [7, 7], [8, 16]),
            2: ([3, 47],),
            3: ([0, 14], [7, 17], [8, 6]),
            4: ([1, 5], [2, 34], [5, 6], [6, 4])}
        expected = np.zeros_like(self.sample.features)
        for value, positions in values_for_positions.items():
            for pos in positions:
                expected[pos[0], pos[1]] = value

        self.assertSequenceEqual(self.sample.features.tolist(), expected.tolist())
