import numpy as np
import os
import unittest
from medaka.features import CountsFeatureEncoder, pileup_counts
from medaka.common import Region

__reads_bam__ = os.path.join(os.path.dirname(__file__), 'data', 'test_reads.bam')
__two_type_bam__ = os.path.join(os.path.dirname(__file__), 'data', 'test_two_type.bam')
__gapped_bam__ = os.path.join(os.path.dirname(__file__), 'data', 'reads_gapped.bam')
__region__ = Region('Consensus_Consensus_Consensus_Consensus_utg000001l', start=50000, end=100000)
__region_start__ = Region('Consensus_Consensus_Consensus_Consensus_utg000001l', start=0, end=200)


class CountsTest(unittest.TestCase):

    def test_001_basic_counting(self):

        # py-style
        kwargs = {'normalise': None}
        encoder = CountsFeatureEncoder(**kwargs)
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
        results = pileup_counts(region, __gapped_bam__)
        self.assertEqual(len(results), 2, 'Number of chunks from gapped alignment')
        for exp_len, chunk in zip(chunk_lengths, results):
            for i in (0, 1):
                # check both pileup and positions
                self.assertEqual(exp_len, len(chunk[i]))

if __name__ == '__main__':
    unittest.main()
