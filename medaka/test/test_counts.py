import copy
import os
import pickle
import tempfile
import unittest

import numpy as np

from .mock_data import simple_data, create_simple_bam
import libmedaka
import medaka.features
from medaka.common import Region, Sample
import medaka.labels

__reads_bam__ = os.path.join(os.path.dirname(__file__), 'data', 'test_reads.bam')
__reads_truth__ = os.path.join(os.path.dirname(__file__), 'data', 'truth_to_ref.bam')
__gapped_bam__ = os.path.join(os.path.dirname(__file__), 'data', 'reads_gapped.bam')
__region__ = Region('utg000001l', start=50000, end=100000)
__region_start__ = Region('utg000001l', start=0, end=200)


class CountsTest(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.expected_width = 86294

    def test_001_basic_counting(self):
        kwargs = {'normalise': None}
        encoder = medaka.features.CountsFeatureEncoder(**kwargs)
        sample = encoder.bam_to_sample(__reads_bam__, __region__)
        self.assertEqual(len(sample), 1)
        sample = sample[0]
        assert tuple(sample.positions.shape) == (self.expected_width,)
        assert tuple(sample.positions[0]) == (50000, 0)
        assert tuple(sample.positions[-1]) == (99999, 1)
        assert sample.features.shape == (self.expected_width, 10)
        # test counts
        np.testing.assert_array_equal(
            sample.features[0],
            np.array([0., 22., 0., 0., 0., 15., 0., 0., 0., 0.]))
        # test mean depth
        np.testing.assert_almost_equal(
            np.mean(np.sum(sample.features, axis=1)),
            18.696468, decimal=6)

    def test_002_raises_on_invalid_norm(self):
        with self.assertRaises(ValueError):
            kwargs = {'normalise': 'nonsense'}
            encoder = medaka.features.CountsFeatureEncoder(**kwargs)

    def test_010_pickleble(self):
        kwargs = {'normalise': None}
        encoder = medaka.features.CountsFeatureEncoder(**kwargs)
        branston = pickle.loads(pickle.dumps(encoder))
        self.assertTrue(hasattr(branston, 'logger'))
        self.assertEqual(branston.normalise, None)

    def test_020_feature_length(self):
        # hardcode value here, if this genuinely needs to change at least
        # the test will make develop think twice about consequences
        kwargs = {'dtypes': ['r9','r10']}
        encoder = medaka.features.CountsFeatureEncoder(**kwargs)
        self.assertEqual(encoder.feature_vector_length, 20)

    def test_030_bams_to_training_samples_simple(self):
        reads_bam = tempfile.NamedTemporaryFile(suffix='.bam').name
        truth_bam = tempfile.NamedTemporaryFile(suffix='.bam').name

        # we had a bug caused by missing qualities and bad indexing...
        data = copy.deepcopy(simple_data['calls'])
        data[0]['quality'] = None

        create_simple_bam(reads_bam, data)
        create_simple_bam(
            truth_bam, [simple_data['truth']])
        encoder = medaka.features.CountsFeatureEncoder(normalise='total')
        label_scheme = medaka.labels.HaploidLabelScheme()
        region = Region('ref', 0, 100)
        result = encoder.bams_to_training_samples(
            truth_bam, reads_bam, region, label_scheme, min_length=0)[0]

        expected = Sample(
            ref_name='ref',
            features=np.array([
                [0.5 , 0.  , 0.  , 0.  , 0.5 , 0.  , 0.  , 0.  , 0.  , 0.  ],
                [0.  , 0.5 , 0.  , 0.  , 0.  , 0.5 , 0.  , 0.  , 0.  , 0.  ],
                [0.5 , 0.  , 0.  , 0.  , 0.5 , 0.  , 0.  , 0.  , 0.  , 0.  ],
                [0.  , 0.25, 0.  , 0.25, 0.  , 0.  , 0.  , 0.25, 0.  , 0.25],
                [0.25, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ],
                [0.  , 0.  , 0.5 , 0.  , 0.  , 0.  , 0.5 , 0.  , 0.  , 0.  ],
                [0.5 , 0.  , 0.  , 0.  , 0.5 , 0.  , 0.  , 0.  , 0.  , 0.  ],
                [0.  , 0.  , 0.  , 0.5 , 0.  , 0.  , 0.  , 0.5 , 0.  , 0.  ],
                [0.  , 0.  , 0.5 , 0.  , 0.  , 0.  , 0.5 , 0.  , 0.  , 0.  ]],
                dtype='float32'),
            # the two insertions with respect to the draft are dropped
            labels=np.array([1, 2, 1, 4, 1, 3, 1, 4, 3]),  # A C A T A G A T C
            ref_seq=None,
            positions=np.array([
                (0, 0), (1, 0), (2, 0), (3, 0), (3, 1), (4, 0), (5, 0), (6, 0), (7, 0)],
                dtype=[('major', '<i8'), ('minor', '<i8')]),
            label_probs=None
        )

        np.testing.assert_equal(result.labels, expected.labels)
        np.testing.assert_equal(result.positions, expected.positions)
        np.testing.assert_equal(result.features, expected.features)


    def test_031_bams_to_training_samples_regression(self):
        encoder = medaka.features.CountsFeatureEncoder(normalise='total')
        label_scheme = medaka.labels.HaploidLabelScheme()
        region = Region(
            'utg000001l',
            149744, 318288)
        result = encoder.bams_to_training_samples(
            __reads_truth__, __reads_bam__, region, label_scheme)[0]

        expected_feature_shape = (177981, 10)
        got_feature_shape = result.features.shape
        self.assertEqual(expected_feature_shape, got_feature_shape)

        expected_label_shape = (177981,)
        got_label_shape = result.labels.shape
        self.assertEqual(expected_label_shape, got_label_shape)



class TrimReadsTest(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        # Create a bam file where we know the alignments
        self.bam = tempfile.NamedTemporaryFile(suffix='.bam').name
        create_simple_bam(self.bam, simple_data['calls'])
        self.reads = [x['seq'] for x in simple_data['calls']]

    def get_reads(self, region):
        reads = list(medaka.features.get_trimmed_reads(region, self.bam))
        reads = [x[1] for x in reads[0][1]]
        reads = reads[1:]  # contains reference first
        return reads

    def test_001_full_region(self):
        region = Region('ref', start=0, end=100000)
        reads = self.get_reads(region)
        self.assertEqual(reads, self.reads)

    def test_002_trim_start(self):
        region = Region('ref', start=0, end=2)
        reads = self.get_reads(region)
        orig = [x[0:2] for x in self.reads]
        self.assertEqual(reads, orig)

    def test_003_trim_end(self):
        region = Region('ref', start=6, end=8)
        reads = self.get_reads(region)
        orig = [x[-2:] for x in self.reads]
        self.assertEqual(reads, orig)

    def test_004_trim_mid(self):
        region = Region('ref', start=1, end=7)
        reads = self.get_reads(region)
        orig = [x[1:-1] for x in self.reads]
        self.assertEqual(reads, orig)

    def test_005_chunked(self):
        region = Region('ref', start=0, end=9)
        n_chunks = len(list(medaka.features.get_trimmed_reads(region, self.bam, region_split=3, chunk_overlap=0)))
        self.assertEqual(n_chunks, 3)

    def test_005_unchunked(self):
        region = Region('ref', start=0, end=9)
        n_chunks = len(list(medaka.features.get_trimmed_reads(region, self.bam, region_split=1000, chunk_overlap=0)))
        self.assertEqual(n_chunks, 1)


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


class SampleGenerator(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.expected_width = 86294

    def test_000_basic_sample_gen(self):
        encoder = medaka.features.CountsFeatureEncoder()
        sample_gen = medaka.features.SampleGenerator(
            __reads_bam__, __region__, encoder, enable_chunking=False)
        # fill features should fill in _source
        sample_gen._fill_features()
        self.assertEqual(len(sample_gen._source), 1)
        self.assertEqual(sample_gen._source[0].positions.shape, (self.expected_width,))

        # reset source to ensure calculated on fly
        sample_gen._source = None
        samples = sample_gen.samples
        self.assertEqual(len(samples), 1, 'Have more than 1 chunk despite chunking disabled.')

    def test_010_chunking(self):
        chunk_len = 1000
        encoder = medaka.features.CountsFeatureEncoder()
        sample_gen = medaka.features.SampleGenerator(
            __reads_bam__, __region__, encoder,
            chunk_len=chunk_len, chunk_overlap=0)
        samples = sample_gen.samples
        self.assertEqual(len(samples), self.expected_width // chunk_len + 1)
        self.assertEqual(len(sample_gen._quarantined), 0, 'Samples were quarantined incorrectly.')

        sample_gen.chunk_len = 1000000
        samples = sample_gen.samples
        self.assertEqual(len(sample_gen._quarantined), 1)

    def test_020_mandatory_label_scheme(self):
        """A SampleGenerator initialised with truth_bam needs a label_scheme"""
        encoder = medaka.features.CountsFeatureEncoder()
        with self.assertRaises(ValueError) as context:
            sample_gen = medaka.features.SampleGenerator(
                __reads_bam__, __region__, encoder, truth_bam=__reads_truth__)


class PileupCounts(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        bam_fname = tempfile.NamedTemporaryFile(suffix='.bam').name
        create_simple_bam(bam_fname, simple_data['calls'])
        cls.region = Region(
            'ref',
            start=0, end=8)
        cls.bam = bam_fname
        cls.expected_counts = np.array(
            [[2, 0, 0, 0, 2, 0, 0, 0, 0, 0],
             [0, 2, 0, 0, 0, 2, 0, 0, 0, 0],
             [2, 0, 0, 0, 2, 0, 0, 0, 0, 0],
             [0, 1, 0, 1, 0, 0, 0, 1, 0, 1],
             [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 2, 0, 0, 0, 2, 0, 0, 0],
             [2, 0, 0, 0, 2, 0, 0, 0, 0, 0],
             [0, 0, 0, 2, 0, 0, 0, 2, 0, 0],
             [0, 0, 2, 0, 0, 0, 2, 0, 0, 0]], dtype=np.uint64)
        cls.expected_positions = np.array(
            [(0, 0), (1, 0), (2, 0), (3, 0), (3, 1),
             (4, 0), (5, 0), (6, 0), (7, 0)],
             dtype=[('major', '<i8'), ('minor', '<i8')])

    def test_000_counts_no_dtypes(self):
        counts, positions = medaka.features.pileup_counts(
            self.region, self.bam)[0]

        self.assertTrue(np.array_equal(counts, self.expected_counts))
        self.assertTrue(np.array_equal(positions, self.expected_positions))

    def test_010_tags(self):
        """Only two reads have tag AA=1."""
        counts, positions = medaka.features.pileup_counts(
            self.region, self.bam, tag_name='AA', tag_value=1,
            keep_missing=False)[0]
        expected_number_reads = {2}
        self.assertEqual(expected_number_reads, set(counts.sum(axis=1)))

    def test_020_tags(self):
        """If keep_missing is True, 3 reads are expected"""
        counts, positions = medaka.features.pileup_counts(
            self.region, self.bam, tag_name='AA', tag_value=1,
            keep_missing=True)[0]
        expected_number_reads = {3}
        self.assertEqual(expected_number_reads, set(counts.sum(axis=1)))

    def test_030_tag_error(self):
        """Tags should be two letters"""
        with self.assertRaises(ValueError) as context:
            counts, positions = medaka.features.pileup_counts(
                self.region, self.bam, tag_name='AAC', tag_value=1,
                keep_missing=True)[0]

        with self.assertRaises(ValueError) as context:
            counts, positions = medaka.features.pileup_counts(
                self.region, self.bam, tag_name='A', tag_value=1,
                keep_missing=True)[0]

    def test_040_dtypes(self):
        counts, positions = medaka.features.pileup_counts(
            self.region, self.bam, dtype_prefixes=['r9', 'r10'])[0]

        expected_shape = (9, 20)
        self.assertEqual(expected_shape, counts.shape)

    def test_050_read_groups(self):
        bam = tempfile.NamedTemporaryFile(suffix='.bam').name
        reads = list()
        region = Region('ref', start=0, end=8)
        for rg in ('first', 'second'):
            rg_reads = copy.deepcopy(simple_data['calls'])
            for read in rg_reads:
                read['tags']['RG'] = rg
                read['query_name'] += '_{}'.format(rg)
            reads.extend(rg_reads)
        create_simple_bam(bam, reads)

        # use everything
        counts, positions = medaka.features.pileup_counts(region, bam)[0]
        self.assertTrue(np.array_equal(counts, 2*self.expected_counts))
        self.assertTrue(np.array_equal(positions, self.expected_positions))

        # use one or other
        for rg in ('first', 'second'):
            counts, positions = medaka.features.pileup_counts(region, bam, read_group=rg)[0]
            self.assertTrue(np.array_equal(counts, self.expected_counts))
            self.assertTrue(np.array_equal(positions, self.expected_positions))

        # use a missing one
        result = medaka.features.pileup_counts(region, bam, read_group='nonsense')
        self.assertTrue(len(result) == 0)


    def test_060_counts_chunked(self):
        """Test that we get the same counts and positions when region_split < region.size

        i.e. the stitching back together of chunks processed in parallel works.
        """
        chunks = medaka.features.pileup_counts(self.region, self.bam, region_split=3)
        self.assertEqual(len(chunks), 1)
        counts, positions = chunks[0]
        self.assertTrue(np.array_equal(counts, self.expected_counts))
        self.assertTrue(np.array_equal(positions, self.expected_positions))


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

        np.testing.assert_equal(self.positions_strat1, self.positions_strat2)

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

    def test_003_check_values(self):
        """Sum of stratified counts should be equal to non-straified counts."""
        strat2_q0 = self.counts_strat2[:,0:10]
        strat2_q1 = self.counts_strat2[:,10:]

        added = strat2_q0 + strat2_q1
        self.assertTrue(np.all(added == self.counts_strat1))


class WeibullSummation(CountsQscoreStratification):

    @classmethod
    def setUpClass(cls):
        temp_file=tempfile.NamedTemporaryFile(suffix='.bam')
        bam_fname = temp_file.name
        region = Region('ref', 0, 100)
        create_simple_bam(bam_fname, simple_data['calls'])
        (counts_strat1, positions_strat1) = medaka.features.pileup_counts(
            region, bam_fname, num_qstrat=1, weibull_summation=True)[0]
        (counts_strat2, positions_strat2) = medaka.features.pileup_counts(
            region, bam_fname, num_qstrat=2, weibull_summation=True)[0]

        cls.counts_strat1 = counts_strat1
        cls.positions_strat1 = positions_strat1
        cls.counts_strat2 = counts_strat2
        cls.positions_strat2 = positions_strat2

    def test_003_check_values(self):
        """Partial array should be the same regardless of stratifications."""
        np.testing.assert_equal(self.counts_strat1, self.counts_strat2[:,:10])
        for i in range(2):
            # for the given parameters we should have non-zero values
            # in both stratifications
            non_zeros = np.sum(self.counts_strat2[:,10*i:10*(i+1)] > 0)
            self.assertTrue(non_zeros > 0)

    def test_004_check_specifics(self):
        temp_file=tempfile.NamedTemporaryFile(suffix='.bam', delete=False)
        bam_fname = temp_file.name
        region = Region('ref', 0, 8)
        create_simple_bam(
            bam_fname, simple_data['calls'][0:1])
        counts, _ = medaka.features.pileup_counts(
            region, bam_fname, num_qstrat=6, weibull_summation=True)[0]

        # should have only one non-zero per row
        non_zeroes = np.nonzero(counts)
        np.testing.assert_equal(non_zeroes[0], np.arange(0,8))
        np.testing.assert_equal(
            non_zeroes[1],
            # acgtACGTdD...
            #         2A  1C  4A  5T  1G  1A  2T  1G
            np.array([14,  5, 34, 47,  6,  4, 17,  6]))


class PileupCountsNormIndices(unittest.TestCase):

    def test_000_single_dtype_no_qstrat(self):
        """If there is a single dtype, no qscore stratification the codes are: acgtACGTdD. """
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
        """Single qtype, but two layers of qscore stratification:
        acgtACGTdDacgtACGTdD
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
        create_simple_bam(RLE_bam, simple_data['calls'])
        sample = encoder.bam_to_sample(RLE_bam, Region('ref', 0, 11))
        cls.bam_fname = RLE_bam
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
        """Check the counts themselves. See above `create_simple_bam` to understand
        the values in the counts. For example, 1 basecalls believe the first
        column of the alignment is `2A`, while two have `2a`. One has `0A` which
        will become a `1A` when corrected.

        Thus:
        [0, 4] (1A) contains 1
        [0,10] (2a) contains 2
        [0,14] (2A) contains 1
        """
        values_for_positions = {
            1.0: [
                (0, 4), (0, 14), (3, 1), (3, 9), (3, 43), (3, 47),
                (4, 0), (7, 7), (7, 17), (8, 6), (8, 16)],
            2.0: [
                (0, 10), (1, 1), (1, 5), (2, 30), (2, 34), (5, 2),
                (5, 6), (6, 0), (6, 4), (7, 13), (8, 2)]}
        expected = np.zeros_like(self.sample.features)
        for value, positions in values_for_positions.items():
            for pos in positions:
                expected[pos[0], pos[1]] = value
        np.testing.assert_equal(self.sample.features, expected)

    def test_005_check_fwd_rev(self):
        """Split normalisation between fwd and rev reads. """
        encoder = medaka.features.HardRLEFeatureEncoder(normalise='fwd_rev')
        region = Region('ref', 0, 11)
        sample = encoder.bam_to_sample(self.bam_fname, region)[0]
        values_for_positions = {
            0.5: [
                (0, 4), (0, 14), (3, 1), (3, 9), (3, 43), (3, 47),
                (4, 0), (7, 7), (7, 17), (8, 6), (8, 16)],
            1.0: [
                (0, 10), (1, 1), (1, 5), (2, 30), (2, 34), (5, 2),
                (5, 6), (6, 0), (6, 4), (7, 13), (8, 2)]}
        expected = np.zeros_like(sample.features)
        for value, positions in values_for_positions.items():
            for pos in positions:
                expected[pos[0], pos[1]] = value
        np.testing.assert_equal(sample.features, expected)

    def test_020_feature_length(self):
        # hardcode value here, if this genuinely needs to change at least
        # the test will make develop think twice about consequences
        kwargs = {'dtypes': ['r9','r10']}

        encoder = medaka.features.HardRLEFeatureEncoder(
            num_qstrat=10, **kwargs)
        self.assertEqual(encoder.feature_vector_length, 200)


class SymHardRLEFeatureEncoder(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        kwargs = {'normalise': None}
        encoder = medaka.features.SymHardRLEFeatureEncoder(**kwargs)

        # Create a bam file where we know the alignments
        RLE_bam = tempfile.NamedTemporaryFile(suffix='.bam').name
        create_simple_bam(RLE_bam, simple_data['calls'])
        sample = encoder.bam_to_sample(RLE_bam, Region('ref', 0, 11))
        cls.bam_fname = RLE_bam
        cls.sample = sample[0]
        cls.num_qstrat = encoder.num_qstrat

    def test_001_check_counts(self):
        """Check the counts themselves. See above `create_simple_bam` to understand
        the values in the counts. For example, 1 basecalls believe the first
        column of the alignment is `2A`, while two have `2a`. One has `0A` which
        will become a `1A` when corrected. Three calls have no insertion, only
        the third basecall has it, so SymHardRLEFeatureEncoder will count 3 indels
        at that position, unlike other FeatureEncoders.

        Thus:
        [4, 8] (1d) will have one counts
        [4, 9] (1D) will have two counts
        """
        values_for_positions = {
            1.0: [
                (0, 4), (0, 14), (3, 1), (3, 9), (3, 43), (3, 47),
                (4, 0), (4, 8), (7, 7), (7, 17), (8, 6), (8, 16)],
            2.0: [
                (0, 10), (1, 1), (1, 5), (2, 30), (2, 34), (4,9),
                (5, 2), (5, 6), (6, 0), (6, 4), (7, 13), (8, 2)]}
        expected = np.zeros_like(self.sample.features)
        for value, positions in values_for_positions.items():
            for pos in positions:
                expected[pos[0], pos[1]] = value
        np.testing.assert_equal(self.sample.features, expected)
