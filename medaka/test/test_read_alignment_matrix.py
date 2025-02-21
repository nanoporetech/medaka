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


class AlignmentMatrixTest(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.expected_width = 86294
        self.expected_depth = 45
        self.base_featlen = 4 # bases, quality, strand, mapping quality

    def test_001_basic_counting(self):
        # check that feature shapes for encoder with out dwells is correct
        encoder = medaka.features.ReadAlignmentFeatureEncoder(include_dwells=False)
        sample = encoder.bam_to_sample(__reads_bam__, __region__)
        self.assertEqual(len(sample), 1)
        sample = sample[0]
        assert tuple(sample.positions.shape) == (self.expected_width, )
        assert tuple(sample.positions[0]) == (50000, 0)
        assert tuple(sample.positions[-1]) == (99999, 1)
        assert sample.features.shape == (self.expected_width, self.expected_depth, self.base_featlen)

    def test_002_totals_match_counts(self):
        encoder = medaka.features.ReadAlignmentFeatureEncoder(include_dwells=False)
        sample = encoder.bam_to_sample(__reads_bam__, __region__)[0]

        counts_encoder = medaka.features.CountsFeatureEncoder(normalise=None, sym_indels=True)
        counts_sample = counts_encoder.bam_to_sample(__reads_bam__, __region__)[0]

        np.testing.assert_array_equal(sample.positions, counts_sample.positions)

        # total reads per position (assuming no-read is encoded as 0)
        NO_READ_VAL = 0
        total_reads = (sample.features[:, :, 0] != NO_READ_VAL).sum(-1)
        counts_total_reads = counts_sample.features.sum(-1)
        np.testing.assert_array_equal(total_reads, counts_total_reads)

        # base counts per position
        base_counts = np.array(
            [(sample.features[:, :, 0] == (i + 1)).sum(-1) for i in range(5)])
        counts_base_counts = np.hstack([
            counts_sample.features[:, :4] + counts_sample.features[:, 4:8],
            counts_sample.features[:, 8][:, None] + counts_sample.features[:, 9][:, None],
        ]).transpose()
        np.testing.assert_array_equal(base_counts, counts_base_counts)

    def test_003_basic_counting_row_per_read(self):
        encoder = medaka.features.ReadAlignmentFeatureEncoder(row_per_read=True, max_reads=1000, include_dwells=False)
        sample = encoder.bam_to_sample(__reads_bam__, __region__)
        self.assertEqual(len(sample), 1)
        sample = sample[0]
        assert tuple(sample.positions.shape) == (self.expected_width, )
        assert tuple(sample.positions[0]) == (50000, 0)
        assert tuple(sample.positions[-1]) == (99999, 1)
        assert sample.features.shape == (self.expected_width, 291, self.base_featlen)

    def test_010_pickleble(self):
        encoder = medaka.features.ReadAlignmentFeatureEncoder()
        branston = pickle.loads(pickle.dumps(encoder))
        self.assertTrue(hasattr(branston, 'logger'))

    def test_020_feature_length(self):
        # hardcode value here, if this genuinely needs to change at least
        # the test will make develop think twice about consequences
        kwargs = {'dtypes': ['r9','r10'], 'include_dwells': False}
        encoder = medaka.features.ReadAlignmentFeatureEncoder(**kwargs)
        self.assertEqual(encoder.feature_vector_length, 5)

    def test_021_feature_length_dwells_no_dtype(self):
        # hardcode value here, if this genuinely needs to change at least
        # the test will make develop think twice about consequences
        kwargs = {'include_dwells': True}
        encoder = medaka.features.ReadAlignmentFeatureEncoder(**kwargs)
        self.assertEqual(encoder.feature_vector_length, 5)

    def test_022_feature_length_dwells_no_dtype(self):
        # hardcode value here, if this genuinely needs to change at least
        # the test will make develop think twice about consequences
        kwargs = {'dtypes': ['r9','r10'], 'include_dwells': True, 'include_haplotype': True}
        encoder = medaka.features.ReadAlignmentFeatureEncoder(**kwargs)
        self.assertEqual(encoder.feature_vector_length, 7)

    def test_030_bams_to_training_samples_simple(self):
        reads_bam = tempfile.NamedTemporaryFile(suffix='.bam').name
        truth_bam = tempfile.NamedTemporaryFile(suffix='.bam').name

        # we had a bug caused by missing qualities and bad indexing...
        data = copy.deepcopy(simple_data['calls'])
        data[0]['quality'] = None

        create_simple_bam(reads_bam, data)
        create_simple_bam(
            truth_bam, [simple_data['truth']])
        encoder = medaka.features.ReadAlignmentFeatureEncoder(include_dwells=False)
        label_scheme = medaka.labels.HaploidLabelScheme()
        region = Region('ref', 0, 100)
        result = encoder.bams_to_training_samples(
            truth_bam, reads_bam, region, label_scheme, min_length=0)[0]

        expected = Sample(
            ref_name='ref',
            features=np.array([[
                # read 1
                    [1, 0, 1, 40],
                    [2, 0, 1, 40],
                    [1, 0, 1, 40],
                    [4, 0, 1, 40],
                    [5, 0, 1, 40],
                    [3, 0, 1, 40],
                    [1, 0, 1, 40],
                    [4, 0, 1, 40],
                    [3, 0, 1, 40],
                ], [
                # read 2
                    [1, 0, 1, 10],
                    [2, 1, 1, 10],
                    [1, 4, 1, 10],
                    [5, 0, 1, 10],
                    [5, 0, 1, 10],
                    [3, 1, 1, 10],
                    [1, 1, 1, 10],
                    [4, 1, 1, 10],
                    [3, 2, 1, 10],
                ], [
                # read 3
                    [1, 2, 0, 16],
                    [2, 1, 0, 16],
                    [1, 4, 0, 16],
                    [4, 5, 0, 16],
                    [1, 1, 0, 16],
                    [3, 1, 0, 16],
                    [1, 1, 0, 16],
                    [4, 2, 0, 16],
                    [3, 1, 0, 16],
                ], [
                # read 4
                    [1,  2, 0, 24],
                    [2,  1, 0, 24],
                    [1,  4, 0, 24],
                    [2,  1, 0, 24],
                    [5,  0, 0, 24],
                    [3,  1, 0, 24],
                    [1,  1, 0, 24],
                    [4,  2, 0, 24],
                    [3,  1, 0, 24],
                ]],
                dtype='int8'
            ).swapaxes(0, 1),
            # the two insertions with respect to the draft are dropped
            labels=np.array([1, 2, 1, 4, 1, 3, 1, 4, 3]),  # A C A T A G A T C
            ref_seq=None,
            positions=np.array([
                (0, 0), (1, 0), (2, 0), (3, 0), (3, 1), (4, 0), (5, 0), (6, 0), (7, 0)],
                dtype=[('major', '<i8'), ('minor', '<i8')]),
            label_probs=None,
            depth=np.array(9 * [10]),
        )

        np.testing.assert_equal(expected.features.shape, (9, 4, 4))
        np.testing.assert_equal(result.features.shape, expected.features.shape)

        np.testing.assert_equal(result.labels, expected.labels)
        np.testing.assert_equal(result.positions, expected.positions)
        np.testing.assert_equal(result.features, expected.features)

    def test_031_bams_to_training_samples_dwells(self):
        reads_bam = tempfile.NamedTemporaryFile(suffix='.bam').name
        truth_bam = tempfile.NamedTemporaryFile(suffix='.bam').name

        # we had a bug caused by missing qualities and bad indexing...
        data = copy.deepcopy(simple_data['calls'])
        data[0]['quality'] = None

        create_simple_bam(reads_bam, data)
        create_simple_bam(
            truth_bam, [simple_data['truth']])
        encoder = medaka.features.ReadAlignmentFeatureEncoder(include_dwells=True)
        label_scheme = medaka.labels.HaploidLabelScheme()
        region = Region('ref', 0, 100)
        result = encoder.bams_to_training_samples(
            truth_bam, reads_bam, region, label_scheme, min_length=0)[0]

        expected = Sample(
            ref_name='ref',
            features=np.array([[
                # read 1
                    [1, 0,  1, 40,  3],
                    [2, 0,  1, 40,  2],
                    [1, 0,  1, 40,  1],
                    [4, 0,  1, 40,  1],
                    [5, 0,  1, 40,  0],
                    [3, 0,  1, 40,  4],
                    [1, 0,  1, 40,  2],
                    [4, 0,  1, 40,  2],
                    [3, 0,  1, 40,  3],
                ], [
                # read 2
                    [1,  0,  1, 10,  1],
                    [2,  1,  1, 10,  3],
                    [1,  4,  1, 10,  1],
                    [5,  0,  1, 10,  0],
                    [5,  0,  1, 10,  0],
                    [3,  1,  1, 10,  4],
                    [1,  1,  1, 10,  2],
                    [4,  1,  1, 10,  1],
                    [3,  2,  1, 10,  3],
                ], [
                # read 3
                    [1,  2, 0, 16, 2],
                    [2,  1, 0, 16, 3],
                    [1,  4, 0, 16, 1],
                    [4,  5, 0, 16, 2],
                    [1,  1, 0, 16, 4],
                    [3,  1, 0, 16, 1],
                    [1,  1, 0, 16, 2],
                    [4,  2, 0, 16, 2],
                    [3,  1, 0, 16, 2],
                ], [
                # read 4, dwells should be blank because of malformed move table
                    [1,  2, 0, 24,  0],
                    [2,  1, 0, 24,  0],
                    [1,  4, 0, 24,  0],
                    [2,  1, 0, 24,  0],
                    [5,  0, 0, 24,  0],
                    [3,  1, 0, 24,  0],
                    [1,  1, 0, 24,  0],
                    [4,  2, 0, 24,  0],
                    [3,  1, 0, 24,  0],
                ]],
                dtype='int8'
            ).swapaxes(0, 1),
            # the two insertions with respect to the draft are dropped
            labels=np.array([1, 2, 1, 4, 1, 3, 1, 4, 3]),  # A C A T A G A T C
            ref_seq=None,
            positions=np.array([
                (0, 0), (1, 0), (2, 0), (3, 0), (3, 1), (4, 0), (5, 0), (6, 0), (7, 0)],
                dtype=[('major', '<i8'), ('minor', '<i8')]),
            label_probs=None,
            depth=np.array(9 * [10]),
        )

        np.testing.assert_equal(expected.features.shape, (9, 4, 5))
        np.testing.assert_equal(result.features.shape, expected.features.shape)

        print(result.features.swapaxes(0, 1))

        np.testing.assert_equal(result.labels, expected.labels)
        np.testing.assert_equal(result.positions, expected.positions)
        np.testing.assert_equal(result.features, expected.features)

    def test_032_bams_to_training_samples_regression(self):
        encoder = medaka.features.ReadAlignmentFeatureEncoder(include_dwells=False)

        label_scheme = medaka.labels.HaploidLabelScheme()
        region = Region(
            'utg000001l',
            149744, 318288)
        result = encoder.bams_to_training_samples(
            __reads_truth__, __reads_bam__, region, label_scheme)[0]

        expected_feature_shape = (177981, 43, 4)
        got_feature_shape = result.features.shape
        self.assertEqual(expected_feature_shape, got_feature_shape)

        expected_label_shape = (177981,)
        got_label_shape = result.labels.shape
        self.assertEqual(expected_label_shape, got_label_shape)


    def test_033_bams_to_training_samples_regression_max_depth(self):
        label_scheme = medaka.labels.HaploidLabelScheme()
        region = Region(
            'utg000001l',
            149744, 318288)
        test_depth = 5

        encoder = medaka.features.ReadAlignmentFeatureEncoder(
            max_reads=test_depth, include_dwells=False)
        result = encoder.bams_to_training_samples(
            __reads_truth__, __reads_bam__, region, label_scheme)[0]
        expected_feature_shape = (177981, test_depth, 4)
        got_feature_shape = result.features.shape
        self.assertEqual(expected_feature_shape, got_feature_shape)


class SampleGenerator(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.expected_width = 86294

    def test_000_basic_sample_gen(self):
        encoder = medaka.features.ReadAlignmentFeatureEncoder()
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
        encoder = medaka.features.ReadAlignmentFeatureEncoder()
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
        encoder = medaka.features.ReadAlignmentFeatureEncoder()
        with self.assertRaises(ValueError) as context:
            sample_gen = medaka.features.SampleGenerator(
                __reads_bam__, __region__, encoder, truth_bam=__reads_truth__)


class ReadMatrixSplitAndJoin(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.counts = np.array([[
            # read 1
                [1,  2,  1, 40, 0],
                [2,  1,  1, 40, 0],
                [1,  4,  1, 40, 0],
                [4,  5,  1, 40, 0],
                [5,  0,  1, 40, 0],
                [3,  1,  1, 40, 0],
                [1,  1,  1, 40, 0],
                [4,  2,  1, 40, 0],
                [3,  1,  1, 40, 0],
            ], [
            # read 2
                [1,  0,  1, 10, 0],
                [2,  1,  1, 10, 0],
                [1,  4,  1, 10, 0],
                [5,  0,  1, 10, 0],
                [5,  0,  1, 10, 0],
                [3,  1,  1, 10, 0],
                [1,  1,  1, 10, 0],
                [4,  1,  1, 10, 0],
                [3,  2,  1, 10, 0],
            ], [
            # read 3
                [1,  2,  0, 16, 0],
                [2,  1,  0, 16, 0],
                [1,  4,  0, 16, 0],
                [4,  5,  0, 16, 0],
                [1,  1,  0, 16, 0],
                [3,  1,  0, 16, 0],
                [1,  1,  0, 16, 0],
                [4,  2,  0, 16, 0],
                [3,  1,  0, 16, 0],
            ], [
            # read 4
                [1,  2,  0, 24, 0],
                [2,  1,  0, 24, 0],
                [1,  4,  0, 24, 0],
                [2,  1,  0, 24, 0],
                [5,  0,  0, 24, 0],
                [3,  1,  0, 24, 0],
                [1,  1,  0, 24, 0],
                [4,  2,  0, 24, 0],
                [3,  1,  0, 24, 0],
            ]],
            dtype='int8'
        ).swapaxes(0, 1)
        cls.positions = np.array(
            [(0, 0), (1, 0), (2, 0), (3, 0), (3, 1),
             (4, 0), (5, 0), (6, 0), (7, 0)],
             dtype=[('major', '<i8'), ('minor', '<i8')])

    def test_001_reorder_matching(self):
        counts_chunks = [self.counts[:5], self.counts[5:]]
        read_ids = [
            (
                np.array(["read_id_1", "read_id_2", "read_id_3", "read_id_4"]), #chunk1 in
                np.array(["read_id_1", "read_id_2", "read_id_3", "read_id_4"])  #chunk1 out
            ), (
                np.array(["read_id_1", "read_id_2", "read_id_3", "read_id_4"]), #chunk2 in
                np.array(["read_id_1", "read_id_2", "read_id_3", "read_id_4"])  #chunk2 out
            )
        ]
        reorder = medaka.features._reorder_reads(counts_chunks, read_ids)

        np.testing.assert_equal(reorder[0], self.counts[:5])
        np.testing.assert_equal(reorder[1], self.counts[5:])

    def test_002_reorder_swapped(self):
        counts_chunks = [self.counts[:5], self.counts[5:, [0, 2, 3, 1], :]]
        read_ids = [
            (
                np.array(["read_id_1", "read_id_2", "read_id_3", "read_id_4"]), #chunk1 in
                np.array(["read_id_1", "read_id_2", "read_id_3", "read_id_4"])  #chunk1 out
            ), (
                np.array(["read_id_1", "read_id_3", "read_id_4", "read_id_2"]), #chunk2 in
                np.array(["read_id_1", "read_id_3", "read_id_4", "read_id_2"])  #chunk2 out
            )
        ]
        reorder = medaka.features._reorder_reads(counts_chunks, read_ids)

        np.testing.assert_equal(reorder[0], self.counts[:5])
        np.testing.assert_equal(reorder[1], self.counts[5:])

    def test_003_reorder_missing_right(self):
        counts_chunks = [self.counts[:5], self.counts[5:, [0, 2, 3], :]]
        read_ids = [
            (
                np.array(["read_id_1", "read_id_2", "read_id_3", "read_id_4"]), #chunk1 in
                np.array(["read_id_1", "blank", "read_id_3", "read_id_4"])  #chunk1 out
            ), (
                np.array(["read_id_1", "read_id_3", "read_id_4"]), #chunk2 in
                np.array(["read_id_1", "read_id_3", "read_id_4"])  #chunk2 out
            )
        ]
        reorder = medaka.features._reorder_reads(counts_chunks, read_ids)
        expected_output = np.concatenate([
                self.counts[5:, [0], :],
                np.zeros_like(self.counts[5:, [0], :]),
                self.counts[5:, 2:, :],
            ],
            axis=1
        )

        np.testing.assert_equal(reorder[0], self.counts[:5])
        np.testing.assert_equal(reorder[1], expected_output)

    def test_004_reorder_missing_left(self):
        counts_chunks = [self.counts[:5, [0, 1, 3]], self.counts[5:]]
        read_ids = [
            (
                np.array(["read_id_1", "read_id_2", "read_id_4"]), #chunk1 in
                np.array(["read_id_1", "read_id_2", "read_id_4"])  #chunk1 out
            ), (
                np.array(["read_id_1", "read_id_2", "read_id_3", "read_id_4"]), #chunk2 in
                np.array(["read_id_1", "read_id_2", "read_id_3", "read_id_4"])  #chunk2 out
            )
        ]
        reorder = medaka.features._reorder_reads(counts_chunks, read_ids)
        expected_output = self.counts[5:, [0, 1, 3, 2], :]

        np.testing.assert_equal(reorder[0], self.counts[:5, [0, 1, 3]])
        np.testing.assert_equal(reorder[1], expected_output)

    def test_005_reorder_mismatch(self):
        counts_chunks = [
            np.concatenate([
                self.counts[:5, :2],
                np.zeros_like(self.counts[:5, [0]]),
                self.counts[:5, 2:],
            ], axis=1),
            np.concatenate([
                self.counts[5:],
                np.zeros_like(self.counts[5:, [0]]),
            ], axis=1),
        ]
        read_ids = [
            (
                np.array(["read_id_1", "read_id_2", "blank_1", "read_id_3", "read_id_4"]), #chunk1 in
                np.array(["read_id_1", "read_id_2", "blank_1", "read_id_3", "read_id_4"]), #chunk1 out
            ), (
                np.array(["read_id_1", "read_id_2", "read_id_3", "read_id_4", "blank_2"]), #chunk2 in
                np.array(["read_id_1", "read_id_2", "read_id_3", "read_id_4", "blank_2"]), #chunk2 in
            )
        ]
        reorder = medaka.features._reorder_reads(counts_chunks, read_ids)
        expected_output1 = counts_chunks[0]
        expected_output2 = np.concatenate([
            # expected: ["read_id_1", "read_id_2", "blank_2", "read_id_3", "read_id_4"]), #chunk2 in
            self.counts[5:, :2],
            np.zeros_like(self.counts[5:, [0]]),
            self.counts[5:, 2:],
        ], axis=1)

        np.testing.assert_equal(reorder[0], expected_output1)
        np.testing.assert_equal(reorder[1], expected_output2)


class FeaturesFromBam(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        bam_fname = tempfile.NamedTemporaryFile(suffix='.bam').name
        create_simple_bam(bam_fname, simple_data['calls'])
        cls.region = Region(
            'ref',
            start=0, end=8)
        cls.bam = bam_fname
        cls.expected_counts = np.array([[
            # read 1
                [1,  2,  1, 40],
                [2,  1,  1, 40],
                [1,  4,  1, 40],
                [4,  5,  1, 40],
                [5,  0,  1, 40],
                [3,  1,  1, 40],
                [1,  1,  1, 40],
                [4,  2,  1, 40],
                [3,  1,  1, 40],
            ], [
            # read 2
                [1,  0,  1, 10],
                [2,  1,  1, 10],
                [1,  4,  1, 10],
                [5,  0,  1, 10],
                [5,  0,  1, 10],
                [3,  1,  1, 10],
                [1,  1,  1, 10],
                [4,  1,  1, 10],
                [3,  2,  1, 10],
            ], [
            # read 3
                [1,  2,  0, 16],
                [2,  1,  0, 16],
                [1,  4,  0, 16],
                [4,  5,  0, 16],
                [1,  1,  0, 16],
                [3,  1,  0, 16],
                [1,  1,  0, 16],
                [4,  2,  0, 16],
                [3,  1,  0, 16],
            ], [
            # read 4
                [1,  2,  0, 24],
                [2,  1,  0, 24],
                [1,  4,  0, 24],
                [2,  1,  0, 24],
                [5,  0,  0, 24],
                [3,  1,  0, 24],
                [1,  1,  0, 24],
                [4,  2,  0, 24],
                [3,  1,  0, 24],
            ]],
            dtype='int8'
        ).swapaxes(0, 1)
        cls.expected_positions = np.array(
            [(0, 0), (1, 0), (2, 0), (3, 0), (3, 1),
             (4, 0), (5, 0), (6, 0), (7, 0)],
             dtype=[('major', '<i8'), ('minor', '<i8')])

    def test_000_counts_no_dtypes(self):
        counts, positions = medaka.features.read_alignment_matrix(
            self.region, self.bam)[0]

        np.testing.assert_equal(counts.shape, self.expected_counts.shape)
        self.assertTrue(np.array_equal(counts, self.expected_counts))
        self.assertTrue(np.array_equal(positions, self.expected_positions))

    def test_010_tags(self):
        """Only two reads have tag AA=1."""
        counts, positions = medaka.features.read_alignment_matrix(
            self.region, self.bam, tag_name='AA', tag_value=1,
            keep_missing=False)[0]
        expected_number_reads = 2
        self.assertEqual(counts.shape[1], expected_number_reads)

    def test_020_tags(self):
        """If keep_missing is True, 3 reads are expected"""
        counts, positions = medaka.features.read_alignment_matrix(
            self.region, self.bam, tag_name='AA', tag_value=1,
            keep_missing=True)[0]
        expected_number_reads = 3
        self.assertEqual(counts.shape[1], expected_number_reads)

    def test_040_dtypes(self):
        counts, positions = medaka.features.read_alignment_matrix(
            self.region, self.bam, dtype_prefixes=['r9', 'r10'])[0]

        expected_shape = (9, 4, 5)
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
        counts, positions = medaka.features.read_alignment_matrix(
            region, bam)[0]
        self.assertEqual(counts.shape[1], 2*self.expected_counts.shape[1])
        self.assertTrue(np.array_equal(counts[:, [0, 1, 4, 5], :], self.expected_counts))
        self.assertTrue(np.array_equal(counts[:, [2, 3, 6, 7], :], self.expected_counts))
        self.assertTrue(np.array_equal(positions, self.expected_positions))

        # use one or other
        for rg in ('first', 'second'):
            counts, positions = medaka.features.read_alignment_matrix(
                region, bam, read_group=rg)[0]
            self.assertTrue(np.array_equal(counts, self.expected_counts))
            self.assertTrue(np.array_equal(positions, self.expected_positions))

        # use a missing one
        result = medaka.features.read_alignment_matrix(
            region, bam, read_group='nonsense')
        self.assertTrue(len(result) == 0)


    def test_060_counts_chunked(self):
        """Test that we get the same counts and positions when region_split < region.size

        i.e. the stitching back together of chunks processed in parallel works.
        """
        chunks = medaka.features.read_alignment_matrix(
            self.region, self.bam, region_split=3)
        self.assertEqual(len(chunks), 1)
        counts, positions = chunks[0]
        self.assertTrue(np.array_equal(counts, self.expected_counts))
        self.assertTrue(np.array_equal(positions, self.expected_positions))
