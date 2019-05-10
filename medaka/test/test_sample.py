from copy import deepcopy
import unittest
import numpy as np
import os

from medaka.common import Sample, OverlapException

root_dir = os.path.abspath(os.path.dirname(__file__))
test_file = os.path.join(root_dir, 'data/test_probs.hdf')


class TestSample(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pos1 = np.array([(0, 0), (0, 1), (1, 0), (2, 0), (2, 1), (2, 2), (3, 0), (4, 0), (4, 1), (4, 2), (4, 3)],
                        dtype=[('major', int), ('minor', int)])
        pos2 = np.array([(4, 1), (4, 2), (4, 3), (4, 5), (5, 0), (6, 0), (6, 1), (7, 0)],
                        dtype=[('major', int), ('minor', int)])
        pos3 = np.array([(5, 0), (6, 0), (6, 1), (7, 0), (7, 1)],
                        dtype=[('major', int), ('minor', int)])
        cls.samples = []
        data_dim = 10
        for pos in pos1, pos2, pos3:
            data = np.random.random_sample(size=data_dim*len(pos)).reshape((len(pos), data_dim))
            cls.samples.append(
                Sample(ref_name='contig1', features=data, ref_seq=None, labels=data, positions=pos, label_probs=data)
            )


    def test_eq(self):
        for sample in self.samples:
           self.assertEqual(sample, deepcopy(sample))


    def test_not_eq(self):
        for sample1, sample2 in zip(self.samples[:-1], self.samples[1:]):
           self.assertFalse(sample1 == sample2)


    def test_first_pos(self):
        for expt, sample in zip([(0, 0), (4, 1)], self.samples):
           self.assertEqual(expt, sample.first_pos)


    def test_last_pos(self):

        for expt, sample in zip([(4, 3), (7, 0)], self.samples):
           self.assertEqual(expt, sample.last_pos)


    def test_span(self):
        for expt, sample in zip([4, 3], self.samples):
           self.assertEqual(expt, sample.span)


    def test_size(self):
        for sample in self.samples:
           self.assertEqual(len(sample.positions), sample.size)


    def test_is_empty(self):
        self.assertFalse(self.samples[0].is_empty)
        empty_sample = Sample('contig', None, None, None, self.samples[0][0:0], None)
        self.assertTrue(empty_sample.is_empty)


    def test_name(self):
        for expt, sample in zip(['contig1:0.0-4.3', 'contig1:4.1-7.0'], self.samples):
           self.assertEqual(expt, sample.name)


    def test_decode_sample_name(self):
        expected = [
            {'ref_name': 'contig1', 'start': '0.0', 'end': '4.3'},
            {'ref_name': 'contig1', 'start': '4.1', 'end': '7.0'},
        ]

        for expt, sample in zip(expected, self.samples):
           self.assertEqual(expt, Sample.decode_sample_name(sample.name))


    def test_slice(self):
        sl = slice(1, 5)
        sliced = self.samples[0].slice(sl)
        self.assertEqual('contig1:0.1-2.1', sliced.name)
        # positions are tuples, need to convert from array to list to use assertSequenceEqual
        self.assertSequenceEqual(list(sliced.positions), list(self.samples[0].positions[sl]))
        for field in ['features', 'labels', 'label_probs']:  # all floats
            np.testing.assert_allclose(getattr(sliced, field), getattr(self.samples[0], field)[sl])
        for field in ['ref_name', 'ref_seq']:  # non-sliceable (ref_seq is None)
            self.assertEqual(getattr(sliced, field), getattr(self.samples[0], field))


    def test_from_samples(self):
        # check we can concat a single sample
        concat1 = Sample.from_samples([self.samples[0]])
        self.assertEqual(concat1, self.samples[0])

        # check we can concat 3 samples
        slices = [slice(0, 5), slice(5, 8), slice(8, 11)]  # these should span entire samples[0]
        sliced = [self.samples[0].slice(sl) for sl in slices]

        concat3 = Sample.from_samples(sliced)
        self.assertEqual(concat3, self.samples[0])

        # also check raises an exception if not forward abutting
        self.assertRaises(ValueError, Sample.from_samples, sliced[::-1])


    def test_relative_position(self):
        # all 9 cases plus major and minor variations of overlapping, abutting and gapped
        #   different_ref_name = 'Samples come from different reference contigs.'
        #   forward_overlap = 'The end of s1 overlaps the start of s2.'
        #   reverse_overlap = 'The end of s2 overlaps the start of s1.'
        #   forward_abutted = 'The end of s1 abuts the start of s2.'
        #   reverse_abutted = 'The end of s2 abuts the start of s1.'
        #   forward_gapped = 's2 follows s1 with a gab inbetween.'
        #   reverse_gapped = 's1 follows s2 with a gab inbetween.'
        #   s2_within_s1 = 's2 is fully contained within s1.'
        #   s1_within_s2 = 's1 is fully contained within s2.'

        slices = [
            slice(0, 3), # (0,0) -> (1,0) in self.samples[0]
            slice(3, 5), # (2,0) -> (2,1) in self.samples[0]
            slice(5, 8), # (2,2) -> (4,0) in self.samples[0]
            slice(8, 11) # (4.1) -> (4.3) in self.samples[0]
        ]
        sliced = [self.samples[0].slice(sl) for sl in slices]
        sample_dict = self.samples[0]._asdict()
        sample_dict['ref_name'] = 'other'
        sample_other = Sample(**sample_dict)

        samples_expt = [
            ([self.samples[0], sample_other], Sample.Relationship.different_ref_name),
            (self.samples[:2], Sample.Relationship.forward_overlap),  # overlap of minor positions
            (self.samples[1:], Sample.Relationship.forward_overlap),  # overlap of major positions
            (self.samples[:2][::-1], Sample.Relationship.reverse_overlap), # overlap of minor positions
            (self.samples[1:][::-1], Sample.Relationship.reverse_overlap),  # overlap of major positions
            (sliced[:2], Sample.Relationship.forward_abutted),  # (1,0) -> (2,0)
            (sliced[1:3], Sample.Relationship.forward_abutted),  # (2,1) -> (2,2)
            (sliced[::-1][2:], Sample.Relationship.reverse_abutted), # (2,0) -> (1,0)
            (sliced[::-1][1:3], Sample.Relationship.reverse_abutted), # (2,2) -> (2,1)
            ([sliced[0], sliced[2]], Sample.Relationship.forward_gapped), # (1,0) -> (2,2)
            ([self.samples[0].slice(slice(0, 4)), sliced[2]], Sample.Relationship.forward_gapped), # (2,0) -> (2,2)
            ([sliced[0], self.samples[0].slice(slice(6, None))], Sample.Relationship.forward_gapped), # (1,0) -> (3,0)
            ([sliced[2], self.samples[0].slice(slice(0, 4))], Sample.Relationship.reverse_gapped), # (2,2) -> (2,0)
            ([self.samples[0].slice(slice(6, None)), sliced[0]], Sample.Relationship.reverse_gapped), # (3,0) -> (1,0)
            ([self.samples[0], sliced[0]], Sample.Relationship.s2_within_s1),
            ([self.samples[0], sliced[1]], Sample.Relationship.s2_within_s1),
            ([self.samples[0], sliced[2]], Sample.Relationship.s2_within_s1),
            ([self.samples[0], sliced[3]], Sample.Relationship.s2_within_s1),
            ([sliced[0], self.samples[0]], Sample.Relationship.s1_within_s2),
            ([sliced[1], self.samples[0]], Sample.Relationship.s1_within_s2),
            ([sliced[2], self.samples[0]], Sample.Relationship.s1_within_s2),
            ([sliced[3], self.samples[0]], Sample.Relationship.s1_within_s2),
        ]
        for samples, expt in samples_expt:
            self.assertIs(Sample.relative_position(*samples), expt)


    def test_overlap_indices(self):
        # check we get right thing for overlap and abutting major and minor
        # and that we get an exception if we give it any other case.
        slices = [
            slice(0, 3), # (0,0) -> (1,0) in self.samples[0]
            slice(3, 5), # (2,0) -> (2,1) in self.samples[0]
            slice(5, 8), # (2,2) -> (4,0) in self.samples[0]
            slice(8, 11) # (4.1) -> (4.3) in self.samples[0]
        ]
        sliced = [self.samples[0].slice(sl) for sl in slices]

        samples_expt = [
            (self.samples[:2], (9, 1)),  # overlap of minor inds with odd number of overlapping positions
            (self.samples[1:], (6, 2)),  # overlap of major inds with even number of overlapping positions
            (sliced[:2], (None, None)),  # abuts
            (sliced[1:3], (None, None)),  # abuts
            (sliced[2:], (None, None)),  # abuts
        ]
        for samples, expt in samples_expt:
            self.assertEqual(Sample.overlap_indices(*samples), expt)

        self.assertRaises(OverlapException, Sample.overlap_indices, samples[1], samples[0])


    def test_chunks(self):
        chunks = self.samples[0].chunks(chunk_len=4, overlap=2)
        slices = [slice(i, i + 4) for i in range(0, len(self.samples[0]) - 4, 2)]
        sliced = [self.samples[0].slice(sl) for sl in slices]
        for expt, got in zip(sliced, chunks):
            self.assertEqual(got, expt)


if __name__ == '__main__':
    unittest.main()
