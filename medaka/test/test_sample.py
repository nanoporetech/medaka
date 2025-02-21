from copy import deepcopy
import unittest
import numpy as np
import os

from medaka.common import Region, Relationship, Sample, OverlapException
import medaka.labels

root_dir = os.path.abspath(os.path.dirname(__file__))
test_file = os.path.join(root_dir, 'data/test_probs.hdf')

def get_test_samples(num_train=1000, num_test=500, chunk_size=1000, read_level_features=False):
    # Adapted from medaka v2.0.1 - the dataloaders now expect a list of samples
    # instead of a tuple of lists of features and labels (these are then
    # converted into a medaka.torch.Batch object which is passed to the model's
    # process_batch or predict_on_batch methods, depending on whether training
    # or doing inference).
    num_samples = num_train + num_test
    train_samples = []
    valid_samples = []
    # insert minor positions with 10% probability
    is_major_position = np.random.random(size=(chunk_size)) < 0.9
    is_major_position[0] = True
    is_minor_position = ~is_major_position
    major_position_index = np.cumsum(is_major_position.astype(int))-1
    indicies = [(0,0)]
    minor_idx = 0
    for i in range(1, chunk_size):
        prev_index = indicies[-1]
        if major_position_index[i] == prev_index[0]:
            minor_idx += 1
        else:
            minor_idx = 0
        indicies.append((major_position_index[i], minor_idx))
    pos = np.array(indicies, dtype=[('major', int), ('minor', int)])

    for i in range(num_samples):
        labels = np.random.randint(0, 5, size=(num_samples, chunk_size))
        # choose random depth between 10 & 20
        depths = np.ones_like(labels)*np.random.randint(10, 20)
        if read_level_features:

            features = np.zeros((num_samples, chunk_size, np.max(depths), 5)
            )
            # bases
            features[:,:,:, 0] = np.random.randint(0, 5, size=(num_samples, chunk_size, np.max(depths)))
            # q scores
            features[:,:,:, 1] = np.random.randint(0, 50, size=(num_samples, chunk_size, np.max(depths)))
        else:
            features = np.random.random(size=(num_samples, chunk_size,10))
            depths = np.ones_like(labels)

        
        sample = medaka.common.Sample(
            ref_name=f"sample_{i}",
            features=features[i],
            labels=labels[i],
            ref_seq="",
            positions=pos,
            label_probs=None,
            depth=depths
        )
        if i < num_train:
            train_samples.append(sample)
        else:
            valid_samples.append(sample)

    return train_samples, valid_samples

class TestRegion(unittest.TestCase):
    def test_self_overlap(self):
        cases = [
            Region('contig1', 0, 100),
            Region('contig1', 0, None),
            Region('contig1', None, None)]
        for c in cases:
            self.assertTrue(c.overlaps(c))

    def test_overlaps(self):
        a = Region('contig1', 50, 100)
        cases = [
            Region('contig1', 49, 150),
            Region('contig1', 50, 150),
            Region('contig1', 51, 150),
            Region('contig1', 50, 99),
            Region('contig1', 50, 100),
            Region('contig1', 50, 101),
        ]
        for c in cases:
            self.assertTrue(a.overlaps(c))
            self.assertTrue(c.overlaps(a))

    def test_no_overlap(self):
        a = Region('contig1', 50, 100)
        cases = [
            Region('contig2', 50, 100),
            Region('contig1', 0, 50),
            Region('contig1', 100, 150),
        ]
        for c in cases:
            self.assertFalse(a.overlaps(c))
            self.assertFalse(c.overlaps(a))

    def test_name(self):
        # tests for Region.from_string are in docs
        # here we test round-tripping
        cases = [
            ['contig1:50-100'] * 2,
            ['contig1:50-'] * 2,
            ['contig1:-100', 'contig1:0-100'],
            ['contig1', 'contig1:0-']
        ]
        for orig, parsed in cases:
            a = Region.from_string(orig)
            self.assertEqual(a.name, parsed)

    def test_size(self):
        a = Region('contig1', 50, 100)
        self.assertEqual(a.size, 50)

    def test_split(self):
        a = Region('contig1', 50, 100)
        regs = a.split(10)
        starts = list(range(50, 100, 10))
        ends = [s + 10 for s in starts]
        self.assertEqual([x.start for x in regs], starts)
        self.assertEqual([x.end for x in regs], ends)

        # case with overlap, and triggering duplicate condition
        regs = a.split(10, 5)
        starts = list(range(50, 95, 5))
        ends = [s + 10 for s in starts]
        self.assertEqual([x.start for x in regs], starts)
        self.assertEqual([x.end for x in regs], ends)

        # case with odd sized remainder - fixed size
        regs = a.split(7)
        starts = [50, 57, 64, 71, 78, 85, 92, 93]
        ends = [57, 64, 71, 78, 85, 92, 99, 100]
        self.assertEqual([x.start for x in regs], starts)
        self.assertEqual([x.end for x in regs], ends)

        # case with odd sized remainder - variable size
        regs = a.split(7, fixed_size=False)
        starts = [50, 57, 64, 71, 78, 85, 92, 99]
        ends = [57, 64, 71, 78, 85, 92, 99, 100]
        self.assertEqual([x.start for x in regs], starts)
        self.assertEqual([x.end for x in regs], ends)


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
            depth = np.array(len(pos) * [10], dtype=int)
            cls.samples.append(
                Sample(ref_name='contig1', features=data, ref_seq=None, labels=data, positions=pos, label_probs=data, depth=depth)
            )

    def test_eq(self):
        for sample in self.samples:
           self.assertEqual(sample, sample)
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
        empty_sample = Sample('contig', None, None, None, self.samples[0][0:0], None, None)
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
        sample_other = self.samples[0].amend(ref_name='other')

        samples_expt = [
            ([self.samples[0], sample_other], Relationship.different_ref_name),
            (self.samples[:2], Relationship.forward_overlap),  # overlap of minor positions
            (self.samples[1:], Relationship.forward_overlap),  # overlap of major positions
            (self.samples[:2][::-1], Relationship.reverse_overlap), # overlap of minor positions
            (self.samples[1:][::-1], Relationship.reverse_overlap),  # overlap of major positions
            (sliced[:2], Relationship.forward_abutted),  # (1,0) -> (2,0)
            (sliced[1:3], Relationship.forward_abutted),  # (2,1) -> (2,2)
            (sliced[::-1][2:], Relationship.reverse_abutted), # (2,0) -> (1,0)
            (sliced[::-1][1:3], Relationship.reverse_abutted), # (2,2) -> (2,1)
            ([sliced[0], sliced[2]], Relationship.forward_gapped), # (1,0) -> (2,2)
            ([self.samples[0].slice(slice(0, 4)), sliced[2]], Relationship.forward_gapped), # (2,0) -> (2,2)
            ([sliced[0], self.samples[0].slice(slice(6, None))], Relationship.forward_gapped), # (1,0) -> (3,0)
            ([sliced[2], self.samples[0].slice(slice(0, 4))], Relationship.reverse_gapped), # (2,2) -> (2,0)
            ([self.samples[0].slice(slice(6, None)), sliced[0]], Relationship.reverse_gapped), # (3,0) -> (1,0)
            ([self.samples[0], sliced[0]], Relationship.s2_within_s1),
            ([self.samples[0], sliced[1]], Relationship.s2_within_s1),
            ([self.samples[0], sliced[2]], Relationship.s2_within_s1),
            ([self.samples[0], sliced[3]], Relationship.s2_within_s1),
            ([sliced[0], self.samples[0]], Relationship.s1_within_s2),
            ([sliced[1], self.samples[0]], Relationship.s1_within_s2),
            ([sliced[2], self.samples[0]], Relationship.s1_within_s2),
            ([sliced[3], self.samples[0]], Relationship.s1_within_s2),
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
            (self.samples[:2], (9, 1, False)),  # overlap of minor inds with odd number of overlapping positions
            (self.samples[1:], (6, 2, False)),  # overlap of major inds with even number of overlapping positions
            (sliced[:2], (None, None, False)),  # abuts
            (sliced[1:3], (None, None, False)),  # abuts
            (sliced[2:], (None, None, False)),  # abuts
        ]
        for samples, expt in samples_expt:
            self.assertEqual(Sample.overlap_indices(*samples), expt)

        self.assertRaises(OverlapException, Sample.overlap_indices, samples[1], samples[0])

    def test_messy_overlap(self):
        dtype = [('major', int), ('minor', int)]
        pos = [
            np.array([
                (0, 0), (1, 0), (2, 0), (2, 1), (2, 2), (3, 0), (4, 0), (4, 1), (4, 2), (4, 3),
                (5, 0), (6, 0), (7, 0), (8, 0), (8, 1), (9, 0),
            ], dtype=dtype),
            np.array([
                (3, 0), (4, 0), (4, 1), (4, 2),  # (4,3) missing
                (5, 0), (6, 0), (7, 0), (8, 0), (8, 1), (9, 0), (10, 0), (10, 1), (10, 2),
            ], dtype=dtype),
            np.array([
                (3, 0), (4, 0), (4, 1), (4, 2),
                (5, 0), (6, 0), (6, 1), (7, 0), (7, 1), (8, 0), (8, 1), (9, 0), (10, 0), (10, 1), (10, 2),
            ], dtype=dtype),
            np.array([
                (3, 0), (4, 0), (4, 1), (4, 2),  # (4,3) missing
                (5, 0), (5, 1), (6, 0), (6, 1), (7, 0), (7, 1), (8, 0), (8, 1), (9, 0), (10, 0), (10, 1), (10, 2),
            ], dtype=dtype),
            np.array([
                (3, 0), (4, 0), (4, 1), (4, 2),  # (4,3) missing
                (5, 0), (5, 1), (6, 0), (6, 1), (7, 0), (7, 1), (8, 0), (9, 0), (10, 0), (10, 1), (10, 2),
            ], dtype=dtype),
        ]

        sample = [Sample(ref_name='contig1', features=None, ref_seq=None, labels=None,
                  positions=p, label_probs=None, depth=None) for p in pos]

        expected = [
            (12, 6),  # (7, 0) is junction
            (10, 4),  # (5, 0) is junction
            (13, 10),  # (8, 0) is junction
            (15, 11),  # (9, 0) is junction
        ]
        for other, exp in enumerate(expected, 1):
            end, start, heuristic = Sample.overlap_indices(sample[0], sample[other])
            self.assertTrue(heuristic)
            self.assertEqual((end, start), exp)
            self.assertEqual(pos[0][exp[0]], pos[other][exp[1]])

    def test_chunks(self):
        chunks = self.samples[0].chunks(chunk_len=4, overlap=2)
        slices = [slice(i, i + 4) for i in range(0, len(self.samples[0]) - 4, 2)]
        sliced = [self.samples[0].slice(sl) for sl in slices]
        for expt, got in zip(sliced, chunks):
            self.assertEqual(got, expt)

    def test_amend(self):
        for s in self.samples:
            new = s.amend()
            self.assertEqual(s, new)

            new = s.amend(ref_name='new_ref')
            self.assertEqual(new.ref_name, 'new_ref')

            new = new.amend(ref_name = s.ref_name)
            self.assertEqual(s, new)

        with self.assertRaises(KeyError):
            s.amend(fake_attr=None)

    def test_depth_filter(self):
        sample = self.samples[0]
        cases = [
            # all high depth, don't do any filtering
            (11 * [10], [slice(None, None)]),
            # single chunk, but trim off first two and last one columns
            ([2, 3, 10, 10, 12, 13, 10, 10, 10, 10, 9], [slice(2, 10)]),
            # split into multiple chunks
            ([2, 3, 10, 10, 4, 13, 10, 8, 10, 9, 1], [slice(2, 4), slice(5,7), slice(8,9)])
        ]
        for i, (depth, expt_slices) in enumerate(cases):
            sample = sample.amend(depth=np.array(depth))
            samples_filt = list(sample.depth_filter(min_depth=10))
            self.assertEqual(len(samples_filt), len(expt_slices))
            for j, (got, sl) in enumerate(zip(samples_filt, expt_slices)):
                expt = sample.slice(sl)
                self.assertEqual(got, expt, msg="Failed for case {} slice {}".format(i, j))

    def test_sample_counts_matrix_and_read_level_features(self):
        samples, _ = get_test_samples(num_train=1, num_test=1, chunk_size=10, read_level_features=True)
        sample = samples[0]
        assert sample.counts_matrix.ndim == 2
        assert sample.features.ndim == 3





class TestTrimSamples(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.ls = medaka.labels.HaploidLabelScheme()

        pos1 = np.array([(0, 0), (0, 1), (1, 0), (2, 0),
                         (2, 1), (2, 2), (3, 0), (4, 0),
                         (4, 1), (4, 2), (4, 3)],
                        dtype=[('major', int), ('minor', int)])
        pos2 = np.array([(4, 1), (4, 2), (4, 3), (4, 5),
                         (5, 0), (6, 0), (6, 1), (7, 0)],
                        dtype=[('major', int), ('minor', int)])
        pos3 = np.array([(6, 0), (6, 1), (7, 0), (7, 1)],
                        dtype=[('major', int), ('minor', int)])
        cls.samples = []
        data_dim = 10
        for pos in pos1, pos2, pos3:
            data = np.random.random_sample(
                size=data_dim*len(pos)).reshape((len(pos), data_dim))
            cls.samples.append(
                medaka.common.Sample(ref_name='contig1', features=data,
                                     ref_seq=None, labels=data, positions=pos,
                                     label_probs=data, depth=None))

        slices = [slice(0, 9),
                  slice(1, 6),
                  slice(1, None)]
        # get the expected sliced chunks
        cls.sliced = [s.slice(sl) for s, sl in zip(cls.samples, slices)]


    def test_works(self):
        # test simple case of 3 chained samples

        trimmed = list(Sample.trim_samples((s for s in self.samples)))

        for i, (expt, (got, is_last_in_contig, heuristic)) in enumerate(
            zip(self.sliced, trimmed)):
            self.assertEqual(got, expt)
            if i == len(self.sliced) - 1:
                self.assertTrue(is_last_in_contig)
            else:
                self.assertFalse(is_last_in_contig)


    def test_gapped(self):
        # check we get is_last_in_contig if we have a gap between samples
        trimmed = list(Sample.trim_samples(
            (s for s in self.samples[::2])))
        for i, (expt, (got, is_last_in_contig, heuristic)) in enumerate(
            zip(self.samples[::2], trimmed)):
            self.assertEqual(got, expt)
            self.assertTrue(is_last_in_contig)


    def test_raises(self):
        # test we get an exception if we e.g. provide samples from different
        # refs or samples out of order
        self.assertRaises(OverlapException,
            lambda x: list(Sample.trim_samples(x)),
                iter(self.samples[::-1]))


    def test_single_sample(self):
        # test that if we provide a single sample, we get the same sample back

        results = list(Sample.trim_samples(iter([self.samples[0]])))
        self.assertEqual(len(results), 1)
        got, is_last_in_contig, heuristic = results[0]
        self.assertEqual(got, self.samples[0])
        self.assertTrue(is_last_in_contig)
