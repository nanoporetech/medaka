import collections
import os
import tempfile
import unittest

import numpy as np

import medaka.common
import medaka.datastore
import medaka.labels
import medaka.variant
from medaka.test.test_labels import haploid_sample_from_labels


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
                                     label_probs=data))

        slices = [slice(0, 9),
                  slice(1, 6),
                  slice(1, None)]
        # get the expected sliced chunks
        cls.sliced = [s.slice(sl) for s, sl in zip(cls.samples, slices)]

    def test_works(self):
        # test simple case of 3 chained samples

        trimmed = list(medaka.variant.trim_samples((s for s in self.samples)))

        for i, (expt, (got, is_last_in_contig)) in enumerate(
            zip(self.sliced, trimmed)):
            self.assertEqual(got, expt)
            if i == len(self.sliced) - 1:
                self.assertTrue(is_last_in_contig)
            else:
                self.assertFalse(is_last_in_contig)


    def test_gapped(self):
        # check we get is_last_in_contig if we have a gap between samples
        trimmed = list(medaka.variant.trim_samples(
            (s for s in self.samples[::2])))
        for i, (expt, (got, is_last_in_contig)) in enumerate(
            zip(self.samples[::2], trimmed)):
            self.assertEqual(got, expt)
            self.assertTrue(is_last_in_contig)


    def test_raises(self):
        # test we get an exception if we e.g. provide samples from different
        # refs or samples out of order
        self.assertRaises(RuntimeError,
            lambda x: list(medaka.variant.trim_samples(x)),
                iter(self.samples[::-1]))


    def test_single_sample(self):
        # test that if we provide a single sample, we get the same sample back

        results = list(medaka.variant.trim_samples(iter([self.samples[0]])))
        self.assertEqual(len(results), 1)
        got, is_last_in_contig = results[0]
        self.assertEqual(got, self.samples[0])
        self.assertTrue(is_last_in_contig)


class TestJoinSamples(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ls = medaka.labels.HaploidLabelScheme()

    def test_not_spanning(self):
        # if variants in a sample do not span sample boundaries confirm we have
        # just rechunking by 1 major position

        # if variants in a sample span sample boundaries confirm we rechunk
        # samples to make sure entire variants are in a single sample

        indel_labels = 'CATGCG****TGCATCG'
        sub_labels =   'CATGCGATACTGCATCG'
        ref_labels =   'CATGCGTCGATGCATCG'
        mix_labels =   'CATGCGAT**TGCATCG'

        # join_samples rolls back to penultimate matching major position
        # since one does not know that a final major does not have minors
        # following in the next chunk

        # CATG, CG****TG, CATCG
        inp_slices = [slice(0, 4), slice(4, 12), slice(12, None)]
        # CAT, GCG****T, GCATCG
        exp_slices = [slice(0, 3), slice(3, 11), slice(11, None)]
        is_last_in_contig = [False, False, True]
        is_last_in_contig_raises = [False, False, False]

        refs_calls = [
            (ref_labels, indel_labels),  # del
            (ref_labels, sub_labels),  # sub
            (indel_labels, ref_labels),  # ins
            (indel_labels, mix_labels),  # ins with some gap calls
        ]

        for i, (ref, call) in enumerate(refs_calls):
            sample, ref_seq = haploid_sample_from_labels(self.ls, ref, call)
            inp_sliced = [sample.slice(sl) for sl in inp_slices]
            exp_sliced = [sample.slice(sl) for sl in exp_slices]
            joined = medaka.variant.join_samples(zip(inp_sliced, is_last_in_contig),
                ref_seq.replace('*', ''), self.ls)

            for j, (expt, got) in enumerate(zip(exp_sliced, joined)):
                self.assertEqual(got.name, expt.name)
                self.assertEqual(got, expt)

            # check we get a ValueError if is_last_in_contig is not False for
            # last chunk
            self.assertRaises(ValueError,
                              lambda x: list(medaka.variant.join_samples(*x)),
                              [zip(inp_sliced, is_last_in_contig_raises),
                               ref_seq.replace('*', ''),
                               self.ls])

    def test_spanning(self):
        # if variants in a sample span sample boundaries confirm we rechunk
        # samples to make sure entire variants are in a single sample

        indel_labels = 'CATGCG****TGCATCG'
        sub_labels =   'CATGCGATACTGCATCG'
        ref_labels =   'CATGCGTCGATGCATCG'
        mix_labels =   'CATGCGAT**TGCATCG'
        # join_samples rolls back to penultimate matching major position

        inp_slices = [slice(0, 8), slice(8, None)]  # CATGCG**, **TGCATCG
        exp_slices = [slice(0, 5), slice(5, None)]  # CATGC, G****TGCATCG

        is_last_in_contig = [False, True]

        refs_calls = [
            (ref_labels, indel_labels),  # del
            (ref_labels, sub_labels),  # sub
            (indel_labels, ref_labels),  # ins
            (indel_labels, mix_labels),  # ins with some gap calls
        ]

        for ref, call in refs_calls:
            sample, ref_seq = haploid_sample_from_labels(self.ls, ref, call)
            inp_sliced = [sample.slice(sl) for sl in inp_slices]
            exp_sliced = [sample.slice(sl) for sl in exp_slices]
            joined = medaka.variant.join_samples(zip(inp_sliced, is_last_in_contig),
                                  ref_seq.replace('*', ''),
                                  self.ls)

            for expt, got in zip(exp_sliced, joined):
                self.assertEqual(got.name, expt.name)
                self.assertEqual(got, expt)


    def test_no_pos_same(self):
        # check we queue up a sample if it is only minor positions or only
        # subs/dels
        indel_labels = 'CATGCG****TGCATCG'
        sub_labels =   'CATGCGATACTGCATCG'
        mix_labels =   'CATGCGAT**TGCATCG'
        ref_labels =   'CATGCGTCGATGCATCG'

        # CATGCG, ****, TGCATCG
        inp_slices = [slice(0, 6), slice(6, 10), slice(10, None)]
        # CATGC, G****TGCATCG
        exp_slices = [slice(0, 5), slice(5, None)]

        is_last_in_contig = [False, False, True]

        refs_calls = [
            (ref_labels, indel_labels),  # del
            (ref_labels, sub_labels),  # sub
            (ref_labels, mix_labels),  # mix
            (indel_labels, ref_labels),  # ins
            (indel_labels, mix_labels),  # ins with some gap calls
        ]

        for ref, call in refs_calls:
            sample, ref_seq = haploid_sample_from_labels(self.ls, ref, call)
            inp_sliced = [sample.slice(sl) for sl in inp_slices]
            exp_sliced = [sample.slice(sl) for sl in exp_slices]
            joined = list(medaka.variant.join_samples(zip(inp_sliced, is_last_in_contig),
                          ref_seq.replace('*', ''), self.ls))

            for expt, got in zip(exp_sliced, joined):
                self.assertEqual(got.name, expt.name)
                self.assertEqual(got, expt)


class TestSamplesToBed(unittest.TestCase):
    def test_samples_to_bed(self):

        samples = [
            ('chr1', 0, 1000),
            ('chr1', 900, 1900),
            ('chr1', 1900, 3000), # should be merged as overlapping
            ('chr1', 3001, 3500),  # should be merged as abutting
            ('chr1', 3502, 4000),  # should not be merged
            ('chr2', 1000, 10000), # single sample in a contig
            ('chr3', 0, 1000),     # check we merge contigs separately
            ('chr3', 2000, 2001),  # short sample
            ('chr3', 3000, 3000),  # sample which starts and ends on same major
        ]
        expected = {
            ('chr1', 0, 3501),
            ('chr1', 3502, 4001),
            ('chr2', 1000, 10001), # single sample in a contig
            ('chr3', 0, 1001),     # check we merge contigs separately
            ('chr3', 2000, 2002),  # short sample
            ('chr3', 3000, 3001),  # sample which starts and ends on same major
        }

        dtype = [('major', int), ('minor', int)]

        _, tmp_hdf = tempfile.mkstemp()
        _, tmp_bed = tempfile.mkstemp()

        with medaka.datastore.DataStore(tmp_hdf, 'w') as ds:
            for contig, start, end in samples:

                pos = np.array([(start, 0), (end, 0)], dtype=dtype)
                s = medaka.common.Sample(contig, None, None, None, pos, None)
                ds.write_sample(s)

        args = collections.namedtuple('args', 'inputs output')
        medaka.variant.samples_to_bed(args([tmp_hdf], tmp_bed))

        intervals = set()
        with open(tmp_bed) as fh:
            for line in fh:
                split = line.split('\t')
                intervals.add((split[0], int(split[1]), int(split[2])))

        os.remove(tmp_hdf)
        os.remove(tmp_bed)

        self.assertEqual(intervals, expected)

