import array
import unittest
from collections import namedtuple
import pickle
import tempfile

import numpy as np

import medaka.common
import medaka.labels
import medaka.rle


def mock_positions_array(ref):

    positions = []
    major, minor = 0, 0
    for l in ref:
        if l == '*':
            minor += 1
            positions.append((major, minor))
            continue
        else:
            major += 1
            minor = 0
            positions.append((major, minor))
    pos = np.array(positions, dtype=[('major', int),
        ('minor', int)])
    pos['major'] -= 1

    return pos


def haploid_sample_from_labels(ls=None,
                               ref=None,
                               pri=None,
                               sec=None,
                               pri_prob=0.6,
                               sec_prob=0.3):
    """Create medaka.common.Sample objects from a specified
    reference string and predicted sequence string for easily
    mocking variant calling scenarios.
    """

    assert len(ref) == len(pri)
    if sec is not None:
        assert len(ref) == len(sec)

    pos = mock_positions_array(ref)

    probs = np.zeros((len(pos), len(ls._decoding)))

    if sec is None:
        pri_prob = pri_prob + sec_prob


    for i, l in enumerate(pri):
        probs[i, ls._encoding[(l,)]] = pri_prob
        if sec is not None:
            assert sec[i] != pri[i]
            probs[i, ls._encoding[(sec[i],)]] = sec_prob
        #set another label to have non-zero prob
        #use the ref if is not in
        #primary or secondary
        other_inds = np.where(probs[i] == 0)[0]
        if ls._encoding[(ref[i],)] in other_inds:
            other_ind = ls._encoding[(ref[i],)]
        else:
            other_ind = other_inds[0]
        probs[i, other_ind] = 1 - np.sum(probs[i])

    s = medaka.common.Sample(ref_name='contig1', features=None,
                             labels=None, ref_seq=None,
                             positions=pos, label_probs=probs)
    return s, ref


class HaploidLabelSchemeTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ls = medaka.labels.HaploidLabelScheme()

    def test_picklable(self):

        with open('data.pickle', 'wb') as f:
            pickle.dump(self.ls, f, pickle.HIGHEST_PROTOCOL)

        with open('data.pickle', 'rb') as f:
            data = pickle.load(f)

    def test_n_element(self):
        self.assertEqual(self.ls.n_elements, 1)

    def test_alignments_start_end(self):
        """Test all alignments start and end in the same position.

        If alignments have a different start/end, `_alignments_to_labels`
        should raise ValueError.
        """
        start, end = np.random.randint(1, 1000, size=2)
        alignment = namedtuple("alignment", ('start', 'end'))

        incorrect_start = [alignment(start, end), alignment(start+1, end)]
        incorrect_end = [alignment(start, end), alignment(start, end+1)]
        incorrect_both = [alignment(start, end), alignment(start+1, end+1)]
        for incorrect in (incorrect_start, incorrect_end, incorrect_both):
            with self.assertRaises(ValueError):
                self.ls._alignments_to_labels(incorrect)

    def test_labels_to_encoded_labels(self):

        dummy = np.array((('A',), ('C',), ('C',), ('T',)))
        expected = np.array([1,2,2,4])
        np.testing.assert_equal(self.ls._labels_to_encoded_labels(dummy),
                                expected)

    def test_encoded_labels_to_training_vectors(self):

        dummy = np.array([1,2,2,4])
        expected = np.array([[1],
                             [2],
                             [2],
                             [4]])
        np.testing.assert_equal(self.ls._encoded_labels_to_training_vectors(dummy),
                                expected)

    def test_encoding(self):

        expected = {('C',): 2, ('G',): 3, ('T',): 4, ('*',): 0, ('A',): 1}
        self.assertEqual(self.ls._encoding, expected)

    def test_prob_to_snp(self):

        self.ls.secondary_threshold = 0.4
        #homozygous ref: ('C',) -> ('C','C')
        var = self.ls._prob_to_snp(np.array([0, 0, 1, 0, 0]),
            10, 'aref', 'C', return_all=True)
        self.assertEqual(var.ref, 'C')
        self.assertEqual(var.alt, ['.'])

        #homozygous alt: ('C',) -> ('A','A'))
        var = self.ls._prob_to_snp(np.array([0, 1, 0, 0, 0]),
            10, 'aref', 'C')
        self.assertEqual(var.ref, 'C')
        self.assertEqual(var.alt, ['A'])
        self.assertEqual(var.sample_dict['GT'], '1/1')

        #heterozygous double: ('C',) -> ('A','T')
        var = self.ls._prob_to_snp(np.array([0, 0.55, 0, 0, 0.45]),
            10, 'aref', 'C')
        self.assertEqual(var.ref, 'C')
        self.assertEqual(var.alt, ['A', 'T'])
        self.assertEqual(var.sample_dict['GT'], '1/2')

        #heterozygous single: ('C',) -> ('C','T')
        var = self.ls._prob_to_snp(np.array([0, 0, 0.55, 0, 0.45]),
            10, 'aref', 'C')
        self.assertEqual(var.ref, 'C')
        self.assertEqual(var.alt, ['T'])
        self.assertEqual(var.sample_dict['GT'], '0/1')

    def test_decode_consensus(self):

        s = unittest.mock.Mock()
        s.label_probs = np.array([[0,1,0,0,0],
                                  [0,0,1,0,0],
                                  [0,0,1,0,0],
                                  [1,0,0,0,0],
                                  [0,0,0,0,1],
                                  [0,0,0,1,0],
                                  [0,0,0,1,0],
                                  [1,0,0,0,0]])
        expected = 'ACCTGG'
        self.assertEqual(self.ls.decode_consensus(s), expected)

    def test_decode_variants(self):

        # test homozygous snp, mnp, sni, mni, snd, mnd and
        # complex mixtures thereof. Ensure we can decode a snp correctly even if
        # the base was not emitted at a major position (i.e. the call was
        # deletion followed by inserted snp).
        # also check we don't get a variant if we have deletion followed by
        # insertion of the deleted base.
        # check also we get the correct decoding if we have a del at the
        # start/end of the reference

        cases = [
            # snp
            ('CATG',          # ref_seq
             'tATG'.upper(),  # call_seq
             0, 'C', 'T'),    # pos, ref, alt
            # snp where SNP is called from deletion followed by insertion
            ('CAT*G',
             'CA*cG'.upper(),
             2, 'T', 'C'),
            # Check we don't get a variant from a deletion followed by insertion
            # of the deleted base
            ('CAT*G',
             'CA*tG'.upper(),
             None, None, None),
            # mnp
            ('CATG',
             'CtgG'.upper(),
             1, 'AT', 'TG'),
            # sni at start of ref
            ('C*ATG',
             'CGATG'.upper(),
             0, 'C', 'CG'),
            # mni at end of ref
            ('CATG**',
             'CATGGT'.upper(),
             3, 'G', 'GGT'),
            # snd at start of ref
            ('CATG',
             '*ATG'.upper(),
             0, 'CA', 'A'),
            # snd at end of ref
            ('CATG',
             'CAT*'.upper(),
             2, 'TG', 'T'),
            # mnd at start of ref
            ('CATG',
             '**TG'.upper(),
             0, 'CAT', 'T'),
            # snd at end of ref
            ('CATG',
             'CA**'.upper(),
             1, 'ATG', 'A'),
            # snp and ins
            ('CA*TG',
             'CgcTG'.upper(),
             1, 'A', 'GC'),
            # snp and ins
            ('CATG',
             'Cg*G'.upper(),
             1, 'AT', 'G'),
            # snp and ins and del
            ('CA*TG',
             'Cgc*G'.upper(),
             1, 'AT', 'GC'),
        ]

        pri_prob = 0.9
        pri_qual = self.ls._phred(1 - pri_prob)
        ref_qual = self.ls._phred(pri_prob)
        qual = pri_qual - ref_qual

        for i, (ref_seq, call, pos, ref, alt) in enumerate(cases):
            msg = 'Failed case {}'.format(i)
            s, ref_seq = haploid_sample_from_labels(self.ls, ref_seq, call,
                pri_prob=pri_prob, sec_prob=0)
            v = self.ls.decode_variants(s, ref_seq=ref_seq.replace('*', ''))

            if pos is None:
                self.assertEqual(len(v), 0)
                continue

            v = v[0]

            # total log likelihood is proportional to number of columns which
            # contribute to the variant call
            n_cols_diff = sum(p != r for p, r in zip(v.info['pred_seq'], v.info['ref_seq']))
            expected_qual = n_cols_diff * qual

            self.assertEqual(v.chrom, s.ref_name, msg=msg)
            self.assertEqual(v.pos, pos, msg=msg)
            self.assertEqual(v.ref, ref, msg=msg)
            self.assertEqual(v.alt, [alt], msg=msg)
            self.assertEqual(v.sample_dict['GT'], '1/1', msg=msg)
            self.assertAlmostEqual(float(v.qual), expected_qual, places=3, msg=msg)
            self.assertAlmostEqual(float(v.sample_dict['GQ']), expected_qual, places=3, msg=msg)


    def test_decode_snps(self):

        # test homozygous and heterozygous where secondary is and is not
        # deletion

        # primary diff, secondary same as ref => heterozygous ref
        # primary diff, secondary diff => heterozygous non-ref
        # primary diff, secondary deletion => homozygous
        # primary diff, secondary below threshold => homozygous
        # primary del, secondary diff  => not a SNP
        # minor position => not a SNP

        ref = 'CATGCGTCGATGCAT*G'
        pri = 'gAgGTGatacT*CATCG'.upper()
        sec = 'Cca***T*c**a**c**'.upper()

        pri_prob = 0.6
        sec_prob = 0.3

        s, ref_seq = haploid_sample_from_labels(self.ls, ref,
            pri, sec, pri_prob, sec_prob)

        qual_homo = self.ls._phred(1 - pri_prob)
        qual_hetero = self.ls._phred(1 - (pri_prob + sec_prob))

        pos_ref_alt_gt = [
            (0, 'C', ['G'], '0/1'),
            (1, 'A', ['C'], '0/1'),
            (2, 'T', ['G', 'A'], '1/2'),
            (4, 'C', ['T'], '1/1'),
            (6, 'T', ['A'], '0/1'),
            (7, 'C', ['T'], '1/1'),
            (8, 'G', ['A', 'C'], '1/2'),
            (9, 'A', ['C'], '1/1'),
            (14, 'T', ['C'], '0/1'),
        ]

        # find loci where the primary differs from the ref
        pos_ref_alt_gt_prim = []
        i = 0
        for r, p in zip(ref, pri):
            if r != p and '*' not in {r, p}:
                pos_ref_alt_gt_prim.append((i, r, [p], '1/1'))
            if ref != '*':
                i += 1

        self.ls.ref_seq = ref.replace('*','')
        self.ls.ref_vcf = None

        # Using a threshold equal to the secondary_prob we should get heterozygous calls
        threshold = sec_prob
        variants = self.ls.decode_snps(s, ref.replace('*', ''),
                                       threshold=threshold)
        variants = sorted(variants, key=lambda x: x.pos)


        # If we increase the threshold, we should only get only homozyous calls
        threshold = 2 * sec_prob
        variants_prim = self.ls.decode_snps(s, ref.replace('*', ''),
                                            threshold = threshold)
        variants_prim = sorted(variants_prim, key=lambda x: x.pos)
        self.assertEqual(len(variants), len(pos_ref_alt_gt))
        self.assertEqual(len(variants_prim), len(pos_ref_alt_gt_prim))

        for v, (pos, ref, alt, gt) in zip(variants + variants_prim, pos_ref_alt_gt + pos_ref_alt_gt_prim):
            self.assertEqual(v.chrom, s.ref_name)
            self.assertEqual(v.pos, pos)
            self.assertEqual(v.ref, ref)
            self.assertEqual(v.alt, alt)
            self.assertEqual(v.sample_dict['GT'], gt)
            gq = qual_homo if len(set(v.gt)) == 1 else qual_hetero
            self.assertAlmostEqual(v.sample_dict['GQ'], gq)


def diploid_sample_from_labels(ls=None,
                               ref=None,
                               hp1=None,
                               hp2=None):
    """Create medaka.common.Sample objects from a specified
    reference string and haplotype strings for easily
    mocking variant calling scenarios.
    """

    assert len(ref) == len(hp1) == len(hp2)

    pos = mock_positions_array(ref)

    # mocking up the network output in terms of 0. and 1.
    # in reality they would be float probabilities

    probs = np.zeros((len(pos), len(ls._decoding)))

    for i in range(len(ref)):

        diploid_label = tuple(sorted((hp1[i], hp2[i])))
        probs[i, ls._encoding[diploid_label]] = 1

    s = medaka.common.Sample(ref_name='contig1', features=None,
                             labels=None, ref_seq=None,
                             positions=pos, label_probs=probs)
    return s, ref


class DiploidLabelSchemeTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ls = medaka.labels.DiploidLabelScheme()

    def test_picklable(self):

        with open('data.pickle', 'wb') as f:
            pickle.dump(self.ls, f, pickle.HIGHEST_PROTOCOL)

        with open('data.pickle', 'rb') as f:
            data = pickle.load(f)

    def test_n_elements(self):

        self.assertEqual(self.ls.n_elements, 2)

    def test_labels_to_encoded_labels(self):

        #expected output from truth alignment to labels
        dummy = np.array((('A','A'), ('C','G'), ('G','C'), ('T','T')))
        expected = np.array([5, 10, 10, 14])
        np.testing.assert_equal(self.ls._labels_to_encoded_labels(dummy), expected)

    def test_encoded_labels_to_training_vectors(self):
        dummy = np.array([5, 10, 10, 14])
        expected = np.array([[5],
                             [10],
                             [10],
                             [14]])
        np.testing.assert_equal(self.ls._encoded_labels_to_training_vectors(dummy),
                                expected)

    def test_encoding(self):

        expected = {('A', 'C'): 6, ('G', 'G'): 12, ('*', '*'): 0,
                    ('C', 'C'): 9, ('G', 'T'): 13, ('*', 'A'): 1,
                    ('T', 'T'): 14, ('A', 'T'): 8, ('*', 'G'): 3,
                    ('A', 'A'): 5, ('C', 'T'): 11, ('C', 'G'): 10,
                    ('*', 'C'): 2, ('*', 'T'): 4, ('A', 'G'): 7}
        self.assertEqual(self.ls._encoding, expected)

    def test_prob_to_snp(self):

        #homozygous ref: ('C',) -> ('C','C')
        var = self.ls._prob_to_snp(np.array([0,0,0,0,0,0,0,0,0,1,0,0,0,0,0]),
            10, 'aref', 'C', return_all=True)
        self.assertEqual(var.ref, 'C')
        self.assertEqual(var.alt, ['.'])

        #homozygous alt: ('C',) -> ('A','A'))
        var = self.ls._prob_to_snp(np.array([0,0,0,0,0,1,0,0,0,0,0,0,0,0,0]),
            10, 'aref', 'C')
        self.assertEqual(var.ref, 'C')
        self.assertEqual(var.alt, ['A'])
        self.assertEqual(var.sample_dict['GT'], '1/1')

        #heterozygous double: ('C',) -> ('A','T')
        var = self.ls._prob_to_snp(np.array([0,0,0,0,0,0,0,0,1,0,0,0,0,0,0]),
            10, 'aref', 'C')
        self.assertEqual(var.ref, 'C')
        self.assertEqual(var.alt, ['A', 'T'])
        self.assertEqual(var.sample_dict['GT'], '1/2')

        #heterozygous single: ('C',) -> ('C','T')
        var = self.ls._prob_to_snp(np.array([0,0,0,0,0,0,0,0,0,0,0,1,0,0,0]),
            10, 'aref', 'C')
        self.assertEqual(var.ref, 'C')
        self.assertEqual(var.alt, ['T'])
        self.assertEqual(var.sample_dict['GT'], '0/1')

    def test_decode_snps(self):

        # minor position => not a SNP

        ref = 'CATGCGTCGATGCAT*G'
        hp1 = 'gAgGTGatacT*CATCG'.upper()
        hp2 = 'Cca***T*c**a**c**'.upper()

        s, ref_seq = diploid_sample_from_labels(self.ls, ref, hp1, hp2)

        primary_prob = 1
        qual = self.ls._phred(1 - primary_prob)

        pos_ref_alt_gt = [
            (0, 'C', ['G'], '0/1'),
            (1, 'A', ['C'], '0/1'),
            (2, 'T', ['A', 'G'], '1/2'),
            (4, 'C', ['T'], '1/1'),
            (6, 'T', ['A'], '0/1'),
            (7, 'C', ['T'], '1/1'),
            (8, 'G', ['A', 'C'], '1/2'),
            (9, 'A', ['C'], '1/1'),
            (11, 'G', ['A'], '1/1'),
            (14, 'T', ['C'], '0/1'),
        ]

        self.ls.ref_seq = ref.replace('*','')
        self.ls.ref_vcf = None

        # Using a threshold equal to the secondary_prob we should get heterozygous calls
        variants = self.ls.decode_snps(s, ref.replace('*', ''))
        variants = sorted(variants, key=lambda x: x.pos)


        self.assertEqual(len(variants), len(pos_ref_alt_gt))

        for v, (pos, ref, alt, gt) in zip(variants, pos_ref_alt_gt):
            self.assertEqual(v.chrom, s.ref_name)
            self.assertEqual(v.pos, pos)
            self.assertEqual(v.ref, ref)
            self.assertEqual(v.alt, alt)
            self.assertEqual(v.sample_dict['GT'], gt)
            self.assertAlmostEqual(v.sample_dict['GQ'], qual)


def diploid_zygosity_sample_from_labels(ls=None,
                                        ref=None,
                                        pri=None,
                                        sec=None,
                                        pri_prob=None,
                                        sec_prob=None,
                                        het=None):
    """Create medaka.common.Sample objects from a specified
    reference string and haplotype strings for easily
    mocking variant calling scenarios.
    """

    assert len(ref) == len(pri) == len(sec)

    pos = mock_positions_array(ref)

    probs = np.zeros((len(pos), len(ls._decoding)))

    for i in range(len(ref)):

        probs[i, ls._unitary_encoding[(pri[i],)]] = pri_prob
        probs[i, ls._unitary_encoding[(sec[i],)]] = sec_prob

        # het is always true
        probs[i,-1] = int(het[i])

    s = medaka.common.Sample(ref_name='contig1', features=None,
                             labels=None, ref_seq=None,
                             positions=pos, label_probs=probs)
    return s, ref


class DiploidZygosityLabelSchemeTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ls = medaka.labels.DiploidZygosityLabelScheme()

    def test_picklable(self):

        with open('data.pickle', 'wb') as f:
            pickle.dump(self.ls, f, pickle.HIGHEST_PROTOCOL)

        with open('data.pickle', 'rb') as f:
            data = pickle.load(f)

    def test_n_elements(self):
        self.assertEqual(self.ls.n_elements, 2)

    def test_labels_to_encoded_labels(self):

        #expected output from truth alignment to labels
        dummy = np.array((('A','A'), ('C','G'), ('G','C'), ('T','T')))
        expected = np.array((((1, 1), (0,)),
                             ((2, 3), (1,)),
                             ((2, 3), (1,)),
                             ((4, 4), (0,))), dtype=object)
        np.testing.assert_equal(self.ls._labels_to_encoded_labels(dummy), expected)

    def test_encoded_labels_to_training_vectors(self):

        dummy = np.array((((1, 1), (0,)),
                          ((2, 3), (1,)),
                          ((2, 3), (1,)),
                          ((4, 4), (0,))), dtype=object)

        expected = np.array([[0,1,0,0,0,0],
                             [0,0,1,1,0,1],
                             [0,0,1,1,0,1],
                             [0,0,0,0,1,0]])
        np.testing.assert_equal(self.ls._encoded_labels_to_training_vectors(dummy),
                                expected)

    def test_encoding(self):

        expected = {('A', 'C'): ((1, 2), (1,)), ('G', 'G'): ((3, 3), (0,)),
                    ('*', '*'): ((0, 0), (0,)), ('C', 'C'): ((2, 2), (0,)),
                    ('G', 'T'): ((3, 4), (1,)), ('*', 'A'): ((0, 1), (1,)),
                    ('T', 'T'): ((4, 4), (0,)), ('A', 'T'): ((1, 4), (1,)),
                    ('*', 'G'): ((0, 3), (1,)), ('A', 'A'): ((1, 1), (0,)),
                    ('C', 'T'): ((2, 4), (1,)), ('C', 'G'): ((2, 3), (1,)),
                    ('*', 'C'): ((0, 2), (1,)), ('*', 'T'): ((0, 4), (1,)),
                    ('A', 'G'): ((1, 3), (1,))}
        self.assertEqual(self.ls._encoding, expected)

    def test_prob_to_snp(self):

        #homozygous ref: ('C',) -> ('C','C')
        var = self.ls._prob_to_snp(np.array([0, 0, 1, 0, 0, 0]),
            10, 'aref', 'C', return_all=True)
        self.assertEqual(var.ref, 'C')
        self.assertEqual(var.alt, ['.'])

        #homozygous alt: ('C',) -> ('A','A'))
        var = self.ls._prob_to_snp(np.array([0, 1, 0, 0, 0, 0]),
            10, 'aref', 'C')
        self.assertEqual(var.ref, 'C')
        self.assertEqual(var.alt, ['A'])
        self.assertEqual(var.sample_dict['GT'], '1/1')

        #heterozygous double: ('C',) -> ('A','T')
        var = self.ls._prob_to_snp(np.array([0, 0.8, 0, 0, 0.6, 1]),
            10, 'aref', 'C')
        self.assertEqual(var.ref, 'C')
        self.assertEqual(var.alt, ['A', 'T'])
        self.assertEqual(var.sample_dict['GT'], '1/2')

        #heterozygous single: ('C',) -> ('C','T')
        var = self.ls._prob_to_snp(np.array([0, 0, 0.8, 0, 0.6, 1]),
            10, 'aref', 'C')
        self.assertEqual(var.ref, 'C')
        self.assertEqual(var.alt, ['T'])
        self.assertEqual(var.sample_dict['GT'], '0/1')

    def test_decode_snps_is_het_always_true(self):

        # minor position => not a SNP

        ref = 'CATGCGTCGATGCAT*G'
        pri = 'gAgGTGatacT*CATCG'.upper()
        sec = 'Cca***T*c**a**c**'.upper()
        # het is always true, and consistent with base probs
        het = '11111111111111111'

        pri_prob = 0.6
        sec_prob = 0.5

        s, ref_seq = diploid_zygosity_sample_from_labels(self.ls,
            ref, pri, sec, pri_prob, sec_prob, het)

        pos_ref_alt_gt = [
            (0, 'C', ['G'], '0/1'),
            (1, 'A', ['C'], '0/1'),
            (2, 'T', ['G', 'A'], '1/2'),
            (4, 'C', ['T'], '1/1'),
            (6, 'T', ['A'], '0/1'),
            (7, 'C', ['T'], '1/1'),
            (8, 'G', ['A', 'C'], '1/2'),
            (9, 'A', ['C'], '1/1'),
            (11, 'G', ['A'], '1/1'),
            (14, 'T', ['C'], '0/1'),
        ]

        self.ls.ref_seq = ref.replace('*','')
        self.ls.ref_vcf = None

        # Using a threshold equal to the secondary_prob we should get heterozygous calls
        variants = self.ls.decode_snps(s, ref.replace('*', ''))
        variants = sorted(variants, key=lambda x: x.pos)

        qual_hom = self.ls._phred(1 - pri_prob)
        qual_het = self.ls._phred(1 - 0.5 * (pri_prob + sec_prob))

        self.assertEqual(len(variants), len(pos_ref_alt_gt))

        for v, (pos, ref, alt, gt) in zip(variants, pos_ref_alt_gt):
            self.assertEqual(v.chrom, s.ref_name)
            self.assertEqual(v.pos, pos)
            self.assertEqual(v.ref, ref)
            self.assertEqual(v.alt, alt)
            self.assertEqual(v.sample_dict['GT'], gt)
            gq = qual_hom if len(set(v.gt)) == 1 else qual_het
            self.assertAlmostEqual(v.sample_dict['GQ'], gq)

    def test_decode_snps_is_het_always_off(self):

        # minor position => not a SNP

        ref = 'CATGCGTCGATGCAT*G'
        pri = 'gAgGTGatacT*CATCG'.upper()
        sec = 'Cca***T*c**a**c**'.upper()
        # het is always true, and consistent with base probs
        het = '00000000000000000'

        pri_prob = 0.6
        sec_prob = 0.5

        s, ref_seq = diploid_zygosity_sample_from_labels(self.ls,
            ref, pri, sec, pri_prob, sec_prob, het)

        pos_ref_alt_gt = [
            (0, 'C', ['G'], '1/1'),
            (2, 'T', ['G'], '1/1'),
            (4, 'C', ['T'], '1/1'),
            (6, 'T', ['A'], '1/1'),
            (7, 'C', ['T'], '1/1'),
            (8, 'G', ['A'], '1/1'),
            (9, 'A', ['C'], '1/1')
        ]

        self.ls.ref_seq = ref.replace('*','')
        self.ls.ref_vcf = None

        # Using a threshold equal to the secondary_prob we should get heterozygous calls
        variants = self.ls.decode_snps(s, ref.replace('*', ''))

        variants = sorted(variants, key=lambda x: x.pos)

        qual_hom = self.ls._phred(1 - pri_prob)
        qual_het = self.ls._phred(1 - 0.5 * (pri_prob + sec_prob))

        self.assertEqual(len(variants), len(pos_ref_alt_gt))

        for v, (pos, ref, alt, gt) in zip(variants, pos_ref_alt_gt):
            self.assertEqual(v.chrom, s.ref_name)
            self.assertEqual(v.pos, pos)
            self.assertEqual(v.ref, ref)
            self.assertEqual(v.alt, alt)
            self.assertEqual(v.sample_dict['GT'], gt)
            gq = qual_hom if len(set(v.gt)) == 1 else qual_het
            self.assertAlmostEqual(v.sample_dict['GQ'], gq)


class RLELabelSchemeTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.scheme = medaka.labels.RLELabelScheme(max_run=3)

    def test_encoding(self):
        """Check some elements of the `_encoding`."""
        expected = {
            ('*', 1): 0, ('A', 1): 1, ('A', 2): 2, ('A', 3): 3,
            ('C', 1): 4, ('C', 2): 5, ('C', 3): 6, ('G', 1): 7,
            ('G', 2): 8, ('G', 3): 9, ('T', 1): 10, ('T', 2): 11,
            ('T', 3): 12}
        encoding = self.scheme._encoding
        self.assertEqual(encoding, expected)

    def test_alignment_to_pairs_001(self):
        """Check output of alignment agrees.

        To fully test the functionality, the alignment contains:
            - insertions
            - deletions
            - a run length larger than `self.max_run`, that will be capped
                 to the max_run
        """
        query_name = 'query'
        reference_id = 1
        reference_start = 10
        query_sequence = 'ACATGATGTAC'
        cigarstring = '3=1I2=1D5='
        flag = 0
        qualities = array.array('B', [2, 1, 4, 5, 1, 1, 2, 16, 2, 3, 4])
        aln = medaka.rle.initialise_alignment(
            query_name, reference_id, reference_start, query_sequence,
            cigarstring, flag, query_qualities=qualities)
        expected = (
            (10, ('A', 2)), (11, ('C', 1)), (12, ('A', 3)), (None, ('T', 3)),
            (13, ('G', 1)), (14, ('A', 1)), (15, ('*', 1)), (16, ('T', 2)),
            (17, ('G', 3)), (18, ('T', 2)), (19, ('A', 3)), (20, ('C', 3)))

        got = tuple(self.scheme._alignment_to_pairs(aln))
        self.assertEqual(got, expected)

    def test_decode_consensus(self):
        """Test the conversion between network outputs and sequence"""
        num_classes = 13  # 3 elements per base * 4 bases + *
        label_probs = np.zeros([6, num_classes])
        label_probs[0, 10] = 0.9   # decodes to (T, 1)
        label_probs[1, 5] = 0.8    # (C, 2)
        label_probs[2, 0] = 0.81   # (*, 1)
        label_probs[3, 3] = 0.95   # (A, 3)
        label_probs[4, 8] = 0.9    # (G, 2)
        label_probs[5, 5] = 0.9    # (C, 2)
        mock = medaka.common.Sample(None, None, None, None, None, label_probs)
        expected = 'TCCAAAGGCC'
        got = self.scheme.decode_consensus(mock)
        self.assertEqual(expected, got)

if __name__ == '__main__':
    unittest.main()
