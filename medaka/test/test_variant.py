import unittest

import numpy as np

import medaka.common
from medaka.variant import yield_trimmed_consensus_chunks, join_chunked_variants
from medaka.variant import SNPDecoder, HaploidVariantDecoder


_label_decoding_ = list(medaka.common._alphabet_ + medaka.common._gap_)
_label_encoding_ = {b:i for i, b in enumerate(_label_decoding_)}
_feature_decoding_ = (
    ('', True, 'A', 1),
    ('', True, 'C', 1),
    ('', True, 'G', 1),
    ('', True, 'T', 1),
    ('', False, 'A', 1),
    ('', False, 'C', 1),
    ('', False, 'G', 1),
    ('', False, 'T', 1),
    ('', True, None, 1),
    ('', False, None, 1)
)

_meta_ = {
    'medaka_label_decoding': _label_decoding_,
    'medaka_feature_decoding': _feature_decoding_,
}


def sample_from_labels(ref_labels='ATCCT*GC', pred_labels='A*CGTTGC',
                       secondary_labels=None, primary_prob=0.6, secondary_prob=0.3):

    assert len(ref_labels) == len(pred_labels)
    if secondary_labels is not None:
        assert len(secondary_labels) == len(ref_labels)

    positions = []
    major, minor = 0, 0
    for l in ref_labels:
        if l == medaka.common._gap_:
            minor += 1
            positions.append((major, minor))
            continue
        else:
            major += 1
            minor = 0
            positions.append((major, minor))
    pos = np.array(positions, dtype=[('major', int), ('minor', int)])
    pos['major'] -= 1
    probs = np.zeros((len(pos), len(_label_decoding_)))

    if secondary_labels is None:
        primary_prob = primary_prob + secondary_prob

    for i, l in enumerate(pred_labels):
        probs[i, _label_encoding_[l]] = primary_prob
        if secondary_labels is not None:
            assert secondary_labels[i] != pred_labels[i]
            probs[i, _label_encoding_[secondary_labels[i]]] = secondary_prob
        # set another label to have non-zero prob - use the ref if is not in
        # primary or secondary
        other_inds = np.where(probs[i] == 0)[0]
        if _label_encoding_[ref_labels[i]] in other_inds:
            other_ind = _label_encoding_[ref_labels[i]]
        else:
            other_ind = other_inds[0]
        probs[i, other_ind] = 1 - np.sum(probs[i])

    labels = np.array([_label_encoding_[b] for b in ref_labels], dtype=int)  # with gaps
    ref_seq_encoded = labels[np.where(labels != _label_encoding_[medaka.common._gap_])]

    s = medaka.common.Sample(ref_name='contig1', features=None, labels=labels, ref_seq=None,
               positions=pos, label_probs=probs)
    return s, ref_seq_encoded


class TestYieldTrimmedConsensusChunks(unittest.TestCase):


    @classmethod
    def setUpClass(cls):
        pos1 = np.array([(0, 0), (0, 1), (1, 0), (2, 0), (2, 1), (2, 2), (3, 0), (4, 0), (4, 1), (4, 2), (4, 3)],
                        dtype=[('major', int), ('minor', int)])
        pos2 = np.array([(4, 1), (4, 2), (4, 3), (4, 5), (5, 0), (6, 0), (6, 1), (7, 0)],
                        dtype=[('major', int), ('minor', int)])
        pos3 = np.array([(6, 0), (6, 1), (7, 0), (7, 1)],
                        dtype=[('major', int), ('minor', int)])
        cls.samples = []
        data_dim = 10
        for pos in pos1, pos2, pos3:
            data = np.random.random_sample(size=data_dim*len(pos)).reshape((len(pos), data_dim))
            cls.samples.append(
                medaka.common.Sample(ref_name='contig1', features=data, ref_seq=None, labels=data, positions=pos, label_probs=data)
            )

        slices = [
            slice(0, 9),
            slice(1, 6),
            slice(1, None),
        ]
        # get the expected sliced chunks
        cls.sliced = [s.slice(sl) for s, sl in zip(cls.samples, slices)]


    def test_works(self):
        # test simple case of 3 chained samples

        trimmed = list(yield_trimmed_consensus_chunks((s for s in self.samples)))

        for i, (expt, (got, is_last_in_contig)) in enumerate(zip(self.sliced, trimmed)):
            self.assertEqual(got, expt)
            if i == len(self.sliced) - 1:
                self.assertTrue(is_last_in_contig)
            else:
                self.assertFalse(is_last_in_contig)


    def test_gapped(self):
        # check we get is_last_in_contig if we have a gap between samples
        trimmed = list(yield_trimmed_consensus_chunks((s for s in self.samples[::2])))
        for i, (expt, (got, is_last_in_contig)) in enumerate(zip(self.samples[::2], trimmed)):
            self.assertEqual(got, expt)
            self.assertTrue(is_last_in_contig)


    def test_raises(self):
        # test we get an exception if we e.g. provide samples from different
        # refs or samples out of order
        self.assertRaises(RuntimeError, lambda x: list(yield_trimmed_consensus_chunks(x)), iter(self.samples[::-1]))


    def test_single_sample(self):
        # test that if we provide a single sample, we get the same sample back

        results = list(yield_trimmed_consensus_chunks(iter([self.samples[0]])))
        self.assertEqual(len(results), 1)
        got, is_last_in_contig = results[0]
        self.assertEqual(got, self.samples[0])
        self.assertTrue(is_last_in_contig)


class TestJoinChunkedVariants(unittest.TestCase):


    def test_not_spanning(self):
        # if variants in a sample do not span sample boundaries confirm we have
        # just rechunking by 1 major position

        # if variants in a sample span sample boundaries confirm we rechunk
        # samples to make sure entire variants are in a single sample

        indel_labels = 'CATGCG****TGCATCG'
        sub_labels =   'CATGCGATACTGCATCG'
        ref_labels =   'CATGCGTCGATGCATCG'
        mix_labels =   'CATGCGAT**TGCATCG'

        # join_chunked_variants rolls back to penultimate matching major position
        # since one does not know that a final major does not have minors
        # following in the next chunk

        inp_slices = [slice(0, 4), slice(4, 12), slice(12, None)]  # CATG, CG****TG, CATCG
        exp_slices = [slice(0, 3), slice(3, 11), slice(11, None)]  # CAT, GCG****T, GCATCG
        is_last_in_contig = [False, False, True]
        is_last_in_contig_raises = [False, False, False]

        refs_calls = [
            (ref_labels, indel_labels),  # del
            (ref_labels, sub_labels),  # sub
            (indel_labels, ref_labels),  # ins
            (indel_labels, mix_labels),  # ins with some gap calls
        ]

        for ref, call in refs_calls:
            sample, ref_seq_encoded = sample_from_labels(ref, call)
            inp_sliced = [sample.slice(sl) for sl in inp_slices]
            exp_sliced = [sample.slice(sl) for sl in exp_slices]
            joined = join_chunked_variants(zip(inp_sliced, is_last_in_contig), ref_seq_encoded, _label_encoding_[medaka.common._gap_])

            for expt, got in zip(exp_sliced, joined):
                self.assertEqual(got.name, expt.name)
                self.assertEqual(got, expt)

            # check we get a ValueError if is_last_in_contig is not False for
            # last chunk
            self.assertRaises(ValueError,
                              lambda x: list(join_chunked_variants(*x)),
                              [zip(inp_sliced, is_last_in_contig_raises), ref_seq_encoded, _label_encoding_[medaka.common._gap_]]
            )


    def test_spanning(self):
        # if variants in a sample span sample boundaries confirm we rechunk
        # samples to make sure entire variants are in a single sample

        indel_labels = 'CATGCG****TGCATCG'
        sub_labels =   'CATGCGATACTGCATCG'
        ref_labels =   'CATGCGTCGATGCATCG'
        mix_labels =   'CATGCGAT**TGCATCG'
        # join_chunked_variants rolls back to penultimate matching major position

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
            sample, ref_seq_encoded = sample_from_labels(ref, call)
            inp_sliced = [sample.slice(sl) for sl in inp_slices]
            exp_sliced = [sample.slice(sl) for sl in exp_slices]
            joined = join_chunked_variants(zip(inp_sliced, is_last_in_contig), ref_seq_encoded, _label_encoding_[medaka.common._gap_])

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

        inp_slices = [slice(0, 6), slice(6, 10), slice(10, None)]  # CATGCG, ****, TGCATCG
        exp_slices = [slice(0, 5), slice(5, None)]  # CATGC, G****TGCATCG

        is_last_in_contig = [False, False, True]

        refs_calls = [
            (ref_labels, indel_labels),  # del
            (ref_labels, sub_labels),  # sub
            (ref_labels, mix_labels),  # mix
            (indel_labels, ref_labels),  # ins
            (indel_labels, mix_labels),  # ins with some gap calls
        ]

        for ref, call in refs_calls:
            sample, ref_seq_encoded = sample_from_labels(ref, call)
            inp_sliced = [sample.slice(sl) for sl in inp_slices]
            exp_sliced = [sample.slice(sl) for sl in exp_slices]
            joined = list(join_chunked_variants(zip(inp_sliced, is_last_in_contig), ref_seq_encoded, _label_encoding_[medaka.common._gap_]))

            for expt, got in zip(exp_sliced, joined):
                self.assertEqual(got.name, expt.name)
                self.assertEqual(got, expt)


class TestSNPDecoder(unittest.TestCase):


    def test_decode_variants(self):

        # test homozygous and heterozygous where secondary is and is not
        # deletion

        # primary diff, secondary same as ref => heterozygous ref
        # primary diff, secondary diff => heterozygous non-ref
        # primary diff, secondary deletion => homozygous
        # primary diff, secondary below threshold => homozygous
        # primary del, secondary diff  => not a SNP
        # minor position => not a SNP

        ref_labels = 'CATGCGTCGATGCAT*G'
        pri_labels = 'gAgGTGatacT*CATCG'.upper()
        sec_labels = 'Cca***T*c**a**c**'.upper()

        primary_prob = 0.6
        secondary_prob = 0.3

        s, ref_seq_encoded = sample_from_labels(ref_labels, pri_labels, secondary_labels=sec_labels,
                               primary_prob=primary_prob, secondary_prob=secondary_prob)

        qual_homo = -10 * np.log10(1 - primary_prob)
        qual_hetero = -10 * np.log10(1 - (primary_prob + secondary_prob))

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
        for ref, prim in zip(ref_labels, pri_labels):
            if ref != prim and medaka.common._gap_ not in {ref, prim}:
                pos_ref_alt_gt_prim.append((i, ref, [prim], '1/1'))
            if ref != medaka.common._gap_:
                i += 1


        # Using a threshold equal to the secondary_prob we should get heterozygous calls
        decoder = SNPDecoder(_meta_, threshold=secondary_prob)
        variants = decoder.decode_variants(iter([s]), ref_seq=ref_labels.replace(medaka.common._gap_, ''))
        variants = sorted(variants, key=lambda x: x.pos)

        # If we increase the threshold, we should only get only homozyous calls
        decoder = SNPDecoder(_meta_, threshold=2*secondary_prob)
        variants_prim = decoder.decode_variants(iter([s]), ref_seq=ref_labels.replace(medaka.common._gap_, ''))
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


class TestHaploidVariantDecoder(unittest.TestCase):


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

        primary_prob = 0.9
        prim_qual = float(-10 * np.log10(1 - primary_prob))
        ref_qual = float(-10 * np.log10(primary_prob))
        qual = prim_qual - ref_qual
        decoder = HaploidVariantDecoder(_meta_)

        for i, (ref_seq, call, pos, ref, alt) in enumerate(cases):
            msg = 'Failed case {}'.format(i)
            s, _ = sample_from_labels(ref_seq, call, primary_prob=primary_prob, secondary_prob=0)
            v = list(decoder.decode_variants(iter([s]),
                                             ref_seq=ref_seq.replace(medaka.common._gap_, ''))
            )
            if pos is None:
                self.assertEqual(len(v), 0)
                continue
            v = v[0]
            # total log likelyhood is proportional to number of columns which
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


if __name__ == '__main__':
    unittest.main()
