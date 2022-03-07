import collections
import itertools
import os
import tempfile
import unittest

import numpy as np
import pysam

import medaka.common
import medaka.datastore
import medaka.labels
import medaka.variant
import medaka.vcf
from medaka.test.test_labels import haploid_sample_from_labels


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
        heuristics = len(is_last_in_contig) * [False]
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
            joined = medaka.variant.join_samples(
                zip(inp_sliced, is_last_in_contig, heuristics),
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
        heuristics = len(is_last_in_contig) * [False]

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
            joined = medaka.variant.join_samples(
                zip(inp_sliced, is_last_in_contig, heuristics),
                ref_seq.replace('*', ''), self.ls)

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
        heuristics = len(is_last_in_contig) * [False]

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
            joined = medaka.variant.join_samples(
                zip(inp_sliced, is_last_in_contig, heuristics),
                ref_seq.replace('*', ''), self.ls)

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
                s = medaka.common.Sample(contig, None, None, None, pos, None,
                                         None)
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


class TestFasta2VCF(unittest.TestCase):

    @staticmethod
    def write_sam(query_seq, ref_seq, rstart, cigar, query_name, ref_name, out_sam_fp):

        header = {
            'HD': {'VN': 1.0},
            'SQ': [{
                'LN': len(ref_seq),
                'SN': ref_name,
            }]
        }
        with pysam.AlignmentFile(out_sam_fp, 'w', header=header) as fh:
            a = pysam.AlignedSegment()
            a.reference_id = 0
            a.query_name = query_name
            a.query_sequence = query_seq
            a.reference_start = rstart
            a.cigarstring = cigar
            a.flag = 0  # not rev_comp
            a.mapping_quality = 60
            a.set_tags([('NM', 1, "i")])  # just so we don't skip this alignment
            fh.write(a)

    def test_decode_variants(self):

        cases = [
            # snp at start
            ('CATG',            # ref_seq
             'tATG',            # call_seq
             ((0, 'C', 'T'),),  # pos, ref, alt
             0, '4M'),          # reference_start, cigar str
            # snp not at start
            ('CATG',
             'CAcG',
             ((2, 'T', 'C',),),
             0, '4M'),
            # snp at end
            ('CATG',
             'CATa',
             ((3, 'G', 'A'),),
             0,'4M'),
            # mnp
            ('CATG',
             'CtgG',
             ((1, 'AT', 'TG',),),
             0,'4M'),
            # sni at start of ref (after first base)
            ('C*ATG',
             'CGATG',
             ((0, 'C', 'CG'),),
             0,'1M1I3M'),
            # sni at start of ref (before first base)
            ('*ATG',
             'GATG',
             ((0, 'A', 'GA'),),
             0,'1I3M'),
            # mni at end of ref
            ('CATG**',
             'CATGGT',
             ((3, 'G', 'GGT',),),
             0,'4M2I'),
            # snd at start of ref
            ('CATG',
             '*ATG',
             ((0, 'CA', 'A'),),
             0,'1D3M'),
            # snd at end of ref
            ('CATG',
             'CAT*',
             ((2, 'TG', 'T'),),
             0,'3M1D'),
            # mnd at start of ref
            ('CATG',
             '**TG',
             ((0, 'CAT', 'T'),),
             0,'2D2M'),
            # mnd at end of ref
            ('CATG',
             'CA**',
             ((1, 'ATG', 'A'),),
             0,'2M2D'),
            # snp and ins
            ('CA*TG',
             'CgcTG',
             ((1, 'A', 'GC'),),
             0, '2M1I2M'),
            # snp and del
            ('CATG',
             'Cg*G',
             ((1, 'AT', 'G'),),
             0, '2M1D1M'),
            # No variant
            ('CATG',
             'CATG',
             (),
             0, '4M'),
            # snp at start, unaligned ref at start
            ('AGCATG',
             '  tATG',
             ((2, 'C', 'T'),),
             2, '4M'),
            # snp not at start, unaligned query at start
            ('  CATG',
             'AACAcG',
             ((2, 'T', 'C'),),
             0, '2S4M'),
            # snp at end, unaligned query at end
            ('CATG  ',
             'CATaAT',
             ((3, 'G', 'A'),),
             0,'4M2S'),
            # No variant
            ('AGCATGAA',
             'TACATGCC',
             (),
             2, '2S4M2S'),
            # more than one variant
            ('CATGTG',
             'CAcGaG',
             ((2, 'T', 'C'), (4, 'T', 'A')),
             0, '6M'),
        ]

        ref_name = 'ref_contig'
        query_name = 'consensus_contig'

        for i, (ref_seq, query_seq, pos_ref_alts, rstart, cigar) in enumerate(cases):
            msg = 'Failed case {}, r {} q {}'.format(i, ref_seq, query_seq)
            ref_seq = ref_seq.replace('*', '').replace(' ', '').upper()
            query_seq = query_seq.replace('*', '').replace(' ', '').upper()
            _, tmp_sam = tempfile.mkstemp()
            self.write_sam(query_seq, ref_seq, rstart, cigar, query_name, ref_name, tmp_sam)
            vs = list(medaka.variant.yield_variants_from_aln(next(pysam.AlignmentFile(tmp_sam)), ref_seq))

            self.assertEqual(len(vs), len(pos_ref_alts), msg=msg)
            if len(pos_ref_alts) == 0:
                continue
            for v, (pos, ref, alt) in zip(vs, pos_ref_alts):
                self.assertEqual(v.chrom, ref_name, msg=msg)
                self.assertEqual(v.pos, pos, msg=msg)
                self.assertEqual(v.ref, ref, msg=msg)
                self.assertEqual(v.alt, [alt], msg=msg)
                self.assertEqual(v.genotype_data['GT'], '1', msg=msg)

            os.remove(tmp_sam)


    @staticmethod
    def grouper(iterable, n, fillvalue=None):
        args = [iter(iterable)] * n
        return itertools.zip_longest(*args, fillvalue=fillvalue)

    @staticmethod
    def mutate(chunk, chunk_start, chrom, max_sub_len=3, max_indel_len=5):
        bases = list('ATGC')
        ref = list(chunk)
        chunk = list(chunk)
        is_sub = bool(np.random.randint(0,2))
        pos = len(chunk) // 2
        if is_sub:
            sub_len = np.random.randint(1, max_sub_len + 1)
            rand_seq = []
            for ref_base in chunk[pos: pos + sub_len]:
                rand_seq.append(np.random.choice([b for b in bases if b != ref_base]))
            chunk[pos: pos + sub_len] = rand_seq

        else:  #indel
            is_del = bool(np.random.randint(0,2))
            indel_len = np.random.randint(1, max_indel_len + 1)
            if is_del:
                chunk[pos: pos + indel_len] = ''
            else:
                rand_seq = ''.join(np.random.choice(bases, size=indel_len))
                chunk[pos] = chunk[pos] + rand_seq
        v = medaka.vcf.Variant(chrom, chunk_start, ''.join(ref), ''.join(chunk), genotype_data={'GT': '1/1'}).trim()

        return ''.join(chunk), v


    def test_vcf_from_fasta(self):
        bases = list('ATGC')
        seq_len = 1000
        ref_seq = np.random.choice(bases, seq_len).tolist()
        tmpdir = tempfile.mkdtemp()
        ref_fasta = os.path.join(tmpdir, 'test_ref_seq.fasta')
        query_fasta = os.path.join(tmpdir, 'test_query_seq.fasta')
        chrom = 'contig1'
        with open(ref_fasta, 'w') as fh:
            fh.write('>{}\n{}\n'.format(chrom, ''.join(ref_seq)))

        chunk_size = 100
        variants = []
        query_seq = ''
        for chunk_num, chunk in enumerate(self.grouper(ref_seq, chunk_size, fillvalue='')):
            chunk_start = chunk_num * chunk_size
            mutated, variant = self.mutate(chunk, chunk_start, chrom)
            variants.append(variant)
            query_seq += mutated

        # chop off query at start and end to test bed creation
        chop = 10
        query_seq = query_seq[chop: -chop]
        with open(query_fasta, 'w') as fh:
            fh.write('>{}\n{}\n'.format(chrom, ''.join(query_seq)))

        MockArgs = collections.namedtuple(
            'args',  ('out_prefix', 'consensus', 'ref_fasta', 'threads', 'bam',
                      'chunk_size', 'pad', 'regions', 'mode'))

        prefix = os.path.join(tmpdir, 'fasta2vcf')
        args = MockArgs(prefix, query_fasta, ref_fasta, 1, None, 100, 100, None, 'HW')

        bed_fp = args.out_prefix + '_coverage.bed'
        gap_bed_fp = args.out_prefix + '_coverage_gaps.bed'
        vcf_fp = args.out_prefix + '.vcf'

        medaka.variant.vcf_from_fasta(args)
        got = list(medaka.vcf.VCFReader(vcf_fp).fetch())
        # there are many ways to represent variants (which is why tools like rtg
        # vcfeval exist) - hence just check we can reconstruct the mutated
        # sequence from the variants.
        query_from_vcf = medaka.variant.apply_variants(got, ref_seq)[chop: -chop]
        self.assertEqual(query_from_vcf, query_seq)

        coverage_intervals = tuple(medaka.common.yield_from_bed(bed_fp))
        self.assertEqual(coverage_intervals, ((chrom, chop, seq_len - chop),))
        gap_intervals = tuple(medaka.common.yield_from_bed(gap_bed_fp))
        self.assertEqual(gap_intervals, ((chrom, 0, chop), (chrom, seq_len - chop, seq_len)))
