from collections import OrderedDict, Counter, namedtuple
import os
import tempfile
import unittest

import intervaltree
import itertools
import numpy as np
np.random.seed(7)
import parasail
import pysam

from medaka.common import yield_from_bed
from medaka.vcf import (VCFWriter, VCFReader, Variant, Haploid2DiploidConverter,
                        split_variants, classify_variant, _merge_variants,
                        MetaInfo, get_padded_haplotypes, align_read_to_haps,
                        align_reads_to_haps, split_mnp)

root_dir = os.path.abspath(os.path.dirname(__file__))
test1_file = os.path.join(root_dir, 'data/test1.vcf')


class TestMetaInfo(unittest.TestCase):
    """Test medaka.vcf.MetaInfo."""

    def test_010_valid_case(self):
        ident = 'ref_prob'
        descr = 'Allele probability'
        for group in ('INFO', 'FILTER', 'FORMAT'):
            for number in ('A', 'R', 'G', '.', 1):
                for typ in ('Integer', 'Float', 'Flag', 'Character', 'String'):
                    try:
                        meta_info = MetaInfo(group, ident, number, typ, descr)
                    except Exception:
                        self.fail('Correct input is raising Exception.')

                    # Check repr() is correctly
                    expected = '{}=<ID={},Number={},Type={},'\
                          'Description="{}">'.format(
                              group, ident, number, typ, descr)

                    got = repr(meta_info)
                    self.assertEqual(expected, got)

    def test_020_invalid_group(self):
        """Any group not in MetaInfo.__valid_groups__ raises ValueError."""
        with self.assertRaises(ValueError):
            meta_info = MetaInfo(
                'nonsense', 'ref_prob', 1, 'Float',
                'Allele probability')

    def test_030_invalid_number(self):
        """Int needs to be an integer or 'A', 'R', 'G' or '.'."""
        with self.assertRaises(ValueError):
            meta_info = MetaInfo(
                'INFO', 'ref_prob', 'T', 'Float',
                'Allele probability')

    def test_040_invalid_type(self):
        """When an invalid data type is entered, ValueError is raised."""
        with self.assertRaises(ValueError):
            meta_info = MetaInfo(
                'INFO', 'ref_prob', 1, 'Double',
                'Allele probability')


class TestVariant(unittest.TestCase):
    """Test medaka.vcf.Variant."""

    @classmethod
    def setUpClass(cls):
        cls.base_parameters = {
            'chrom': 'chr1', 'pos': 14369, 'ref': 'G', 'alt': ['A'],
            'ident': 'rs6054257', 'qual': 29, 'filt': 'PASS',
            'info': {'NS': 3, 'DP': 14, 'AF': 0.5, 'DB': True, 'H2': True},
            'genotype_data': OrderedDict(
                [('GT', '1|0'), ('DP', '8'), ('GQ', '48'), ('HQ', '51,51')])}
        cls.variant = Variant(**cls.base_parameters)

    def test_010_initialisation(self):
        """Test attributes are correct after initialising an instance."""

        expected = self.base_parameters
        for key, exp in expected.items():
            got = getattr(self.variant, key)
            self.assertEqual(got, exp)

    def test_020_inequalities(self):
        """Check equality of two variants."""

        # Create an alternative parameters with all values different
        alternative_parameters = {
            'chrom': 'chr2', 'pos': 1, 'ref': 'T', 'alt': ['C'],
            'ident': 'rt', 'qual': 1, 'filt': '.',
            'info': {'NS': 1, 'DP': 1, 'AF': 0.1, 'DB': False, 'H2': False},
            'genotype_data': OrderedDict(
                [('GT', '0|0'), ('DP', '7'), ('GQ', '12'), ('HQ', '5,5')])
        }

        # A variant created with the same parameters should be equal,
        # but not the same object
        variant1 = Variant(**self.base_parameters)
        self.assertTrue(id(self.variant) != id(variant1))
        self.assertTrue(self.variant == variant1)

        # Changing one thing at a time
        for changing_key in alternative_parameters.keys():
            new_parameters = self.base_parameters.copy()
            new_parameters[changing_key] = alternative_parameters[changing_key]
            variant1 = Variant(**new_parameters)
            self.assertTrue(variant1 != self.variant)

    def test_030_genotype_info(self):
        # self.variant has genotype information, variant2 will not have any
        variant2 =  Variant(
            'chr1', 14369, 'G', alt=['A'], ident='rs6054257', qual=29)

        expected_genotype_keys = 'GT:DP:GQ:HQ'  # Resorted keys
        self.assertEqual(self.variant.genotype_keys, expected_genotype_keys)

        expected_genotype_values = '1|0:8:48:51,51'
        self.assertEqual(self.variant.genotype_values, expected_genotype_values)

        expected_gt = (1, 0)
        self.assertEqual(self.variant.gt, expected_gt)

        # If no genotype info is present, expected None
        expected_gt = None
        self.assertEqual(variant2.gt, expected_gt)


class TestReader(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        try:
            cls.vcf_reader = VCFReader(test1_file)
        except Exception as e:
            cls.fail('setUpClass raised {} unexpectedly'.format(e))


    def test_010_check_meta(self):
        expected = ('fileformat=VCFv4.0',
                    'fileDate=20090805',
                    'source=myImputationProgramV3.1',
                    'reference=1000GenomesPilot-NCBI36',
                    'phasing=partial',
                    'INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">',
                    'INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
                    'INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">',
                    'INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">',
                    'INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">',
                    'INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">',
                    'FILTER=<ID=q10,Description="Quality below 10">',
                    'FILTER=<ID=s50,Description="Less than 50% of samples have data">',
                    'FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                    'FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
                    'FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
                    'FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">')
        self.assertSequenceEqual(self.vcf_reader.meta, expected, "Meta-data incorrect.")


    def test_020_check_header(self):
        expected = ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
                    'FORMAT', 'SAMPLE')
        self.assertSequenceEqual(self.vcf_reader.header, expected, "Header incorrect.")


    def test_030_check_default_lines(self):
        expected = 5
        result = sum([1 for _ in self.vcf_reader.fetch()])
        self.assertEqual(expected, result)


    def test_040_check_restricted_lines(self):
        start, end = 15000, 1231000
        result = sum([1 for _ in self.vcf_reader.fetch(start=start, end=end)])
        expected = 3
        self.assertEqual(result, expected)


    def test_050_check_variant_contents(self):
        expected = [
                Variant('chr1', 14369, 'G', alt=['A'], ident='rs6054257', qual=29, filt='PASS',
                    info={'NS': 3, 'DP': 14, 'AF': 0.5, 'DB': True, 'H2':True},
                    genotype_data=OrderedDict([('GT', '1|0'), ('GQ', '48'), ('DP', '8'), ('HQ', '51,51')]),
                ),
                Variant('chr2', 17329, 'T', alt=['A'], ident='.', qual=3, filt='q10',
                    info={'NS': 3, 'DP': 11, 'AF': 0.017},
                    genotype_data=OrderedDict([('GT', '0|0'), ('GQ', '49'), ('DP', '3'), ('HQ', '58,50')]),
                ),
                Variant('chr10', 1110695, 'A', alt=['G', 'T'], ident='rs6040355', qual=67, filt='PASS',
                    info={'NS': 2, 'DP': 10, 'AF': [0.333,0.667], 'AA': 'T', 'DB': True},
                    genotype_data=OrderedDict([('GT', '1|2'), ('GQ', '21'), ('DP', '6'), ('HQ', '23,27')]),
                ),
                Variant('chr20', 1230236, 'T', alt=['.'], ident='.', qual=47, filt='PASS',
                    info={'NS': 3, 'DP': 13, 'AA': 'T'},
                    genotype_data=OrderedDict([('GT', '0|0'), ('GQ', '54'), ('DP', '7'), ('HQ', '56,60')]),
                ),
                Variant('chrX', 1234566, 'GTCT', alt=['G','GTACT'], ident='microsat1', qual=50, filt='PASS',
                    info={'NS': 3, 'DP': 9, 'AA': 'G'},
                    genotype_data=OrderedDict([('GT', '1/1'), ('GQ', '40'), ('DP', '3')]),
                )
        ]
        result = list(self.vcf_reader.fetch())
        self.assertSequenceEqual(result, expected)

    def test_060_using_cache_does_not_affect_results(self):
        uncached = list(VCFReader(test1_file, cache=False).fetch())
        cached = list(VCFReader(test1_file, cache=True).fetch())
        self.assertSequenceEqual(uncached, cached, 'Caching does not affect parsing.')


class TestWriter(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        try:
            cls.vcf_reader = VCFReader(test1_file)
        except Exception as e:
            cls.fail('setUpClass raised {} unexpectedly'.format(e))


    def test_010_check_output(self):
        with tempfile.NamedTemporaryFile() as test_file:
            # Write contents
            with VCFWriter(test_file.name) as vcf_writer:
                for variant in self.vcf_reader.fetch():
                    vcf_writer.write_variant(variant)


            # Read them back and compare with original file
            vcf_test = VCFReader(test_file.name).fetch()
            for variant_original, variant_test in zip(self.vcf_reader.fetch(), vcf_test):
                for key in  ('chrom', 'pos', 'ident', 'ref', 'alt', 'qual', 'filt', 'info', 'genotype_keys', 'genotype_values'):
                    expected = getattr(variant_original, key)
                    result = getattr(variant_test, key)
                    self.assertEqual(expected, result, 'Round trip failed for {}.'.format(key))


class TestTrim(unittest.TestCase):

    def test_001_check_trim_start(self):
        v_orig = Variant('20', 14369, 'GGC', alt=['GGA'])
        v_expt = Variant('20', 14371, 'C', alt=['A'])
        v_trim = v_orig.trim()
        self.assertEqual(v_expt, v_trim, 'Trimming failed for {}.'.format(v_expt))


    def test_002_check_trim_end(self):
        v_orig = Variant('20', 14369, 'CGG', alt=['AGG'])
        v_expt = Variant('20', 14369, 'C', alt=['A'])
        v_trim = v_orig.trim()
        self.assertEqual(v_expt, v_trim, 'Trimming failed for {}.'.format(v_expt))


    def test_003_check_trim_both(self):
        v_orig = Variant('20', 14369, 'ATCGG', alt=['ATAGG'])
        v_expt = Variant('20', 14371, 'C', alt=['A'])
        v_trim = v_orig.trim()
        self.assertEqual(v_expt, v_trim, 'Trimming failed for {}.'.format(v_expt))


    def test_004_check_trim_ref(self):
        v_orig = Variant('20', 14369, 'ATCGG', alt=['ATCGG'])
        v_expt = Variant('20', 14369, 'A', alt=['A'])
        v_trim = v_orig.trim()
        self.assertEqual(v_expt, v_trim, 'Trimming failed for {}.'.format(v_expt))


    def test_005_check_trim_ref(self):
        v_orig = Variant('20', 14369, 'ATCGG', alt=['ATAGG', 'ATGGG'])
        v_expt = Variant('20', 14371, 'C', alt=['A', 'G'])
        v_trim = v_orig.trim()
        self.assertEqual(v_expt, v_trim, 'Trimming failed for {}.'.format(v_expt))


    def test_006_check_trim_ref(self):
        v_orig = Variant('20', 14369, 'CCTG', alt=['C'])
        v_expt = v_orig
        v_trim = v_orig.trim()
        self.assertEqual(v_expt, v_trim, 'Trimming failed for {}.'.format(v_expt))


    def test_007_check_trim_ref(self):
        # if ref and alt are the same, make sure we don't completely remove ref
        # and alt
        v_orig = Variant('20', 14369, ref='TAGTCACAG', alt=['TCACAG'])
        v_expt = Variant('20', 14369, ref='TAGT', alt=['T'])
        v_trim = v_orig.trim()
        self.assertEqual(v_expt, v_trim, 'Trimming failed for {}.'.format(v_expt))


class TestSplitHaplotypes(unittest.TestCase):

    def test_001_check_hetero(self):
        v_orig = Variant('20', 14369, 'G', alt=['A', 'C'], qual=10,
                         genotype_data=OrderedDict([('GT', '1/2'), ('GQ', 10.0)]))
        genotype_data = v_orig.genotype_data.copy()
        genotype_data['GT'] = '1/1'
        expt = tuple([
            (1, Variant('20', 14369, 'G', alt=['A'], qual=10, genotype_data=genotype_data)),
            (2, Variant('20', 14369, 'G', alt=['C'], qual=10, genotype_data=genotype_data)),
        ])
        got = v_orig.split_haplotypes()
        self.assertEqual(expt, got, 'Splitting haplotypes failed for {}'.format(v_orig))


    def test_002_check_hetero(self):
        v_orig = Variant('20', 14369, 'G', alt=['A'], qual=10,
                         genotype_data=OrderedDict([('GT', '1/0'), ('GQ', 10.0)]))
        genotype_data = v_orig.genotype_data.copy()
        genotype_data['GT'] = '1/1'
        expt = tuple([
            (1, Variant('20', 14369, 'G', alt=['A'], qual=10, genotype_data=genotype_data)),
            (2, None),
        ])
        got = v_orig.split_haplotypes()
        self.assertEqual(expt, got, 'Splitting haplotypes failed for {}'.format(v_orig))


    def test_003_check_hetero(self):
        v_orig = Variant('20', 14369, 'G', alt=['T', 'A'], qual=10,
                         genotype_data=OrderedDict([('GT', '0/1'), ('GQ', 10.0)]))
        genotype_data = v_orig.genotype_data.copy()
        genotype_data['GT'] = '1/1'
        expt = tuple([
            (1, None),
            (2, Variant('20', 14369, 'G', alt=['T'], qual=10, genotype_data=genotype_data)),
        ])
        got = v_orig.split_haplotypes()
        self.assertEqual(expt, got, 'Splitting haplotypes failed for {}'.format(v_orig))


    def test_004_check_homo(self):
        v_orig = Variant('20', 14369, 'G', alt=['T', 'A'], qual=10,
                         genotype_data=OrderedDict([('GT', '1/1'), ('GQ', 10.0)]))
        genotype_data = v_orig.genotype_data.copy()
        genotype_data['GT'] = '1/1'
        expt = tuple([
            (1, Variant('20', 14369, 'G', alt=['T'], qual=10, genotype_data=genotype_data)),
            (2, Variant('20', 14369, 'G', alt=['T'], qual=10, genotype_data=genotype_data)),
        ])
        got = v_orig.split_haplotypes()
        self.assertEqual(expt, got, 'Splitting haplotypes failed for {}'.format(v_orig))


    def test_005_check_homo(self):
        v_orig = Variant('20', 14369, 'G', alt=['T', 'A'], qual=10,
                         genotype_data=OrderedDict([('GT', '2/2'), ('GQ', 10.0)]))
        genotype_data = v_orig.genotype_data.copy()
        genotype_data['GT'] = '1/1'
        expt = tuple([
            (1, Variant('20', 14369, 'G', alt=['A'], qual=10, genotype_data=genotype_data)),
            (2, Variant('20', 14369, 'G', alt=['A'], qual=10, genotype_data=genotype_data)),
        ])
        got = v_orig.split_haplotypes()
        self.assertEqual(expt, got, 'Splitting haplotypes failed for {}'.format(v_orig))


    def test_006_check_homo(self):
        v_orig = Variant('20', 14369, 'G', alt=['T', 'A'], qual=10,
                         genotype_data=OrderedDict([('GT', '0/0'), ('GQ', 10.0)]))
        genotype_data = v_orig.genotype_data.copy()
        expt = tuple([
            (1, None),
            (2, None),
        ])
        got = v_orig.split_haplotypes()
        self.assertEqual(expt, got, 'Splitting haplotypes failed for {}'.format(v_orig))


class TestMergeAndSplitVCFs(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.vcf_merged = os.path.join(root_dir, 'data/test_merged.vcf')
        cls.ref_fasta = os.path.join(root_dir, 'data/test_ref.fasta')
        cls.vcf1 = os.path.join(root_dir, 'data/test_hap1.vcf')
        cls.vcf2 = os.path.join(root_dir, 'data/test_hap2.vcf')
        cls.split_fps = []

    @classmethod
    def tearDownClass(cls):
        for f in cls.split_fps:
            if os.path.isfile(f):
                os.remove(f)


    @staticmethod
    def intervaltree_prep(h1, h2, ref_seq):
        """Mock up trees as would be done by the `VCFReader` class so we can test medaka.vcf._merge_variants

        :param h1, h2: iterable of variants in first and second haplotype, respectively.
        :param ref_seq: str, reference sequence

        :returns: (`intervaltree.Interval` containing interable of all variants,
                   [`intervaltree.IntervalTree` for each haplotype])
        """
        trees = []
        for variants in h1, h2:
            trees.append(intervaltree.IntervalTree())
            for v in variants:
                trees[-1].add(intervaltree.Interval(v.pos, v.pos + len(v.ref), data=v))
        only_overlapping=True
        comb_tree = intervaltree.IntervalTree(trees[0].all_intervals.union(trees[1].all_intervals))
        # if strict, merge only overlapping intervals (not adjacent ones)
        comb_tree.merge_overlaps(strict=only_overlapping, data_initializer=list(), data_reducer=lambda x,y: x + [y])
        comb_interval = list(comb_tree.all_intervals)[0]
        return comb_interval, trees


    def test_001_check_merge(self):
        converter = Haploid2DiploidConverter(self.vcf1, self.vcf2, self.ref_fasta,
                                    only_overlapping=True, discard_phase=False,
                                    detailed_info=True)
        merged = converter.variants()
        for expt, found in zip(VCFReader(self.vcf_merged, cache=False).fetch(), merged):
            for key in  ('chrom', 'pos', 'ref', 'alt', 'info_string', 'gt', 'phased'):
                expected = getattr(expt, key)
                result = getattr(found, key)
                self.assertEqual(expected, result, 'Merging failed for {}:{} {}.'.format(expt.chrom, expt.pos+1, key))
        self.assertEqual(len(converter.meta_info), 10)


    def test_002_check_merge_snps(self):
        ref_seq = 'ATGGTATGCGATTGACC'
        chrom='chrom1'

        h1 = [Variant(chrom, 0, 'A', alt='C', qual=5, genotype_data={'GT':'1/1'})]
        h2 = [Variant(chrom, 0, 'A', alt='T', qual=10, genotype_data={'GT':'1/1'})]
        expt = Variant(chrom, 0, 'A', alt=['C', 'T'], qual=7.5, genotype_data={'GT': '1|2'})

        comb_interval, trees = self.intervaltree_prep(h1, h2, ref_seq)
        # preserve phase otherwise alts could be switched around
        got = _merge_variants(comb_interval, trees, ref_seq, discard_phase=False)
#        self.assertEqual(expt, got)

        for key in  ('chrom', 'pos', 'ref', 'qual', 'alt', 'gt', 'phased'):
            expected = getattr(expt, key)
            result = getattr(got, key)
            self.assertEqual(expected, result, 'Merging failed for {}:{} {}.'.format(expt.chrom, expt.pos+1, key))


    def test_003_check_merge_indels(self):
        ref_seq = 'ATGGTATGCGATTGACC'
        chrom='chrom1'

        h1 = [Variant(chrom, 0, 'ATG', alt='G', qual=5, genotype_data={'GT':'1/1'})]
        h2 = [Variant(chrom, 1, 'T', alt='TT', qual=10, genotype_data={'GT':'1/1'})]
        expt = Variant(chrom, 0, 'ATG', alt=['G', 'ATTG'], qual=7.5, genotype_data={'GT': '1|2'})

        comb_interval, trees = self.intervaltree_prep(h1, h2, ref_seq)
        # preserve phase otherwise alts could be switched around
        got = _merge_variants(comb_interval, trees, ref_seq, discard_phase=False)
        for key in  ('chrom', 'pos', 'ref', 'qual', 'alt', 'gt', 'phased'):
            expected = getattr(expt, key)
            result = getattr(got, key)
            self.assertEqual(expected, result, 'Merging failed for {}:{} {}.'.format(expt.chrom, expt.pos+1, key))


    def test_003_check_merge_multi(self):
        ref_seq = 'ATGGTATGCGATTGACC'
        chrom='chrom1'

        h1 = [Variant(chrom, 0, 'ATG', alt='G', qual=5, genotype_data={'GT':'1/1'}),
              Variant(chrom, 4, 'T', alt='G', qual=5, genotype_data={'GT':'1/1'}),
              Variant(chrom, 7, 'G', alt='GG', qual=5, genotype_data={'GT':'1/1'}),
              Variant(chrom, 9, ref_seq[9], alt=ref_seq[9] + 'T', qual=5, genotype_data={'GT':'1/1'}),
              ]
        h2 = [Variant(chrom, 1, 'T', alt='TT', qual=10, genotype_data={'GT':'1/1'}),
              Variant(chrom, 2, ref_seq[2:10], alt=ref_seq[2], qual=10, genotype_data={'GT':'1/1'}),
              ]

        # POS  0    1   2   3   4   5   6   7   8   9   10
        # REF  A    T   G   G   T   A   T   G   C   G   A
        # H1   -    -   G   G   g   A   T   Gg  C   Gt  A
        # H2   A    Tt  G   -   -   -   -   -   -   -   A

        # expected merged variants
        ref_expt =  "ATGGTATGCG"
        alt1_expt = "GGGATGGCGT"
        alt2_expt = "ATTG"

        expt = Variant(chrom, 0, ref_expt, alt=[alt1_expt, alt2_expt], qual=7.5, genotype_data={'GT': '1|2'})

        comb_interval, trees = self.intervaltree_prep(h1, h2, ref_seq)
        # preserve phase otherwise alts could be switched around
        got = _merge_variants(comb_interval, trees, ref_seq, discard_phase=False)
        for key in  ('chrom', 'pos', 'ref', 'qual', 'alt', 'gt', 'phased'):
            expected = getattr(expt, key)
            result = getattr(got, key)
            self.assertEqual(expected, result, 'Merging failed for {}:{} {}.'.format(expt.chrom, expt.pos+1, key))

    def test_004_check_merge_multi_bug(self):
        # if we have two indels on one haplotype that cancel each other out
        # (e.g. insertion of a T followed by a deletion of a T)
        # check we don't have an alt that is the same as the ref.
        ref_seq = 'TTTTTTTTTT'
        chrom='chrom1'

        h1 = [Variant(chrom, 0, 'TTTTT', alt='T', qual=5, genotype_data={'GT':'1/1'})]
        h2 = [Variant(chrom, 1, 'T', alt='TT', qual=10, genotype_data={'GT':'1/1'}),
              Variant(chrom, 3, 'TT', alt='T', qual=10, genotype_data={'GT':'1/1'})]
        expt = Variant(chrom, 0, 'TTTTT', alt='T', qual=5, genotype_data={'GT':'1|0'})

        comb_interval, trees = self.intervaltree_prep(h1, h2, ref_seq)
        # preserve phase otherwise alts could be switched around
        got = _merge_variants(comb_interval, trees, ref_seq, discard_phase=False)
        for key in  ('chrom', 'pos', 'ref', 'qual', 'alt', 'gt', 'phased'):
            expected = getattr(expt, key)
            result = getattr(got, key)
            self.assertEqual(expected, result, 'Merging failed for {}:{} {}.'.format(expt.chrom, expt.pos+1, key))

    def test_002_check_split(self):
        self.split_fps.extend(split_variants(self.vcf_merged))

        # 1-based positions on either position to exclude
        # these are where two variants have been merged and cannot be easily
        # separated without alignment, or where an indel has been rolled
        # forwards due to the merge and splitting apart process.
        expt_excluded = [
            {675, 677, 1582, 1734, 1775},  # hap 1
            {370, 1194},  # hap 2
        ]

        for expt_vcf, got_vcf, excluded in zip([self.vcf1, self.vcf2], self.split_fps, expt_excluded):
            expt_vcfr = VCFReader(expt_vcf)
            got_vcfr = VCFReader(got_vcf)

            for expt in expt_vcfr.fetch():
                if expt.pos + 1 in excluded:
                    continue
                got = list(got_vcfr.fetch(expt.chrom, expt.pos, expt.pos + len(expt.ref) + 1))
                self.assertEqual(len(got), 1, 'Could not find split variant for {}:{}.'.format(expt.chrom, expt.pos+1))
                got = got[0]
                for key in  ('chrom', 'pos', 'ref', 'alt'):
                    expected = getattr(expt, key)
                    result = getattr(got, key)
                    self.assertEqual(expected, result, 'Splitting failed for {}:{} {}.'.format(expt.chrom, expt.pos+1, key))


class TestSplitMNP(unittest.TestCase):
    chrom = 'chr20'
    qual = 5
    info = {'my_meta': 10}

    def _make_variant(self, pos, ref, alt, gt):
        return Variant(self.chrom, pos, ref, alt,
                       genotype_data={'GT': '{}|{}'.format(*gt),
                                      'GQ': self.qual},
                       info=self.info)

    def test_classify_variant(self):
        cases = [
            # initial pos, ref, alt, gt
            # followed by split pos, ref, alt, gt
            ((60795, 'GG', ['AC', 'AG'], (1, 2)),
             (60795, 'G', ['A'], (1, 1)),
             (60796, 'G', ['C'], (1, 0))),
            ((60863, 'AA',['CC', 'AC'], (1, 2)),
             (60863, 'A', ['C'], (1, 0)),
             (60864, 'A', ['C'], (1, 1))),
            ((175688, 'GG', ['CT'], (1, 0)),
             (175688, 'G', ['C'], (1, 0)),
             (175689, 'G', ['T'], (1, 0))),
            ((359953, 'TGC', ['CAT'], (0, 1)),
             (359953, 'T', ['C'], (0, 1)),
             (359954, 'G', ['A'], (0, 1)),
             (359955, 'C', ['T'], (0, 1))),
            ((1000, 'TGC', ['C'], (0, 1)),  # non-MNP should come back unchanged
             (1000, 'TGC', ['C'], (0, 1))),

        ]
        for mnp, *split_vars in cases:
            got = split_mnp(self._make_variant(*mnp))
            expt = [self._make_variant(*i) for i in split_vars]
            self.assertEqual(len(expt), len(got))
            for g, e in zip(got, expt):
                msg = 'Attr {} failed for {}'
                for attr in 'chrom', 'pos', 'alt', 'info', 'gt':
                    self.assertEqual(getattr(g, attr), getattr(e, attr),
                                     msg=msg.format(attr, mnp))


class TestClassifyVariant(unittest.TestCase):

    def test_classify_variant(self):
        cases = [
            ('snp', 'G', ['A']),
            ('mnp', 'GG', ['AT']),
            ('snd', 'GA', ['G']),
            ('snd', 'GA', ['A']),
            ('mnd', 'GAT', ['G']),
            ('mnd', 'GAA', ['A']),
            ('sni', 'G', ['GT']),
            ('sni', 'G', ['AG']),
            ('mni', 'G', ['GTC']),
            ('mni', 'G', ['ATG']),
            ('other', 'G', ['ATA']),
            ('other', 'GAC', ['T']),
            # Classed as sub + ins, but with alignment one could class it as a mnd
            ('other', 'G', ['TGC']),
            # Classed as sub + del, but with alignment one could class it as a mni
            ('other', 'GAT', ['A']),
            ('indel', 'GG', ['G', 'GGA']),
            ('mnd', 'GGG', ['G', 'GG']),
            ('mni', 'G', ['GG', 'GGG']),
            ('mnp', 'GA', ['CT', 'CA']),
        ]
        for klass, ref, alts in cases:
            var = Variant('20', 14369, ref, alt=alts)
            self.assertEqual(klass, classify_variant(var), 'Classification failed for {} {} {}'.format(ref, alts, klass))


class TestGetPaddedHaplotypes(unittest.TestCase):

    def test_get_padded_haplotypes(self):
        chrom = 'my_chrom'
        ref_seq = 'ATGCTACTGC'
        # (pos, ref, alt), pad, padded ref, padded alt, start, end
        cases = [
            ((4, 'T', 'G'), 2, 'GCTAC', 'GCGAC', 2, 7),  #  sub
            ((4, 'T', 'TA'), 2, 'GCTAC', 'GCTAAC', 2, 7), #  ins
            ((4, 'T', 'GA'), 2, 'GCTAC', 'GCGAAC', 2, 7), #  sub ins
            ((4, 'TA', 'T'), 2, 'GCTACT', 'GCTCT', 2, 8), #  del
            ((4, 'TA', 'G'), 2, 'GCTACT', 'GCGCT', 2, 8), #  sub del
            # test what happens for variant at start and end of chrom
            ((0, 'A', 'G'), 2, 'ATG', 'GTG', 0, 3), #  sub at start
            ((0, 'A', 'AG'), 2, 'ATG', 'AGTG', 0, 3), #  ins at start
            ((0, 'AT', 'T'), 2, 'ATGC', 'TGC', 0, 4), #  del at start
            ((9, 'C', 'G'), 2, 'TGC', 'TGG', 7, 10), #  sub at end
            ((9, 'C', 'CG'), 2, 'TGC', 'TGCG', 7, 10), #  ins at end
            ((8, 'GC', 'G'), 2, 'CTGC', 'CTG', 6, 10), #  del at end
        ]
        for ((pos, ref, alt), pad, pad_ref, pad_alt, start, end) in cases:
            var = Variant(chrom, pos, ref, alt)
            padded, region = get_padded_haplotypes(var, ref_seq, pad)
            self.assertEqual(pad_ref, padded[0])
            self.assertEqual(pad_alt, padded[1])
            self.assertEqual(region.start, start)
            self.assertEqual(region.end, end)

    def test_raises(self):
        chrom = 'my_chrom'
        ref_seq = 'ATGCTACTGC'
        var = Variant(chrom, 2, 'GT', 'G')  # ref should be GC
        with self.assertRaises(ValueError):
            get_padded_haplotypes(var, ref_seq, 2)


def strip(r):
    """Return r.upper().replace('*'. '')"""
    return r.upper().replace('*', '')


class TestAlignReadToHaps(unittest.TestCase):

    def test_align_read_to_haps(self):
        g_open = 5
        g_ext = 3
        matrix = parasail.dnafull
        match = matrix.matrix[0, 0]
        mismatch = matrix.matrix[0, 1]
        counts = dict(zip(*np.unique(matrix.matrix[:4, :4],
                                     return_counts=True)))
        diag = np.unique(matrix.matrix.diagonal()[:4])[0]
        msg = 'Parasail matrix is not symmetric.'
        self.assertEqual(counts, {mismatch: 12, match: 4}, msg=msg)
        self.assertEqual(diag, match, msg=msg)

        read = 'ATGCTTTTTGCTAC'
        haps_scores = [
            ('ATGCTTTTTGCTAC',  len(read) * match),
            ('ATGCTTaTTGCTAC',  (len(read) - 1) * match + mismatch),
            ('ATGCTTTT*GCTAC',   (len(read) - 1) * match - g_open),
            ('ATGCTTT**GCTAC',  (len(read) - 2) * match - g_open - g_ext),
        ]
        scores = align_read_to_haps(read, [strip(h[0]) for h in haps_scores],
                                    g_open, g_ext, matrix)
        for (hap, exp), got in zip(haps_scores, scores):
            self.assertEqual(exp, got, msg='Failed for hap {}'.format(hap))

class TestAlignReadToHaps(unittest.TestCase):

    def test_align_read_to_haps(self):
        haps = [
            'ATGCTTTTT*GCTAC',  # ref
            'ATGCTTaTT*GCTAC',  # alt 1
            'ATGCTTTTTTGCTAC',  # alt 2
        ]
        reads = [
            'AaGCTTTTT*GCcAC',  # ref with a couple of sub errors
            'ATcCTTaTT*GCTgC',  # alt 1 with a couple of sub errors
            'ATGgTTTTTTGCcAC',  # alt 2 with a couple of sub errors
            'ATGCTTgTT*GCTAC',  # neither ref nor alt 1
        ]

        # remove * and make upper
        haps = [strip(h) for h in haps]
        reads = [strip(r) for r in reads]
        def align(r, h, go=5, ge=3, matrix=parasail.dnafull):
            return parasail.sw_trace_striped_32(r, h, go, ge, matrix).score
        # expected count for each read against first two haps
        for is_rev in False, True:
            # trimmed reads are always represented as fwd strand even if they
            # were reverse complement alignments
            reads_ws = [(is_rev, r) for r in reads]
            exp_cnts_pr = [  # expected result for each read
                Counter({(is_rev, 0): 1}),  # reads[0] aligns best to haps[0]
                Counter({(is_rev, 1): 1}),  # reads[1] aligns best to haps[1]
                Counter({(is_rev, 0): 1}),  # reads[2] aligns best to haps[1]
                Counter({(is_rev, None): 1}),  # reads[3] ambig wrt haps[:2]
            ]
            exp_scores_pr = [
                Counter({(is_rev, i): align(r[1], h)
                         for i, h in enumerate(haps[:2])})
                for r in reads_ws
            ]
            # one read at a time
            for read, cnts, scrs in zip(reads_ws, exp_cnts_pr, exp_scores_pr):
                cnts_got, scrs_got = align_reads_to_haps([read], haps[:2])
                self.assertEqual(cnts_got, cnts)
                self.assertEqual(scrs_got, scrs)
            # all reads together
            exp_cnts, exp_scrs = Counter(), Counter()
            for cnts, scrs in zip(exp_cnts_pr, exp_scores_pr):
                exp_cnts.update(cnts)
                exp_scrs.update(scrs)
            cnts_got, scrs_got = align_reads_to_haps(reads_ws, haps[:2])
            self.assertEqual(exp_cnts, cnts_got)
            self.assertEqual(exp_scrs, scrs_got)
        # counts for all three haps, where reads[3] should match haps[2]
        reads_ws = [(False, r) for r in reads]
        exp = Counter({(False, 0): 2, (False, 1): 1, (False, 2): 1})
        self.assertAlmostEqual(align_reads_to_haps(reads_ws, haps)[0], exp)
