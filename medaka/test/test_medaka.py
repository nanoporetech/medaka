import argparse
import contextlib
import dataclasses
import gzip
import io
import os
import tempfile
import unittest

import medaka.medaka
import medaka.options


class ParseDictArgTest(unittest.TestCase):

    def test_001_basic_counting(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('--dict', action=medaka.medaka.StoreDict, nargs='+')

        args = (
            'none=None str=string strs=r9,r10 an_int=1 a_float=1.0 '
            'numbers=1,10.0 a_bool=False bools=True,true,TRUE,False,false,FALSE').split()
        expected = {
            'none':None, 'str':'string', 'strs': ['r9', 'r10'], 'an_int': 1, 'a_float': 1.0,
            'numbers': [1, 10.0], 'a_bool': False, 'bools': [True, True, True, False, False, False]}
        parsed = parser.parse_args(['--dict'] + args)
        self.assertEqual(vars(parsed)['dict'], expected)


class TestTandemRegionParser(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root_dir = os.path.abspath(os.path.dirname(__file__))
        cls.test_bam = os.path.join(cls.root_dir, 'data', 'tandem', 'test.bam')
        cls.test_ref = os.path.join(
            cls.root_dir, 'data', 'tandem', 'chr20.fa.gz')

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        medaka.common.remove_directory(self.temp_dir)

    def _parse_tandem(self, regions):
        parser = medaka.medaka.medaka_parser()
        return parser.parse_args([
            'tandem',
            self.test_bam,
            self.test_ref,
            regions,
            'male',
            self.temp_dir,
            '--ignore_read_groups',
        ])

    def _assert_parser_error(self, regions, msg):
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            with self.assertRaises(SystemExit) as raised:
                self._parse_tandem(regions)
        self.assertEqual(raised.exception.code, 2)
        self.assertIn(msg, stderr.getvalue())

    def test_missing_bed_path_fails_fast(self):
        missing_bed = os.path.join(self.temp_dir, 'missing_regions.bed')
        self._assert_parser_error(
            missing_bed, "BED file for 'regions' does not exist")

    def test_empty_bed_fails_fast(self):
        empty_bed = os.path.join(self.temp_dir, 'empty_regions.bed')
        with open(empty_bed, 'w'):
            pass
        self._assert_parser_error(
            empty_bed, "No valid tandem regions were parsed from BED input")

    def test_non_bed_region_string_fails_fast(self):
        self._assert_parser_error(
            'chr20:10-50', "BED file for 'regions' does not exist")

    def test_non_positive_region_fails_fast(self):
        bad_bed = os.path.join(self.temp_dir, 'bad_regions.bed')
        with open(bad_bed, 'w') as fh:
            fh.write("chr20\t10\t10\n")
        self._assert_parser_error(
            bad_bed, "Tandem regions require end > start")

    def test_valid_bed_parses(self):
        good_bed = os.path.join(self.temp_dir, 'good_regions.bed')
        with open(good_bed, 'w') as fh:
            fh.write("chr20\t100\t140\n")
        args = self._parse_tandem(good_bed)
        self.assertEqual(
            args.regions, [medaka.common.Region('chr20', 100, 140)])

    def test_valid_gzipped_bed_parses(self):
        good_bed = os.path.join(self.temp_dir, 'good_regions.bed.gz')
        with gzip.open(good_bed, 'wt') as fh:
            fh.write("chr20\t100\t140\n")
        args = self._parse_tandem(good_bed)
        self.assertEqual(
            args.regions, [medaka.common.Region('chr20', 100, 140)])


class TestCheckCompatible(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.root_dir = os.path.abspath(os.path.dirname(__file__))
        cls.no_dwell_bam = os.path.join(cls.root_dir, 'data', 'test_reads.bam')
        cls.dwell_bam = os.path.join(cls.root_dir, 'data', 'test_reads_dwells.bam') # same bam with dummy move table 
        cls.no_dwell_fasta = os.path.join(cls.root_dir, 'data', 'test_reads.fastq')
        cls.dwell_fasta = os.path.join(cls.root_dir, 'data', 'test_reads_dwells.fastq')
        cls.no_dwells_model = medaka.models.resolve_model(medaka.options.default_models['consensus'])
        cls.dwells_model = medaka.models.resolve_model('r1041_e82_400bps_hac_v5.0.0_rl_lstm384_dwells')

    def test_001_check_bam_compatible_dwells(self):
        self.assertFalse(medaka.common.check_bam_for_dwells(self.no_dwell_bam))
        self.assertTrue(medaka.common.check_bam_for_dwells(self.dwell_bam))

    def test_002_check_fastx_compatible_dwells(self):
        self.assertFalse(medaka.common.check_fastx_for_dwells(self.no_dwell_fasta))
        self.assertTrue(medaka.common.check_fastx_for_dwells(self.dwell_fasta))

    def test_003_check_compatible(self):
        @dataclasses.dataclass
        class DummyArgs:
            model: str
            data: str

        # check non-existant model names fail
        with self.assertRaises(FileNotFoundError):
            medaka.medaka.check_compatible(
                DummyArgs(model='NotAModel', data=self.no_dwell_bam)
            )

        # check non dwells model is compatible with non dwells bam
        medaka.medaka.check_compatible(
            DummyArgs(model=self.no_dwells_model, data=self.no_dwell_bam)
        )
        # check non-dwells model is compatible with non-dwells fastq
        medaka.medaka.check_compatible(
            DummyArgs(model=self.no_dwells_model, data=self.no_dwell_fasta)
        )

        # check non-dwells model is compatible with dwells fastq
        medaka.medaka.check_compatible(
            DummyArgs(model=self.no_dwells_model, data=self.dwell_fasta)
        )
        
        # check dwells model is not compatible with non dwells bam
        with self.assertRaises(ValueError):
            medaka.medaka.check_compatible(
                DummyArgs(model=self.dwells_model, data=self.no_dwell_bam)
            )
                                        
