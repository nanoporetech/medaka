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


@dataclasses.dataclass
class DummyArgs:
    model: str
    data: str


class TestCheckModelNeedsDwells(unittest.TestCase):

    def test_001_check_no_dwells_model(self):
        # check non dwells model needs neither dwells nor haplotypes
        model = medaka.models.resolve_model('r1041_e82_400bps_hac_v5.0.0_rl_lstm384_no_dwells')
        needs_dwells, needs_haplotypes = medaka.medaka.model_needs_dwells_and_haplotype(model)
        self.assertFalse(needs_dwells)
        self.assertFalse(needs_haplotypes)

    def test_002_check_dwells_model(self):
        # check dwells model needs dwells, not haplotypes
        model = medaka.models.resolve_model('r1041_e82_400bps_hac_v5.0.0_rl_lstm384_dwells')
        needs_dwells, needs_haplotypes = medaka.medaka.model_needs_dwells_and_haplotype(model)
        self.assertTrue(needs_dwells)
        self.assertFalse(needs_haplotypes)

    def test_003_check_pileup_from_dict_model(self):
        # check counts model created with model_from_dict, needs neither
        model = medaka.models.resolve_model('r1041_e82_400bps_hac_v6.0.0')
        needs_dwells, needs_haplotypes = medaka.medaka.model_needs_dwells_and_haplotype(model)
        self.assertFalse(needs_dwells)
        self.assertFalse(needs_haplotypes)

    def test_004_check_pileup_torch_model(self):
        # check counts model created with build_torch_model, needs neither
        model = medaka.models.resolve_model('r1041_e82_400bps_hac_v5.0.0')
        needs_dwells, needs_haplotypes = medaka.medaka.model_needs_dwells_and_haplotype(model)
        self.assertFalse(needs_dwells)
        self.assertFalse(needs_haplotypes)


class TestCheckCompatible(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.root_dir = os.path.abspath(os.path.dirname(__file__))
        cls.no_dwell_bam = os.path.join(cls.root_dir, 'data', 'test_reads.bam')
        cls.dwell_bam = os.path.join(cls.root_dir, 'data', 'test_reads_dwells.bam') # same bam with dummy move table
        cls.no_dwell_fasta = os.path.join(cls.root_dir, 'data', 'test_reads.fastq')
        cls.dwell_fasta = os.path.join(cls.root_dir, 'data', 'test_reads_dwells.fastq')
        cls.no_dwells_model = medaka.models.resolve_model('r1041_e82_400bps_hac_v5.0.0_rl_lstm384_no_dwells')
        cls.dwells_model = medaka.models.resolve_model('r1041_e82_400bps_hac_v5.0.0_rl_lstm384_dwells')
        cls.counts_model = medaka.models.resolve_model('r1041_e82_400bps_hac_v5.0.0')

    def test_001_check_bam_compatible_dwells(self):
        self.assertFalse(medaka.common.check_bam_for_dwells(self.no_dwell_bam))
        self.assertTrue(medaka.common.check_bam_for_dwells(self.dwell_bam))

    def test_002_check_fastx_compatible_dwells(self):
        self.assertFalse(medaka.common.check_fastx_for_dwells(self.no_dwell_fasta))
        self.assertTrue(medaka.common.check_fastx_for_dwells(self.dwell_fasta))

    def test_003_check_compatible_no_file(self):
        # check non-existant model names fail
        with self.assertRaises(FileNotFoundError):
            medaka.medaka.check_compatible(
                DummyArgs(model='NotAModel', data=self.no_dwell_bam)
            )

    def test_004_check_compatible_no_dwells_model(self):
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

    def test_005_check_compatible_dwells_model(self):
        # check dwells model is not compatible with non dwells bam
        with self.assertRaises(ValueError):
            medaka.medaka.check_compatible(
                DummyArgs(model=self.dwells_model, data=self.no_dwell_bam)
            )

        # check dwells model is compatible with dwells bam
        medaka.medaka.check_compatible(
            DummyArgs(model=self.dwells_model, data=self.dwell_bam)
        )
        # check dwells model is compatible with dwells fastq
        medaka.medaka.check_compatible(
            DummyArgs(model=self.dwells_model, data=self.dwell_fasta)
        )

    def test_006_check_compatible_counts_model(self):
        # check counts model is always compatible
        medaka.medaka.check_compatible(
            DummyArgs(model=self.counts_model, data=self.dwell_bam)
        )
        medaka.medaka.check_compatible(
            DummyArgs(model=self.counts_model, data=self.no_dwell_bam)
        )
