import argparse
import dataclasses
import os
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
                                        

