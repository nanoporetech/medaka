import argparse
import unittest

import medaka.medaka


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
