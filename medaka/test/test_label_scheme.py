import unittest
from collections import Counter
import numpy as np

import medaka.labels


class HaploidLabelSchemeTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.counts = Counter(
            {
                # ploidy * (base, length) tuples
                ((0, 1),): 15,   # *
                ((1, 1),): 10,   # A
                ((2, 1),): 20,   # C
                ((3, 1),): 5,    # G
                ((4, 1),): 18,   # T
                ((1, 2),): 1,    # AA
                ((1, 4),): 1,    # AAAA
            }
        )
        cls.ls = medaka.labels.HaploidLabelScheme(cls.counts, max_label_len=2)


    def test_truncation(self):
        """Test that labels are correctly truncated to max_label_len"""

        expt = Counter(
            {
                ((0, 1),): 15,
                ((1, 1),): 10,
                ((1, 2),): 2,
                ((2, 1),): 20,
                ((3, 1),): 5,
                ((4, 1),): 18,
            }
        )
        self.assertEqual(self.ls._truncated_input_counts_, expt)


    def test_label_description(self):

        self.assertEqual(self.ls.label_description,
            (
                '* haploid',
                'A haploid',
                'AA haploid',
                'C haploid',
                'G haploid',
                'T haploid'
            )
        )


    def test_label_encoding(self):

        expt = {
            ((0, 1),): (0,),
            ((1, 1),): (1,),
            ((1, 2),): (2,),
            ((2, 1),): (3,),
            ((3, 1),): (4,),
            ((4, 1),): (5,),
        }
        self.assertEqual(self.ls.label_encoding, expt)


    def test_label_decoding(self):

        self.assertEqual(self.ls.label_decoding, ('*', 'A', 'AA', 'C', 'G', 'T'))


    def test_ploidy(self):

        self.assertEqual(self.ls.ploidy, 1)


class FactoredBaseZygosityLabelScheme(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.counts = Counter(
            {
                # ploidy * (base, length) tuples
                ((0, 1), (0, 1)): 15,   # *, *
                ((1, 1), (1, 1)): 10,   # A, A
                ((2, 1), (3, 1)): 20,   # C, G
                ((4, 1), (4, 1)): 18,   # T, T
                ((1, 2), (1, 2)): 1,    # AA, AA
                ((1, 1), (1, 2)): 1,    # A, AA
                ((1, 4), (1, 5)): 1,    # AAAA, AAAAA
            }
        )
        cls.ls = medaka.labels.FactoredBaseZygosityLabelScheme(cls.counts, max_label_len=2)


    def test_truncation(self):
        """Test that labels are correctly truncated to max_label_len"""

        expt = Counter(
            {
                ((0, 1), (0, 1)): 15,
                ((1, 1), (1, 1)): 10,
                ((1, 1), (1, 2)): 1,
                ((1, 2), (1, 2)): 2,
                ((2, 1), (3, 1)): 20,
                ((4, 1), (4, 1)): 18
            }
        )
        self.assertEqual(self.ls._truncated_input_counts_, expt)


    def test_label_description(self):

        self.assertEqual(self.ls.label_description,
            (
                '* multi-label factored base & zygosity',
                'A multi-label factored base & zygosity',
                'AA multi-label factored base & zygosity',
                'C multi-label factored base & zygosity',
                'G multi-label factored base & zygosity',
                'T multi-label factored base & zygosity',
                'homozygous multi-label factored base & zygosity',
                'heterozygous multi-label factored base & zygosity',
            )
        )


    def test_label_encoding(self):

        expt = {
            ((0, 1), (0, 1 )): (0, 6),
            ((1, 1), (1, 1)): (1, 6),
            ((1, 1), (1, 2)): (1, 2, 7),
            ((1, 2), (1, 2)): (2, 6),
            ((2, 1), (3, 1)): (3, 4, 7),
            ((4, 1), (4, 1)): (5, 6)
        }
        self.assertEqual(self.ls.label_encoding, expt)


    def test_label_decoding(self):

        self.assertEqual(self.ls.label_decoding, ('*', 'A', 'AA', 'C', 'G', 'T', 'homozygous', 'heterozygous'))


    def test_ploidy(self):

        self.assertEqual(self.ls.ploidy, 2)


if __name__ == '__main__':
    unittest.main()
