import os
import types
import unittest

from medaka import smolecule

root_dir = os.path.abspath(os.path.dirname(__file__))
test_fasta = os.path.join(root_dir, 'data/smolecule.fasta')
test_mfasta = os.path.join(root_dir, 'data/smolecule_multi.fasta')


class TestRead(unittest.TestCase):


    def test_00_read_single(self):
        read = smolecule.Read.from_fastx(test_fasta)
        self.assertIsInstance(read, smolecule.Read)
        self.assertFalse(read._initialized, 'Read is uninitialized')
        self.assertEqual(read.nseqs, 9, 'Read has correct number of subreads.')


    def test_01_read_multi(self):
        reads = smolecule.Read.multi_from_fastx(test_mfasta)
        self.assertIsInstance(reads, types.GeneratorType, 'multi-read gives generator.')
        reads = [x for x in reads]
        n_subreads = [9, 22]
        self.assertEqual(len(reads), len(n_subreads), 'Retrieved correct number of subreads.')
        for i, (nsr, read) in enumerate(zip(n_subreads, reads)):
            self.assertEqual(read.nseqs, nsr, 'Read {} has correct number of subreads.'.format(i))


    def test_10_initialize(self):
        read = smolecule.Read.from_fastx(test_fasta)
        self.assertIsInstance(read, smolecule.Read)
        read.initialize()
        self.assertTrue(read._initialized, 'Read is initialized after .initialize().')
        self.assertFalse(read._alignments is None, '.alignments is not None after .initialize().')
        self.assertTrue(read._alignments_valid, '.alignments_valid is True after .initialize().')
        self.assertFalse(read._orient is None, '.orients is not None after .initialize().')


    def test_20_basic_consensus(self):
        read = smolecule.Read.from_fastx(test_fasta)
        cons = read.poa_consensus()
        self.assertTrue(read._initialized, 'Read is initialized after poa.')
        self.assertEqual(cons, read.consensus, 'Returned sequence is self.consensus.')
        self.assertFalse(read._alignments_valid, '.alignments_valid is False after .poa_consensus.()')
