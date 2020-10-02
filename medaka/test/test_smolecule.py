import io
import os
import types
import unittest
import warnings

import pytest

from medaka import smolecule
import medaka.common

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

    def test_012_interleave(self):
        read = smolecule.Read.from_fastx(test_fasta)
        orient, subreads = read.interleaved_subreads
        exp_orient = (True, False, True, False, True, False, True, False, True)
        self.assertEqual(orient, exp_orient, "Orientations interleaved")
        # note, the sort is stable so this test is only useful
        # because some sorting needs to be done on the test set
        orig_order = [r.name for r in subreads]
        new_order = [r.name for r in read.subreads]
        self.assertNotEqual(orig_order, new_order, "Reads are reordered.")

    @pytest.mark.skipif("CITEST" in os.environ, reason="CI instruction-set issue")
    def test_20_basic_consensus(self):
        read = smolecule.Read.from_fastx(test_fasta)
        cons = read.poa_consensus()
        self.assertTrue(read._initialized, 'Read is initialized after poa.')
        self.assertEqual(cons, read.consensus, 'Returned sequence is self.consensus.')
        self.assertFalse(read._alignments_valid, '.alignments_valid is False after .poa_consensus.()')

    @pytest.mark.skipif("CITEST" in os.environ, reason="CI instruction-set issue")
    def test_25_basic_consensus_racon(self):
        read = smolecule.Read.from_fastx(test_fasta)
        cons = read.poa_consensus(method='racon')
        self.assertTrue(read._initialized, 'Read is initialized after poa.')
        self.assertEqual(cons, read.consensus, 'Returned sequence is self.consensus.')
        self.assertFalse(read._alignments_valid, '.alignments_valid is False after .poa_consensus.()')

    def test_30_parasail_align(self):
        revcom = medaka.common.reverse_complement
        seq_mult = 100
        seq = 'ACGACTACGACTACGACT' * seq_mult
        sub_reads = [
            (smolecule.Subread('read_0', seq), 0),
            (smolecule.Subread('read_1', seq), 0),
            (smolecule.Subread('read_2', revcom(seq)), 16)]
        read = smolecule.Read('test', [s[0] for s in sub_reads])

        expected = [
            smolecule.Alignment(
                'test', sr.name, flag, 0,
                sr.seq if flag == 0 else revcom(sr.seq),
                '{}='.format(len(seq)))
            for i, (sr, flag) in enumerate(sub_reads)]

        for aligner in ('align_to_template', 'mappy_to_template'):
            func = getattr(read, aligner)
            alignments = func(sub_reads[0][0].seq, 'test')
            self.assertEqual(len(alignments), len(expected))
            for aln, exp in zip(alignments, expected):
                for attr, exp_attr in zip(aln, exp):
                    # mappy doesn't report equality, just match
                    if aligner == 'mappy_to_template' and exp_attr == '{}='.format(len(seq)):
                        exp_attr = '{}M'.format(len(seq))
                    self.assertEqual(attr, exp_attr)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            read.mappy_to_template(sub_reads[0][0].seq, 'test', align=False)
            assert issubclass(w[-1].category, DeprecationWarning)

