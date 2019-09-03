"""Testing for medaka.rle module."""
import array
from collections import namedtuple
import tempfile
import unittest

import pysam

import medaka.rle


class RLE(unittest.TestCase):
    """Test medaka.rle.rle function."""

    def test_rle(self):
        """Test conversion of a basecall to RLE."""
        basecall = 'ACCGTTTA'
        expected = [
            (1, 0, 'A'), (2, 1, 'C'), (1, 3, 'G'), (3, 4, 'T'), (1, 7, 'A')]
        got = medaka.rle.rle(basecall).tolist()
        self.assertTrue(expected, got)

    def test_low_mem_vs_high_mem(self):
        """Test consistent output whether low_mem is True or False."""
        basecall = 'ACCGTTTA'
        low_mem = medaka.rle.rle(basecall, low_mem=True)
        high_mem = medaka.rle.rle(basecall, low_mem=False)
        self.assertTrue(all(low_mem == high_mem))


class RLEConversion(unittest.TestCase):
    """Test medaka.rle.RLEConverter class."""

    def test_compression(self):
        """Check recovery of original basecall from RLE version."""
        basecall = 'CCCTAGGTTA'
        rle_object = medaka.rle.RLEConverter(basecall)
        bases = rle_object.compact_basecall
        lengths = rle_object.homop_length
        got = ''.join(b * l for b, l in zip(bases, lengths))
        self.assertEqual(basecall, got)

    def test_alignment_001(self):
        """Check coordinate conversion between full and compact indices."""
        basecall = 'CCCTAGGTTA'
        rle_object = medaka.rle.RLEConverter(basecall)
        full_indices = range(len(basecall))

        # each base in a homopolymer maps to a base in the compressed sequence
        expected = (0, 0, 0) + (1, 2) + (3, 3) + (4, 4) + (5,)
        got = [
            rle_object.coord_full_to_compact(index) for index in full_indices]
        self.assertSequenceEqual(expected, got)

    def test_alignment_002(self):
        """Check coordinate conversion from compact to full indices."""
        basecall = 'CCCTAGGTTA'
        rle_object = medaka.rle.RLEConverter(basecall)
        compact_indices = range(len(rle_object.compact_basecall))

        expected = (0, 3, 4, 5, 7, 9)
        got = [
            rle_object.coord_compact_to_full(index)
            for index in compact_indices]
        self.assertSequenceEqual(expected, got)


class ParasailAlignment(unittest.TestCase):
    """Check medaka.rle.parasail_alignment function."""

    ref = 'AGCATGTTAGATAAGATA'

    def test_alignment_001(self):
        """Test alignment of read starting mid-reference and clip at end.

        AGCATGTTAGATAA**GATA
        TTAGATAAAGGATActg
        """
        seq = 'TTAGATAAAGGATACTG'
        rstart, cigar = medaka.rle.parasail_alignment(seq, self.ref)
        expected = (6, '8=2I4=3S')
        self.assertEqual(cigar, expected[1], 'cigars differ!')
        self.assertEqual(rstart, expected[0], 'rstart wrong!')

    def test_alignment_002(self):
        """Test alignment of read starting mid-reference with clip at begining.

        AGCATGTTAGATAA*GATA
        aaaAGATAAGGATA
        """
        seq = 'AAAAGATAAGGATA'
        rstart, cigar = medaka.rle.parasail_alignment(seq, self.ref)
        expected = (8, '3S6=1I4=')
        self.assertEqual(cigar, expected[1], 'cigars differ!')
        self.assertEqual(rstart, expected[0], 'rstart wrong!')

    def test_alignment_003(self):
        """Test alignment with clipping in both ends.

        AGCATGTTAGATAAGATA
        ccCATGTTAGATAAGcc
        """
        seq = 'CCCATGTTAGATAAGCC'
        rstart, cigar = medaka.rle.parasail_alignment(seq, self.ref)
        expected = (2, '2S13=2S')
        self.assertEqual(cigar, expected[1], 'cigars differ!')
        self.assertEqual(rstart, expected[0], 'rstart wrong!')

    def test_alignment_004(self):
        """Test alignment of read starting before the reference.

        AGCATGTTAGATAAGATA
        ctgAGCATGT
        """
        seq = 'CTGAGCATGT'
        rstart, cigar = medaka.rle.parasail_alignment(seq, self.ref)
        expected = (0, '3S7=')
        self.assertEqual(cigar, expected[1], 'cigars differ!')
        self.assertEqual(rstart, expected[0], 'rstart wrong!')


class Clipping(unittest.TestCase):
    """Test soft clipping of parasail cigar strings."""

    cigar1 = '4S2=1D40=1X29=2D12=1I10=3S'
    cigar2 = '4S2=1D40=1X29=2D12=1I10='
    cigar3 = '2=1D40=1X29=2D12=1I10=3S'
    cigar4 = '2=1D40=1X29=2D12=1I10='
    cigars = (cigar1, cigar2, cigar3, cigar4)

    def test_no_added_clipping(self):
        """Do not add anything, cigar should return untouched."""
        for cigar in self.cigars:
            got = medaka.rle.add_extra_clipping(cigar, 0, 0)
            self.assertEqual(cigar, got)

    def test_add_start(self):
        """Add clipping only at the beginning of the cigar."""
        start_clipped = 4
        expected = (
            '8S2=1D40=1X29=2D12=1I10=3S',
            '8S2=1D40=1X29=2D12=1I10=',
            '4S2=1D40=1X29=2D12=1I10=3S',
            '4S2=1D40=1X29=2D12=1I10=')
        got = tuple(
            medaka.rle.add_extra_clipping(cigar, start_clipped, 0)
            for cigar in self.cigars)
        self.assertEqual(expected, got)

    def test_add_end(self):
        """Add clipping at the end of the cigar."""
        end_clipped = 8
        expected = (
            '4S2=1D40=1X29=2D12=1I10=11S',
            '4S2=1D40=1X29=2D12=1I10=8S',
            '2=1D40=1X29=2D12=1I10=11S',
            '2=1D40=1X29=2D12=1I10=8S')
        got = tuple(
            medaka.rle.add_extra_clipping(cigar, 0, end_clipped)
            for cigar in self.cigars)
        self.assertEqual(expected, got)

    def test_add_both(self):
        """Add clipping both ends of the cigar."""
        start_clipped = 6
        end_clipped = 9
        expected = (
            '10S2=1D40=1X29=2D12=1I10=12S',
            '10S2=1D40=1X29=2D12=1I10=9S',
            '6S2=1D40=1X29=2D12=1I10=12S',
            '6S2=1D40=1X29=2D12=1I10=9S')
        got = tuple(
            medaka.rle.add_extra_clipping(cigar, start_clipped, end_clipped)
            for cigar in self.cigars)
        self.assertEqual(expected, got)


class InitialiseAlignment(unittest.TestCase):
    """Test medaka.rle.initialise_alignment."""

    input_kwargs = {
        'query_name': 'test',
        'reference_id': 0,
        'reference_start': 2,
        'query_sequence': 'GCCCTGTTGATCTT',
        'cigarstring': "1S3=1D8=2S",
        'flag': 0,
        'mapping_quality': 60}
    len_query = len(input_kwargs['query_sequence'])
    input_kwargs['query_qualities'] = array.array('B', [20] * len_query)

    def test_inputs(self):
        """Test inputs are correctly passed to alignment."""
        alignment = medaka.rle.initialise_alignment(**self.input_kwargs)
        for key, expected in self.input_kwargs.items():
            got = getattr(alignment, key)
            self.assertEqual(expected, got)

    def test_derived(self):
        """Test arguments derived from inputs."""
        alignment = medaka.rle.initialise_alignment(**self.input_kwargs)
        expected_kwargs = {
            'query_alignment_start': 1,
            'query_alignment_end': 12,
            'query_alignment_sequence': 'CCCTGTTGATC'}

        for key, expected in expected_kwargs.items():
            got = getattr(alignment, key)
            self.assertEqual(got, expected)


class CompressAlignment(unittest.TestCase):
    """Test medaka.rle._compress_alignment function."""

    alignment_kwargs = {
        'query_name': 'test',
        'reference_id': 0,
        'reference_start': 2,
        'query_sequence': 'GCCCAGTTGATCTT',
        'cigarstring': "1S3=1D8=2S",
        'flag': 0,
        'mapping_quality': 60}
    ref = 'TACCCATGTTGATCG'

    def test_compression(self):
        """Compress alignment.

        ref: TACCCATGTTGATCG  --> TACATGTGATCG
        seq:  gCCCA*GTTGATCtt -->  gCA*GTGATCt
        cigar:    1S4=1D7=2S  -->  1S2=1D6=1S
        """
        alignment = medaka.rle.initialise_alignment(**self.alignment_kwargs)
        ref_rle = medaka.rle.RLEConverter(self.ref)
        compressed_alignment = medaka.rle._compress_alignment(
            alignment, ref_rle)
        real_outputs = {
            'cigarstring': '1S2=1D6=1S',
            'query_sequence': 'GCAGTGATCT',
            'query_alignment_start': 1,
            'query_alignment_end': 9,
            'query_alignment_sequence': 'CAGTGATC',
            'reference_start': 2,
            'reference_end': 11}
        for key, expected in real_outputs.items():
            got = getattr(compressed_alignment, key)
            self.assertEqual(got, expected)


class CompressSeq(unittest.TestCase):
    """Test medaka.rle.compress_seq function."""

    def test_compress_read(self):
        """Given an input record, check the returned RLE compressed version."""
        read = pysam.FastxRecord(
            name='test',
            comment='runid=b81',
            sequence='ACCGTTTAC')
        compressed_read = medaka.rle.compress_seq(read)
        true_output = {
            'name': read.name,
            'comment': read.comment,
            'sequence': 'ACGTAC',
            'quality': '"#"$""'}

        for key, expected in true_output.items():
            got = getattr(compressed_read, key)
            self.assertEqual(got, expected)


class CompressBasecallsTest(unittest.TestCase):
    """Test medaka.rle.compress_basecalls function."""

    @classmethod
    def setUpClass(cls):
        """Prepare for tests.

        Create an input fasta, prepare an args object with correct attributes.
        """
        tmp_input = tempfile.NamedTemporaryFile(suffix='.fa').name
        with open(tmp_input, 'w') as input_handle:
            input_handle.write('>{}\n{}\n'.format('input1', 'AGTTTGGCTCCCCCA'))
            input_handle.write('>{}\n{}\n'.format('input2', 'AAACTTTCCCC'))

        tmp_output = tempfile.NamedTemporaryFile(suffix='.fq').name

        args = namedtuple('args', ['input', 'output', 'threads'])
        cls.args = args(tmp_input, tmp_output, 2)

    def test_compress_runs(self):
        """Test the code runs without raising exceptions."""
        try:
            medaka.rle.compress_basecalls(self.args)
        except Exception:
            self.fail('medaka.rle.compress_basecalls() raised an exception.')

    def test_output(self):
        """Compare output fastq to expected sequence and scores."""
        expected = {
            'input1': ('AGTGCTCA', (1, 1, 3, 2, 1, 1, 5, 1)),
            'input2': ('ACTC', (3, 1, 3, 4))}
        got = {}
        with pysam.FastqFile(self.args.output) as output_handle:
            for read in output_handle:
                lengths = tuple(read.get_quality_array())
                got[read.name] = (read.sequence, lengths)
        print('Javi: {}'.format(got))
        self.assertEqual(expected, got)


if __name__ == '__main__':
    unittest.main()
