"""Testing for medaka.rle module."""
import array
from collections import namedtuple
import os
import tempfile
import unittest

import h5py
import numpy as np
import pysam

import medaka.rle


def mock_fast5_file(read_info):
    """Create a fast5 file with fake shape/scale rle data"""
    fast5_file = tempfile.NamedTemporaryFile(suffix='.fast5').name
    data_path = (
        'read_{}/Analyses/Basecall_1D_000/'
        'BaseCalled_template/RunlengthBasecall')

    with h5py.File(fast5_file, 'w') as h:
        for read_id, (bases, shapes, scales) in read_info.items():
            arr = np.fromiter(
                [(b, s, sc) for (b, s, sc) in zip(bases, shapes, scales)],
                dtype=[('base', 'S1'), ('shape', '>f4'), ('scale', '>f4')])
            h.create_dataset(data_path.format(read_id), data=arr)

    return fast5_file


def mock_summary_file(read_ids, fast5_fnames):
    """Create a summary file to link read_id and fast5 filename"""
    summary_file = tempfile.NamedTemporaryFile(suffix='.txt').name
    with open(summary_file, 'w') as output:
        output.write('read_id\tfilename\n')
        for read_id, filename in zip(read_ids, fast5_fnames):
            output.write('{}\t{}\n'.format(read_id, filename))

    return summary_file

class RLE(unittest.TestCase):
    """Test medaka.rle.rle function."""

    def test_rle(self):
        """Test conversion of a basecall to RLE."""
        basecalls = ['ACCGTTTA', 'AA', 'AT', 'A']
        expected = [
            [(1, 0, 'A'), (2, 1, 'C'), (1, 3, 'G'), (3, 4, 'T'), (1, 7, 'A')],
            [(2, 0, 'A'),],
            [(1, 0, 'A'), (1, 1, 'T')],
            [(1, 0, 'A'),]]
        for call, exp in zip(basecalls, expected):
            got = medaka.rle.rle(call).tolist()
            self.assertTrue(exp, got)

    def test_input_must_be_1D(self):
        """Test that passing in a 2D array raises TypeError"""
        invalid = np.array([
            ['A', 'C', 'A'],
            ['C', 'T', 'G']], dtype='U1')
        with self.assertRaises(TypeError):
            medaka.rle.rle(invalid)


class RLEConversion(unittest.TestCase):
    """Test medaka.rle.RLEConverter class."""

    def test_compression(self):
        """Check recovery of original basecall from RLE version."""
        basecalls = ['CCCTAGGTTA', 'AA', 'AT', 'A']
        for basecall in basecalls:
            rle_object = medaka.rle.RLEConverter(basecall)
            bases = rle_object.compact_basecall
            lengths = rle_object.homop_length
            got = ''.join(b * l for b, l in zip(bases, lengths))
            self.assertEqual(basecall, got)

    def test_alignment_001(self):
        """Check coordinate slicing between full and compact indices."""
        basecall = 'CCCTAGGTTA'
        rle_seq = 'CTAGTA'
        rle_object = medaka.rle.RLEConverter(basecall)
        self.assertEqual(rle_object.compact_basecall, rle_seq)

        slcs = [(0, 10), (1, 9), (2, 8), (3, 7), (4, 6)]
        exp =  [(0,  6), (0, 5), (0, 5), (1, 4), (2, 4)]

        for (in_s, in_e), (exp_s, exp_e) in zip(slcs, exp):
            s, e = rle_object.transform_coords(in_s, in_e)
            trimmed_input = basecall[in_s:in_e]
            trimmed_input_compressed = medaka.rle.RLEConverter(trimmed_input).compact_basecall
            trimmed_compressed = rle_object.trimmed_compact(in_s, in_e)

            self.assertSequenceEqual(
                (exp_s, exp_e), (s, e),
                'Slice co-ordinates incorrect.')
            self.assertEqual(
                trimmed_compressed, rle_seq[s:e],
                "Direct trimming did not give same result.")
            self.assertEqual(
                trimmed_input_compressed, trimmed_compressed,
                "Trimming and compressing did not commute.")


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

    def test_unmapped_return_None(self):
        """Unmapped or secondary reads are skipped"""
        expected = None
        alignment = pysam.AlignedSegment()
        alignment.is_unmapped = True
        got = medaka.rle._compress_alignment(alignment, None)
        self.assertEqual(expected, got)


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


args_class = namedtuple(
    'args', [
        'bam_input', 'bam_output', 'ref_fname', 'threads', 'regions',
        'use_fast5_info'])


class CompressBamTest(unittest.TestCase):
    """Test medaka.rle.compress_bam"""

    @classmethod
    def setUpClass(cls):
        """Create temporary files and bam file

        Ref     T  T  A  A    C  T  T  T  G
        Read1         A  A    C  T  T  T  G
        Read2      T  A  A  A C  T  T  T  G
        """
        cls.bam_input = tempfile.NamedTemporaryFile(suffix='.bam').name
        cls.bam_output = tempfile.NamedTemporaryFile(suffix='.bam').name
        cls.ref_fname = tempfile.NamedTemporaryFile(suffix='.fasta').name

        with open(cls.ref_fname, 'w') as fasta:
            fasta.write('>ref\n')
            fasta.write('TTAACTTTG\n')

        header = {
            'HD': {'VN': '1.0'},
            'SQ': [{'LN': 9, 'SN': 'ref'}, ]}

        basecalls = {
            'read1': {
                'query_name': 'read1',
                'reference_id': 0,
                'reference_start': 2,
                'query_sequence': 'AACTTTG',
                'cigarstring': '7=',
                'flag': 0,
                'mapping_quality': 50},
            'read2': {
                'query_name': 'read2',
                'reference_id': 0,
                'reference_start': 1,
                'query_sequence': 'TAAACTTTG',
                'cigarstring': '3=1I5=',
                'flag': 0,
                'mapping_quality': 50}}

        tmp_file = '{}.tmp'.format(cls.bam_input)
        with pysam.AlignmentFile(tmp_file, 'wb', header=header) as bam:
            for basecall in basecalls.values():
                record = medaka.rle.initialise_alignment(**basecall)
                bam.write(record)

        pysam.sort("-o", cls.bam_input, tmp_file)
        os.remove(tmp_file)
        pysam.index(cls.bam_input)

    def test_output_rle_bam(self):
        expected = {
            ('read1', 1, 'ACTG',  '4=', (2, 1, 3, 1)),
            ('read2', 0, 'TACTG', '5=', (1, 3, 1, 3, 1))}

        # None vs full region, no difference
        for regions in (None, ['ref:0-10']):
            args = args_class(
                self.bam_input, self.bam_output, self.ref_fname,
                2, regions, None)
            medaka.rle.compress_bam(args)
            got = set()
            for read in pysam.AlignmentFile(self.bam_output):
                data = (
                    read.query_name, read.reference_start, read.query_sequence,
                    read.cigarstring, tuple(read.get_forward_qualities()))
                got.add(data)
            self.assertEqual(got, expected)

    def test_bam_compression_with_RLE_parameters(self):
        read_info = {
            'read1': (
                ['A', 'C', 'T', 'G'],
                [0.1, 0.2, 0.1, 0.5],
                [1., 2., 3., 4.]),
            'read2': (
                ['T', 'A', 'C', 'T', 'G'],
                [0.3, 0.1, 0.7, 0.1, 0.2],
                [5., 6., 7., 8., 9.])}

        # Create a mock fast5 file with only an RLE table
        fast5_path = mock_fast5_file(read_info)
        fast5_dir = os.path.dirname(fast5_path)
        fast5_fname = os.path.basename(fast5_path)

        # Create a mock summary file with read_id and filename
        read_ids = list(read_info.keys())
        summary_file = mock_summary_file(read_ids, [fast5_fname] * len(read_ids))

        regions = ['ref:0-10']
        args = args_class(
            self.bam_input, self.bam_output, self.ref_fname,
            2, regions, (fast5_dir, summary_file))
        medaka.rle.compress_bam(args)

        for read in pysam.AlignmentFile(self.bam_output):
            expected = read_info[read.query_name][1]
            got = list(read.get_tag('WL'))
            np.testing.assert_allclose(got, expected)

            expected = read_info[read.query_name][2]
            got = list(read.get_tag('WK'))
            np.testing.assert_allclose(got, expected)
