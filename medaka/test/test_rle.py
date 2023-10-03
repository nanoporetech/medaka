"""Testing for medaka.rle module."""
from collections import namedtuple
import os
import tempfile
import unittest

import numpy as np
import parasail
import pysam

from medaka import common
from .mock_data import create_simple_bam, mock_summary_file, \
    mock_fast5_file, simple_data
import medaka.rle


class RLE(unittest.TestCase):
    """Test medaka.common.rle function."""

    def test_rle(self):
        """Test conversion of a basecall to RLE."""
        basecalls = ['ACCGTTTA', 'AA', 'AT', 'A']
        expected = [
            [(1, 0, 'A'), (2, 1, 'C'), (1, 3, 'G'), (3, 4, 'T'), (1, 7, 'A')],
            [(2, 0, 'A'), ],
            [(1, 0, 'A'), (1, 1, 'T')],
            [(1, 0, 'A'), ]]
        for call, exp in zip(basecalls, expected):
            got = medaka.common.rle(call).tolist()
            self.assertTrue(exp, got)

    def test_input_must_be_1D(self):
        """Test that passing in a 2D array raises TypeError"""
        invalid = np.array([
            ['A', 'C', 'A'],
            ['C', 'T', 'G']], dtype='U1')
        with self.assertRaises(TypeError):
            medaka.common.rle(invalid)


class RLEConversion(unittest.TestCase):
    """Test medaka.rle.RLEConverter class."""

    @unittest.skipIf(parasail.dnafull is None, "Using fake parasail.")
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
        exp = [(0,  6), (0, 5), (0, 5), (1, 4), (2, 4)]

        for (in_s, in_e), (exp_s, exp_e) in zip(slcs, exp):
            s, e = rle_object.transform_coords(in_s, in_e)
            trimmed_input = basecall[in_s:in_e]
            trimmed_input_compressed = medaka.rle.RLEConverter(
                trimmed_input).compact_basecall
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

    @unittest.skipIf(parasail.dnafull is None, "Using fake parasail.")
    def test_compression(self):
        """Compress alignment.

        ref: TACCCATGTTGATCG  --> TACATGTGATCG
        seq:  gCCCA*GTTGATCtt -->  gCA*GTGATCt
        cigar:    1S4=1D7=2S  -->  1S2=1D6=1S
        """
        alignment = common.initialise_alignment(**self.alignment_kwargs)
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
                record = common.initialise_alignment(**basecall)
                bam.write(record)

        pysam.sort("-o", cls.bam_input, tmp_file)
        os.remove(tmp_file)
        pysam.index(cls.bam_input)

    @unittest.skipIf(parasail.dnafull is None, "Using fake parasail.")
    def test_output_rle_bam(self):
        expected = {
            ('read1', 1, 'ACTG',  '4=', (2, 1, 3, 1)),
            ('read2', 0, 'TACTG', '5=', (1, 3, 1, 3, 1))}

        # None vs full region, no difference
        for regions in (None, [medaka.common.Region('ref', 0, 10)]):
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


@unittest.skipIf(parasail.dnafull is None, "Using fake parasail.")
class RLEParamsBam(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.bam_input = tempfile.NamedTemporaryFile(suffix='.bam').name
        create_simple_bam(cls.bam_input, simple_data['calls'])

        cls.bam_output = tempfile.NamedTemporaryFile(suffix='.bam').name
        cls.ref_fname = tempfile.NamedTemporaryFile(suffix='.fasta').name

        with open(cls.ref_fname, 'w') as fasta:
            fasta.write('>ref\n')
            fasta.write('{}\n'.format(simple_data['ref']))

        # The mock fast5 file contains the RLE table
        fast5_path = mock_fast5_file()
        cls.fast5_dir = os.path.dirname(fast5_path)
        cls.fast5_fname = os.path.basename(fast5_path)

        # Create a mock summary file with read_id and filename
        cls.summary_file = mock_summary_file(cls.fast5_fname)

        regions = [medaka.common.Region('ref', 0, 10)]
        args = args_class(
            cls.bam_input, cls.bam_output, cls.ref_fname,
            2, regions, (cls.fast5_dir, cls.summary_file))
        medaka.rle.compress_bam(args)

    def test_bam_compression_with_RLE_parameters(self):
        """Test RLE parameter extraction from fast5 file."""

        for read in pysam.AlignmentFile(self.bam_output):
            for tag_name in ['WL', 'WK']:
                expected = [
                    x['tags'][tag_name] for x in simple_data['calls']
                    if x['query_name'] == read.query_name][0]

            got = read.get_tag(tag_name)
            np.testing.assert_allclose(got, expected)
