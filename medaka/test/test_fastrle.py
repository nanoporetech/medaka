"""Testing for libmedaka.fastrle program."""
import subprocess
import tempfile
import unittest 

import pysam


class RLE(unittest.TestCase):
    """Test libmedaka.fastrle function."""

    input_fasta = tempfile.NamedTemporaryFile(suffix='.fasta').name
    output_fastqrle = tempfile.NamedTemporaryFile(suffix='.fastqrle').name
    basecalls = ['ACCCCCCCGTTTA', 'CCCCCCCC']

    @classmethod
    def setUpClass(cls):
        """Create input fasta file."""
        with open(cls.input_fasta, 'w') as f:
            for index, basecall in enumerate(cls.basecalls):
                f.write('>read_{}\n{}\n'.format(index, basecall))

    def test_rle(self):
        """Test the conversion of basecalls into fastqrle file."""

        block_size = 3
        with open(self.output_fastqrle, 'w') as f:
            subprocess.call(['medaka', 'fastrle', self.input_fasta, '--block_size', str(block_size)], stdout=f)

        expected_results = (
            [('A', 1), ('C', 3), ('C', 3), ('C', 1), ('G', 1), ('T', 3), ('A', 1)], 
            [('C', 3), ('C', 3), ('C', 2)])
  
        with pysam.FastqFile(self.output_fastqrle) as f:
            for index, entry in enumerate(f):
                bases = entry.sequence
                qualities = entry.get_quality_array()
                got = list(zip(bases, qualities))
                expected = expected_results[index]
                self.assertEqual(expected, got, "Expected and got differ: ({} != {})".format(expected, got))
