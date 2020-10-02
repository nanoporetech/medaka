import os
import tempfile
import unittest

import pysam

import medaka.align
from medaka import medaka as medaka_parsers
import medaka.wrappers


class HaploidVariant(unittest.TestCase):
    """HaploidVariant integration test."""

    __ref_fasta__ = os.path.join(os.path.dirname(__file__), 'data', 'draft_ref.fasta')
    __reads_bam__ = os.path.join(os.path.dirname(__file__), 'data', 'test_reads.bam')
    __reads_fasta__ = tempfile.NamedTemporaryFile(suffix='.fasta')
    __output_dir__ = tempfile.TemporaryDirectory()

    @classmethod
    def setUpClass(cls):
        # convert bam to reads
        pysam.fasta(cls.__reads_bam__, save_stdout=cls.__reads_fasta__.name)


    def test_run(self):
        """Test that the pipeline runs and that we get a vcf file containing at least one variant.
        """
        from tensorflow.python.eager import context
        context._context = None
        #context._create_context()

        parser = medaka_parsers._haploid_variant_argparser()
        args = [self.__reads_fasta__.name, self.__ref_fasta__, '-o', self.__output_dir__.name]
        args = parser.parse_args(args)
        medaka_parsers._validate_common_args(args)
        medaka.wrappers.haploid_variant(args)
        vcf_fp = os.path.join(args.output_dir, 'consensus_to_ref.vcf')
        vcf = medaka.vcf.VCFReader(vcf_fp, cache=False)
        variant = next(vcf.fetch())
        self.assertEqual(variant.chrom, 'utg000001l')


