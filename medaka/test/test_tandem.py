import unittest
import os
import subprocess
import tempfile
import medaka.common
import medaka.tandem.tandem
import medaka.medaka
import pysam
import gzip

def decompress_gz_file(input_path, output_path):
    """
    Decompress a gzipped file to a specified output file using the gzip command.

    Args:
        input_path (str): Path to the gzipped file.
        output_path (str): Path to save the decompressed file.

    Returns:
        None
    """
    try:
        with gzip.open(input_path, 'rb') as gz_in, open(output_path, 'wb') as out_f:
          out_f.write(gz_in.read())

        print(f"Decompressed file saved to {output_path}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while decompressing: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

class TestTandem(unittest.TestCase):
   @classmethod
   def setUpClass(cls):
      cls.root_dir = os.path.abspath(os.path.dirname(__file__))
      cls.test_bam = os.path.join(cls.root_dir, 'data', 'tandem', 'test.bam')
      cls.test_fasta_compressed = os.path.join(cls.root_dir, 'data', 'tandem', 'chr20.fa.gz')
      cls.test_bed = os.path.join(cls.root_dir, 'data', 'tandem', 'test.bed')
      cls.test_dup_bed = os.path.join(cls.root_dir, 'data', 'tandem', 'test_duplicates.bed')


   def setUp(self):
      self.temp_folder = os.path.join(self.root_dir, 'temp', tempfile.mkdtemp())
      os.makedirs(self.temp_folder, exist_ok=True)

      self.test_fasta = os.path.join(self.temp_folder, 'chr20.fa')
      decompress_gz_file(self.test_fasta_compressed, self.test_fasta)

   def run_tandem(self, args, bed=None):
      bed = self.test_bed if bed is None else bed
      args += [
         self.test_bam,
         self.test_fasta,
         bed,
         "male",
         self.temp_folder,
         "--ignore_read_groups"
      ]
      parser = medaka.medaka.medaka_parser()
      medaka.tandem.tandem.main(parser.parse_args(args))

   def tearDown(self):
      # Clean up temporary files after all tests
      medaka.common.remove_directory(self.temp_folder)

   @staticmethod
   def parse_fasta(file_path):
      fasta_dict = {}
      with open(file_path, 'r') as f:
         header = None
         for line in f:
            line = line.strip()
            if line.startswith(">"):
               header = line[1:]  # Remove the '>' character
               fasta_dict[header] = ""
            elif header:
               fasta_dict[header] += line  # Append sequence lines
      return fasta_dict

   def assertFastaEqual(self, truth_fasta, test_fasta):
      truth_dict = self.parse_fasta(truth_fasta)
      test_dict = self.parse_fasta(test_fasta)
      self.assertEqual(truth_dict.keys(), test_dict.keys(), f"Fasta files are not equal {truth_fasta}")
      for key in truth_dict.keys():
         self.assertEqual(truth_dict[key], test_dict[key], f"Fasta files are not equal {truth_fasta}") 

   def assertFilesEqual(self,file1, file2, sort_files=False):
      with open(file1, 'r') as f1, open(file2, 'r') as f2:
         lines1 = f1.readlines()
         lines2 = f2.readlines()

         if sort_files:
               lines1 = sorted(lines1)
               lines2 = sorted(lines2)

         try:
            assert len(lines1) == len(lines2), f"Files {file1} and {file2} are not equal"
            for line1, line2 in zip(lines1, lines2):
               self.assertEqual(line1.strip(), line2.strip())
         except ValueError:
            self.fail(f"Files {file1} and {file2} are not equal")
            
   def assertBAMEqual(self, bam1, bam2):
         """
         Compare two BAM files for header, read count, and read content equality.
         """
         # Open BAM files
         with pysam.AlignmentFile(bam1, "rb") as bam1_file, pysam.AlignmentFile(bam2, "rb") as bam2_file:
            header1 = bam1_file.header.to_dict()
            header2 = bam2_file.header.to_dict()

            # Remove @PG tags if they exist
            header1.pop("PG", None)
            header2.pop("PG", None)


            self.assertEqual(header1.keys(), header2.keys(), f"Header keys are not equal in {bam1} and {bam2}")
            for key in header1.keys():
                  items_1 = header1[key]
                  items_2 = header2[key]
                  if key == 'SQ':
                     items_1 = sorted([(item['SN'], item['LN']) for item in items_1])
                     items_2 = sorted([(item['SN'], item['LN']) for item in items_2])
                  self.assertEqual(items_1, items_2, f"Header values are not equal for key: {key} in {bam1} and {bam2}")

            # Extract and sort reads by query name
            reads1 = sorted(list(bam1_file), key=lambda read: read.query_name)
            reads2 = sorted(list(bam2_file), key=lambda read: read.query_name)

            # Compare read counts
            self.assertEqual(len(reads1), len(reads2), "Read counts are not equal")

            # Compare each read
            for read1, read2 in zip(reads1, reads2):
                  self.assertEqual(
                     read1.to_dict(), 
                     read2.to_dict(), 
                     f"Reads do not match: {read1} vs {read2}"
                  )
   def check_files(self, truth_folder):
      fastas = ["consensus.fasta", "poa.fasta"]
      for fasta in fastas:
         truth_fasta = os.path.join(os.path.dirname(__file__), 'data', 'tandem', 'truth', truth_folder, fasta)
         test_fasta = os.path.join(self.temp_folder, fasta)
         self.assertFastaEqual(truth_fasta, test_fasta)

      beds = ["skipped.bed", "skipped_large.bed"]
      for bed in beds:
         truth_bed = os.path.join(os.path.dirname(__file__), 'data', 'tandem', 'truth', truth_folder, bed)
         test_bed = os.path.join(self.temp_folder, bed)
         self.assertFilesEqual(truth_bed, test_bed, sort_files=True)

      truth_vcf = os.path.join(os.path.dirname(__file__), 'data', 'tandem', 'truth', truth_folder, 'medaka_to_ref.TR.vcf')
      test_vcf = os.path.join(self.temp_folder, 'medaka_to_ref.TR.vcf')
      self.assertFilesEqual(truth_vcf, test_vcf)

      metrics = ["abpoa_region_metrics.txt", "prephased_region_metrics.txt", "unphased_region_metrics.txt"]
      for metric in metrics:
         truth_metrics = os.path.join(os.path.dirname(__file__), 'data', 'tandem', 'truth', truth_folder, metric)
         test_metrics = os.path.join(self.temp_folder, metric)
         self.assertFilesEqual(truth_metrics, test_metrics, sort_files=True)

      bams = ["medaka_to_ref.bam", "trimmed_reads_to_poa.bam"]
      for bam in bams:
         truth_bam = os.path.join(os.path.dirname(__file__), 'data', 'tandem', 'truth', truth_folder, bam)
         test_bam = os.path.join(self.temp_folder, bam)
         self.assertBAMEqual(truth_bam, test_bam)

   def test_tandem_single_thread(self):
      truth = "replace"  # Test name
      args = ["tandem", "--model", "r1041_e82_400bps_sup_v5.0.0", "--workers", "1"]
      self.run_tandem(args)
      self.check_files(truth)

   def test_tandem_hybrid_replace_style(self):
      truth = "replace"  # Test name
      args = ["tandem", "--model", "r1041_e82_400bps_sup_v5.0.0", "--workers", "4"]
      self.run_tandem(args)
      self.check_files(truth)

   def test_tandem_hybrid_decompose_style(self):
      truth = "decompose"
      args = ["tandem", "--model", "r1041_e82_400bps_sup_v5.0.0", "--workers", "4", "--decompose"]
      self.run_tandem(args)
      self.check_files(truth)

   def test_tandem_hybrid_decompose_style_single_thread(self):
      truth = "decompose"
      args = ["tandem", "--model", "r1041_e82_400bps_sup_v5.0.0", "--workers", "1", "--decompose"]
      self.run_tandem(args)
      self.check_files(truth)

   def test_tandem_hybrid_unphased_with_readnames_single_thread(self):
      truth = "unphased"
      args = ["tandem", "--model", "r1041_e82_400bps_sup_v5.0.0", "--workers", "1", "--add_read_names", "--phasing", "unphased"]
      self.run_tandem(args)
      self.check_files(truth)

   def test_tandem_hybrid_unphased_with_readnames(self):
      truth = "unphased"
      args = ["tandem", "--model", "r1041_e82_400bps_sup_v5.0.0", "--workers", "4", "--add_read_names", "--phasing", "unphased"]
      self.run_tandem(args)
      self.check_files(truth)

   def test_duplicated_bed_file(self):
      truth = "replace"  # Test name
      args = ["tandem", "--model", "r1041_e82_400bps_sup_v5.0.0", "--workers", "4"]
      self.run_tandem(args, bed=self.test_dup_bed)
      self.check_files(truth)

class TestCheckReadLevelModel(unittest.TestCase):
   def setUp(self):
      self.valid_models = [
         "r1041_e82_400bps_sup_v5.0.0_rl_lstm384_dwells_model_pt.tar.gz",
         "r1041_e82_400bps_sup_v5.0.0_rl_lstm384_no_dwells_model_pt.tar.gz",
         "/fake/path/to/r1041_e82_400bps_sup_v5.0.0_rl_lstm384_dwells_model_pt.tar.gz",
         "/fake/path/to/r1041_e82_400bps_sup_v5.0.0_rl_lstm384_no_dwells_model_pt.tar.gz"
      ]
      self.invalid_models = [
         "random_string_model.tar.gz",
         "",
         None,
         "/fake/path/to/r1041_e82_400bps_hac_v5.0.0_model_pt.tar.gz",
         "/fake/path/to/r1041_e82_400bps_sup_v5.0.0_model_pt.tar.gz"
      ]

   def test_check_read_level_model(self):
      for model in self.valid_models:
         with self.subTest(model=model):
            self.assertTrue(medaka.tandem.tandem.check_read_level_model(model))

      for model in self.invalid_models:
         with self.subTest(model=model):
            self.assertFalse(medaka.tandem.tandem.check_read_level_model(model))

   