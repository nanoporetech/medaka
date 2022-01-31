import logging
import os
import random
import tempfile
import unittest

import pysam

from medaka.align import initialise_alignment
from medaka.common import get_bam_regions, Region
from medaka.features import BAMHandler, CountsFeatureEncoder
from medaka.prediction import DataLoader


class TestDataLoader(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # create a trivial bam one read, no insertions
        # this means we know exactly how many chunks to expect
        read = ''.join(random.choices("ATCG", k=5000))
        cls.outdir = tempfile.TemporaryDirectory()
        cls.bam = os.path.join(cls.outdir.name, "reads.bam")
        with pysam.AlignmentFile(cls.bam, 'wb', reference_names=['ref'], reference_lengths=[5000]) as bam:
            bam.write(initialise_alignment('read', 0, 0, read, '5000M', 0))
        pysam.index(cls.bam)
        logging.getLogger('medaka').setLevel(logging.CRITICAL)

    def _run_one(self, batch_size, chunk_len, chunk_overlap, exp_batches, exp_samples, exp_remains=0, regions=None):
        if regions is None:
            regions = get_bam_regions(self.bam)
        bam = BAMHandler(self.bam, size=1)
        loader = DataLoader(
            bam, regions, batch_size,
            batch_cache_size=4, bam_workers=4,
            feature_encoder=CountsFeatureEncoder(),
            chunk_len=chunk_len, chunk_overlap=chunk_overlap,
            enable_chunking=True)
        batches = list(x[1] for x in loader)
        self.assertEqual(len(loader.remainders), exp_remains)
        self.assertEqual(sum(len(b) for b in batches), exp_samples)
        self.assertEqual(len(batches), exp_batches)

    def test_010_singleton(self):
        self._run_one(200, 5000, 100, exp_batches=1, exp_samples=1)

    def test_011_singleton_remainder(self):
        self._run_one(200, 10000, 100, exp_batches=0, exp_samples=0, exp_remains=1)

    def test_012_multi_remain(self):
        regions = [get_bam_regions(self.bam)[0]] * 5
        self._run_one(200, 10000, 100, exp_batches=0, exp_samples=0, exp_remains=5, regions=regions)

    def test_020_half(self):
        self._run_one(200, 2500, 0, exp_batches=1, exp_samples=2)

    def test_021_half_overlap(self):
        self._run_one(200, 2500, 100, exp_batches=1, exp_samples=3)

    def test_030_half_overlap_batch(self):
        self._run_one(2, 2500, 100, exp_batches=2, exp_samples=3)

    def test_040_lots(self):
        self._run_one(1, 5, 0, exp_batches=1000, exp_samples=1000)

    def test_050_lots_one_batch(self):
        self._run_one(5000, 1, 0, exp_batches=1, exp_samples=5000)

    def test_051_lots_nearly_one_batch(self):
        self._run_one(4999, 1, 0, exp_batches=2, exp_samples=5000)

    def test_100_convoluted(self):
        # the chunk generator will give a larger overlap in the last chunk
        # to give 4 chunks per region and no remainders for this:
        regions = [get_bam_regions(self.bam)[0]] * 5
        # but augment with 7 small subregions that will fall through
        # as remainders:
        regions += [Region.from_string("ref:0-1000")] * 7
        self._run_one(2, 1300, 0, exp_batches=10, exp_samples=20, exp_remains=7, regions=regions)

    def test_101_spam(self):
        regions = [get_bam_regions(self.bam)[0]] * 100
        regions += [Region.from_string("ref:0-100")] * 1000  # spam
        self._run_one(19, 250, 0, exp_batches=106, exp_samples=2000, exp_remains=1000, regions=regions)

