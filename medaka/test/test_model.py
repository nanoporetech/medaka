from concurrent.futures import ThreadPoolExecutor
import os
import tempfile
import unittest

import numpy as np
import torch

from medaka import models
from medaka.common import Sample
from medaka.features import BaseFeatureEncoder
from medaka.labels import BaseLabelScheme
import medaka.options
import medaka.models


class TestModelFiles(unittest.TestCase):

    def test_000_total_bundled_size(self):
        total = 0
        for name in medaka.options.current_models:
            model_file = models.resolve_model(name)
            total += os.path.getsize(model_file)
        self.assertLess(total / 1024 / 1024, 45, "Bundled model file size too large")

    def test_001_default_models(self):
        for name in medaka.options.default_models.values():
            if name not in medaka.options.current_models:
                self.fail(f'Default Model {name} not in current_models')

    def test_010_failed_download(self):
        name = 'garbage'
        medaka.options.known_models.append(name)
        with self.assertRaises(medaka.models.DownloadError):
            models.resolve_model(name)
        medaka.options.known_models.pop()

    def test_011_success_download(self):
        name = 'r1041_e82_400bps_hac_v5.0.0'
        # Temporarily add in old suffix format so we can test
        # download that already exists in lfs
        medaka.models.model_suffixes += ["_model.tar.gz",]
        model_file = models.resolve_model(name)
        tmp_file = "{}.tmp".format(model_file)
        os.rename(model_file, tmp_file)
        new_file = models.resolve_model(name)
        self.assertTrue(os.path.isfile(new_file))
        os.remove(new_file)
        os.rename(tmp_file, model_file)

    def test_020_basecaller_model(self):
        name = next(iter(medaka.options.basecaller_models.keys()))
        for variety in ('consensus', 'variant'):
            try:
                model = models.resolve_model(f"{name}:{variety}")
            except medaka.options.DeprecationError:
                continue

    def test_030_resolve_deprecated(self):
        name = 'r941_min_high_g351'
        with self.assertRaises(medaka.options.DeprecationError):
            models.resolve_model(name)

    def test_999_load_all_models(self):
        with ThreadPoolExecutor(max_workers=4) as executor:
            for name in medaka.options.allowed_models:
                model_file = models.resolve_model(name)
                executor.submit(self._load_one, args=[self, model_file])

    def _load_one(self, model_file):
        print(f"loading: {model_file}")
        with medaka.models.open_model(model_file) as ds:
            model = ds.load_model()
            self.assertIsInstance(model, torch.nn.Module)
            feature_encoder = ds.get_meta('feature_encoder')
            self.assertIsInstance(feature_encoder, BaseFeatureEncoder)
            label_scheme = ds.get_meta('label_scheme')
            self.assertIsInstance(label_scheme, BaseLabelScheme)


class TestScrapBasecaller(unittest.TestCase):

    root_dir = os.path.abspath(os.path.dirname(__file__))
    bam = os.path.join(root_dir, 'data/bc_model_scrape.bam')
    fastq = os.path.join(root_dir, 'data/bc_model_scrape.fastq.gz')
    fastq_minknow = os.path.join(root_dir, 'data/bc_model_scrape_minknow.fastq.gz')

    def test_000_from_bam_consensus(self):
        model = models.model_from_basecaller(self.bam, variant=False)
        self.assertEqual(model, "r1041_e82_400bps_hac_v4.2.0")

    def test_001_from_bam_variant(self):
        model = models.model_from_basecaller(self.bam, variant=True)
        self.assertEqual(model, "r1041_e82_400bps_hac_variant_v4.2.0")

    def test_010_from_fastq_consensus(self):
        model = models.model_from_basecaller(self.fastq, variant=False)
        self.assertEqual(model, "r1041_e82_400bps_hac_v4.2.0")

    def test_011_from_fastq_variant(self):
        model = models.model_from_basecaller(self.fastq, variant=True)
        self.assertEqual(model, "r1041_e82_400bps_hac_variant_v4.2.0")

    def test_020_from_fastq_minknow(self):
        model = models.model_from_basecaller(self.fastq_minknow, variant=False)
        self.assertEqual(model, "r1041_e82_400bps_sup_v4.2.0")

    def test_021_from_fastq_minknow_variant(self):
        model = models.model_from_basecaller(self.fastq_minknow, variant=True)
        self.assertEqual(model, "r1041_e82_400bps_sup_variant_v4.2.0")

    def test_030_from_bam_consensus_bacteria(self):
        model = models.model_from_basecaller(self.bam, variant=False,
                                             bacteria=True)
        self.assertEqual(model, "r1041_e82_400bps_bacterial_methylation")

    def test_031_from_fastq_consensus_bacteria(self):
        model = models.model_from_basecaller(self.fastq, variant=False,
                                             bacteria=True)
        self.assertEqual(model, "r1041_e82_400bps_bacterial_methylation")
        

class TestBuildModel(unittest.TestCase):

    def test_000_build_all_models(self):
        num_classes, time_steps, feat_len = 5, 5, 5
        for name, func in models.model_builders.items():
            model = func(feat_len, num_classes, time_steps=time_steps)


class TestMajorityModel(unittest.TestCase):

    def test_000_initialise_majority_model(self):
        in_features, out_labels = 10, 5
        majority_model = models.build_majority(in_features, out_labels)
        self.assertIsInstance(majority_model, torch.nn.Module)

        batch_size = 4
        positions = 100
        test_data = torch.rand((batch_size, positions, in_features), dtype=float)

        result = majority_model.predict_on_batch(test_data)
        expected_size = (batch_size, positions, out_labels)
        self.assertSequenceEqual(list(result.shape), expected_size)
