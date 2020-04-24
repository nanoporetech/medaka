import os
import tempfile
import unittest

import numpy as np
import tensorflow

import medaka.medaka as _medaka
from medaka import models
from medaka.common import Sample
from medaka.datastore import DataStore
from medaka.features import BaseFeatureEncoder
from medaka.labels import BaseLabelScheme


class TestModels(unittest.TestCase):

    def test_000_load_all_models(self):
        # only check models in package data
        # (don't check models stored in ~/.medaka)
        model_dir = _medaka.model_stores[0]
        print(model_dir)
        for name, model_file in _medaka.model_dict.items():
            model_fp = os.path.join(model_dir, model_file)
            if not os.path.exists(model_fp):
                continue
            model = models.load_model(model_fp)
            self.assertIsInstance(model, tensorflow.keras.models.Model)
            # Check we can get necessary functions for inference
            with DataStore(model_fp) as ds:
                feature_encoder = ds.get_meta('feature_encoder')
                self.assertIsInstance(feature_encoder, BaseFeatureEncoder)
                label_scheme = ds.get_meta('label_scheme')
                self.assertIsInstance(label_scheme, BaseLabelScheme)

    def test_001_default_models(self):
        for model_file in (_medaka.default_consensus_model, _medaka.default_snp_model, _medaka.default_variant_model):
            if model_file not in _medaka.model_dict:
                self.fail('Model {} not in model_dict'.format(model_file))


    def test_001_build_all_models(self):
        num_classes, time_steps, feat_len = 5, 5, 5
        for name, func in models.model_builders.items():
            model = func(feat_len, num_classes, time_steps=time_steps)


class TestMajorityModel(unittest.TestCase):

    def test_000_initialise_majority_model(self):
        majority_model = models.build_majority(10, 5)
        self.assertIsInstance(majority_model, tensorflow.keras.models.Model)
