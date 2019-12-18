import tempfile
import unittest

import numpy as np
import tensorflow

from medaka import models
from medaka.medaka import model_dict
from medaka.common import Sample
from medaka.datastore import DataStore
from medaka.features import BaseFeatureEncoder
from medaka.labels import BaseLabelScheme


class TestModels(unittest.TestCase):

    def test_000_load_all_models(self):
        for name, model_file in model_dict.items():
            model = models.load_model(model_file)
            self.assertIsInstance(model, tensorflow.keras.models.Model)
            # Check we can get necessary functions for inference
            with DataStore(model_file) as ds:
                feature_encoder = ds.get_meta('feature_encoder')
                self.assertIsInstance(feature_encoder, BaseFeatureEncoder)
                label_scheme = ds.get_meta('label_scheme')
                self.assertIsInstance(label_scheme, BaseLabelScheme)


class TestMajorityModel(unittest.TestCase):

    def test_000_initialise_majority_model(self):
        majority_model = models.build_majority(10, 5)
        self.assertIsInstance(majority_model, tensorflow.keras.models.Model)
