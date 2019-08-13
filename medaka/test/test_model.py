import tempfile
import unittest

import numpy as np
import tensorflow

from medaka.medaka import model_dict
from medaka.common import Sample
from medaka import models


class TestModels(unittest.TestCase):

    def test_000_load_all_models(self):
        for name, model_file in model_dict.items():
            model = models.load_model(model_file)
            self.assertIsInstance(model, tensorflow.keras.models.Model)
