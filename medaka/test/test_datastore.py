import tempfile
import unittest


import numpy as np
from medaka.common import Sample
from medaka import datastore


class TestStore(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pos = np.array([(0, 0), (0, 1), (1, 0), (2, 0), (2, 1), (2, 2), (3, 0), (4, 0), (4, 1), (4, 2), (4, 3)],
                        dtype=[('major', int), ('minor', int)])
        data_dim = 10
        data = np.zeros(shape=(len(pos), data_dim))
        cls.sample = Sample(
            ref_name='contig1', features=data, ref_seq=None,
            labels=data, positions=pos, label_probs=data)
        cls.file = tempfile.NamedTemporaryFile()
        with datastore.DataStore(cls.file.name, 'w') as store:
            store.write_sample(cls.sample)


    def test_000_round_trip(self):
        with datastore.DataStore(self.file.name, 'r') as store:
            grp = store.fh['samples/data']
            name = self.sample.name
            self.assertListEqual(list(grp.keys()), [name], "Stored one sample")
            sample = store.load_sample(self.sample.name)
            # Check name (a simple string) and features (numpy data)
            self.assertEqual(sample.name, self.sample.name, "Loaded sample has correct name.")
            self.assertIsInstance(sample.features, np.ndarray, "Sample features is a numpy array.")
            self.assertSequenceEqual(sample.features.shape, self.sample.features.shape, "Sample features were loaded.")


    def test_001_round_trip_through_index(self):
        index = datastore.DataIndex([self.file.name])
        samples = list(index.yield_from_feature_files())
        self.assertEqual(len(samples), 1, "Retrieve one sample")
        self.assertEqual(samples[0].name, self.sample.name, "Loaded sample has correct name")
        self.assertIsInstance(samples[0].features, np.ndarray, "Sample features is a numpy array.")
        self.assertSequenceEqual(samples[0].features.shape, self.sample.features.shape, "Sample features were loaded.")
