import os
import tempfile
import unittest

import numpy as np

from medaka.common import Region, Sample
from medaka import datastore


class TestStore(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pos = np.array([(0, 0), (0, 1), (1, 0), (2, 0), (2, 1), (2, 2), (3, 0), (4, 0), (4, 1), (4, 2), (4, 3)],
                        dtype=[('major', int), ('minor', int)])
        data_dim = 10
        depth = np.array(len(pos) * [10], dtype=int)
        data = np.zeros(shape=(len(pos), data_dim))
        cls.sample = Sample(
            ref_name='contig1', features=data, ref_seq=None,
            labels=data, positions=pos, label_probs=data, depth=depth)
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
            self.assertSequenceEqual(list(sample.depth), list(self.sample.depth), "Sample depth was loaded.")


    def test_001_round_trip_through_index(self):
        index = datastore.DataIndex([self.file.name])
        samples = list(index.yield_from_feature_files())
        self.assertEqual(len(samples), 1, "Retrieve one sample")
        self.assertEqual(samples[0].name, self.sample.name, "Loaded sample has correct name")
        self.assertIsInstance(samples[0].features, np.ndarray, "Sample features is a numpy array.")
        self.assertSequenceEqual(samples[0].features.shape, self.sample.features.shape, "Sample features were loaded.")
        self.assertSequenceEqual(list(samples[0].depth), list(self.sample.depth), "Sample depth was loaded.")


    def test_002_test_filtered_yield(self):
        index = datastore.DataIndex([self.file.name])
        region_specs = [
            (None, 1),
            ([Region('contig1', 0, 5)], 1),
            ([Region('contig1', 5, 10)], 0),
            ([Region('contig2', None, None)], 0),
        ]
        for regs, exp_len in region_specs:
            samples = list(index.yield_from_feature_files(regions=regs))
            self.assertEqual(len(samples), exp_len)
