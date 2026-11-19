"""Regression tests for background I/O and inference failure handling."""
import concurrent.futures
import os
import queue
import tempfile
import threading
import time
import unittest
from unittest import mock

import numpy as np

from medaka.common import Region, Sample
from medaka.datastore import DataIndex, DataStore
from medaka.prediction import DataLoader, run_prediction
from medaka.torch_ext import Batch
from medaka.training import run_training


def sample():
    """A single position with a unanimous A vote."""
    return Sample(
        ref_name='ref',
        labels=None, ref_seq=None, label_probs=None, depth=None,
        positions=np.array([(0, 0)], dtype=[('major', int), ('minor', int)]),
        features=np.array([[1., 0., 0., 0., 0., 0., 0., 0., 0., 0.]]))


class TestRuntimeSafety(unittest.TestCase):
    def test_write_failure_is_reported_and_file_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            store = DataStore(os.path.join(directory, 'samples.hdf'), 'w')
            with mock.patch.object(
                    store, '_write_dataset', side_effect=OSError('disk full')):
                with self.assertRaisesRegex(OSError, 'disk full'):
                    with store:
                        store.write_sample(sample())
            self.assertFalse(store.fh.id.valid)

    def test_direct_close_waits_for_writes(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'samples.hdf')
            store = DataStore(path, 'w')
            store.write_sample(sample())
            store.close()
            with DataStore(path) as reader:
                np.testing.assert_array_equal(
                    reader.load_sample(sample().name).features,
                    sample().features)

    def test_writer_applies_backpressure(self):
        with tempfile.TemporaryDirectory() as directory:
            with DataStore(
                    os.path.join(directory, 'samples.hdf'), 'w',
                    max_pending_writes=2) as store:
                entered, release = threading.Event(), threading.Event()
                original_write = store._write_dataset

                def slow_write(*args):
                    entered.set()
                    release.wait(5)
                    original_write(*args)

                def produce():
                    for i in range(10):
                        store.write_sample(sample().amend(ref_name=f'ref{i}'))

                with mock.patch.object(store, '_write_dataset', slow_write):
                    with concurrent.futures.ThreadPoolExecutor(1) as executor:
                        future = executor.submit(produce)
                        try:
                            self.assertTrue(entered.wait(2))
                            with self.assertRaises(
                                    concurrent.futures.TimeoutError):
                                future.result(timeout=0.05)
                            self.assertLessEqual(len(store.write_futures), 2)
                        finally:
                            release.set()
                        future.result(timeout=5)
            with DataStore(store.filename) as reader:
                self.assertEqual(reader.n_samples, 10)

    def test_index_rejects_unreadable_later_file(self):
        with tempfile.TemporaryDirectory() as directory:
            valid = os.path.join(directory, 'valid.hdf')
            invalid = os.path.join(directory, 'invalid.hdf')
            with DataStore(valid, 'w') as store:
                store.write_sample(sample())
            with open(invalid, 'w') as handle:
                handle.write('not an HDF file')
            with self.assertRaisesRegex(RuntimeError, 'invalid.hdf'):
                DataIndex([valid, invalid])

    def test_collation_preserves_majority_vote_baseline(self):
        batch = Batch.collate([sample()], counts_matrix=True)
        self.assertIsNotNone(batch.majority_vote_probs)
        np.testing.assert_array_equal(
            batch.majority_vote_probs.numpy(), [[[0., 1., 0., 0., 0.]]])

    def test_inference_skips_training_baseline(self):
        batch = Batch.collate([sample()])
        self.assertIsNone(batch.majority_vote_probs)

    def test_training_without_model_selects_default(self):
        batcher = mock.Mock(feature_shape=(1, 10))
        with mock.patch(
                'medaka.models.model_from_dict',
                side_effect=RuntimeError('default selected')):
            with self.assertRaisesRegex(RuntimeError, 'default selected'):
                run_training('unused', batcher)

    def assert_loader_failure(self, loader, message):
        """Ensure a failure reaches the consumer within a bounded time."""
        result = queue.Queue()

        def consume():
            try:
                list(loader)
            except Exception as exc:
                result.put(exc)
            else:
                result.put(None)

        consumer = threading.Thread(target=consume, daemon=True)
        consumer.start()
        try:
            try:
                error = result.get(timeout=2)
            except queue.Empty:
                self.fail('inference consumer hung after worker failure')
            self.assertIsInstance(error, ValueError)
            self.assertIn(message, str(error))
        finally:
            if hasattr(loader, 'close'):
                loader.close()
            else:
                # Unblock the original implementation during the red run.
                loader._samples.put_nowait(StopIteration)
                loader._batches.put_nowait(StopIteration)
            consumer.join(timeout=2)

    def test_region_failure_reaches_consumer(self):
        with mock.patch.object(
                DataLoader, '_run_region',
                side_effect=ValueError('bad region')):
            loader = DataLoader(None, ['region'], 1, bam_workers=1)
            self.assert_loader_failure(loader, 'bad region')

    def test_collation_failure_reaches_consumer(self):
        with mock.patch.object(
                DataLoader, '_run_region', return_value=([sample()], [])), \
                mock.patch.object(
                    Batch, 'collate', side_effect=ValueError('bad batch')):
            loader = DataLoader(None, ['region'], 1, bam_workers=1)
            self.assert_loader_failure(loader, 'bad batch')

    def test_close_releases_workers_with_full_queues(self):
        with mock.patch.object(
                DataLoader, '_run_region',
                return_value=([sample()] * 30, [])):
            loader = DataLoader(
                None, ['region'], 1, batch_cache_size=1, bam_workers=2)
            try:
                deadline = time.monotonic() + 2
                while not (loader._batches.full() and loader._samples.full()):
                    if time.monotonic() >= deadline:
                        self.fail('loader did not fill its bounded queues')
                    time.sleep(0.01)
            finally:
                loader.close()
            self.assertFalse(loader._bthread.is_alive())
            self.assertTrue(all(
                not worker.is_alive() for worker in loader._sample_workers))

    def test_prediction_failure_closes_loader(self):
        with tempfile.TemporaryDirectory() as directory, mock.patch.object(
                DataLoader, '_run_region',
                return_value=([sample()] * 30, [])):
            loader = DataLoader(
                None, ['region'], 1, batch_cache_size=1, bam_workers=1)
            model = mock.Mock()
            model.predict_on_batch.side_effect = ValueError(
                'prediction failed')
            try:
                with mock.patch(
                        'medaka.prediction.DataLoader', return_value=loader):
                    with self.assertRaisesRegex(
                            ValueError, 'prediction failed'):
                        run_prediction(
                            os.path.join(directory, 'output.hdf'), None,
                            [Region('ref', 0, 1)], model, None, 1, 0,
                            batch_size=1)
                self.assertFalse(loader._bthread.is_alive())
                self.assertFalse(loader._sample_workers[0].is_alive())
            finally:
                loader.close()

    def test_exhausted_loader_stays_exhausted(self):
        loader = DataLoader(None, [], 1, bam_workers=2)
        self.assertEqual(list(loader), [])
        result = queue.Queue()
        consumer = threading.Thread(
            target=lambda: result.put(list(loader)), daemon=True)
        consumer.start()
        try:
            self.assertEqual(result.get(timeout=2), [])
        finally:
            if hasattr(loader, 'close'):
                loader.close()
            else:
                loader._batches.put_nowait(StopIteration)
            consumer.join(timeout=2)
