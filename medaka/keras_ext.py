from collections import Counter, deque
from concurrent.futures import Future
import functools
import queue
import threading
import time
from timeit import default_timer as now

import numpy as np
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.utils import Sequence

import medaka.common
import medaka.datastore

# define subclassess here to avoid top-level keras import

class ModelMetaCheckpoint(ModelCheckpoint):

    def __init__(self, medaka_meta, *args, **kwargs):
        """Custom ModelCheckpoint to add medaka-specific metadata to model files."""
        super(ModelMetaCheckpoint, self).__init__(*args, **kwargs)
        self.medaka_meta = medaka_meta

    def on_epoch_end(self, epoch, logs=None):
        super(ModelMetaCheckpoint, self).on_epoch_end(epoch, logs)
        filepath = self.filepath.format(epoch=epoch + 1, **logs)
        with medaka.datastore.DataStore(filepath, 'a') as ds:
            ds.meta.update(self.medaka_meta)


class SequenceBatcher(Sequence):
    def __init__(self, batcher, dataset='train', mini_epochs=1, seed=None):
        """Interface for keras to a `TrainBatcher` for training and validation
        batches.

        :param batcher: a `medaka.inference.TrainBatcher` instance.
        :param dataset: one of 'train' or 'validation'.
        :param mini_epochs: factor by which to rescale the number of batches
            in an epoch (useful to output checkpoints more frequently).
        :param seed: random seed for shuffling data.

        """
        self.batcher = batcher
        self.dataset = dataset
        self.mini_epochs = mini_epochs
        self.batch_size = self.batcher.batch_size
        if seed is not None:
            np.random.seed(seed)
        self.epoch = 0

        if dataset == 'train':
            self.data = batcher.train_samples
        elif dataset == 'validation':
            self.data = batcher.valid_samples
            if mini_epochs != 1:
                raise ValueError("'mini_epochs' must be equal to 1 for validation data.")
        else:
            raise ValueError("'dataset' should be 'train' or 'validation'.")

        original_size = len(self.data)
        self.n_batches = len(self.data) // self.batch_size
        self.data = self.data[:self.n_batches*self.batch_size]
        np.random.shuffle(self.data)
        self.logger = medaka.common.get_named_logger('{}Batcher'.format(dataset.capitalize()))
        self.logger.info(
            '{} batches of {} samples ({}), from {} original.'.format(
            self.n_batches, self.batch_size, len(self.data), original_size
        ))


    def __len__(self):
        # report our length (to keras) as being the size of a mini_epoch batch
        return self.n_batches // self.mini_epochs


    def __getitem__(self, idx):
        t0 = now()
        bs = self.batch_size
        # offset into where miniepoch starts: miniepoch * samples in miniepoch
        mb_offset = (self.epoch % self.mini_epochs) * (len(self.data) // self.mini_epochs)
        start, stop = mb_offset + idx * bs, mb_offset + (idx + 1) * bs
        if self.dataset == 'validation':
            self.logger.debug("Request for batch {}: [{}:{}], {}".format(idx, start, stop, len(self.data)))
        samples = self.data[start:stop]
        batch = self.batcher.samples_to_batch(samples)
        if self.dataset == 'validation':
            self.logger.debug("Took {:5.3}s to load batch {}. (epoch {})".format(now() - t0, idx, self.epoch))
        return batch


    def on_epoch_end(self):
        # shuffle data at end of miniepoch
        self.epoch += 1
        if self.dataset == 'train' and self.epoch % self.mini_epochs == 0:
            self.logger.debug("Shuffling data")
            np.random.shuffle(self.data)


