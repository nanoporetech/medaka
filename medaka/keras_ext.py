from timeit import default_timer as now

import numpy as np
from keras.callbacks import ModelCheckpoint 
from keras.utils import Sequence

from medaka.common import get_named_logger
from medaka.datastore import DataStore

# define subclassess here to avoid top-level keras import

class ModelMetaCheckpoint(ModelCheckpoint):
    """Custom ModelCheckpoint to add medaka-specific metadata to model files"""
    def __init__(self, medaka_meta, *args, **kwargs):
        super(ModelMetaCheckpoint, self).__init__(*args, **kwargs)
        self.medaka_meta = medaka_meta

    def on_epoch_end(self, epoch, logs=None):
        super(ModelMetaCheckpoint, self).on_epoch_end(epoch, logs)
        filepath = self.filepath.format(epoch=epoch + 1, **logs)
        with DataStore(filepath, 'a') as ds:
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
        if seed is not None:
            np.random.seed(seed)
        self.epoch = 1
        if dataset == 'train':
            self.data = batcher.train_samples
            self.n_batches = self.batcher.n_train_batches
        elif dataset == 'validation':
            self.data = batcher.valid_samples
            self.n_batches = self.batcher.n_valid_batches
            if mini_epochs != 1:
                raise ValueError("'mini_epochs' must be equal to 1 for validation data.")
        else:
            raise ValueError("'dataset' should be 'train' or 'validation'.")
        self.logger = get_named_logger('{}Batcher'.format(dataset.capitalize()))


    def __len__(self):
        return self.n_batches // self.mini_epochs


    def __getitem__(self, idx):
        t0 = now()
        bs = self.batcher.batch_size
        samples = self.data[idx * bs:(idx + 1) * bs]
        batch = self.batcher.samples_to_batch(samples)
        self.logger.debug("Took {:5.3}s to load batch {}. (epoch {})".format(now() - t0, idx, self.epoch))
        return batch


    def on_epoch_end(self):
        # shuffle data (keras only shuffles batches)
        # TODO: respect mini_epochs to cycle through all data
        self.epoch += 1
        np.random.shuffle(self.data)
