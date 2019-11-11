"""Extensions to keras API for medaka."""
import os
from timeit import default_timer as now

import numpy as np
import tensorflow as tf
from tensorflow.keras.callbacks import ModelCheckpoint, TensorBoard
from tensorflow.keras.utils import Sequence

import medaka.common
import medaka.datastore

# define subclassess here to avoid top-level keras import


class ModelMetaCheckpoint(ModelCheckpoint):
    """Custom ModelCheckpoint to add medaka-specific metadata."""

    def __init__(self, medaka_meta, *args, **kwargs):
        """Initialize checkpointing.

        :param medaka_meta: dictionary of meta data to store in checkpoint
            files.
        :param args: positional arguments for baseclass.
        :param kwargs: keyword arguments for baseclass.

        """
        required_meta = set(
            ('model_function', 'label_scheme', 'feature_encoder'))
        super(ModelMetaCheckpoint, self).__init__(*args, **kwargs)
        self.medaka_meta = medaka_meta
        if not set(medaka_meta.keys()).issubset(required_meta):
            raise KeyError(
                '`medaka_meta may only contain: {}'.format(required_meta))

    def on_epoch_end(self, epoch, logs=None):
        """Perform actions at the end of an epoch."""
        super(ModelMetaCheckpoint, self).on_epoch_end(epoch, logs)
        filepath = self.filepath.format(epoch=epoch + 1, **logs)
        with medaka.datastore.DataStore(filepath, 'a') as ds:
            for k, v in self.medaka_meta.items():
                ds.set_meta(v, k)


class TrainValTensorBoard(TensorBoard):
    """Modification of tensorboard to plot test and validation together."""

    def __init__(self, log_dir='./logs', **kwargs):
        """Initialise log writing."""
        # the strategy here is to log training and validation to different
        # subdirectories and rename validation metrics to be the same as
        # the training metrics (remove their 'val' prefix)
        training_log_dir = os.path.join(log_dir, 'training')
        super(TrainValTensorBoard, self).__init__(training_log_dir, **kwargs)
        self.val_log_dir = os.path.join(log_dir, 'validation')

    def set_model(self, model):
        """Set writer for validation metrics."""
        self.val_writer = tf.summary.FileWriter(self.val_log_dir)
        super(TrainValTensorBoard, self).set_model(model)

    def on_epoch_end(self, epoch, logs=None):
        """Write logs on epoch end."""
        # take validation logs, save separately renaming keys.
        # epoch is added as this is what `on_epoch_end` does
        logs = logs or {}
        val_logs = {
            k.replace('val_', 'epoch_'): v
            for k, v in logs.items() if k.startswith('val_')}
        for name, value in val_logs.items():
            summary = tf.Summary()
            summary_value = summary.value.add()
            summary_value.simple_value = value.item()
            summary_value.tag = name
            self.val_writer.add_summary(summary, epoch)
        self.val_writer.flush()

        # take remaining logs and handle normally
        logs = {k: v for k, v in logs.items() if not k.startswith('val_')}
        super(TrainValTensorBoard, self).on_epoch_end(epoch, logs)

    def on_train_end(self, logs=None):
        """Close the validation writer on exit."""
        super(TrainValTensorBoard, self).on_train_end(logs)
        self.val_writer.close()


class SequenceBatcher(Sequence):
    """Interface for keras to a `TrainBatcher` for training and validation."""

    def __init__(self, batcher, dataset='train', mini_epochs=1, seed=None):
        """Initialize batching for training.

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
                raise ValueError(
                    "'mini_epochs' must be equal to 1 for validation data.")
        else:
            raise ValueError("'dataset' should be 'train' or 'validation'.")

        original_size = len(self.data)
        self.n_batches = len(self.data) // self.batch_size
        self.data = self.data[:self.n_batches*self.batch_size]
        np.random.shuffle(self.data)
        self.logger = medaka.common.get_named_logger(
            '{}Batcher'.format(dataset.capitalize()))
        self.logger.info(
            '{} batches of {} samples ({}), from {} original.'.format(
                self.n_batches, self.batch_size, len(self.data),
                original_size))

    def __len__(self):
        """Return the number of batches."""
        # report our length (to keras) as being the size of a mini_epoch batch
        return self.n_batches // self.mini_epochs

    def __getitem__(self, idx):
        """Return the ith batch."""
        t0 = now()
        bs = self.batch_size
        # offset into where miniepoch starts: miniepoch * samples in miniepoch
        mb_offset = (self.epoch % self.mini_epochs) * \
            (len(self.data) // self.mini_epochs)
        start, stop = mb_offset + idx * bs, mb_offset + (idx + 1) * bs
        if self.dataset == 'validation':
            self.logger.debug("Request for batch {}: [{}:{}], {}".format(
                idx, start, stop, len(self.data)))
        samples = self.data[start:stop]
        batch = self.batcher.samples_to_batch(samples)
        if self.dataset == 'validation':
            self.logger.debug(
                "Took {:5.3}s to load batch {}. (epoch {})".format(
                    now() - t0, idx, self.epoch))
        return batch

    def on_epoch_end(self):
        """Perform actions at the end of an epoch."""
        # shuffle data at end of miniepoch
        self.epoch += 1
        if self.dataset == 'train' and self.epoch % self.mini_epochs == 0:
            self.logger.debug("Shuffling data")
            np.random.shuffle(self.data)
