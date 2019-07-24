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
        self.epoch = 1

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
        return self.n_batches // self.mini_epochs


    def __getitem__(self, idx):
        t0 = now()
        bs = self.batch_size
        start, stop = idx * bs, (idx + 1) * bs
        if self.dataset == 'validation':
            self.logger.debug("Request for batch {}: [{}:{}], {}".format(idx, start, stop, len(self.data)))
        samples = self.data[start:stop]
        batch = self.batcher.samples_to_batch(samples)
        if self.dataset == 'validation':
            self.logger.debug("Took {:5.3}s to load batch {}. (epoch {})".format(now() - t0, idx, self.epoch))
        return batch


    def on_epoch_end(self):
        # shuffle data (keras only shuffles batches)
        # TODO: respect mini_epochs to cycle through all data
        self.epoch += 1
        np.random.shuffle(self.data)


class BatchQueue(object):
    def  __init__(self, samples, prep_func, batch_size, executor, seed=None, name='Train', maxsize=100):
        """Load and queue training samples into batches from `.hdf` files.

        :param samples: tuples of (filename, hdf sample key).
        :param prep_func: function to transform a sample to x,y data.
        :param batch_size: group samples by this number.
        :param executor: `ThreadPoolExecutor` instance.
        :param seed: seed for shuffling.
        :param name: str, name for logger.
        :param maxsize: int, maximum queue size.

        Once initialized batches can be retrieved using batch_q._queue.get().

        """
        self.samples = samples
        self.prep_func = prep_func
        self.batch_size = batch_size

        if seed is not None:
            np.random.seed(seed)

        self.name = name
        self.logger = medaka.common.get_named_logger('{}Batcher'.format(name.capitalize()))
        self.maxsize = maxsize
        self._queue = queue.Queue(maxsize=self.maxsize)
        self.executor = executor
        self.stopped = threading.Event()
        self.qthread = threading.Thread(target=self._fill_queue_batch)
        self.qthread.daemon = True

        original_size = len(self.samples)
        self.n_batches = len(self.samples) // self.batch_size
        self.samples = self.samples[:self.n_batches * self.batch_size]
        self.logger.info(
            '{} batches of {} samples ({}), from {} original.'.format(
            self.n_batches, self.batch_size, len(self.samples), original_size
        ))
        if self.n_batches == 0:
            raise ValueError("Number of batches is zero.")

        self.qthread.start()
        time.sleep(2)
        self.logger.info("Started reading samples from files with queue size {}".format(maxsize))


    def stop(self, timeout=5):
        self.logger.info("About to stop.")
        self.stopped.set()
        self.logger.info("Waiting for read thread.")
        self.qthread.join(timeout)
        if self.qthread.is_alive:
            self.logger.critical("Read thread did not terminate after {}s.".format(timeout))


    @staticmethod
    def samples_to_batch(samples, prep_func, name, batch, epoch):
        t0 = now()
        items = [prep_func(s) for s in samples]
        xs, ys = zip(*items)
        x, y = np.stack(xs), np.stack(ys)
        medaka.common.get_named_logger(name).debug("Took {:5.3}s to load batch {} (epoch {})".format(now()-t0, batch, epoch))
        return x, y


    def _fill_queue_batch(self):
        epoch = 0
        self.loaded_batches = 0
        self.submitted_batches = 0
        self.taken_batches = 0
        while not self.stopped.is_set():
            batch = 0
            np.random.shuffle(self.samples)
            for samples in medaka.common.grouper(iter(self.samples), batch_size=self.batch_size):
                if self.stopped.is_set(): # the loop is potentially long-running
                    self.logger.info("Batching stopped.")
                    return
                res = self.executor.submit(self.samples_to_batch, samples, self.prep_func, self.name, batch, epoch)
                res.add_done_callback(self._count_finished)
                while True: # keep an eye on the stopped flag
                    try:
                        self._queue.put(res, timeout=1)
                    except queue.Full:
                        #self.logger.debug("Queue is full ({}), cannot put, trying again.".format(self._queue.qsize()))
                        if self.stopped.is_set():
                            self.logger.info("Batching stopped.")
                            return
                    else:
                        #self.logger.debug("Successfully put sample (future).")
                        self.submitted_batches += 1
                        break
                batch += 1
            epoch += 1
        self.logger.info("Batching stopped.")
        return

    def _count_finished(self, future):
        self.loaded_batches += 1


    #@medaka.common.threadsafe_generator
    def yield_batches(self):
        time_between = deque(maxlen=50)
        get_time = deque(maxlen=50)
        t0 = now()
        try:
            while True:
                t0, t1 = now(), t0
                ta = now()
                res = self._queue.get()
                if isinstance(res, Future):
                    res =  res.result()
                get_time.append(now() - ta)
                time_between.append(t0 - t1)
                get_rate = np.mean(get_time)
                req_rate = np.mean(time_between) - get_rate
                self.logger.debug(
                    "Request every: {:5.3}s. Fetch time: {:5.3}.".format(
                    np.mean(time_between), np.mean(get_time)
                ))
                self.logger.debug("Queue state: {}/{} ready.".format(
                    self.loaded_batches - self.taken_batches,
                    self.submitted_batches - self.taken_batches
                ))
                self.taken_batches += 1
                yield res
        except Exception as e:
            self.logger.critical("Exception caught why yielding batches: {}".format(e))
            self.stop()
            raise e


