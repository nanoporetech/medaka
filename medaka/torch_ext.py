"""Extensions to torch API for medaka."""
from abc import ABC, abstractmethod
import math
import os
import pickle
import shutil
import tarfile
from time import perf_counter

import numpy as np
import torch
from tqdm import tqdm

import medaka.common
import medaka.datastore

# define subclassess here to avoid top-level torch import


# Subclass to save model and metadata

class ModelMetaCheckpoint:
    """Custom ModelCheckpoint to add medaka-specific metadata to hdf5 files."""

    def __init__(self, medaka_meta, outdir, *args, **kwargs):
        """Initialize checkpointing.

        :param medaka_meta: dict of meta data to store in checkpoint files.
        :param outdir: directory in which to output saved models
        :param args: positional arguments for baseclass.
        :param kwargs: keyword arguments for baseclass.

        """
        required_meta = set(
            ('model_function', 'label_scheme', 'feature_encoder'))
        self.medaka_meta = medaka_meta
        if not set(medaka_meta.keys()).issubset(required_meta):
            raise KeyError(
                '`medaka_meta may only contain: {}'.format(required_meta))
        self.filepath = os.path.join(outdir, "model-{epoch}")
        self.logger = medaka.common.get_named_logger("CheckPnt")

    def on_epoch_end(self, epoch, model):
        """Perform actions at the end of an epoch."""
        self.epoch_fp = self.filepath.format(epoch=epoch)
        os.makedirs(self.epoch_fp, exist_ok=True)
        if os.path.exists(self.epoch_fp):
            model_path = os.path.join(self.epoch_fp, "weights.pt")
            self.logger.debug(f"Saving model weights to {model_path}")
            torch.save(model.state_dict(), model_path)
            self.pack_meta(clean=True)

    def pack_meta(self, clean=True):
        """Write meta to pickled file in model directory and zip."""
        with open(os.path.join(self.epoch_fp, 'meta.pkl'), 'wb') as handle:
            self.logger.debug(f"Saving model meta to {handle.name}")
            pickle.dump(self.medaka_meta, handle)
        with tarfile.open(self.epoch_fp + '.tar.gz', "w:gz") as tar:
            self.logger.debug(f"Compressing model to {tar.fileobj.name}")
            tar.add(
                self.epoch_fp,
                arcname=medaka.datastore.ModelStoreTGZ.top_level_dir)
        if clean:
            shutil.rmtree(self.epoch_fp)


# Training datasets and other useful classes

class Sequence(torch.utils.data.Dataset):
    """Interface for a torch Dataset to a sample array in a `TrainBatcher`."""

    def __init__(
            self, data, sampler_func, dataset="train", mini_epochs=1,
            seed=None):
        """Initialize batching for training.

        :param data: array of samples from a `TrainBatcher`.
        :param sampler_func: function to get a sample given an index in data.
        :param dataset: one of 'train' or 'validation'.
        :param mini_epochs: factor by which to rescale the number of batches
            in an epoch (useful to output checkpoints more frequently).
        :param seed: random seed for shuffling data.

        """
        self.data = data
        self.sampler_func = sampler_func
        self.dataset = dataset
        self.mini_epochs = mini_epochs
        if seed is not None:
            np.random.seed(seed)

        self.logger = medaka.common.get_named_logger(
            '{}Batcher'.format(self.dataset.capitalize()))

        if dataset == 'train':
            if self.mini_epochs == 1:
                self.logger.info(
                    "Not using mini_epochs, an epoch is a full traversal "
                    "of the training data")
            else:
                self.logger.info(
                    "Using mini_epochs, an epoch is a traversal of 1/{} "
                    "of the training data".format(self.mini_epochs))
        elif dataset == 'validation':
            if mini_epochs != 1:
                raise ValueError(
                    "'mini_epochs' must be equal to 1 for validation data.")
        else:
            raise ValueError("'dataset' should be 'train' or 'validation'.")

        np.random.shuffle(self.data)
        self.samples_per_epoch = len(self.data) // self.mini_epochs

        self.mini_epoch_counter = 0
        self.epoch_start = 0

    def __len__(self):
        """Return the number of samples in one epoch."""
        return self.samples_per_epoch

    def __getitem__(self, index: int):
        """Get the sample of data at a given index position."""
        return self.sampler_func(self.data[index + self.epoch_start])

    def on_epoch_start(self):
        """Perform actions at the start of an epoch."""
        self.epoch_start = self.mini_epoch_counter * self.samples_per_epoch
        self.logger.debug(
            f"Mini-epoch {self.mini_epoch_counter + 1} of {self.mini_epochs}")
        self.logger.debug(
            f"Sample start index: {self.epoch_start}, "
            f"end index: {self.epoch_start + self.samples_per_epoch}")

    def on_epoch_end(self):
        """Perform actions at the end of an epoch."""
        if self.dataset == 'train':
            self.mini_epoch_counter += 1
            if self.mini_epoch_counter == self.mini_epochs:
                self.logger.debug("Shuffling data")
                np.random.shuffle(self.data)
                self.mini_epoch_counter = 0


class SequenceBatcher(torch.utils.data.DataLoader):
    """Dataloader operating on a Sequence dataset."""

    def __init__(self, sequence, batch_size=1, threads=1):
        """Initialize Dataloader.

        :param sequence: a `medaka.torch_ext.Sequence` instance.
        :param batch_size: number of samples per batch.
        :param threads: number of parallel loader threads.
        """
        super().__init__(
            sequence,
            batch_size=batch_size,
            shuffle=False,
            num_workers=threads,
            drop_last=True
        )
        self.dataset.logger.info(
            '{} batches of {} samples ({}) per epoch.'.format(
                len(self.dataset) // self.batch_size, self.batch_size,
                self.batch_size * (len(self.dataset) // self.batch_size),
            )
        )


class ClipGrad:
    """Gradient clipping by quantile."""

    def __init__(self, quantile=0.5, factor=2.0, buffer_size=100):
        """Initialise gradient clipper.

        :param quantile: quantile to use for setting gradient threshold.
        :param factor: factor multiplied by the quantile to calculate
            gradient threshold.
        :param buffer_size: size of gradient history to store.
        """
        self.buffer = np.full(buffer_size, fill_value=1e6)
        self.quantile = quantile
        self.factor = factor
        self.i = 0

    def append(self, grad_norm):
        """Add grad_norm to buffer."""
        self.buffer[self.i] = grad_norm
        self.i = (self.i + 1) % len(self.buffer)

    def __call__(self, parameters):
        """Clip gradient norms.

        :param parameters: model parameters.

        :return: clipped gradient norms.
        """
        max_norm = self.factor * np.quantile(self.buffer, self.quantile)
        grad_norm = torch.nn.utils.clip_grad_norm_(
            parameters, max_norm=max_norm).item()
        if not math.isnan(grad_norm):
            self.append(grad_norm)
        return grad_norm


# Training loop functions

class ProgressBar(tqdm):
    """Container for a progress bar to display epoch progress."""

    def __init__(self, total, ncols=100):
        """Initialise progress bar.

        :param total: expected number of samples per epoch.
        :param ncols: number of columns to display progress bar.
        """
        super().__init__(
            total=total, desc='[0/{}]'.format(total),
            ascii=True, leave=True, ncols=ncols,
            bar_format='{l_bar}{bar}| [{elapsed}{postfix}]'
        )

    def on_batch(self, samples, loss):
        """Update with batch results."""
        self.set_postfix(loss='%.4f' % loss)
        self.set_description(
            "[{}/{}]".format(samples, self.total))
        self.update()


def train_one_epoch(
        model, dataloader, loss_fn, optimizer, scaler=None, clip_grad=None,
        lr_scheduler=None, loss_log=None):
    """Run training for one epoch.

    :param model
    :param dataloader: SequenceBatcher, iterable of training batches.
    :param loss_fn: loss function.
    :param optimizer: torch.optim.Optimizer object.
    :param scaler: optional torch.cuda.amp.GradScaler, gradient scaler.
    :param clip_grad: optional, gradient clipping function.
    :param lr_scheduler: optional, learning rate scheduler.
    :param loss_log: optional CSVLogger, log loss per batch to file.

    :returns: mean training loss per sample.
    """
    model.train()
    device = model.device()
    dataloader.dataset.on_epoch_start()
    train_loss = 0
    samples = 0
    t0 = perf_counter()
    model.normalise = False

    progress_bar = ProgressBar(len(dataloader))
    for X, targets in dataloader:
        optimizer.zero_grad()

        samples += X.shape[0]
        X, targets = X.to(device), targets.to(device)
        try:
            pred = model(X)
            # Model output order is (N,T,C) but torch loss functions
            # require (N,C,...), so apply swapaxes here.
            loss = loss_fn(pred.swapaxes(-1, -2), targets.squeeze())

            if lr_scheduler is not None:
                lr_scheduler.step()
            if scaler is not None:
                scaler.scale(loss).backward()
                scaler.unscale_(optimizer)
            else:
                loss.backward()
            if clip_grad is not None:
                grad_norm = clip_grad(model.parameters())
            if scaler is not None:
                scaler.step(optimizer)
                scaler.update()
            else:
                optimizer.step()
        except RuntimeError as ex:
            dump_dir = os.path.dirname(loss_log.filename)
            dump_state(model, (X, targets), optimizer, dump_dir)
            raise ex

        train_loss += loss.item()

        if loss_log is not None:
            lr = (
                lr_scheduler.get_last_lr()
                if lr_scheduler is not None
                else [pg["lr"] for pg in optimizer.param_groups]
            )
            if len(lr) == 1:
                lr = lr[0]
            loss_log.append({
                'samples': samples,
                'time': perf_counter() - t0,
                'grad_norm': grad_norm,
                'lr': lr,
                'loss': loss.item()
            })

        progress_bar.on_batch(samples, loss)

    dataloader.dataset.on_epoch_end()

    return train_loss / (samples)


def validate_one_epoch(model, dataloader, loss_fn, metrics={}):
    """Run validation for one epoch.

    :param model
    :param dataloader: SequenceBatcher, iterable of training batches.
    :param loss_fn: loss function.
    :param metrics: dict of secondary evaluation metrics to compute.

    :returns: (mean validation loss per sample, mean of each metric per sample)
    """
    model.eval()
    device = model.device()
    dataloader.dataset.on_epoch_start()
    test_loss = 0.0
    samples = 0
    metric_values = {k: 0 for k in metrics.keys()}
    model.normalise = False
    progress_bar = ProgressBar(len(dataloader))
    with torch.no_grad():
        for X, targets in dataloader:
            samples += X.shape[0]
            X, targets = X.to(device), targets.to(device)
            pred = model(X)
            loss = loss_fn(pred.swapaxes(-1, -2), targets.squeeze())
            test_loss += loss.item()
            for k, f in metrics.items():
                v = f(targets, pred)
                if v.ndim > 1:
                    v = v.mean(dim=-1)
                metric_values[k] += v.sum().item()

            progress_bar.on_batch(samples, loss)
    dataloader.dataset.on_epoch_end()

    return (
        test_loss / samples,
        {k: m / samples for k, m in metric_values.items()}
    )


# Score metrics

class Metric(ABC):
    """Abstract class for accuracy metrics."""

    @abstractmethod
    def __init__(self, *args, **kwargs):
        """Initialise an accuracy metric."""
        raise NotImplementedError

    @abstractmethod
    def __call__(self, y_true: torch.Tensor, y_pred: torch.Tensor):
        """Calculate accuracy.

        :param y_true: tensor of true class labels.
        :param y_pred: class output scores from network.

        :returns: tensor, binary whether each entry is correct.
        """
        raise NotImplementedError


class SparseCategoricalAccuracy(Metric):
    """Sparse categorical accuracy."""

    def __init__(self):
        """Initialise."""
        pass

    def __call__(self, y_true, y_pred):
        """Loss function for sparse_categorical_accuracy.

        :param y_true: tensor of true integer class labels.
        :param y_pred: tensor of class output scores from network.

        :returns: tensor, binary whether each entry is correct.
        """
        # reshape in case it's in shape (num_samples, 1) instead of
        # (num_samples,)
        true_type = y_true.dtype
        def_float_type = torch.get_default_dtype()
        if len(y_true.shape) == len(y_pred.shape):
            y_true = torch.squeeze(y_true, -1)
        # convert dense predictions to labels
        y_pred_labels = torch.argmax(y_pred, dim=-1)
        y_pred_labels = y_pred_labels.to(true_type)
        return (y_true == y_pred_labels).to(def_float_type)


class CategoricalAccuracy(SparseCategoricalAccuracy):
    """Categorical accuracy."""

    def __init__(self):
        """Initialise categorical accuracy."""
        super().__init__()

    def __call__(self, y_true, y_pred):
        """Loss function for binary accuracy.

        :param y_true: tensor of one-hot true class labels.
        :param y_pred: tensor of class output scores from network.

        :returns: tensor, binary whether each entry is correct.
        """
        return super().__call__(torch.argmax(y_true, dim=-1), y_pred)


class QScore(Metric):
    """Accuracy converted to q-score."""

    def __init__(
        self,
        eps: float = 1e-10,
        accuracy: Metric = CategoricalAccuracy,
        accuracy_kwargs: dict = {}
    ):
        """Initialize q-score.

        accuracy specifies the Metric function that will be called to get
        the accuracy which will be converted to Qscore.
        """
        self.eps = eps
        self.accuracy = accuracy(**accuracy_kwargs)

    def __call__(self, y_true, y_pred):
        """Calculate accuracy and convert to q-score."""
        error = 1 - self.accuracy(y_true, y_pred)
        error = error.mean(dim=-1).clip(self.eps, 1.)
        return -10.0 * 0.434294481 * torch.log(error)


def dump_state(model, batch, optimizer, dirname):
    """Dump current batch and state tensors to files, for debugging."""
    torch.save(model, os.path.join(dirname, "nan_model.pt"))
    torch.save(optimizer, os.path.join(dirname, "nan_optim.pt"))
    features, targets = batch
    torch.save(features, os.path.join(dirname, "nan_batch_features.pt"))
    torch.save(targets, os.path.join(dirname, "nan_batch_targets.pt"))


# Learning rate scheduler (taken from bonito.schedule)

def linear_schedule(y0, y1):
    """Linear Scheduler."""
    return lambda t: y0 + (y1 - y0) * t


def cosine_decay_schedule(y0, y1):
    """Cosine Decay Scheduler."""
    return lambda t: y1 + 0.5 * (y0 - y1) * (np.cos(t * np.pi) + 1.0)


def piecewise_schedule(knots, funcs):
    """Piecewise Scheduler."""
    def f(t):
        i = np.searchsorted(knots, t)
        t0 = 0.0 if i == 0 else knots[i - 1]
        t1 = 1.0 if i == len(knots) else knots[i]
        return funcs[i]((t - t0) / (t1 - t0))
    return f


def func_scheduler(
        optimizer, func, total_steps, warmup_steps=None, warmup_ratio=0.1,
        start_step=0):
    """Learning Rate Scheduler."""
    if warmup_steps:
        y0 = func(0.0)
        func = piecewise_schedule(
            [warmup_steps / total_steps],
            [linear_schedule(warmup_ratio * y0, y0), func]
        )
    return torch.optim.lr_scheduler.LambdaLR(
        optimizer, (lambda step: func((step + start_step) / total_steps)))


def linear_warmup_cosine_decay(end_ratio=0.01, warmup_steps=500, **kwargs):
    """Linear warmup, cosine decay scheduler."""
    return lambda optimizer, train_loader, epochs, last_epoch: func_scheduler(
        optimizer=optimizer,
        func=cosine_decay_schedule(1.0, end_ratio),
        total_steps=epochs * len(train_loader),
        warmup_steps=warmup_steps,
        start_step=last_epoch * len(train_loader),
    )


def no_schedule(warmup_steps=None, **kwargs):
    """Linear warmup, cosine decay scheduler."""
    return lambda optimizer, train_loader, epochs, last_epoch: func_scheduler(
        optimizer=optimizer,
        func=lambda x: 1.0,
        total_steps=epochs * len(train_loader),
        warmup_steps=warmup_steps,
        start_step=last_epoch * len(train_loader),
    )
