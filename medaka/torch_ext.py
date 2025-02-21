"""Extensions to torch API for medaka."""
from collections import defaultdict
from dataclasses import dataclass
import math
import os
import pickle
import shutil
import tarfile
from time import perf_counter

import numpy as np
import toml
import torch
from tqdm import tqdm

import medaka.common
import medaka.datastore

EXPORT_CONFIG_VERSION = 3


# Subclass to save model and metadata
class ModelMetaCheckpoint:
    """Custom ModelCheckpoint to add medaka-specific metadata to hdf5 files."""

    def __init__(self, medaka_meta, outdir):
        """Initialize checkpointing.

        :param medaka_meta: dict of meta data to store in checkpoint files.
        :param outdir: directory in which to output saved models
        """
        required_meta = {'model_function', 'label_scheme', 'feature_encoder'}
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
            self, data, sampler_func, dataset="train",
            seed=None):
        """Initialize batching for training.

        :param data: array of samples from a `TrainBatcher`.
        :param sampler_func: function to get a sample given an index in data.
        :param seed: random seed for shuffling data.

        """
        self.data = data
        self.sampler_func = sampler_func
        self.dataset = dataset
        if seed is not None:
            np.random.seed(seed)

        if dataset not in ['train', 'validation']:
            raise ValueError(
                "dataset must be one of 'train' or 'validation'")

        self.logger = medaka.common.get_named_logger(
            '{}Batcher'.format(self.dataset.capitalize()))

        np.random.shuffle(self.data)

    def __len__(self):
        """Return the number of samples in one epoch."""
        return len(self.data)

    def __getitem__(self, index: int):
        """Get the sample of data at a given index position."""
        return self.sampler_func(self.data[index])


@dataclass
class Batch:
    """Batch of samples, used for both training and inference."""

    read_level_features: torch.Tensor = None
    counts_matrix: torch.Tensor = None
    labels: torch.Tensor = None
    majority_vote_probs: torch.Tensor = None

    @classmethod
    def collate(cls, samples, counts_matrix=False):
        """Construct batch from a list of samples.

        :param samples: List of `medaka.common.Sample` objects.
        :param counts_matrix: bool: Whether to calculate the counts matrix
            if features are read-level. Disabled at inference time for speed,
            enabled at training so models can be compared with argmax. Note,
            this does not affect features already saved as counts matrices.
        :returns: torch_ext.Batch
        """
        batch_size = len(samples)
        feature_shape = samples[0].features.shape
        features = [torch.from_numpy(s.features) for s in samples]
        batch_dict = {}

        def pad_to_max_depth(read_level_features):
            depths = [rlf.shape[1] for rlf in read_level_features]
            npos, _, nfeats = feature_shape
            depths = [rlf.shape[1] for rlf in read_level_features]
            max_depth = max(depths)
            rlf_shape = (batch_size, npos, max_depth, nfeats)
            feature_matrix = np.zeros(rlf_shape, dtype=np.uint8)
            for i, rlf in enumerate(read_level_features):
                feature_matrix[i, :, : depths[i], :] = rlf
            return torch.from_numpy(feature_matrix).to(torch.uint8)

        if features[0].ndim == 3:
            # 3-dimensional features are read level
            batch_dict['read_level_features'] = pad_to_max_depth(features)
            if counts_matrix:
                batch_dict['counts_matrix'] = torch.stack(
                    [
                        torch.from_numpy(s.counts_matrix)
                        for s in samples
                    ]).float()
        elif features[0].ndim == 2:
            batch_dict['counts_matrix'] = torch.stack(features).float()
        else:
            raise ValueError(
                f"Unknown feature dimension {features[0].ndim}. Expect 3 for"
                "read level features or 2 for counts matrices.")

        if getattr(batch_dict, 'counts_matrix', None) is not None:
            # also get the majority vote probs for comparison of model
            # with argmax
            batch_dict['majority_vote_probs'] = torch.stack(
                [torch.from_numpy(s.majority_vote_probs) for s in samples]
            )

        if samples[0].labels is not None:
            batch_dict['labels'] = torch.stack([
                torch.from_numpy(s.labels) for s in samples])

        # now make batch from the tensors
        return cls(**batch_dict)

    @property
    def features(self):
        """Return the features tensor."""
        if self.read_level_features is None:
            return self.counts_matrix
        return self.read_level_features


class SequenceBatcher(torch.utils.data.DataLoader):
    """Dataloader operating on a Sequence dataset."""

    def __init__(
            self, sequence, batch_size=1, threads=1,
            shuffle=False, collate_fn=Batch.collate):
        """Initialize Dataloader.

        :param sequence: a `medaka.torch_ext.Sequence` instance.
        :param batch_size: number of samples per batch.
        :param threads: number of parallel loader threads.
        :param collate_fn: function to collate samples into a batch.
        :param shuffle: shuffle the data at the start of each epoch.
        """
        super().__init__(
            sequence,
            batch_size=batch_size,
            num_workers=threads,
            shuffle=shuffle,
            collate_fn=collate_fn,
            drop_last=True,
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

    def on_batch(self, n_samples, loss, model_acc, argmax_acc, batch_size):
        """Update with batch results."""

        def qacc(x):
            return -10 * np.log10(np.max([1e-6, 1 - x]))

        model_qacc = qacc(model_acc)
        argmax_qacc = qacc(argmax_acc)

        self.set_postfix(
            loss="%.4f" % loss,
            model_q="%.2f" % model_qacc,
            argmax_q="%.2f" % argmax_qacc,
            delta="%.2f" % (model_qacc - argmax_qacc),
        )
        self.set_description("[{}/{}]".format(n_samples, self.total))
        self.update(batch_size)


def run_epoch(
        model, dataloader, loss_fn, optimizer, scaler=None, clip_grad=None,
        lr_scheduler=None, loss_log=None,
        is_training_epoch=False, total_num_samples=None):
    """Run training for one epoch.

    :param model
    :param dataloader: SequenceBatcher, iterable of training batches.
    :param loss_fn: loss function.
    :param optimizer: torch.optim.Optimizer object.
    :param scaler: optional torch.cuda.amp.GradScaler, gradient scaler.
    :param clip_grad: optional, gradient clipping function.
    :param lr_scheduler: optional, learning rate scheduler.
    :param loss_log: optional CSVLogger, log loss per batch to file.
    :param is_training_epoch (bool), whether this is a training epoch.
    :param total_num_samples: optional int, number of samples in the epoch.

    :returns: mean training loss per sample.
    """
    model.train() if is_training_epoch else model.eval()
    sum_epoch_loss = 0
    n_samples = 0
    t0 = perf_counter()
    model.normalise = not isinstance(loss_fn, torch.nn.CrossEntropyLoss)
    total_epoch_metrics = defaultdict(float)

    # use the entire dataset if total_num_samples is not provided
    if total_num_samples is None:
        total_num_samples = len(dataloader.dataset)

    total_num_batches = total_num_samples // dataloader.batch_size

    progress_bar = ProgressBar(total_num_samples)
    for batch_idx, batch in enumerate(dataloader):
        n_samples += batch.labels.shape[0]
        if optimizer is not None:
            optimizer.zero_grad()

        with torch.set_grad_enabled(is_training_epoch):
            loss, batch_metrics = model.process_batch(batch, loss_fn)

        if is_training_epoch:
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            if clip_grad is not None:
                grad_norm = clip_grad(model.parameters())
            else:
                grad_norm = 0
            if scaler is not None:
                scaler.step(optimizer)
                scaler.update()
            else:
                optimizer.step()

        sum_epoch_loss += loss.item()

        batch_model_acc = (
            batch_metrics["n_model_correct"]
            / batch_metrics["n_positions"]
        )
        if "n_argmax_correct" in batch_metrics:
            batch_argmax_acc = (
                batch_metrics["n_argmax_correct"]
                / batch_metrics["n_positions"]
            )
        else:
            batch_argmax_acc = 0
        for key, count in batch_metrics.items():
            total_epoch_metrics[key] += count

        # record batch-level losses
        if loss_log is not None:
            lr = (
                lr_scheduler.get_last_lr()
                if lr_scheduler is not None
                else [pg["lr"] for pg in optimizer.param_groups]
            )
            if len(lr) == 1:
                lr = lr[0]
            loss_log.append(
                {
                    "samples": n_samples,
                    "time": perf_counter() - t0,
                    "grad_norm": grad_norm,
                    "lr": lr,
                    "loss": loss.item(),
                    "model_correct": batch_model_acc,
                    "argmax_correct": batch_argmax_acc,
                    "n_positions": batch_metrics["n_positions"],
                }
            )

        if lr_scheduler is not None:
            lr_scheduler.step()

        progress_bar.on_batch(
            n_samples,
            loss,
            batch_model_acc,
            batch_argmax_acc,
            dataloader.batch_size,
        )

        if batch_idx >= total_num_batches - 1:
            break

    epoch_model_acc = (
        total_epoch_metrics["n_model_correct"]
        / total_epoch_metrics["n_positions"]
    )
    epoch_argmax_acc = (
        total_epoch_metrics["n_argmax_correct"]
        / total_epoch_metrics["n_positions"]
    )

    # calculate metrics to return
    epoch_mean_batch_loss = sum_epoch_loss / (
        n_samples / dataloader.batch_size
    )  # divide sum loss by number of batches
    epoch_metrics = {
        "model_tot_acc": epoch_model_acc,
        "argmax_tot_acc": epoch_argmax_acc,
    }

    for key, count in total_epoch_metrics.items():
        epoch_metrics[key] = count

    return epoch_mean_batch_loss, epoch_metrics


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


def export_model(args):
    """Export a model to a torchscript model and config file."""
    logger = medaka.common.get_named_logger('ModelExport')
    model_fp = args.model
    output_fp = args.output
    force = args.force
    script = args.script
    supported_basecallers = args.supported_basecallers

    if not os.path.exists(model_fp):
        raise FileNotFoundError(f"Model file not found: {model_fp}")

    if output_fp is None:
        output_fp = os.path.basename(model_fp).replace(".tar.gz", "_export")

    logger.info(f"Exporting model from {model_fp} to {output_fp}.tar.gz")
    logger.info("Supported basecallers: {}".format(supported_basecallers))
    if force:
        os.makedirs(output_fp, exist_ok=True)
    else:
        os.mkdir(output_fp)

    msg = "Output file path cannot be the same as the model file path"
    assert model_fp != output_fp, msg

    model_store = medaka.models.open_model(model_fp)
    model = model_store.load_model()
    label_scheme = model_store.get_meta("label_scheme")
    feature_encoder = model_store.get_meta("feature_encoder")

    config = {
        'config_version': EXPORT_CONFIG_VERSION,
        'model': model.to_dict(),
        'feature_encoder': feature_encoder.to_dict(),
        'supported_basecallers': supported_basecallers,
        'label_scheme': label_scheme.to_dict()}

    with open(os.path.join(output_fp, 'config.toml'), 'w') as f:
        toml.dump(config, f)

    if script:
        for name, module in model.named_modules():
            if (
                isinstance(module, torch.nn.LSTM) or
                isinstance(module, torch.nn.GRU)
            ):
                module.flatten_parameters()
        scripted_model = torch.jit.script(model)
        scripted_model.save(os.path.join(output_fp, 'model.pt'))

    # save state dict
    w = {k: v.cpu() for k, v in model.state_dict().items()}
    torch.save(w, os.path.join(output_fp, 'weights.pt'))

    # now add to tarfile
    output_tarfile = output_fp + ".tar.gz"
    with tarfile.open(output_tarfile, "w:gz") as tar:
        tar.add(output_fp, arcname='model')

    shutil.rmtree(output_fp)
