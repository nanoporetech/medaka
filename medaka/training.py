"""Training program and ancillary functions."""
import csv
import functools
import os

import numpy as np
import toml
import torch

import medaka.common
import medaka.datastore
import medaka.labels
import medaka.models
import medaka.torch_ext


def train(args):
    """Training program."""
    train_name = args.train_name
    medaka.common.mkdir_p(train_name, info='Results will be overwritten.')

    logger = medaka.common.get_named_logger('Training')
    logger.debug("Loading datasets:\n{}".format('\n'.join(args.features)))
    logger.info(f"Starting training with seed: {args.seed}")

    if args.validation_features is None:
        args.validation = args.validation_split
    else:
        args.validation = args.validation_features

    batcher = TrainBatcher(
        args.features, args.validation,
        args.seed, args.batch_size,
        args.max_samples, args.max_valid_samples,
        threads=args.threads_io)

    with torch.cuda.device("cuda:{}".format(args.device)):
        if args.validate_only:
            run_validation(
                train_name,
                batcher,
                model_path=args.model,
                threads_io=args.threads_io,
                loss_args=args.loss_args,
                amp=args.amp,
            )
        else:
            run_training(
                train_name, batcher, model_fp=args.model, epochs=args.epochs,
                samples_per_training_epoch=args.samples_per_training_epoch,
                threads_io=args.threads_io, optimizer=args.optimizer,
                optim_args=args.optim_args,
                loss_args=args.loss_args,
                amp=args.amp, use_lr_schedule=args.use_lr_schedule,)

    # stop batching threads
    logger.info("Training finished.")


def run_training(
        train_name, batcher, model_fp=None,
        epochs=10, threads_io=1,
        samples_per_training_epoch=None, optimizer="adam", optim_args=None,
        loss_args=None, quantile_grad_clip=True,
        amp=False, use_lr_schedule=True, save_optim_every=5):
    """Run training."""
    from medaka.torch_ext import ModelMetaCheckpoint

    logger = medaka.common.get_named_logger('RunTraining')

    time_steps, feat_dim = batcher.feature_shape[0:2]
    if model_fp.endswith('.gz') or model_fp.endswith('.hdf'):
        msg = "Loading pretrained model from model store: {}"
        logger.info(msg.format(model_fp))
        model_store = medaka.models.open_model(model_fp)
        partial_model_function = model_store.get_meta('model_function')
    elif model_fp is None or model_fp.endswith('toml'):
        if model_fp is None:
            model_dict = medaka.models.DEFAULT_MODEL_DICT
            msg = "No model file given, using default."
        else:
            msg = f"Loading model config from {model_fp}"
            model_dict = toml.load(model_fp)
        partial_model_function = functools.partial(
            medaka.models.model_from_dict, model_dict
        )
    else:
        raise ValueError(f"Unknown model file type: {model_fp}")

    # now build model
    model = partial_model_function().to('cuda')
    logger.info("Model loaded: {}".format(model))

    model_metadata = {
        'model_function': partial_model_function,
        'label_scheme': batcher.label_scheme,
        'feature_encoder': batcher.feature_encoder}

    model.check_feature_encoder_compatibility(batcher.feature_encoder)
    metacheckpoint = ModelMetaCheckpoint(model_metadata, train_name)

    metrics = ['val_model_tot_acc']
    loss = torch.nn.CrossEntropyLoss(
        **(loss_args if loss_args is not None else {}))
    logger.info("Using {} loss function".format(loss))

    if optimizer == 'nadam':
        if optim_args is None:
            optim_args = {
                'lr': 0.002,
                'betas': (0.9, 0.99),
                'eps': 1e-07,
                # 'momentum_decay': 0.0
            }
        optimizer = functools.partial(torch.optim.NAdam, **optim_args)
    elif optimizer == 'adam':
        if optim_args is None:
            optim_args = {
                'lr': 0.0001,
                'betas': (0.9, 0.99),
                'eps': 1e-07,
                # 'momentum_decay': 0.0
            }
        optimizer = functools.partial(torch.optim.Adam, **optim_args)
    elif optimizer == 'rmsprop':
        if optim_args is None:
            optim_args = {
                'lr': 0.001,
                'alpha': 0.9,
                'eps': 1e-07,
                'momentum': 0.0
            }
        optimizer = functools.partial(torch.optim.RMSprop, **optim_args)
    elif optimizer == 'sgd':
        if optim_args is None:
            optim_args = {
                'lr': 0.001,
            }
        optimizer = functools.partial(torch.optim.SGD, **optim_args)
    else:
        raise ValueError('Unknown optimizer: {}'.format(optimizer))
    optimizer = optimizer(model.parameters())
    logger.info("Optimizer: {}".format(optimizer))

    model.half_precision = amp
    scaler = torch.amp.GradScaler("cuda")

    if quantile_grad_clip:
        clip_grad = medaka.torch_ext.ClipGrad()
    else:
        def clip_grad_fn(params, max_norm=2.0):
            return torch.nn.utils.clip_grad_norm_(
                params, max_norm=max_norm).item()
        clip_grad = clip_grad_fn

    train_loader = batcher.train_loader(threads=threads_io)
    valid_loader = batcher.valid_loader(threads=threads_io)

    if use_lr_schedule:
        lr_scheduler = medaka.torch_ext.linear_warmup_cosine_decay()
        logger.info("Using cosine learning rate decay.")
    else:
        lr_scheduler = medaka.torch_ext.no_schedule(warmup_steps=None)
        logger.info("Using constant learning rate.")
    lr_scheduler = lr_scheduler(optimizer, train_loader, epochs, 0)

    best_val_loss, best_val_loss_epoch = np.inf, -1
    best_metrics = {k: 0 for k in metrics}

    if os.environ.get('MEDAKA_DETECT_ANOMALY') is not None:
        logger.info("Activating anomaly detection")
        torch.autograd.set_detect_anomaly(True, check_nan=True)

    with CSVLogger(
            os.path.join(train_name, 'training.csv'), write_cache=0) \
            as training_log:
        for n in range(epochs):
            with CSVLogger(
                    os.path.join(train_name, 'losses_{}.csv'.format(n))) \
                    as loss_log:
                train_loss, train_metrics = medaka.torch_ext.run_epoch(
                    model, train_loader, loss, optimizer, scaler, clip_grad,
                    is_training_epoch=True,
                    total_num_samples=samples_per_training_epoch,
                    lr_scheduler=lr_scheduler, loss_log=loss_log,
                )

            # save optimizer state dict
            if n % save_optim_every == 0:
                with open("optim_{}.pt".format(n), "wb") as f:
                    torch.save(optimizer.state_dict(), f)

            val_loss, val_metrics = medaka.torch_ext.run_epoch(
                model, valid_loader, loss, optimizer, scaler, clip_grad,
                is_training_epoch=False, lr_scheduler=lr_scheduler,
                loss_log=None)

            logger.info(
                f"epoch {n + 1}/{epochs}: "
                f"val_loss={val_loss}, train_loss={train_loss}"
            )

            def qacc(x):
                return -10 * np.log10(np.max([1e-6, 1 - x]))

            train_metrics = {"train_" + k: v for k, v in train_metrics.items()}
            val_metrics = {"val_" + k: v for k, v in val_metrics.items()}
            all_metrics = {**train_metrics, **val_metrics}

            for k, v in all_metrics.items():
                logger.info(f"{k}={v:.5f} ({qacc(v):.4f} Q)")

            # log epoch results
            training_log.append({
                "epoch": n,
                "train_loss": train_loss,
                "val_loss": val_loss,
                **all_metrics
            })

            # save model checkpoints
            metacheckpoint.on_epoch_end(n, model)

            for k, v in val_metrics.items():
                if k in best_metrics and v > best_metrics[k]:
                    metacheckpoint.on_epoch_end("best_{}".format(k), model)
                    best_metrics[k] = v

            if val_loss < best_val_loss:
                metacheckpoint.on_epoch_end("best_val_loss", model)
                best_val_loss = val_loss
                best_val_loss_epoch = n

            # early stopping
            if n >= best_val_loss_epoch + 20:
                break


def run_validation(
        train_name,
        batcher,
        model_path,
        threads_io,
        loss_args,
        amp=False,
):
    """Run single validation epoch."""
    logger = medaka.common.get_named_logger("RunValidation")

    time_steps, feat_dim = batcher.feature_shape[0:2]
    if model_path.endswith(".hdf5") or model_path.endswith(".gz"):
        logger.info("Loading pretrained model from {}".format(model_path))
        model_store = medaka.models.open_model(model_path)
        model = model_store.load_model(time_steps=time_steps)

    else:
        raise ValueError("Not a valid trained model: {}".format(model_path))

    logger.info("Number of trainable parameters: {}".format(
        model.count_parameters()
    ))
    model = model.to("cuda")

    loss = torch.nn.CrossEntropyLoss(
        **(loss_args if loss_args is not None else {})
    )
    logger.info("Using {} loss function".format(loss))

    valid_loader = batcher.valid_loader(threads=threads_io)

    if os.environ.get("MEDAKA_DETECT_ANOMALY") is not None:
        logger.info("Activating anomaly detection")
        torch.autograd.set_detect_anomaly(True, check_nan=True)

    val_loss, val_metrics = medaka.torch_ext.run_epoch(
        model,
        valid_loader,
        loss_fn=loss,
        optimizer=None,
        is_training_epoch=False,
        loss_log=None,
    )

    logger.info(f"val_loss={val_loss}")

    def qacc(x):
        return -10 * np.log10(1 - x + 1e-6)

    val_metrics = {"val_" + k: v for k, v in val_metrics.items()}

    for k, v in val_metrics.items():
        logger.info(f"{k}={v:.5f} ({qacc(v):.4f} Q)")


class TrainBatcher:
    """Batching of training and validation samples."""

    def __init__(
            self, features, validation=0.2, seed=0, batch_size=500,
            max_samples=None, max_valid_samples=None, threads=1):
        """Serve up batches of training or validation data.

        :param features: iterable of str, training feature files.
        :param validation: float, fraction of batches to use for validation, or
            iterable of str, validation feature files.
        :param seed: int, random seed for separation of batches into
            training/validation.
        :param batch_size: int, number of samples per batch.
        :param max_samples: int, number of samples to use for training, `None`
            to use the full dataset (default=`None`).
        :param max_valid_samples: int, number of samples to use for validation,
            `None` to use the full dataset (default=`None`).
        :param threads: int, number of threads to use for indexing samples.
        """
        self.logger = medaka.common.get_named_logger('TrainBatcher')
        self.features = features
        self.validation = validation
        self.seed = seed
        self.batch_size = batch_size

        di = medaka.datastore.DataIndex(self.features, threads=threads)
        self.samples = di.samples.copy()

        self.label_scheme = di.metadata['label_scheme']
        self.feature_encoder = di.metadata['feature_encoder']

        # check sample size using first batch
        test_sample, test_fname = self.samples[0]
        with medaka.datastore.DataStore(test_fname) as ds:
            # TODO: this should come from feature_encoder
            self.feature_shape = ds.load_sample(test_sample).features.shape
            assert len(self.feature_shape) in (
                2,
                3,
            ), "Bad feature shape. Must be 2D (counts) or 3D (full alignment)"
            self.feature_type = (
                "pileup" if len(self.feature_shape) == 2 else "full_alignment"
            )
        self.logger.info(
            "Sample features have shape {}".format(self.feature_shape))

        generator = np.random.default_rng(self.seed)
        if isinstance(self.validation, float):
            generator.shuffle(self.samples)
            n_sample_train = int((1 - self.validation) * len(self.samples))
            self.train_samples = self.samples[:n_sample_train]
            self.valid_samples = self.samples[n_sample_train:]
            self.logger.info(
                'Randomly selected {} ({:3.2%}) of features for '
                'validation (seed {})'.format(
                    len(self.valid_samples), self.validation, self.seed))
        else:
            self.train_samples = self.samples
            self.valid_samples = medaka.datastore.DataIndex(
                self.validation).samples.copy()
            msg = 'Found {} validation samples, to {:3.2%} of all the data'
            fraction = len(self.valid_samples) / \
                (len(self.valid_samples) + len(self.train_samples))
            self.logger.info(msg.format(len(self.valid_samples), fraction))

        msg1 = 'Got {} samples ({} labels) for {}'
        msg2 = '{} set, using {} samples for {}'

        self.logger.info(msg1.format(
            len(self.train_samples),
            len(self.train_samples) * self.feature_shape[0],
            'training'))
        if max_samples is not None and max_samples < len(self.train_samples):
            generator = np.random.default_rng(self.seed)
            generator.shuffle(self.train_samples)
            self.train_samples = self.train_samples[:max_samples]
            self.logger.info(
                msg2.format(
                    "max_samples", len(self.train_samples), "training"))

        self.logger.info(msg1.format(
            len(self.valid_samples),
            len(self.valid_samples) * self.feature_shape[0],
            'validation'))
        if (
            max_valid_samples is not None and
            max_valid_samples < len(self.valid_samples)
        ):
            generator = np.random.default_rng(self.seed)
            generator.shuffle(self.valid_samples)
            self.valid_samples = self.valid_samples[:max_valid_samples]
            self.logger.info(
                msg2.format(
                    "max_valid_samples", len(self.valid_samples), "validation"
                ))

    def train_loader(self, threads=1, **kwargs):
        """Create a DataLoader for the training samples array.

        :param threads: int, number of parallel threads for loading samples.

        :returns: SequenceBatcher
        """
        return medaka.torch_ext.SequenceBatcher(
            medaka.torch_ext.Sequence(
                self.train_samples,
                sampler_func=functools.partial(
                    self.load_sample_worker,
                    label_scheme=self.label_scheme),
                seed=self.seed,
                **kwargs
            ),
            batch_size=self.batch_size,
            threads=threads,
            collate_fn=functools.partial(
                medaka.torch_ext.Batch.collate,
                counts_matrix=True),
            shuffle=True
        )

    def valid_loader(self, threads=1, **kwargs):
        """Create a DataLoader for the validation samples array.

        :param threads: int, number of parallel threads for loading samples.

        :returns: SequenceBatcher
        """
        return medaka.torch_ext.SequenceBatcher(
            medaka.torch_ext.Sequence(
                self.valid_samples,
                sampler_func=functools.partial(
                    self.load_sample_worker, label_scheme=self.label_scheme),
                seed=self.seed,
                **kwargs
            ),
            batch_size=self.batch_size,
            threads=threads,
            collate_fn=functools.partial(
                medaka.torch_ext.Batch.collate,
                counts_matrix=True
            ),
            shuffle=False
        )

    @staticmethod
    def load_sample_worker(sample, label_scheme):
        """Load a `common.Sample` object.

        :param sample: (sample key, filename).
        :param label_scheme: `LabelScheme` obj.

        :returns: (np.ndarray of inputs, np.ndarray of labels)

        """
        sample_key, sample_file = sample

        with medaka.datastore.DataStore(sample_file) as ds:
            s = ds.load_sample(sample_key)
            # precompute and cache the majority_vote_probs
            # while still in worker thread - minor speed up
            s.majority_vote_probs

        return s


def get_argmax_preds_from_pileup(pileup):
    """Calculate majority vote and support from pileup features.

    :param pileup: np.ndarray, shape (n, 10), pileup (trad medaka) matrix.

    :returns: (np.ndarray of shape n, np.ndarray of shape n).
    """
    b2i = medaka.common.base2index
    # slice pileup and sum for reverse and forward base counts
    b = pileup[:, b2i['a']:b2i['t']+1] + pileup[:, b2i['A']:b2i['G']+1]
    # sum deletion counts (indexing in this way retains correct shape)
    d = pileup[:, b2i['d']:b2i['d']+1] + pileup[:, b2i['D']:b2i['D']+1]
    out = np.concatenate(
        [d, b], axis=-1
    )  # d first - it is the deletion class and labelled 0 in training data
    out[:, 0] += 1 - out.sum(axis=-1)
    return np.argmax(out, axis=-1), np.max(out, axis=-1)


class CSVLogger:
    """Common interface for writing logs to file.

    Based on bonito.io.CSVLogger.
    """

    def __init__(self, filename, sep=',', write_cache=100):
        """Initialise file for output.

        :param filename: output file path.
        :param sep: delimiter (default=',')
        :param write_cache: number of lines to cache before writing.
        """
        self.filename = str(filename)
        if os.path.exists(self.filename):
            with open(self.filename) as f:
                self.columns = csv.DictReader(f).fieldnames
        else:
            self.columns = None
        self.fh = open(self.filename, 'a', newline='')
        self.csvwriter = csv.writer(self.fh, delimiter=sep)
        self.count = 0
        self.write_cache = write_cache

    def set_columns(self, columns):
        """Set list of column headers."""
        if self.columns:
            raise Exception('Columns already set')
        self.columns = list(columns)
        self.csvwriter.writerow(self.columns)

    def append(self, row):
        """Append a row entry to the writer queue."""
        if self.columns is None:
            self.set_columns(row.keys())
        self.csvwriter.writerow([row.get(k, '-') for k in self.columns])
        self.count += 1
        if self.count > self.write_cache:
            self.count = 0
            self.fh.flush()

    def close(self):
        """Close output file handle."""
        self.fh.close()

    def __enter__(self):
        """Context manager entry for opening output file."""
        return self

    def __exit__(self, *args):
        """Context manager exit for closing output file."""
        self.close()
