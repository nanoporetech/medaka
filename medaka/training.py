"""Training program and ancillary functions."""
import csv
import functools
import os

import numpy as np
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

    if args.validation_features is None:
        args.validation = args.validation_split
    else:
        args.validation = args.validation_features

    batcher = TrainBatcher(
        args.features, args.validation,
        args.seed, args.batch_size,
        args.max_samples, args.max_valid_samples,
        threads=args.threads_io)

    with torch.cuda.device('cuda:{}'.format(args.device)):
        run_training(
            train_name, batcher, model_fp=args.model, epochs=args.epochs,
            n_mini_epochs=args.mini_epochs,
            threads_io=args.threads_io, optimizer=args.optimizer,
            optim_args=args.optim_args,
            loss_args=args.loss_args,
            lr_schedule=args.lr_schedule)

    # stop batching threads
    logger.info("Training finished.")


def run_training(
        train_name, batcher, model_fp=None,
        epochs=10, n_mini_epochs=1, threads_io=1,
        optimizer='rmsprop', optim_args=None,
        loss_args=None, quantile_grad_clip=True,
        lr_schedule='none'):
    """Run training."""
    import torch.optim as optimizers
    from medaka.torch_ext import ModelMetaCheckpoint

    logger = medaka.common.get_named_logger('RunTraining')

    time_steps, feat_dim = batcher.feature_shape

    if model_fp is not None:
        model_store = medaka.models.open_model(model_fp)
        partial_model_function = model_store.get_meta('model_function')
        model = model_store.load_model(time_steps=time_steps)
    else:
        num_classes = batcher.label_scheme.num_classes
        model_name = "two_layer_bidirectional_TorchGRU"
        # model_name = medaka.models.default_model
        model_function = medaka.models.model_builders[model_name]
        partial_model_function = functools.partial(
            model_function, feat_dim, num_classes)
        model = partial_model_function(time_steps=time_steps)
    model = model.to('cuda')

    model_metadata = {
        'model_function': partial_model_function,
        'label_scheme': batcher.label_scheme,
        'feature_encoder': batcher.feature_encoder}

    if isinstance(
            batcher.label_scheme, medaka.labels.DiploidLabelScheme):
        logger.debug("Training with DiploidLabelScheme is untested")

    accuracy = medaka.torch_ext.SparseCategoricalAccuracy
    metrics = {
        "val_cat_acc": accuracy(),
        "val_qscore": medaka.torch_ext.QScore(eps=1e-6, accuracy=accuracy,)
    }
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
        optimizer = functools.partial(optimizers.NAdam, **optim_args)
    elif optimizer == 'rmsprop':
        if optim_args is None:
            optim_args = {
                'lr': 0.001,
                'alpha': 0.9,
                'eps': 1e-07,
                'momentum': 0.0
            }
        optimizer = functools.partial(optimizers.RMSprop, **optim_args)
    elif optimizer == 'sgd':
        if optim_args is None:
            optim_args = {
                'lr': 0.001,
            }
        optimizer = functools.partial(optimizers.SGD, **optim_args)
    else:
        raise ValueError('Unknown optimizer: {}'.format(optimizer))
    optimizer = optimizer(model.parameters())
    logger.info("Optimizer: {}".format(optimizer))

    scaler = torch.cuda.amp.GradScaler()

    if quantile_grad_clip:
        clip_grad = medaka.torch_ext.ClipGrad()
    else:
        def clip_grad_fn(params, max_norm=2.0):
            return torch.nn.utils.clip_grad_norm_(
                params, max_norm=max_norm).item()
        clip_grad = clip_grad_fn

    metacheckpoint = ModelMetaCheckpoint(model_metadata, train_name)

    train_loader = batcher.train_loader(n_mini_epochs, threads=threads_io)
    valid_loader = batcher.valid_loader(threads=threads_io)

    total_mini_epochs = epochs * n_mini_epochs

    if lr_schedule == 'none':
        lr_scheduler = medaka.torch_ext.no_schedule(warmup_steps=None)
        logger.info("Using constant learning rate.")
    elif lr_schedule == 'cosine':
        lr_scheduler = medaka.torch_ext.linear_warmup_cosine_decay()
        logger.info("Using cosine learning rate decay.")
    else:
        raise ValueError(f"Unknown learning rate schedule: {lr_schedule}.")
    lr_scheduler = lr_scheduler(optimizer, train_loader, total_mini_epochs, 0)

    best_val_loss, best_val_loss_epoch = np.inf, -1
    best_metrics = {k: 0 for k in metrics.keys()}

    if os.environ.get('MEDAKA_DETECT_ANOMALY') is not None:
        logger.info("Activating anomaly detection")
        torch.autograd.set_detect_anomaly(True, check_nan=True)

    with CSVLogger(
            os.path.join(train_name, 'training.csv'), write_cache=0) \
            as training_log:
        for n in range(total_mini_epochs):
            with CSVLogger(
                    os.path.join(train_name, 'losses_{}.csv'.format(n))) \
                    as loss_log:
                train_loss = medaka.torch_ext.train_one_epoch(
                    model, train_loader, loss, optimizer, scaler, clip_grad,
                    lr_scheduler=lr_scheduler, loss_log=loss_log)

            val_loss, val_metrics = \
                medaka.torch_ext.validate_one_epoch(
                    model, valid_loader, loss, metrics=metrics)
            logger.info(
                f"epoch {n + 1}/{total_mini_epochs}: val_loss={val_loss}, "
                ", ".join(f'{k}={v:.4f}' for k, v in val_metrics.items())
            )

            # log epoch results
            training_log.append({
                "epoch": n,
                "train_loss": train_loss,
                "val_loss": val_loss,
                **val_metrics
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
        self.logger.info(
            "Sample features have shape {}".format(self.feature_shape))

        if isinstance(self.validation, float):
            np.random.seed(self.seed)
            np.random.shuffle(self.samples)
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
            np.random.shuffle(self.train_samples)
            self.train_samples = self.train_samples[:max_samples]
            self.logger.info(
                msg2.format(
                    "max_samples", len(self.training_samples), "training"))

        self.logger.info(msg1.format(
            len(self.valid_samples),
            len(self.valid_samples) * self.feature_shape[0],
            'validation'))
        if (
            max_valid_samples is not None and
            max_valid_samples < len(self.valid_samples)
        ):
            np.random.shuffle(self.valid_samples)
            self.valid_samples = self.valid_samples[:max_valid_samples]
            self.logger.info(
                msg2.format(
                    "max_valid_samples", len(self.valid_samples), "validation"
                ))

    def train_loader(self, mini_epochs=1, threads=1, **kwargs):
        """Create a DataLoader for the training samples array.

        :param mini_epochs: int, number of mini_epochs per full epoch.
        :param threads: int, number of parallel threads for loading samples.

        :returns: SequenceBatcher
        """
        return medaka.torch_ext.SequenceBatcher(
            medaka.torch_ext.Sequence(
                self.train_samples,
                sampler_func=functools.partial(
                        self.sample_to_x_y_bq_worker,
                        label_scheme=self.label_scheme),
                dataset="train",
                mini_epochs=mini_epochs,
                seed=self.seed,
                **kwargs
            ),
            batch_size=self.batch_size,
            threads=threads,
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
                        self.sample_to_x_y_bq_worker,
                        label_scheme=self.label_scheme),
                dataset="validation",
                mini_epochs=1,
                seed=self.seed,
                **kwargs
            ),
            batch_size=self.batch_size,
            threads=threads,
        )

    def sample_to_x_y(self, sample):
        """Convert a `common.Sample` object into a training x, y tuple.

        This method is synonymous to the static method
        sample_to_x_y_bq_worker.

        :param sample: (filename, sample key)

        :returns: (np.ndarray of inputs, np.ndarray of labels)

        """
        return self.sample_to_x_y_bq_worker(sample, self.label_scheme)

    def samples_to_batch(self, samples):
        """Convert a set of `common.Sample` objects into a training X, Y tuple.

        The function wraps `.sample_to_x_y` and stacks the outputs.

        :param samples: (filename, sample key) tuples

        :returns: (np.ndarray of inputs, np.ndarray of labels)

        """
        items = [self.sample_to_x_y(s) for s in samples]
        xs, ys = zip(*items)
        x, y = np.stack(xs), np.stack(ys)
        return x, y

    @staticmethod
    def sample_to_x_y_bq_worker(sample, label_scheme):
        """Convert a `common.Sample` object into a training x, y tuple.

        :param sample: (sample key, filename).
        :param label_scheme: `LabelScheme` obj.

        :returns: (np.ndarray of inputs, np.ndarray of labels)

        """
        sample_key, sample_file = sample

        with medaka.datastore.DataStore(sample_file) as ds:
            s = ds.load_sample(sample_key)
        if s.labels is None:
            raise ValueError("Sample {} in {} has no labels.".format(
                sample_key, sample_file))
        x = s.features
        y = label_scheme.encoded_labels_to_training_vectors(s.labels)
        return x, y


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
