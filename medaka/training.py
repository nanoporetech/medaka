"""Training program and ancillary functions."""
import functools
import os

import numpy as np

import medaka.common
import medaka.datastore
import medaka.labels
import medaka.models


def qscore(y_true, y_pred):
    """Keras metric function for calculating scaled error.

    :param y_true: tensor of true class labels.
    :param y_pred: class output scores from network.

    :returns: class error expressed as a phred score.
    """
    from tensorflow.keras import backend as K
    error = K.cast(K.not_equal(
        K.max(y_true, axis=-1), K.cast(K.argmax(y_pred, axis=-1), K.floatx())),
        K.floatx()
    )
    error = K.sum(error) / K.sum(K.ones_like(error))
    return -10.0 * 0.434294481 * K.log(error)


def cat_acc(y_true, y_pred):
    """Keras loss function for sparse_categorical_accuracy.

    :param y_true: tensor of true class labels.
    :param y_pred: class output scores from network.

    :returns: categorical accuracy.
    """
    # sparse_categorical_accuracy is broken in keras 2.2.4
    #   https://github.com/keras-team/keras/issues/11348#issuecomment-439969957
    # this is taken from e59570ae
    from tensorflow.keras import backend as K
    # reshape in case it's in shape (num_samples, 1) instead of (num_samples,)
    if K.ndim(y_true) == K.ndim(y_pred):
        y_true = K.squeeze(y_true, -1)
    # convert dense predictions to labels
    y_pred_labels = K.argmax(y_pred, axis=-1)
    y_pred_labels = K.cast(y_pred_labels, K.floatx())
    return K.cast(K.equal(y_true, y_pred_labels), K.floatx())


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
        args.seed, args.batch_size, threads=args.threads_io)

    import tensorflow as tf
    with tf.device('/gpu:{}'.format(args.device)):
        run_training(
            train_name, batcher, model_fp=args.model, epochs=args.epochs,
            n_mini_epochs=args.mini_epochs,
            threads_io=args.threads_io, optimizer=args.optimizer,
            optim_args=args.optim_args)

    # stop batching threads
    logger.info("Training finished.")


def run_training(
        train_name, batcher, model_fp=None,
        epochs=5000, class_weight=None, n_mini_epochs=1, threads_io=1,
        optimizer='rmsprop', optim_args=None, allow_cudnn=True):
    """Run training."""
    from tensorflow.keras.callbacks import \
        CSVLogger, EarlyStopping, TerminateOnNaN
    from tensorflow.keras import optimizers
    from medaka.keras_ext import \
        ModelMetaCheckpoint, SequenceBatcher

    logger = medaka.common.get_named_logger('RunTraining')

    time_steps, feat_dim = batcher.feature_shape

    if model_fp is not None:
        with medaka.datastore.DataStore(model_fp) as ds:
            partial_model_function = ds.get_meta('model_function')
            model = partial_model_function(time_steps=time_steps,
                                           allow_cudnn=allow_cudnn)
        try:
            model.load_weights(model_fp)
            logger.info("Loading weights from {}".format(model_fp))
        except Exception:
            logger.info("Could not load weights from {}".format(model_fp))

    else:
        num_classes = batcher.label_scheme.num_classes
        model_name = medaka.models.default_model
        model_function = medaka.models.model_builders[model_name]
        partial_model_function = functools.partial(
            model_function, feat_dim, num_classes)
        model = partial_model_function(time_steps=time_steps,
                                       allow_cudnn=allow_cudnn)

    model_metadata = {'model_function': partial_model_function,
                      'label_scheme': batcher.label_scheme,
                      'feature_encoder': batcher.feature_encoder}

    opts = dict(verbose=1, save_best_only=True, mode='max')

    if isinstance(batcher.label_scheme,
                  medaka.labels.DiploidZygosityLabelScheme):

        metrics = ['binary_accuracy']
        call_back_metrics = metrics
        loss = 'binary_crossentropy'
        logger.info(
            "Using {} loss function for multi-label training".format(loss))
    else:
        metrics = [cat_acc, qscore]
        call_back_metrics = {'cat_acc': cat_acc}
        loss = 'sparse_categorical_crossentropy'
        logger.info("Using {} loss function".format(loss))

    if optimizer == 'nadam':
        if optim_args is None:
            # defaults from docs as of 01/09/2019
            optim_args = {
                'lr': 0.002, 'beta_1': 0.9, 'beta_2': 0.999,
                'epsilon': 1e-07, 'schedule_decay': 0.004}
        optimizer = optimizers.Nadam(**optim_args)
    elif optimizer == 'rmsprop':
        if optim_args is None:
            optim_args = {
                'lr': 0.001, 'rho': 0.9, 'epsilon': 1e-07, 'decay': 0.0}
        optimizer = optimizers.RMSprop(**optim_args)
    else:
        raise ValueError('Unknown optimizer: {}'.format(optimizer))
    model.compile(
       loss=loss,
       optimizer=optimizer,
       metrics=metrics,
    )

    logger.info('Model metrics: {}'.format(model.metrics_names))

    callbacks = []
    for metric in call_back_metrics:
        for m in metric, 'val_{}'.format(metric):
            best_fn = 'model.best.{}.hdf5'.format(m)
            improv_fn = 'model-' + metric + '-improvement-{epoch:02d}-{' \
                + metric + ':.2f}.hdf5'
            for fn in best_fn, improv_fn:
                callbacks.append(ModelMetaCheckpoint(
                    model_metadata, os.path.join(train_name, fn),
                    monitor=m, **opts))
    callbacks.extend([
        # Stop when no improvement
        EarlyStopping(monitor='val_loss', patience=20),
        # Log of epoch stats
        CSVLogger(os.path.join(train_name, 'training.log'), separator='\t'),
        # Allow us to run tensorboard to see how things are going
        # TrainValTensorBoard(
        #    log_dir=os.path.join(train_name, 'logs'),
        #    histogram_freq=5, batch_size=100, write_graph=True,
        #    write_grads=True, write_images=True)
        # terminate training when a Nan loss is encountered
        TerminateOnNaN()
    ])

    if n_mini_epochs == 1:
        logger.info(
            "Not using mini_epochs, an epoch is a full traversal "
            "of the training data")
    else:
        logger.info(
            "Using mini_epochs, an epoch is a traversal of 1/{} "
            "of the training data".format(n_mini_epochs))

    model.fit_generator(
        SequenceBatcher(batcher, mini_epochs=n_mini_epochs),
        validation_data=SequenceBatcher(batcher, 'validation'),
        max_queue_size=2*threads_io, workers=threads_io,
        use_multiprocessing=True, epochs=epochs, callbacks=callbacks,
        class_weight=class_weight)


class TrainBatcher():
    """Batching of training and validation samples."""

    def __init__(
            self, features, validation=0.2, seed=0,
            batch_size=500, threads=1):
        """Serve up batches of training or validation data.

        :param features: iterable of str, training feature files.
        :param validation: float, fraction of batches to use for validation, or
            iterable of str, validation feature files.
        :param seed: int, random seed for separation of batches into
            training/validation.
        :param batch_size: int, number of samples per batch.
        :param threads: int, number of threads to use for preparing batches.

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

        msg = 'Got {} samples in {} batches ({} labels) for {}'
        self.logger.info(msg.format(
            len(self.train_samples),
            len(self.train_samples) // batch_size,
            len(self.train_samples) * self.feature_shape[0],
            'training'))
        self.logger.info(msg.format(
            len(self.valid_samples),
            len(self.valid_samples) // batch_size,
            len(self.valid_samples) * self.feature_shape[0],
            'validation'))

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
