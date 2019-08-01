from concurrent.futures import ProcessPoolExecutor
import functools
import inspect
import logging
import os
from timeit import default_timer as now

import numpy as np

import medaka.common
import medaka.datastore
import medaka.features
import medaka.labels
import medaka.models


def qscore(y_true, y_pred):
    from tensorflow.keras import backend as K
    error = K.cast(K.not_equal(
        K.max(y_true, axis=-1), K.cast(K.argmax(y_pred, axis=-1), K.floatx())),
        K.floatx()
    )
    error = K.sum(error) / K.sum(K.ones_like(error))
    return -10.0 * 0.434294481 * K.log(error)


def cat_acc(y_true, y_pred):
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


def run_training(train_name, batcher, model_fp=None,
                 epochs=5000, class_weight=None, n_mini_epochs=1, threads_io=1, multi_label=False):
    """Run training."""
    from tensorflow.keras.callbacks import CSVLogger, TensorBoard, EarlyStopping, ReduceLROnPlateau
    from medaka.keras_ext import ModelMetaCheckpoint, SequenceBatcher

    logger = medaka.common.get_named_logger('RunTraining')

    if model_fp is None:
        model_name = medaka.models.default_model
        model_kwargs = {
            k:v.default for (k,v) in
            inspect.signature(medaka.models.model_builders[model_name]).parameters.items()
            if v.default is not inspect.Parameter.empty
        }
    else:
        with medaka.datastore.DataStore(model_fp) as ds:
            model_name = ds.meta['medaka_model_name']
            model_kwargs = ds.meta['medaka_model_kwargs']

    opt_str = '\n'.join(['{}: {}'.format(k,v) for k, v in model_kwargs.items()])
    logger.info('Building {} model with: \n{}'.format(model_name, opt_str))
    num_classes = len(batcher.label_scheme.label_decoding)
    timesteps, feat_dim = batcher.feature_shape
    model = medaka.models.model_builders[model_name](timesteps, feat_dim, num_classes, **model_kwargs)

    if model_fp is not None:
        try:
            model.load_weights(model_fp)
            logger.info("Loading weights from {}".format(model_fp))
        except:
            logger.info("Could not load weights from {}".format(model_fp))

    msg = "feat_dim: {}, timesteps: {}, num_classes: {}"
    logger.info(msg.format(feat_dim, timesteps, num_classes))
    model.summary()

    model_details = batcher.meta.copy()

    model_details['medaka_model_name'] = model_name
    model_details['medaka_model_kwargs'] = model_kwargs
    model_details['medaka_multi_label'] = multi_label
    model_details['medaka_label_decoding'] = batcher.label_scheme.label_decoding
    model_details['medaka_label_description'] = batcher.label_scheme.label_description
    model_details['medaka_label_counts'] = batcher.label_scheme.label_counts
    model_details['medaka_label_scheme'] = batcher.label_scheme.__class__.__name__

    opts = dict(verbose=1, save_best_only=True, mode='max')

    if multi_label:
        metrics = ['binary_accuracy']
        call_back_metrics = metrics
        loss = 'binary_crossentropy'
        logger.info("Using {} loss function for multi-label training".format(loss))
    else:
        metrics=[cat_acc, qscore]
        call_back_metrics = {'cat_acc': cat_acc}
        if class_weight is not None:
            loss = weighted_categorical_crossentropy(class_weight)
            logger.info("Using weighted_categorical_crossentropy loss function")
        else:
            loss = 'sparse_categorical_crossentropy'
            logger.info("Using {} loss function".format(loss))

    model.compile(
       loss=loss,
       optimizer='nadam',
       metrics=metrics,
    )

    logger.info('Model metrics: {}'.format(model.metrics_names))

    callbacks = []
    for metric in call_back_metrics:
        for m in metric, 'val_{}'.format(metric):
            best_fn = 'model.best.{}.hdf5'.format(m)
            improv_fn = 'model-' + metric + '-improvement-{epoch:02d}-{' + metric + ':.2f}.hdf5'
            for fn in best_fn, improv_fn:
                callbacks.append(ModelMetaCheckpoint(model_details, os.path.join(train_name, fn), monitor=m, **opts))
    callbacks.extend([
        ## Reduce learning rate when no improvement
        #ReduceLROnPlateau(monitor='val_loss', factor=0.1, patience=5,
        #    verbose=1, min_delta=0.0001, cooldown=0, min_lr=0),
        # Stop when no improvement
        EarlyStopping(monitor='val_loss', patience=20),
        # Log of epoch stats
        CSVLogger(os.path.join(train_name, 'training.log')),
        # Allow us to run tensorboard to see how things are going. Some
        #   features require validation data, not clear why.
        #TensorBoard(log_dir=os.path.join(train_name, 'logs'),
        #            histogram_freq=5, batch_size=100, write_graph=True,
        #            write_grads=True, write_images=True)
    ])


    if n_mini_epochs == 1:
        logger.info("Not using mini_epochs, an epoch is a full traversal of the training data")
    else:
        logger.info("Using mini_epochs, an epoch is a traversal of 1/{} of the training data".format(n_mini_epochs))

    model.fit_generator(
        SequenceBatcher(batcher, mini_epochs=n_mini_epochs),
        validation_data=SequenceBatcher(batcher, 'validation'),
        max_queue_size=2*threads_io, workers=threads_io, use_multiprocessing=True,
        epochs=epochs,
        callbacks=callbacks,
        class_weight=class_weight,
    )


class TrainBatcher():
    def __init__(self, features, label_scheme_cls, max_label_len, validation=0.2, seed=0, batch_size=500, threads=1):
        """
        Class to server up batches of training / validation data.

        :param features: iterable of str, training feature files.
        :param label_scheme_cls, LabellingScheme class.
        :param max_label_len: int, maximum label length, longer labels will be truncated.
        :param validation: float, fraction of batches to use for validation, or
                iterable of str, validation feature files.
        :param seed: int, random seed for separation of batches into training/validation.
        :param batch_size: int, number of samples per batch.
        :param threads: int, number of threads to use for preparing batches.
        """
        self.logger = medaka.common.get_named_logger('TrainBatcher')

        self.features = features
        self.max_label_len = max_label_len
        self.validation = validation
        self.seed = seed
        self.sparse_labels = label_scheme_cls.sparse_labels
        self.batch_size = batch_size

        di = medaka.datastore.DataIndex(self.features, threads=threads)
        self.samples = di.samples.copy()
        self.meta = di.meta.copy()
        self.label_counts = self.meta['medaka_label_counts']
        self.label_scheme = label_scheme_cls(self.label_counts, max_label_len=self.max_label_len)

        # check sample size using first batch
        test_sample, test_fname = self.samples[0]
        with medaka.datastore.DataStore(test_fname) as ds:
            self.feature_shape = ds.load_sample(test_sample).features.shape
        self.logger.info("Sample features have shape {}".format(self.feature_shape))

        if isinstance(self.validation, float):
            np.random.seed(self.seed)
            np.random.shuffle(self.samples)
            n_sample_train = int((1 - self.validation) * len(self.samples))
            self.train_samples = self.samples[:n_sample_train]
            self.valid_samples = self.samples[n_sample_train:]
            msg = 'Randomly selected {} ({:3.2%}) of features for validation (seed {})'
            self.logger.info(msg.format(len(self.valid_samples), self.validation, self.seed))
        else:
            self.train_samples = self.samples
            self.valid_samples = medaka.datastore.DataIndex(self.validation).samples.copy()
            msg = 'Found {} validation samples equivalent to {:3.2%} of all the data'
            fraction = len(self.valid_samples) / len(self.valid_samples) + len(self.train_samples)
            self.logger.info(msg.format(len(self.valid_samples), fraction))

        msg = 'Got {} samples in {} batches ({} labels) for {}'
        self.logger.info(msg.format(len(self.train_samples),
                                    len(self.train_samples) // batch_size,
                                    len(self.train_samples) * self.feature_shape[0],
                                    'training'))
        self.logger.info(msg.format(len(self.valid_samples),
                                    len(self.valid_samples) // batch_size,
                                    len(self.valid_samples) * self.feature_shape[0],
                                    'validation'))


    def sample_to_x_y(self, sample):
        """Convert a `medaka.common.Sample` object into an x,y tuple for training.

        :param sample: (filename, sample key)

        :returns: (np.ndarray of inputs, np.ndarray of labels)

        """
        return self.sample_to_x_y_bq_worker(sample, self.label_scheme)


    def samples_to_batch(self, samples):
        """Convert a set of `medaka.common.Sample` objects into an X, Y tuple for training.

        :param samples: (filename, sample key) tuples

        :returns: (np.ndarray of inputs, np.ndarray of labels)

        """
        t0 = now()
        items = [self.sample_to_x_y(s) for s in samples]
        xs, ys = zip(*items)
        x, y = np.stack(xs), np.stack(ys)
        return x, y


    @staticmethod
    def sample_to_x_y_bq_worker(sample, label_scheme):
        """Convert a `medaka.common.Sample` object into an x,y tuple for training.

        :param sample: (filename, sample key)
        :param label_scheme: `LabellingScheme` obj
        :returns: (np.ndarray of inputs, np.ndarray of labels)

        """
        sample_key, sample_file = sample

        with medaka.datastore.DataStore(sample_file) as ds:
            s = ds.load_sample(sample_key)
        if s.labels is None:
            raise ValueError("Sample {} in {} has no labels.".format(sample_key, sample_file))
        x = s.features

        # trim label lengths to max_label_len
        s.labels['run_length'] = np.minimum(s.labels['run_length'], label_scheme.max_label_len, out=s.labels['run_length'])

        if label_scheme.sparse_labels:
            # label encoding values are tuples of int labels, but for sparse labels
            # there will only be one, so take first one.
            # need to use tolist() to convert numpy dtypes to python dtypes
            int_tups = (label_scheme.label_encoding[tuple(l)][0] for l in s.labels.tolist())
            y = np.fromiter(int_tups, dtype=int, count=len(s.labels)).reshape(-1, 1)

        else:
            y = np.zeros(shape=(len(s.labels), len(label_scheme.label_decoding)), dtype=int)
            encoded = (label_scheme.label_encoding[tuple(l)] for l in s.labels.tolist())
            # TODO can this be implemented more efficiently up, e.g. converting
            # this to sparse array in scipy / numpy then to non-sparse array?
            for row_ind, col_inds in enumerate(encoded):
                for col_ind in col_inds:
                    y[row_ind, col_ind] = 1

        return x, y


def run_prediction(output, bam, regions, model, model_file,
        chunk_len, chunk_ovlp, batch_size=200,
        save_features=False, tag_name=None, tag_value=None,
        tag_keep_missing=False, enable_chunking=True):
    """Inference worker."""

    logger = medaka.common.get_named_logger('PWorker')

    remainder_regions = list()
    def sample_gen():
        # chain all samples whilst dispensing with generators when done
        #   (they hold the feature vector in memory until they die)
        for region in regions:
            data_gen = medaka.features.SampleGenerator(
                bam, region, model_file,
                chunk_len=chunk_len, chunk_overlap=chunk_ovlp,
                tag_name=tag_name, tag_value=tag_value,
                tag_keep_missing=tag_keep_missing,
                enable_chunking=enable_chunking)
            yield from data_gen.samples
            remainder_regions.extend(data_gen._quarantined)
    batches = medaka.common.background_generator(
        medaka.common.grouper(sample_gen(), batch_size), 10
    )

    total_region_mbases = sum(r.size for r in regions) / 1e6
    logger.info("Running inference for {:.1f}M draft bases.".format(total_region_mbases))

    with medaka.datastore.DataStore(output, 'a', verify_on_close=False) as ds:
        mbases_done = 0

        t0 = now()
        tlast = t0
        for data in batches:
            x_data = np.stack([x.features for x in data])
            class_probs = model.predict_on_batch(x_data)
            # calculate bases done taking into account overlap
            new_bases = 0
            for x in data:
                if chunk_ovlp < x.size:
                    new_bases += x.last_pos[0] - x._get_pos(chunk_ovlp)[0]
                else:
                    new_bases += x.span
            mbases_done += new_bases / 1e6
            mbases_done = min(mbases_done, total_region_mbases)  # just to avoid funny log msg
            t1 = now()
            if t1 - tlast > 10:
                tlast = t1
                msg = '{:.1%} Done ({:.1f}/{:.1f} Mbases) in {:.1f}s'
                logger.info(msg.format(mbases_done / total_region_mbases, mbases_done, total_region_mbases, t1 - t0))

            best = np.argmax(class_probs, -1)
            for sample, prob, pred, feat in zip(data, class_probs, best, x_data):
                # write out positions and predictions for later analysis
                sample_d = sample._asdict()
                sample_d['label_probs'] = prob
                sample_d['features'] = feat if save_features else None
                ds.write_sample(medaka.common.Sample(**sample_d))

    logger.info("All done, {} remainder regions.".format(len(remainder_regions)))
    return remainder_regions


def predict(args):
    """Inference program."""
    logger_level = logging.getLogger(__package__).level
    if logger_level > logging.DEBUG:
        os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"

    import tensorflow as tf
    from tensorflow.keras.models import load_model
    from tensorflow.keras import backend as K

    args.regions = medaka.common.get_regions(args.bam, region_strs=args.regions)
    logger = medaka.common.get_named_logger('Predict')
    logger.info('Processing region(s): {}'.format(' '.join(str(r) for r in args.regions)))

    # write class names to output
    with medaka.datastore.DataStore(args.model) as ds:
        meta = ds.meta
    with medaka.datastore.DataStore(args.output, 'w', verify_on_close=False) as ds:
        ds.update_meta(meta)

    logger.info("Setting tensorflow threads to {}.".format(args.threads))
    tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
    K.set_session(tf.Session(
        config=tf.ConfigProto(
            intra_op_parallelism_threads=args.threads,
            inter_op_parallelism_threads=args.threads)
    ))

    # Split overly long regions to maximum size so as to not create
    #   massive feature matrices
    MAX_REGION_SIZE = int(1e6)  # 1Mb
    regions = []
    for region in args.regions:
        if region.size > MAX_REGION_SIZE:
            regs = region.split(MAX_REGION_SIZE, args.chunk_ovlp)
        else:
            regs = [region]
        regions.extend(regs)

    logger.info("Processing {} long region(s) with batching.".format(len(regions)))
    model = medaka.models.load_model(args.model, time_steps=args.chunk_len)
    # the returned regions are those where the pileup width is smaller than chunk_len
    remainder_regions = run_prediction(
        args.output, args.bam, regions, model, args.model,
        args.chunk_len, args.chunk_ovlp,
        batch_size=args.batch_size, save_features=args.save_features,
        tag_name=args.tag_name, tag_value=args.tag_value, tag_keep_missing=args.tag_keep_missing
    )

    # short/remainder regions: just do things without chunking. We can do this
    # here because we now have the size of all pileups (and know they are small).
    # TODO: can we avoid calculating pileups twice whilst controlling memory?
    if len(remainder_regions) > 0:
        logger.info("Processing {} short region(s).".format(len(remainder_regions)))
        model = medaka.models.load_model(args.model, time_steps=None)
        for region in remainder_regions:
            new_remainders = run_prediction(
                args.output, args.bam, [region[0]], model, args.model,
                args.chunk_len, args.chunk_ovlp, # these won't be used
                batch_size=args.batch_size, save_features=args.save_features,
                tag_name=args.tag_name, tag_value=args.tag_value, tag_keep_missing=args.tag_keep_missing,
                enable_chunking=False
            )
            if len(new_remainders) > 0:
                # shouldn't get here
                ignored = [x[0] for x in new_remainders]
                n_ignored = len(ignored)
                logger.warning("{} regions were not processed: {}.".format(n_ignored, ignored))

    logger.info("Finished processing all regions.")

    if args.check_output:
        logger.info("Validating and finalising output data.")
        with medaka.datastore.DataStore(args.output, 'a') as ds:
            pass


def train(args):
    """Training program."""
    train_name = args.train_name
    medaka.common.mkdir_p(train_name, info='Results will be overwritten.')

    logger = medaka.common.get_named_logger('Training')
    logger.debug("Loading datasets:\n{}".format('\n'.join(args.features)))

    args.validation = args.validation_features if args.validation_features is not None else args.validation_split

    label_scheme_cls = medaka.labels.label_schemes[args.label_scheme]
    batcher = TrainBatcher(args.features, label_scheme_cls, args.max_label_len, args.validation,
                           args.seed, args.batch_size, threads=args.threads_io)


    import tensorflow as tf
    with tf.device('/gpu:{}'.format(args.device)):
        run_training(
            train_name, batcher, model_fp=args.model, epochs=args.epochs,
            n_mini_epochs=args.mini_epochs,
            threads_io=args.threads_io, multi_label=args.multi_label)

    # stop batching threads
    logger.info("Training finished.")
