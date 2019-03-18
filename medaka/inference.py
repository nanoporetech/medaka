from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import functools
import inspect
import itertools
import os
from timeit import default_timer as now

import numpy as np
import pysam

from medaka import vcf
from medaka.datastore import DataStore, DataIndex
from medaka.common import (get_regions, decoding, grouper, mkdir_p, Sample,
                           _gap_, background_generator,
                           get_named_logger)

from medaka.features import SampleGenerator
import medaka.models


def weighted_categorical_crossentropy(weights):
    """
    A weighted version of keras.objectives.categorical_crossentropy
    @url: https://gist.github.com/wassname/ce364fddfc8a025bfab4348cf5de852d
    @author: wassname

    Variables:
        weights: numpy array of shape (C,) where C is the number of classes

    Usage:
        weights = np.array([0.5,2,10]) # Class one at 0.5, class 2 twice the normal weights, class 3 10x.
        loss = weighted_categorical_crossentropy(weights)
        model.compile(loss=loss,optimizer='adam')
    """
    from keras import backend as K
    weights = K.variable(weights)

    def loss(y_true, y_pred):
        # scale predictions so that the class probas of each sample sum to 1
        y_pred /= K.sum(y_pred, axis=-1, keepdims=True)
        # clip to prevent NaN's and Inf's
        y_pred = K.clip(y_pred, K.epsilon(), 1 - K.epsilon())
        # calc
        loss = y_true * K.log(y_pred) * weights
        loss = -K.sum(loss, -1)
        return loss

    return loss


def qscore(y_true, y_pred):
    from keras import backend as K
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
    from keras import backend as K
    # reshape in case it's in shape (num_samples, 1) instead of (num_samples,)
    if K.ndim(y_true) == K.ndim(y_pred):
        y_true = K.squeeze(y_true, -1)
    # convert dense predictions to labels
    y_pred_labels = K.argmax(y_pred, axis=-1)
    y_pred_labels = K.cast(y_pred_labels, K.floatx())
    return K.cast(K.equal(y_true, y_pred_labels), K.floatx())


def run_training(train_name, batcher, model_fp=None,
                 epochs=5000, class_weight=None, n_mini_epochs=1, threads_io=1):
    """Run training."""
    from keras.callbacks import CSVLogger, TensorBoard, EarlyStopping, ReduceLROnPlateau
    from medaka.keras_ext import ModelMetaCheckpoint, SequenceBatcher, BatchQueue

    logger = get_named_logger('RunTraining')

    if model_fp is None:
        model_name = medaka.models.default_model
        model_kwargs = {
            k:v.default for (k,v) in
            inspect.signature(medaka.models.model_builders[model_name]).parameters.items()
            if v.default is not inspect.Parameter.empty
        }
    else:
        with DataStore(model_fp) as ds:
            model_name = ds.meta['medaka_model_name']
            model_kwargs = ds.meta['medaka_model_kwargs']

    opt_str = '\n'.join(['{}: {}'.format(k,v) for k, v in model_kwargs.items()])
    logger.info('Building {} model with: \n{}'.format(model_name, opt_str))
    num_classes = len(batcher.label_counts)
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
    model_details['medaka_label_decoding'] = batcher.label_decoding

    opts = dict(verbose=1, save_best_only=True, mode='max')

    callbacks = [
        # Best model according to training set accuracy
        ModelMetaCheckpoint(model_details, os.path.join(train_name, 'model.best.hdf5'),
                            monitor='cat_acc', **opts),
        # Best model according to validation set accuracy
        ModelMetaCheckpoint(model_details, os.path.join(train_name, 'model.best.val.hdf5'),
                        monitor='val_cat_acc', **opts),
        # Best model according to validation set qscore
        ModelMetaCheckpoint(model_details, os.path.join(train_name, 'model.best.val.qscore.hdf5'),
                        monitor='val_qscore', **opts),
        # Checkpoints when training set accuracy improves
        ModelMetaCheckpoint(model_details, os.path.join(train_name, 'model-acc-improvement-{epoch:02d}-{cat_acc:.2f}.hdf5'),
                        monitor='cat_acc', **opts),
        ModelMetaCheckpoint(model_details, os.path.join(train_name, 'model-val_acc-improvement-{epoch:02d}-{val_cat_acc:.2f}.hdf5'),
                        monitor='val_cat_acc', **opts),
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
    ]

    if class_weight is not None:
        loss = weighted_categorical_crossentropy(class_weight)
        logger.info("Using weighted_categorical_crossentropy loss function")
    else:
        loss = 'sparse_categorical_crossentropy'
        logger.info("Using {} loss function".format(loss))

    model.compile(
       loss=loss,
       optimizer='nadam',
       metrics=[cat_acc, qscore],
    )

    if n_mini_epochs == 1:
        logger.info("Not using mini_epochs, an epoch is a full traversal of the training data")
    else:
        logger.info("Using mini_epochs, an epoch is a traversal of 1/{} of the training data".format(n_mini_epochs))


    with ProcessPoolExecutor(threads_io) as executor:
        logger.info("Starting data queues.")
        prep_function = functools.partial(
            batcher.sample_to_x_y_bq_worker,
            max_label_len=batcher.max_label_len,
            label_encoding=batcher.label_encoding,
            sparse_labels=batcher.sparse_labels,
            n_classes=batcher.n_classes
        )
        # TODO: should take mini_epochs into account here
        train_queue = BatchQueue(
            batcher.train_samples, prep_function, batcher.batch_size, executor,
            seed=batcher.seed, name='Train', maxsize=100
        )
        valid_queue = BatchQueue(
            batcher.valid_samples, prep_function, batcher.batch_size, executor,
            seed=batcher.seed, name='Valid', maxsize=100
        )

        # run training
        logger.info("Starting training.")
        model.fit_generator(
            generator=train_queue.yield_batches(), steps_per_epoch=train_queue.n_batches // n_mini_epochs,
            validation_data=valid_queue.yield_batches(), validation_steps=valid_queue.n_batches,
            max_queue_size=2*threads_io, workers=1, use_multiprocessing=False,
            epochs=epochs,
            callbacks=callbacks,
            class_weight=class_weight,
        )
        logger.info("Training finished.")
        train_queue.stop()
        valid_queue.stop()

    #TODO: understand why this is buggy (occasionally hangs during validation)
    #model.fit_generator(
    #    SequenceBatcher(batcher, mini_epochs=n_mini_epochs),
    #    validation_data=SequenceBatcher(batcher, 'validation'),
    #    max_queue_size=2*threads_io, workers=threads_io, use_multiprocessing=True,
    #    epochs=epochs,
    #    callbacks=callbacks,
    #    class_weight=class_weight,
    #)


class TrainBatcher():
    def __init__(self, features, max_label_len, validation=0.2, seed=0, sparse_labels=True, batch_size=500, threads=1):
        """
        Class to server up batches of training / validation data.

        :param features: iterable of str, training feature files.
        :param max_label_len: int, maximum label length, longer labels will be truncated.
        :param validation: float, fraction of batches to use for validation, or
                iterable of str, validation feature files.
        :param seed: int, random seed for separation of batches into training/validation.
        :param sparse_labels: bool, create sparse labels.

        """
        self.logger = get_named_logger('TrainBatcher')

        self.features = features
        self.max_label_len = max_label_len
        self.validation = validation
        self.seed = seed
        self.sparse_labels = sparse_labels
        self.batch_size = batch_size

        di = DataIndex(self.features, threads=threads)
        self.samples = di.samples.copy()
        self.meta = di.meta.copy()
        self.label_counts = self.meta['medaka_label_counts']

        # check sample size using first batch
        test_sample, test_fname = self.samples[0]
        with DataStore(test_fname) as ds:
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
            self.valid_samples = DataIndex(self.validation).samples.copy()
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

        self.n_classes = len(self.label_counts)

        # get label encoding, given max_label_len
        self.logger.info("Max label length: {}".format(self.max_label_len if self.max_label_len is not None else 'inf'))
        self.label_encoding, self.label_decoding, self.label_counts = process_labels(self.label_counts, max_label_len=self.max_label_len)


    def sample_to_x_y(self, sample):
        """Convert a `Sample` object into an x,y tuple for training.

        :param sample: (filename, sample key)

        :returns: (np.ndarray of inputs, np.ndarray of labels)

        """
        return self.sample_to_x_y_bq_worker(
            sample, self.max_label_len, self.label_encoding,
            self.sparse_labels, self.n_classes)


    def samples_to_batch(self, samples):
        """Convert a set of `Sample` objects into an X, Y tuple for training.

        :param samples: (filename, sample key) tuples

        :returns: (np.ndarray of inputs, np.ndarray of labels)

        """
        t0 = now()
        items = [self.sample_to_x_y(s) for s in samples]
        xs, ys = zip(*items)
        x, y = np.stack(xs), np.stack(ys)
        return x, y


    @staticmethod
    def sample_to_x_y_bq_worker(sample, max_label_len, label_encoding, sparse_labels, n_classes):
        """Convert a `Sample` object into an x,y tuple for training.

        :param sample: (filename, sample key)
        :param max_label_len: int, maximum label length, longer labels will be truncated.
        :param label_encoding: {label: int encoded label}.
        :param sparse_labels: bool, create sparse labels.
        :param n_classes: int, number of label classes.

        :returns: (np.ndarray of inputs, np.ndarray of labels)

        """
        sample_key, sample_file = sample

        with DataStore(sample_file) as ds:
            s = ds.load_sample(sample_key)
        if s.labels is None:
            raise ValueError("Sample {} in {} has no labels.".format(sample_key, sample_file))
        x = s.features
        # labels can either be unicode strings or (base, length) integer tuples
        if isinstance(s.labels[0], np.unicode):
            # TODO: is this ever used now we have dispensed with tview code?
            y = np.fromiter((label_encoding[l[:min(max_label_len, len(l))]]
                               for l in s.labels), dtype=int, count=len(s.labels))
        else:
            y = np.fromiter((label_encoding[tuple((l['base'], min(max_label_len, l['run_length'])))]
                             for l in s.labels), dtype=int, count=len(s.labels))
        y = y.reshape(y.shape + (1,))
        if not sparse_labels:
            from keras.utils.np_utils import to_categorical
            y = to_categorical(y, num_classes=n_classes)
        return x, y


class VarQueue(list):

    @property
    def last_pos(self):
        if len(self) == 0:
            return None
        else:
            return self[-1].pos

    def write(self, vcf_fh):
        if len(self) > 1:
            are_dels = all(len(x.ref) == 2 for x in self)
            are_same_ref = len(set(x.chrom for x in self)) == 1
            if are_dels and are_same_ref:
                name = self[0].chrom
                pos = self[0].pos
                ref = ''.join((x.ref[0] for x in self))
                ref += self[-1].ref[-1]
                alt = ref[0]

                merged_var = vcf.Variant(name, pos, ref, alt, info=info)
                vcf_fh.write_variant(merged_var)
            else:
                raise ValueError('Cannot merge variants: {}.'.format(self))
        elif len(self) == 1:
            vcf_fh.write_variant(self[0])
        del self[:]


class VCFChunkWriter(object):
    def __init__(self, fname, chrom, start, end, reference_fasta, label_decoding):
        vcf_region_str = '{}:{}-{}'.format(chrom, start, end) #is this correct?
        self.label_decoding = label_decoding
        self.logger = get_named_logger('VCFWriter')
        self.logger.info("Writing variants for {}".format(vcf_region_str))

        vcf_meta = ['region={}'.format(vcf_region_str)]
        self.writer = vcf.VCFWriter(fname, meta_info=vcf_meta)
        self.ref_fasta = pysam.FastaFile(reference_fasta)

    def __enter__(self):
        self.writer.__enter__()
        self.ref_fasta.__enter__()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.writer.__exit__(exc_type, exc_val, exc_tb)
        self.ref_fasta.__exit__(exc_type, exc_val, exc_tb)

    def add_chunk(self, sample, pred):
        # Write consensus alts to vcf
        cursor = 0
        var_queue = list()
        ref_seq = self.ref_fasta.fetch(sample.ref_name)
        for pos, grp in itertools.groupby(sample.positions['major']):
            end = cursor + len(list(grp))
            alt = ''.join(self.label_decoding[x] for x in pred[cursor:end]).replace(_gap_, '')
            # For simple insertions and deletions in which either
            #   the REF or one of the ALT alleles would otherwise be
            #   null/empty, the REF and ALT Strings must include the
            #   base before the event (which must be reflected in
            #   the POS field), unless the event occurs at position
            #   1 on the contig in which case it must include the
            #   base after the event
            if alt == '':
                # deletion
                if pos == 0:
                    # the "unless case"
                    ref = ref_seq[1]
                    alt = ref_seq[1]
                else:
                    # the usual case
                    pos = pos - 1
                    ref = ref_seq[pos:pos+2]
                    alt = ref_seq[pos]
            else:
                ref = ref_seq[pos]

            # Merging of variants produced by considering major.{minor} positions
            # These are of the form:
            #    X -> Y          - subs
            #    prev.X -> prev  - deletion
            #    X -> Xyy..      - insertion
            # In the second case we may need to merge variants from consecutive
            # major positions.
            if alt == ref:
                self.write(var_queue)
                var_queue = list()
            else:
                var = vcf.Variant(sample.ref_name, pos, ref, alt)
                if len(var_queue) == 0 or pos - var_queue[-1].pos == 1:
                    var_queue.append(var)
                else:
                    self.write(var_queue)
                    var_queue = [var]
            cursor = end
        self.write(var_queue)


    def write(self, var_queue):
        if len(var_queue) > 1:
            are_dels = all(len(x.ref) == 2 for x in var_queue)
            are_same_ref = len(set(x.chrom for x in var_queue)) == 1
            if are_dels and are_same_ref:
                name = var_queue[0].chrom
                pos = var_queue[0].pos
                ref = ''.join((x.ref[0] for x in var_queue))
                ref += var_queue[-1].ref[-1]
                alt = ref[0]

                merged_var = vcf.Variant(name, pos, ref, alt)
                self.writer.write_variant(merged_var)
            else:
                raise ValueError('Cannot merge variants: {}.'.format(var_queue))
        elif len(var_queue) == 1:
            self.writer.write_variant(var_queue[0])


def run_prediction(output, bam, regions, model, model_file, rle_ref,
        read_fraction, chunk_len, chunk_ovlp, batch_size=200,
        save_features=False, tag_name=None, tag_value=None,
        tag_keep_missing=False, enable_chunking=True):
    """Inference worker."""

    logger = get_named_logger('PWorker')

    remainder_regions = list()
    def sample_gen():
        # chain all samples whilst dispensing with generators when done
        #   (they hold the feature vector in memory until they die)
        for region in regions:
            data_gen = SampleGenerator(
                bam, region, model_file, rle_ref, read_fraction,
                chunk_len=chunk_len, chunk_overlap=chunk_ovlp,
                tag_name=tag_name, tag_value=tag_value,
                tag_keep_missing=tag_keep_missing,
                enable_chunking=enable_chunking)
            yield from data_gen.samples
            remainder_regions.extend(data_gen._quarantined)
    batches = background_generator(
        grouper(sample_gen(), batch_size), 10
    )

    total_region_mbases = sum(r.size for r in regions) / 1e6
    logger.info("Running inference for {:.1f}M draft bases.".format(total_region_mbases))

    with DataStore(output, 'a') as ds:
        mbases_done = 0

        t0 = now()
        tlast = t0
        for data in batches:
            x_data = np.stack([x.features for x in data])
            class_probs = model.predict_on_batch(x_data)
            mbases_done += sum(x.span for x in data) / 1e6
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
                ds.write_sample(Sample(**sample_d))

    logger.info("All done, {} remainder regions.".format(len(remainder_regions)))
    return remainder_regions


def predict(args):
    """Inference program."""
    os.environ["TF_CPP_MIN_LOG_LEVEL"]="2"
    from keras.models import load_model
    from keras import backend as K

    args.regions = get_regions(args.bam, region_strs=args.regions)
    logger = get_named_logger('Predict')
    logger.info('Processing region(s): {}'.format(' '.join(str(r) for r in args.regions)))

    # write class names to output
    with DataStore(args.model) as ds:
        meta = ds.meta
    with DataStore(args.output, 'w') as ds:
        ds.update_meta(meta)

    logger.info("Setting tensorflow threads to {}.".format(args.threads))
    K.tf.logging.set_verbosity(K.tf.logging.ERROR)
    K.set_session(K.tf.Session(
        config=K.tf.ConfigProto(
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
        args.output, args.bam, regions, model, args.model, args.rle_ref, args.read_fraction,
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
                args.output, args.bam, [region[0]], model, args.model, args.rle_ref, args.read_fraction,
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


def process_labels(label_counts, max_label_len=10):
    """Create map from full labels to (encoded) truncated labels.

    :param label_counrs: `Counter` obj of label counts.
    :param max_label_len: int, maximum label length, longer labels will be truncated.
    :returns:
    :param label_encoding: {label: int encoded label}.
    :param sparse_labels: bool, create sparse labels.
    :param n_classes: int, number of label classes.
    :returns: ({label: int encoding}, [label decodings], `Counter` of truncated counts).
    """
    logger = get_named_logger('Labelling')

    old_labels = [k for k in label_counts.keys()]
    if type(old_labels[0]) == tuple:
        new_labels = (l[1] * decoding[l[0]].upper() for l in old_labels)
    else:
        new_labels = [l for l in old_labels]

    if max_label_len < np.inf:
        new_labels = [l[:max_label_len] for l in new_labels]

    old_to_new = dict(zip(old_labels, new_labels))
    label_decoding = list(sorted(set(new_labels)))
    label_encoding = { l: label_decoding.index(old_to_new[l]) for l in old_labels}
    logger.info("Label encoding dict is:\n{}".format('\n'.join(
        '{}: {}'.format(k, v) for k, v in label_encoding.items()
    )))

    new_counts = Counter()
    for l in old_labels:
        new_counts[label_encoding[l]] += label_counts[l]
    logger.info("New label counts {}".format(new_counts))

    return label_encoding, label_decoding, new_counts


def train(args):
    """Training program."""
    train_name = args.train_name
    mkdir_p(train_name, info='Results will be overwritten.')

    logger = get_named_logger('Training')
    logger.debug("Loading datasets:\n{}".format('\n'.join(args.features)))

    sparse_labels = not args.balanced_weights

    args.validation = args.validation_features if args.validation_features is not None else args.validation_split

    batcher = TrainBatcher(args.features, args.max_label_len, args.validation,
                           args.seed, sparse_labels, args.batch_size, threads=args.threads_io)

    if args.balanced_weights:
        n_labels = sum(batcher.label_counts.values())
        n_classes = len(batcher.label_counts)
        class_weight = {k: float(n_labels)/(n_classes * count) for (k, count) in batcher.label_counts.items()}
        class_weight = np.array([class_weight[c] for c in sorted(class_weight.keys())])
    else:
        class_weight = None

    h = lambda d, i: d[i] if d is not None else 1
    logger.info("Label statistics are:\n{}".format('\n'.join(
        '{} ({}) {} (w. {:9.6f})'.format(i, l, batcher.label_counts[i], h(class_weight, i))
            for i, l in enumerate(batcher.label_decoding)
    )))

    import tensorflow as tf
    with tf.device('/gpu:{}'.format(args.device)):
        run_training(
            train_name, batcher, model_fp=args.model, epochs=args.epochs,
            class_weight=class_weight, n_mini_epochs=args.mini_epochs,
            threads_io=args.threads_io)

    # stop batching threads
    logger.info("Training finished.")
