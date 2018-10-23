import glob
import h5py
import inspect
import itertools
import numpy as np
import os
import pysam
from collections import Counter
from functools import partial
from timeit import default_timer as now

from keras import backend as K
from keras.models import Sequential, load_model
from keras.layers import Dense, GRU, Dropout
from keras.layers.wrappers import Bidirectional
from keras.callbacks import ModelCheckpoint, CSVLogger, TensorBoard, EarlyStopping

# set some Tensorflow session options
"""
from keras.backend.tensorflow_backend import set_session
import tensorflow as tf
config = tf.ConfigProto()
config.gpu_options.per_process_gpu_memory_fraction = 0.3
set_session(tf.Session(config=config))
"""

import logging
logger = logging.getLogger(__name__)


from medaka import features
from medaka.common import (encode_sample_name, _label_counts_path_,
                           _label_decod_path_, chain_thread_safe, decoding,
                           get_sample_index_from_files, grouper,
                           load_sample_from_hdf, load_yaml_data,
                           _model_opt_path_, mkdir_p, Sample, sample_to_x_y,
                           threadsafe_generator, write_samples_to_hdf,
                           write_sample_to_hdf, write_yaml_data,
                           yield_from_feature_files, load_from_feature_files,
                           gen_train_batch, _label_batches_path_, _feature_batches_path_,
                           yield_batches_from_hdfs)




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


def build_model(chunk_size, feature_len, num_classes, gru_size=128, input_dropout=0.0,
                inter_layer_dropout=0.0, recurrent_dropout=0.0):
    """Builds a bidirectional GRU model
    :param chunk_size: int, number of pileup columns in a sample.
    :param feature_len: int, number of features for each pileup column.
    :param num_classes: int, number of output class labels.
    :param gru_size: int, size of each GRU layer.
    :param input_dropout: float, fraction of the input feature-units to drop.
    :param inter_layer_dropout: float, fraction of units to drop between layers.
    :param recurrent_dropout: float, fraction of units to drop within the recurrent state.
    :returns: `keras.models.Sequential` object.
    """

    model = Sequential()

    gru1 = GRU(gru_size, activation='tanh', return_sequences=True, name='gru1',
               dropout=input_dropout, recurrent_dropout=recurrent_dropout)
    gru2 = GRU(gru_size, activation='tanh', return_sequences=True, name='gru2',
               dropout=inter_layer_dropout, recurrent_dropout=recurrent_dropout)

    # Bidirectional wrapper takes a copy of the first argument and reverses
    #   the direction. Weights are independent between components.
    model.add(Bidirectional(gru1, input_shape=(chunk_size, feature_len)))

    model.add(Bidirectional(gru2, input_shape=(chunk_size, feature_len)))

    if inter_layer_dropout > 0:
        model.add(Dropout(inter_layer_dropout))

    model.add(Dense(num_classes, activation='softmax', name='classify'))

    return model


def qscore(y_true, y_pred):
    error = K.cast(K.not_equal(
        K.max(y_true, axis=-1), K.cast(K.argmax(y_pred, axis=-1), K.floatx())),
        K.floatx()
    )
    error = K.sum(error) / K.sum(K.ones_like(error))
    return -10.0 * 0.434294481 * K.log(error)


def run_training(train_name, sample_gen, valid_gen, n_batch_train, n_batch_valid, label_decoding,
                 timesteps, feat_dim, model_fp=None, epochs=5000, batch_size=100, class_weight=None):
    """Run training."""

    logging.info("Got {} batches for training, {} for validation.".format(n_batch_train, n_batch_valid))

    num_classes = len(label_decoding)

    if model_fp is None:
        model_kwargs = { k:v.default for (k,v) in inspect.signature(build_model).parameters.items()
                         if v.default is not inspect.Parameter.empty}
    else:
        model_kwargs = load_yaml_data(model_fp, _model_opt_path_)
        assert model_kwargs is not None

    opt_str = '\n'.join(['{}: {}'.format(k,v) for k, v in model_kwargs.items()])
    logging.info('Building model with: \n{}'.format(opt_str))
    model = build_model(timesteps, feat_dim, num_classes, **model_kwargs)

    if model_fp is not None and os.path.splitext(model_fp)[-1] != '.yml':
        old_model = load_model(model_fp, custom_objects={'qscore': qscore})
        old_model_feat_dim = old_model.get_input_shape_at(0)[2]
        assert old_model_feat_dim == feat_dim
        logging.info("Loading weights from {}".format(model_fp))
        model.load_weights(model_fp)

    logging.info("feat_dim: {}, timesteps: {}, num_classes: {}".format(feat_dim, timesteps, num_classes))
    logging.info("\n{}".format(model.summary()))

    model_details = {_label_decod_path_: label_decoding,
                     _model_opt_path_: model_kwargs}
    write_yaml_data(os.path.join(train_name, 'training_config.yml'), model_details)

    opts = dict(verbose=1, save_best_only=True, mode='max')
    callbacks = [
        # Best model according to training set accuracy
        ModelCheckpoint(os.path.join(train_name, 'model.best.hdf5'),
                        monitor='acc', **opts),
        # Best model according to validation set accuracy
        ModelCheckpoint(os.path.join(train_name, 'model.best.val.hdf5'),
                        monitor='val_acc', **opts),
        # Best model according to validation set qscore
        ModelCheckpoint(os.path.join(train_name, 'model.best.val.qscore.hdf5'),
                        monitor='val_qscore', **opts),
        # Checkpoints when training set accuracy improves
        ModelCheckpoint(os.path.join(train_name, 'model-improvement-{epoch:02d}-{acc:.2f}.hdf5'),
                        monitor='acc', **opts),
        # Stop when no improvement, patience is number of epochs to allow no improvement
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
        logging.info("Using weighted_categorical_crossentropy loss function")
    else:
        loss = 'sparse_categorical_crossentropy'
        logging.info("Using {} loss function".format(loss))

    model.compile(
       loss=loss,
       optimizer='rmsprop',
       metrics=['accuracy', qscore],
    )
    # fit generator
    model.fit_generator(
        sample_gen, steps_per_epoch=n_batch_train,
        validation_data=valid_gen, validation_steps=n_batch_valid,
        max_queue_size=8, workers=8, use_multiprocessing=False,
        epochs=epochs,
        callbacks=callbacks,
        class_weight=class_weight,
    )

    # append label_decoding and model building options to the hdf files
    for hd5_model_path in glob.glob(os.path.join(train_name, '*.hdf5')):
            write_yaml_data(hd5_model_path, model_details)


def run_prediction(sample_gen, model, label_decoding, reference_fasta, output_prefix='consensus_chunks',
                   batch_size=200, predictions_file=None, n_samples=None):
    """Run inference."""
    batches = grouper(sample_gen, batch_size)

    if predictions_file is not None:
        pred_h5 = h5py.File(predictions_file, 'w')

    first_batch = True
    n_samples_done = 0
    with open('{}.fa'.format(output_prefix), 'w') as fasta_out, \
        vcf.VCFWriter('{}.vcf'.format(output_prefix)) as vcf_out \
        pysam.FastaFile(reference_fasta) as ref_fasta:
        
        t0 = now()
        tlast = t0
        for data in batches:
            x_data = np.stack((x.features for x in data))

            if first_batch:
                logging.info('Updating model state with first batch')
                class_probs = model.predict(x_data, batch_size=batch_size, verbose=1)
                first_batch = False

            class_probs = model.predict(x_data, batch_size=batch_size, verbose=1)
            n_samples_done += x_data.shape[0]
            t1 = now()
            if t1 - tlast > 10:
                tlast = t1
                if n_samples is not None:  # we know how many samples we will be processing
                    msg = '{:.1%} Done ({}/{} samples) in {:.1f}s'
                    logging.info(msg.format(n_samples_done / n_samples, n_samples_done, n_samples, t1 - t0))
                else:
                    logging.info('Done {} samples in {:.1f}s'.format(n_samples_done, t1 - t0, x_data.shape))
             best = np.argmax(class_probs, -1)
 
             count = 0
             best = np.argmax(class_probs, -1)
             for sample, prob, pred in zip(data, class_probs, best):
                 # write consensus to fasta TODO: remove
                 seq = ''.join(label_decoding[x] for x in pred).replace(_gap_, '')
                 key = encode_sample_name(sample)
                 fasta_out.write(">{}\n{}\n".format(key, seq))

                 # write out positions and predictions for later analysis
                 if predictions_file is not None:
                     sample_d = sample._asdict()
                     sample_d['label_probs'] = prob
                     sample_d['features'] = None  # to keep file sizes down
                     write_sample_to_hdf(Sample(**sample_d), pred_h5)
                 count += 1

                 # Write consensus alts to vcf
                 ref_seq = ref_fasta.fetch(sample.ref_name)
                 cursor = 0
                 for pos, grp in groupby(sample.positions['major']):
                     end = cursor + len(grp)
                     alt = ''.join(label_decoding[x] for x in pred[cursor:end]).replace(_gap_, '')
                     ref = ref_seq[pos]
                     vcf_out.write(vcf.Variant(sample.ref_name), pos, ref, alt)
                     cursor = end
 
    if predictions_file is not None:
        pred_h5.close()
        write_yaml_data(predictions_file, {_label_decod_path_: label_decoding})

    logging.info('All done')


def process_labels(label_counts, max_label_len=10):
    """
    Create map from full labels to (encoded) truncated labels.
    """

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
    logging.info("Label encoding dict is:\n{}".format('\n'.join(
        '{}: {}'.format(k, v) for k, v in label_encoding.items()
    )))

    new_counts = Counter()
    for l in old_labels:
        new_counts[label_encoding[l]] += label_counts[l]
    logging.info("New label counts {}".format(new_counts))

    return label_encoding, label_decoding, new_counts


def train(args):
    """Training program."""
    train_name = args.train_name
    mkdir_p(train_name, info='Results will be overwritten.')

    logging.info("Loading datasets:\n{}".format('\n'.join(args.features)))

    # get counts of labels in training samples
    label_counts = Counter()
    for f in args.features:
        label_counts.update(load_yaml_data(f, _label_counts_path_))

    logging.info("Total labels {}".format(sum(label_counts.values())))


    is_batched = True
    handles = {fname: h5py.File(fname, 'r') for fname in args.features}
    for f, h5 in handles.items():
        has_batches = _label_batches_path_ in h5 and _feature_batches_path_ in h5
        msg = 'Found batches in {}.' if has_batches else 'No batches in {}.'
        logging.info(msg.format(f))
        is_batched = is_batched and has_batches

    if is_batched:

        batches = [ (fname, k) for (fname, fh) in handles.items() for k in fh[_feature_batches_path_]]

        logging.info("Got {} batches.".format(len(batches)))
        # check batch size using first batch
        test_f, test_batch = batches[0]
        h5 = handles[test_f]
        batch_shape = h5['{}/{}'.format(_feature_batches_path_, test_batch)].shape
        label_shape = h5['{}/{}'.format(_label_batches_path_, test_batch)].shape
        logging.info("Got {} batches with feat shape {}, label shape {}".format(len(batches), batch_shape, label_shape))
        batch_size, timesteps, feat_dim = batch_shape
        if not sum(label_counts.values()) != len(batches) * timesteps:
            raise ValueError('Label counts not consistent with number of batches')

        n_batch_train = int((1 - args.validation_split) * len(batches))
        train_batches = batches[:n_batch_train]
        valid_batches = batches[n_batch_train:]
        n_batch_valid = len(valid_batches)

        if args.balanced_weights:
            sparse_labels = False
        else:
            sparse_labels = True

        gen_train = yield_batches_from_hdfs(handles, train_batches, sparse_labels=sparse_labels,
                                            n_classes=len(label_counts))
        gen_valid = yield_batches_from_hdfs(handles, valid_batches, sparse_labels=sparse_labels,
                                            n_classes=len(label_counts))
        label_decoding = load_yaml_data(args.features[0], _label_decod_path_)

    else:
        sample_index = get_sample_index_from_files(args.features, max_samples=args.max_samples)
        refs = [k for k in sample_index.keys()]
        logging.info("Got the following references for training:\n{}".format('\n'.join(refs)))
        n_samples = sum([len(sample_index[k]) for k in sample_index])
        logging.info("Got {} pileup chunks for training.".format(n_samples))
        # get label encoding, given max_label_len
        logging.info("Max label length: {}".format(args.max_label_len if args.max_label_len is not None else 'inf'))
        label_encoding, label_decoding, label_counts = process_labels(label_counts, max_label_len=args.max_label_len)

        # create seperate generators of x,y for training and validation
        # shuffle samples before making split so each ref represented in training
        # and validation and batches are not biased to a particular ref.
        samples = [(d['key'], d['filename']) for ref in sample_index.values() for d in ref]
        np.random.shuffle(samples)  # shuffle in place
        n_samples_train = int((1 - args.validation_split) * len(samples))
        samples_train = samples[:n_samples_train]
        samples_valid = samples[n_samples_train:]
        # all batches need to be the same size, so pad samples_train and samples_valid
        n_extra_train = args.batch_size - len(samples_train) % args.batch_size
        n_extra_valid = args.batch_size - len(samples_valid) % args.batch_size
        samples_train += [samples_train[i] for i in np.random.choice(len(samples_train), n_extra_train)]
        samples_valid += [samples_valid[i] for i in np.random.choice(len(samples_valid), n_extra_valid)]

        msg = '{} training samples padded to {}, {} validation samples padded to {}'
        logging.info(msg.format(len(set(samples_train)), len(samples_train),
                                len(set(samples_valid)), len(samples_valid)))

        # load one sample to figure out timesteps and data dim
        key, fname = samples_train[0]
        with h5py.File(fname) as h5:
            first_sample = load_sample_from_hdf(key, h5)
        timesteps = first_sample.features.shape[0]
        feat_dim = first_sample.features.shape[1]

        if not sum(label_counts.values()) // timesteps == n_samples:
            raise ValueError('Label counts not consistent with number of samples')

        gen_train_samples = yield_from_feature_files(args.features, samples=itertools.cycle(samples_train))
        gen_valid_samples = yield_from_feature_files(args.features, samples=itertools.cycle(samples_valid))
        s2xy = partial(sample_to_x_y, encoding=label_encoding)
        gen_train = gen_train_batch(map(s2xy, gen_train_samples), args.batch_size, name='training')
        gen_valid = gen_train_batch(map(s2xy, gen_valid_samples), args.batch_size, name='validation')
        n_batch_train = len(samples_train) // args.batch_size
        n_batch_valid = len(samples_valid) // args.batch_size

    if args.balanced_weights:
        n_samples = sum(label_counts.values())
        n_classes = len(label_counts)
        class_weight = {k: float(n_samples)/(n_classes * count) for (k, count) in label_counts.items()}
        class_weight = np.array([class_weight[c] for c in sorted(class_weight.keys())])
    else:
        class_weight = None

    h = lambda d, i: d[i] if d is not None else 1
    logging.info("Label encoding is:\n{}".format('\n'.join(
        '{} ({}, {:9.6f}): {}'.format(i, label_counts[i], h(class_weight, i), l) for i, l in enumerate(label_decoding)
    )))

    run_training(train_name, gen_train, gen_valid, n_batch_train, n_batch_valid, label_decoding,
                 timesteps, feat_dim, model_fp=args.model, epochs=args.epochs, batch_size=args.batch_size,
                 class_weight=class_weight)

    for fh in handles.values():
        fh.close()


def predict(args):
    """Inference program."""
    n_samples = None
    if args.features:
        logging.info("Loading features from file for refs {}.".format(args.ref_names))
        index = get_sample_index_from_files(args.features)
        refs_n_samples = {r: len(l) for (r, l) in index.items()}
        ref_names = args.ref_names if args.ref_names is not None else refs_n_samples.keys()
        msg = "\n".join(["{}: {}".format(r, refs_n_samples[r]) for r in ref_names])
        logging.info("Number of samples per reference:\n{}\n".format(msg))
        data = yield_from_feature_files(args.features, ref_names=args.ref_names)
        n_samples = sum(refs_n_samples.values())
    else:
        logging.info("Generating features from bam.")
        raise ValueError("Unsupported args.features value.")

    # take a sneak peak at the first sample
    first_sample = next(data)
    data = chain_thread_safe([first_sample], data)
    timesteps = first_sample.features.shape[0]
    feat_dim = first_sample.features.shape[1]

    # re-build model with desired number of sequence steps per sample
    # derived from the features so this works wherever we got features from.
    model_data = {}
    for path in (_model_opt_path_, _label_decod_path_):
        opt = load_yaml_data(args.model, path)
        if opt is None and args.model_yml is not None:
            opt = load_yaml_data(args.model_yml, path)
        if opt is None:
            msg = '{} was not present in the input model, use the --model_yml to provide this data.'
            raise KeyError(msg.format(path))
        model_data[path] = opt

    num_classes = len(model_data[_label_decod_path_])

    opt_str = '\n'.join(['{}: {}'.format(k,v) for k, v in model_data[_model_opt_path_].items()])
    logging.info('Building model with: \n{}'.format(opt_str))
    model = build_model(timesteps, feat_dim, num_classes, **model_data[_model_opt_path_])
    logging.info("Loading weights from {}".format(args.model))
    model.load_weights(args.model)
    # check new model and old are consistent in size
    old_model = load_model(args.model, custom_objects={'qscore': qscore})
    get_feat_dim = lambda m: m.get_input_shape_at(0)[2]
    get_label_dim = lambda m: m.get_output_shape_at(-1)[-1]
    if not get_feat_dim(model) == get_feat_dim(old_model):
        msg = 'Incorrect feature dimension: got {}, model expects {}'
        raise ValueError(msg.format(feat_dim, get_feat_dim(old_model)))
    if not get_label_dim(model) == get_label_dim(old_model):
        msg = 'Incorrect label dimension: got {}, model expects {}'
        raise ValueError(msg.format(num_classes, get_label_dim(old_model)))

    logging.info("Label decoding is:\n{}".format('\n'.join(
        '{}: {}'.format(i, x) for i, x in enumerate(model_data[_label_decod_path_])
    )))
    logging.info('\n{}'.format(model.summary()))

    run_prediction(
        data, model, label_decoding=model_data[_label_decod_path_], args.reference
        output_file=args.output_fasta,
        predictions_file=args.output_probs,
        batch_size=args.batch_size, n_samples=n_samples,
    )
