import glob
import h5py
import itertools
import json
import numpy as np
import os
import sys
from collections import OrderedDict
from timeit import default_timer as now

from keras import backend as K
from keras.models import Sequential, load_model, save_model
from keras.layers import Dense, GRU
from keras.layers.wrappers import Bidirectional, TimeDistributed
from keras.callbacks import ModelCheckpoint, CSVLogger, TensorBoard, EarlyStopping

import logging
logger = logging.getLogger(__name__)


from medaka.tview import _ref_gap_, _gap_, load_pileup, generate_pileup_chunks, rechunk
from medaka import features
from medaka.common import mkdir_p

_encod_path_ = 'medaka_label_encoding'


def generate_features(pileup_gen, coverage_filter=features.coverage_filter,
                      feature_func=features.counts, training=True):
    """Generate training features and labels from a pileup generator

    :param pileup_gen: generator of (pileups, labels) tuples;
        pileups: list of `Pileup` objects
        labels: np.array of labels or None.
    :param feature_func: function to generate features
    :param training: bool, whether features suitable for training should be
        generated. This requires a labelled pileup.

    :returns dict of (features, labels) If `training` is False,
        `labels` will be None.
    """

    all_data = OrderedDict()  # so chunks remain in order
    t0 = now()
    for pileups, labels in coverage_filter(pileup_gen):
        pos = pileups[0].positions
        ref_name = pileups[0].ref_name
        grp = '{}:{}-{}'.format(ref_name, pos['major'][0], pos['major'][-1] + 1)
        if training and labels is None:
            raise ValueError("Cannot train without labels")
        all_data[grp] = (feature_func(pileups), labels)
    t1 = now()
    logging.info("Generating features took {:.3f}s.".format(t1 - t0))
    return all_data


def build_model(chunk_size, feature_len, num_classes, gru_size=128):
    """Builds a bidirectional GRU model"""

    model = Sequential()
    # Bidirectional wrapper takes a copy of the first argument and reverses
    #   the direction. Weights are independent between components.
    model.add(Bidirectional(
        GRU(gru_size, activation='tanh', return_sequences=True, name='gru1'), 
        input_shape=(chunk_size, feature_len))
    )
    model.add(Bidirectional(
        GRU(gru_size, activation='tanh', return_sequences=True, name='gru2'))
    )
    model.add(Dense(num_classes, activation='softmax', name='classify'))
    return model


def save_encoding(fname, encoding):
    """Save label encodings to json."""
    with open(fname, 'w') as json_file:
        json_file.write(json.dumps({'encoding': encoding}))


def load_encoding(fname):
    with open(fname) as fh:
        encoding = json.load(fh)['encoding']
    return encoding


def write_encoding_to_hdf(fname, encoding):
    """Write model encoding to model with label encoding."""
    with h5py.File(fname) as h5:
        if _encod_path_ in h5:
            del h5[_encod_path_]
        h5[_encod_path_] = [s.encode() for s in encoding]


def save_model_hdf(fname, model, encoding):
    """Save model with label encoding."""
    # TODO: save feature function and coverage_filter function with model ?
    save_model(model, fname)
    write_encoding_to_hdf(fname, encoding)


def load_model_hdf(model_path, encoding_json=None, need_encoding=True):
    """Load a model from a .hdf file. If label encodings are not present,
    try to load them from enncoding_json."""
    # try to get structure + encoding from hdf, else use the
    # one provided (or raise an exception if it is not provided).
    m = load_model(model_path, custom_objects={'qscore': qscore})
    encoding = None
    with h5py.File(model_path, 'r') as h5:
        if _encod_path_ in h5:
            encoding = [s.decode() for s in h5[_encod_path_][()]]
            logging.info("Loaded encoding from {}.".format(model_path))

    if encoding is None and encoding_json is not None:
        encoding = load_encoding(encoding_json)
        logging.info("Loaded encoding from {}.".format(encoding_json))

    if encoding is None and need_encoding:
        raise KeyError("Could not find label encodings in the model, please provide an encoding json")

    return m, encoding


def save_feature_file(fname, data):
    """Save feature dictionary."""
    np.save(fname, data)


def load_feature_file(fname):
    """Load the result of `save_feature_file` back to the original
    representation.
    """
    data = OrderedDict()  # so chunks stay in order
    src = np.load(fname)
    for fname in src[()].keys():
        data[fname] = src[()][fname]
    return data


def qscore(y_true, y_pred):
    error = K.cast(K.not_equal(
        K.max(y_true, axis=-1), K.cast(K.argmax(y_pred, axis=-1), K.floatx())),
        K.floatx()
    )
    error = K.sum(error) / K.sum(K.ones_like(error))
    return -10.0 * 0.434294481 * K.log(error)


def run_training(train_name, x_train, y_train, encoding, model_data=None):
    """Run training."""
    data_dim = x_train.shape[2]
    timesteps = x_train.shape[1]
    num_classes = len(encoding)

    if model_data is None:
        model = build_model(timesteps, data_dim, num_classes)
    else:
        model, old_encoding = load_model_hdf(model_data, need_encoding=False)
        if old_encoding is not None:
            logging.info("Old label encoding was:\n{}".format('\n'.join(
                '{}: {}'.format(i, x) for i, x in enumerate(old_encoding)
            )))
        # TODO: should check model data dimensions match data dimensions

    logging.info("data_dim: {}, timesteps: {}, num_classes: {}".format(data_dim, timesteps, num_classes))
    logging.info("\n{}".format(model.summary()))

    save_model_hdf(os.path.join(train_name, 'model_structure.hdf'), model, encoding)
    encoding_name = os.path.join(train_name, '{}_label_encodings.json'.format(train_name))
    save_encoding(encoding_name, encoding)

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
        TensorBoard(log_dir=os.path.join(train_name, 'logs'),
                    histogram_freq=5, batch_size=100, write_graph=True,
                    write_grads=True, write_images=True)
    ]

    model.compile(
       loss='sparse_categorical_crossentropy',
       optimizer='rmsprop',
       metrics=['accuracy', qscore],
    )

    # maybe possible to increase batch_size for faster processing
    model.fit(
        x_train, y_train,
        batch_size=100, epochs=5000,
        validation_split=0.2,
        callbacks=callbacks,
    )

    # append encoding to the hdf files
    for hd5_model_path in glob.glob(os.path.join(train_name, '*.hdf5')):
        write_encoding_to_hdf(hd5_model_path, encoding)


def to_sparse_categorical(data):
    """Create a sparse categorical matrix from a list of label vectors

    :param data: list of label vectors.

    :returns: tuple (sparse categorical array, encoding list)
    """
    all_labels = set()
    for d in data:
        all_labels |= set(d)
    # sort to get reproducible ordering for identical datasets
    encoding = sorted(list(all_labels))
    decoding = {a: i for i, a in enumerate(encoding)}

    matrix = np.stack([
        [[decoding[x], ] for x in sample]
        for sample in data
    ])
    return matrix, encoding


def run_prediction(data, model_hdf, encoding_json=None, output_file='basecalls.fasta', batch_size=1500):
    """Run inference."""
    model, encoding = load_model_hdf(model_hdf, encoding_json, need_encoding=True)
    logging.info("Label encoding is:\n{}".format('\n'.join(
        '{}: {}'.format(i, x) for i, x in enumerate(encoding)
    )))
    logging.info('\n{}'.format(model.summary()))

    x_data = np.stack((x[0] for x in itertools.islice(data.values(), batch_size)))

    t0 = now()
    class_probs = model.predict(x_data, batch_size=batch_size, verbose=1)
    t1 = now()
    logging.info('Running network took {}s for data of shape {}'.format(t1 - t0, x_data.shape))

    t0 = now()
    count = 0
    # write out contig name and position in fasta, ref_name:start-end
    with open(output_file, 'w') as fasta:
        best = np.argmax(class_probs, -1)
        for key, seq in zip(data.keys(), (''.join(encoding[x] for x in sample) for sample in best)):
            for gap_sym in (_gap_, _ref_gap_):
                seq = seq.replace(gap_sym, '')
            fasta.write(">{}\n{}\n".format(key, seq))
            count += 1
        t1 = now()
        logging.info('Decoding took {}s for {} chunks.'.format(t1 - t0, count))


def train(args):
    """Training program."""
    train_name = args.train_name
    mkdir_p(train_name, info='Results will be overwritten and may use pregenerated features.')
    dataset_name = os.path.join(train_name, '{}_squiggles.npy'.format(train_name))
    logging.info("Using {} for feature storage/reading.".format(dataset_name))

    if not os.path.isfile(dataset_name):
        logging.info("Creating dataset. This may take a while.")
        chunks = rechunk(load_pileup(args.pileupdata))
        data = generate_features(chunks)
        save_feature_file(dataset_name, data)
        if args.features:
            logging.info("Stopping as only feature generation requested.")
            sys.exit(0)
    else:
        logging.info("Loading dataset from file.")
        data = load_feature_file(dataset_name)
    logging.info("Got {} pileup chunks for training.".format(len(data)))

    x_data, y_labels = zip(*data.values())
    # stack the individual samples into one big tensor
    x_data = np.stack(x_data)
    # find unique labels and create a sparse encoding
    y_labels, encoding = to_sparse_categorical(y_labels)

    logging.info("Label encoding is:\n{}".format('\n'.join(
        '{}: {}'.format(i, x) for i, x in enumerate(encoding)
    )))

    run_training(train_name, x_data, y_labels, encoding, model_data=args.model)


def predict(args):
    """Inference program."""
    if args.pileupdata:
        # TODO make use of start and end options when loading pileup from hdf ?
        logging.info("Loading pileup from file.")
        chunks = rechunk(load_pileup(args.pileupdata))
        data = generate_features(chunks, training=False)
    elif args.feature_file:
        logging.info("Loading features from file.")
        data = load_feature_file(args.feature_file)
    else:
        logging.info("Generating features from bam.")
        bam, ref_fasta, ref_name = args.alignments
        chunks = rechunk(generate_pileup_chunks((bam,), ref_fasta, ref_name=ref_name, start=args.start, end=args.end))
        data = generate_features(chunks, training=False)

    logging.info("Got {} pileup chunks for inference.".format(len(data)))
    run_prediction(data, args.model, encoding_json=args.encoding, output_file=args.output_fasta)
