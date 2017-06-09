"""
Train a neural network error model for correcting sequences.
"""

import argparse
import h5py

import numpy as np
from keras.models import Sequential
from keras.callbacks import CSVLogger, ModelCheckpoint
from keras.utils import to_categorical
from keras.layers import Dense, LSTM, Dropout
from medaka.util.generators import (serve_sample, serve_sample_batch,
                                    serve_data, serve_data_batch)
from medaka.util.sequences import get_benchmark, write_test_sequences
from medaka.util.reshape import trim_to_step_multiple


def preprocess_tvt_data(data, label, lims, batch_size, window_size):
    """Split data into train, validate and test sets.

    :param data: entire feature array
    :param label: entire label sequence
    :param lims: indices separating data into train, validate and test subsets
    :param batch_size: int batch size
    :param window_size: int (odd) window of positions flanking predicted base
    :returns: tuple (X_train, y_train, train_steps, X_validate, y_validate,
              validate_steps, X_test, y_test, test_steps)

              - a data set, encoded label and exact number of steps
                required to consume the data in multipes of batch_size
              - labels are converted to one-hot encoding with `to_categorical`
    """
    params = (batch_size, window_size)
    X_train, train_steps = trim_to_step_multiple(data[lims[0]:lims[1]], *params)
    y_train, _ = trim_to_step_multiple(
        to_categorical(label[lims[0]:lims[1]], num_classes=6), *params)
    X_validate, validate_steps = trim_to_step_multiple(data[lims[1]:lims[2]], *params)
    y_validate, _ = trim_to_step_multiple(
        to_categorical(label[lims[1]:lims[2]], num_classes=6), *params)
    X_test, test_steps = trim_to_step_multiple(data[lims[2]:lims[3]], *params)
    y_test, _ = trim_to_step_multiple(
        to_categorical(label[lims[2]:lims[3]], num_classes=6), *params)
    return (X_train, y_train, train_steps,
            X_validate, y_validate, validate_steps,
            X_test, y_test, test_steps)


def build_lstm(batch_shape):
    """Define model architecture using keras Sequential Model API

    :returns: keras model object

    .. note::

        see `keras <https://keras.io>`_ for full documentation
    """
    model = Sequential()
    model.add(LSTM(32, return_sequences=True, batch_input_shape=batch_shape))
    model.add(Dropout(0.5))
    model.add(LSTM(32))
    model.add(Dropout(0.5))
    model.add(Dense(6, activation='softmax'))
    model.compile(optimizer='rmsprop',
                  loss='categorical_crossentropy',
                  metrics=['accuracy'])
    return model


def calculate_tvt_limits(data_size, tvt_ratio, data_pp):
    """
    Calculate indices separating training, validate and test set

    :param data_size: int number of rows in feature array
    :param tvt_ratio: list ratio of [training, validatation, test] set sizes
    :param data_pp: proportion of data to use (e.g. 0.5 for first half of data)
    :returns: 4-element numpy array of indices
    """
    tvt_ratio = np.array(tvt_ratio, 'float')
    lims = [0] +  [int(l) for l in (np.cumsum(tvt_ratio) / np.sum(tvt_ratio)) *
                   (data_size * data_pp)]
    return lims


def train_lstm(data, label, tvt_ratio, data_pp, batch_size, window_size,
               out_prefix, epochs):
    """
    Train LSTM model.
    """
    # indices marking limits of train, validate, test sets.
    limits = calculate_tvt_limits(data.shape[0], tvt_ratio, data_pp)

    params = (batch_size, window_size)

    # get data split into train, validate and test subsets
    (X_train, y_train, train_steps,
     X_validate, y_validate, validate_steps,
     X_test, y_test, test_steps) = preprocess_tvt_data(data,
                                                       label,
                                                       limits,
                                                       *params)

    # establish the baseline accuracy of the reference
    train_benchmark = get_benchmark(X_train, y_train,
                                    train_steps, *params)
    validate_benchmark = get_benchmark(X_validate, y_validate,
                                       validate_steps, *params)
    test_benchmark = get_benchmark(X_test, y_test,
                                   test_steps, *params)
    print('benchmark accuracy on training set {}, '
          '(reference gives correct base).'.format(train_benchmark))
    print('benchmark accuracy on validation set {}, '
          '(reference gives correct base).'.format(validate_benchmark))
    print('benchmark accuracy on test set {}, '
          '(reference gives correct base).'.format(test_benchmark))

    np.random.seed(0)

    # setup model
    feature_length = X_train.shape[1]
    batch_shape = (batch_size, window_size, feature_length)
    model = build_lstm(batch_shape)
    model.summary()

    # setup data generators
    train_generator = serve_sample_batch(X_train, y_train, *params)
    validate_generator = serve_sample_batch(X_validate, y_validate, *params)
    test_generator = serve_sample_batch(X_test, y_test, *params)

    # setup callbacks
    training_outfile = '_'.join([out_prefix, 'training_history.txt'])
    csv_logger = CSVLogger(training_outfile, separator='\t')
    checkpoint_out = '_'.join([out_prefix, 'checkpoint_model.epoch{epoch:02d}.h5'])
    model_checkpoint = ModelCheckpoint(checkpoint_out)
    callbacks = [csv_logger, model_checkpoint]

    # train model
    model.fit_generator(train_generator, steps_per_epoch=train_steps,
                        epochs=epochs, validation_data=validate_generator,
                        validation_steps=validate_steps, callbacks=callbacks)

    # report final results on test data
    loss_and_metrics = model.evaluate_generator(test_generator, test_steps)
    print('test set results {}'.format(loss_and_metrics))

    model_outfile = '_'.join([out_prefix, 'neural_network.h5'])
    model.save(model_outfile)

    write_test_sequences(X_test, model, test_steps, *params)


def get_parser():
    parser = argparse.ArgumentParser(
        description="""Train model using preprocessed training data.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--out_prefix', type=str, default='',
                        help='name for outputs.')
    parser.add_argument('--epochs', type=int, default=10,
                        help='number of training epochs.')
    parser.add_argument('--window_size', type=int, default=3,
                        help='width of pileup window fed to trainer.')
    parser.add_argument('--batch_size', type=int, default=256,
                        help='number of training samples per batch.')
    parser.add_argument('--data_pp', type=float, default=1.0,
                        help='proportion of data to process (from start).')
    parser.add_argument('--tvt_ratio', nargs=3, default=[64, 16, 20],
                        metavar=('train','validate','test'),
                        help='select train:validate:test ratio.')
    parser.add_argument('datafile', help='.h5 generated by medaka_prepare')
    return parser


def load_data(datafile, data_pp):
    with h5py.File(datafile) as f:
        limit = int(f['data'].shape[0] * data_pp)
        data = f['data'][:limit]
        label = f['label'][:limit]
    return data, label


def main():
    args = get_parser().parse_args()
    data, label = load_data(args.datafile, args.data_pp)
    train_lstm(
        data,
        label,
        args.tvt_ratio,
        args.data_pp,
        args.batch_size,
        args.window_size,
        args.out_prefix,
        args.epochs
    )


if __name__ == '__main__':
    main()
