import functools
import os
import tempfile
import unittest

import numpy as np
import tensorflow as tf

import medaka.models
import medaka.keras_ext
import medaka.datastore

def get_test_data(num_train=1000, num_test=500, input_shape=(10,),
                  output_shape=(2,),
                  classification=True, num_classes=2):
    # Taken from keras.utils.test_utils
    samples = num_train + num_test
    if classification:
        y = np.random.randint(0, num_classes, size=(samples,))
        X = np.zeros((samples,) + input_shape, dtype=np.float32)
        for i in range(samples):
            X[i] = np.random.normal(loc=y[i], scale=0.7, size=input_shape)
    else:
        y_loc = np.random.random((samples,))
        X = np.zeros((samples,) + input_shape, dtype=np.float32)
        y = np.zeros((samples,) + output_shape, dtype=np.float32)
        for i in range(samples):
            X[i] = np.random.normal(loc=y_loc[i], scale=0.7, size=input_shape)
            y[i] = np.random.normal(loc=y_loc[i], scale=0.7, size=output_shape)

    return (X[:num_train], y[:num_train]), (X[num_train:], y[num_train:])


class ModelAndData(unittest.TestCase):
    input_dim = 2
    num_hidden = 4
    num_classes = 2
    batch_size = 5
    num_train = 20
    num_test = 20

    model_function = functools.partial(medaka.models.build_model,
        feature_len=10, num_classes=6, gru_size=128,
        classify_activation='softmax', time_steps=None,
        allow_cudnn=None)
    model_meta = {
        'model_function': model_function,
        'label_scheme': None,
        'feature_encoder': None}

    callback_opts = dict(verbose=0, save_best_only=True, mode='max')
    metrics = ['binary_accuracy']

    @classmethod
    def setUpClass(self):
        (X_train, y_train), (X_test, y_test) = self._get_data_callbacks(self)
        y_train = tf.keras.utils.to_categorical(y_train)
        y_test = tf.keras.utils.to_categorical(y_test)

        self.val_model_file = tempfile.NamedTemporaryFile()
        self.val_model_fname = self.val_model_file.name
        self.train_model_file = tempfile.NamedTemporaryFile()
        self.train_model_fname = self.train_model_file.name
        self.tb_dir = tempfile.mkdtemp()
        callbacks = [
            medaka.keras_ext.ModelMetaCheckpoint(
                self.model_meta, self.val_model_fname,
                monitor='val_{}'.format(self.metrics[0]),
                **self.callback_opts),
            medaka.keras_ext.ModelMetaCheckpoint(
                self.model_meta, self.train_model_fname,
                monitor=self.metrics[0],
                **self.callback_opts),
            medaka.keras_ext.TrainValTensorBoard(log_dir=self.tb_dir)
        ]
        model = self._get_model(self)
        model.fit(
            X_train, y_train, validation_data=(X_test, y_test),
            batch_size=2, epochs=4, callbacks=callbacks)

    def _get_model(self):
        dense = tf.keras.layers.Dense
        layers = [
            dense(10, activation='relu', input_dim=self.input_dim),
            dense(self.num_classes, activation='softmax')
        ]
        model = tf.keras.models.Sequential(layers=layers)
        model.compile(optimizer='adam', loss='binary_crossentropy', metrics=self.metrics)
        return model


    def _get_data_callbacks(self):
        return get_test_data(
            num_train=self.num_train,
            num_test=self.num_test,
            input_shape=(self.input_dim,),
            classification=True,
            num_classes=self.num_classes)


class TestCallbacks(ModelAndData):

    def test_000_checkpoint_saves_meta(self):
        with medaka.datastore.DataStore(self.val_model_fname, 'r') as ds:
           model_func = ds.get_meta('model_function')
           # we can not simply assert equality of partial functions
           # https://bugs.python.org/issue3564
           # we need to test .func, .args and. keywords separately
           self.assertEqual(model_func.func, self.model_meta['model_function'].func)
           self.assertEqual(model_func.args, self.model_meta['model_function'].args)
           self.assertEqual(model_func.keywords, self.model_meta['model_function'].keywords)

    def test_001_checkpoint_raise_with_bad_meta(self):
        meta = {'thing1': None}
        model_file = tempfile.NamedTemporaryFile()
        model_fname = model_file.name
        with self.assertRaises(KeyError):
            medaka.keras_ext.ModelMetaCheckpoint(meta, model_fname)

    def test_010_tensorboard(self):
        # check that both training and validation logs contain correctly names metrics
        for group in ('training', 'validation'):
            grp_path = os.path.join(self.tb_dir, group)
            fname = [x for x in os.listdir(grp_path) if x != 'plugins' and not x.endswith('profile-empty')]
            self.assertTrue(len(fname) == 1)
            fname = fname[0]
            found = False
            with open(os.path.join(grp_path, fname), 'rb') as fh:
                found = any(b'epoch_binary_accuracy' in line for line in fh)
                self.assertTrue(found, "Failed to find correctly named metric for {}".format(group))


class TestBatcher(ModelAndData):
    pass
