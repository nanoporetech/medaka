import tempfile
import unittest

import numpy as np
import tensorflow as tf

from medaka import keras_ext
from medaka.datastore import DataStore

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


class ModelAndData(object):
    input_dim = 2
    num_hidden = 4
    num_classes = 2
    batch_size = 5
    num_train = 20
    num_test = 20
    # note DataStore only reads certain groups as meta
    model_meta = dict(medaka_model_name='dumbmodel')
    callback_opts = dict(verbose=0, save_best_only=True, mode='max')
    metrics = ['binary_accuracy']

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


class TestCheckpoint(unittest.TestCase, ModelAndData):

    def test_000_checkpoint_saves_meta(self):
        (X_train, y_train), (X_test, y_test) = self._get_data_callbacks()
        y_train = tf.keras.utils.to_categorical(y_train)
        y_test = tf.keras.utils.to_categorical(y_test)
        model_file = tempfile.NamedTemporaryFile()
        model_fname = model_file.name
        callbacks = [
            keras_ext.ModelMetaCheckpoint(
                self.model_meta, model_fname, monitor='val_{}'.format(self.metrics[0]), **self.callback_opts)]
        model = self._get_model()
        model.fit(
            X_train, y_train, validation_data=(X_test, y_test),
            batch_size=2, epochs=5, callbacks=callbacks)

        with DataStore(model_fname, 'r') as ds:
           meta = ds.meta
           self.assertDictEqual(meta, self.model_meta)


class TestBatcher(unittest.TestCase, ModelAndData):
    pass
