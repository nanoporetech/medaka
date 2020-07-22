import copy
import functools
import os
import tempfile
import unittest

import numpy as np
import tensorflow as tf

import medaka.models
import medaka.keras_ext
import medaka.datastore
import medaka.training

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
        classify_activation='softmax', time_steps=None)
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

        self.val_model_file = tempfile.NamedTemporaryFile(suffix='.h5')
        self.val_model_fname = self.val_model_file.name
        self.train_model_file = tempfile.NamedTemporaryFile(suffix='.h5')
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
                **self.callback_opts)
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


class ModelAndDataTF(unittest.TestCase):
    input_dim = 2
    num_hidden = 4
    num_classes = 2
    batch_size = 5
    num_train = 20
    num_test = 20

    model_function = functools.partial(medaka.models.build_model,
        feature_len=10, num_classes=6, gru_size=128,
        classify_activation='softmax', time_steps=None)
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

        #self.val_model_file = tempfile.TemporaryDirectory()
        #self.val_model_fname = self.val_model_file.name
        self.val_model_fname = tempfile.mkdtemp() 
        #self.train_model_file = tempfile.TemporaryDirectory()
        #self.train_model_fname = self.train_model_file.name
        self.train_model_fname = tempfile.mkdtemp() 
        
        self.tb_dir = tempfile.mkdtemp()
        callbacks = [
            medaka.keras_ext.ModelMetaCheckpointTF(
                self.model_meta, self.val_model_fname,
                monitor='val_{}'.format(self.metrics[0]),
                **self.callback_opts),
            medaka.keras_ext.ModelMetaCheckpointTF(
                self.model_meta, self.train_model_fname,
                monitor=self.metrics[0],
                **self.callback_opts)
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


class TestCallbacksTF(ModelAndDataTF):

    def test_000_checkpoint_saves_meta(self):
        fname = self.val_model_fname + '.tar.gz'
        with medaka.datastore.ModelStoreTF(fname) as ds:
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
            medaka.keras_ext.ModelMetaCheckpointTF(meta, model_fname)

training_features = os.path.join(os.path.dirname(__file__), 'data', 'training_features.hdf5')


class TestSequenceBatcher(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.batch_size = 5

        self.tb = medaka.training.TrainBatcher(
            [training_features], batch_size=self.batch_size)

        self.n_training_samples = len(self.tb.train_samples)
        self.n_valid_samples = len(self.tb.valid_samples)

        self.sb_train = medaka.keras_ext.SequenceBatcher(
            self.tb, dataset='train')

        self.sb_valid = medaka.keras_ext.SequenceBatcher(
            self.tb, dataset='validation')

        self.mini_epochs = 3

        self.sb_mini = medaka.keras_ext.SequenceBatcher(
            self.tb, mini_epochs=self.mini_epochs)

    def test_000_correct_batch_size(self):
        first_batch = self.sb_train.__getitem__(0)
        self.assertEqual(len(first_batch[0]), self.batch_size) # features
        self.assertEqual(len(first_batch[1]), self.batch_size) # labels

    def test_001_correct_num_train_batches(self):
        self.assertEqual(len(self.sb_train),
            self.n_training_samples // self.batch_size)

    def test_002_correct_num_validation_batches(self):
        self.assertEqual(len(self.sb_valid),
            self.n_valid_samples // self.batch_size)

    def test_003_shuffle_occurs_on_epoch_end(self):
        pre_shuffle = copy.deepcopy(self.sb_train.data)
        self.sb_train.on_epoch_end()
        post_shuffle = self.sb_train.data
        self.assertNotEqual(pre_shuffle, post_shuffle)
        self.assertEqual(set(pre_shuffle), set(post_shuffle))

    def test_004_correct_num_batches_with_miniepochs(self):
        self.assertEqual(len(self.sb_mini),
            self.n_training_samples // self.batch_size // self.mini_epochs)

    def test_005_validation_init_with_miniepochs_gt_1_fails(self):
        with self.assertRaises(ValueError) as context:
            sb = medaka.keras_ext.SequenceBatcher(
                self.tb, dataset='validation', mini_epochs=3)

    def test_006_sb_init_with_incorrect_dataset_name_fails(self):
        with self.assertRaises(ValueError) as context:
            sb = medaka.keras_ext.SequenceBatcher(
                self.tb, dataset='sasquatch') # only 'train' or 'validation' permitted

    def test_007_sb_with_random_seed(self):
        sb = medaka.keras_ext.SequenceBatcher(self.tb, seed=2)
