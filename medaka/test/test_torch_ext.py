import copy
import functools
import os
import tempfile
import unittest

import numpy as np
import torch

import medaka.datastore
import medaka.training
import medaka.torch_ext
import medaka.models


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


class TestNamespace:

    class ModelAndData:
        input_dim = 2
        num_hidden = 4
        num_classes = 2
        batch_size = 2
        num_train = 20
        num_test = 20

        model_function = functools.partial(medaka.models.build_model_torch,
            feature_len=10, num_classes=6, gru_size=128,
            classify_activation='softmax', time_steps=None)
        model_meta = {
            'model_function': model_function,
            'label_scheme': None,
            'feature_encoder': None}

        class MockModel(medaka.models.TorchModel):
            def __init__(self, input_dim, num_hidden, num_classes):
                super().__init__()
                self.linear1 = torch.nn.Linear(input_dim, num_hidden)
                self.linear2 = torch.nn.Linear(num_hidden, num_classes)

            def forward(self, x):
                x = self.linear1(x)
                x = self.linear2(x)
                return x

        @classmethod
        def setUpClass(self):
            data_train, data_test = self._get_data(self)
            data_train, data_test = list(zip(*data_train)), list(zip(*data_test))

            self.dataloader_train = medaka.torch_ext.SequenceBatcher(
                medaka.torch_ext.Sequence(
                    data_train,
                    sampler_func=lambda x: x,
                    dataset="train",
                ),
                batch_size=self.batch_size,
            )
            self.dataloader_test = medaka.torch_ext.SequenceBatcher(
                medaka.torch_ext.Sequence(
                    data_test,
                    sampler_func=lambda x: x,
                    dataset="validation",
                ),
                batch_size=self.batch_size,
            )

            model_dict = self._get_model(self)
            _ = medaka.torch_ext.train_one_epoch(
                dataloader=self.dataloader_train,
                **model_dict
            )
            train_checkpoint = self.CheckpointType(self.model_meta, self.train_model_fname)
            train_checkpoint.on_epoch_end(0, model_dict['model'])

            _ = medaka.torch_ext.validate_one_epoch(
                model=model_dict["model"],
                dataloader=self.dataloader_test,
                loss_fn=model_dict["loss_fn"],
                metrics={"val_cat_acc": medaka.torch_ext.SparseCategoricalAccuracy(),}
            )
            valid_checkpoint = self.CheckpointType(self.model_meta, self.val_model_fname)
            valid_checkpoint.on_epoch_end(0, model_dict['model'])

        def _get_model(self):
            model = self.MockModel(self.input_dim, self.num_hidden, self.num_classes)
            loss = torch.nn.CrossEntropyLoss()
            optimizer = torch.optim.RMSprop(model.parameters())

            return {
                'model': model,
                'loss_fn': loss,
                'optimizer': optimizer,
                'scaler': None,
                'clip_grad': None,
            }

        def _get_data(self):
            return get_test_data(
                num_train=self.num_train,
                num_test=self.num_test,
                input_shape=(self.input_dim,),
                classification=True,
                num_classes=self.num_classes)

        def test_000_checkpoint_saves_meta(self):
            with self.ModelType(self.val_model_load_fname) as model_store:
               model_func = model_store.get_meta('model_function')
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
                self.CheckpointType(meta, model_fname)

        def test_010_no_scheduler(self):
            model_dict = self._get_model()
            scheduler_func = medaka.torch_ext.no_schedule()
            lr_scheduler = scheduler_func(
                model_dict["optimizer"], self.dataloader_train, 1, 0)
            expected_lrs = [
                0.01, # default lr for RMSProp
            ] * (self.num_train // self.batch_size)
            self.assertIsInstance(lr_scheduler, torch.optim.lr_scheduler.LambdaLR)
            for n, _ in enumerate(iter(self.dataloader_train)):
                self.assertListEqual(lr_scheduler.get_last_lr(), [expected_lrs[n],])
                lr_scheduler.step()

        def test_011_warmup_scheduler(self):
            model_dict = self._get_model()
            warmup_steps = 3
            scheduler_func = medaka.torch_ext.no_schedule(warmup_steps=warmup_steps)
            lr_scheduler = scheduler_func(
                model_dict["optimizer"], self.dataloader_train, 1, 0)
            expected_lrs = [
                0.001,
                0.004,
                0.007,
            ] + [0.01, ] * ((self.num_train // self.batch_size) - warmup_steps)
            for n, _ in enumerate(iter(self.dataloader_train)):
                self.assertAlmostEqual(lr_scheduler.get_last_lr()[0], expected_lrs[n])
                lr_scheduler.step()

        def test_012_cosine_scheduler(self):
            model_dict = self._get_model()
            scheduler_func = medaka.torch_ext.linear_warmup_cosine_decay(
                end_ratio=0.01,
                warmup_steps=0)
            lr_scheduler = scheduler_func(
                model_dict["optimizer"], self.dataloader_train, 1, 0)
            expected_lrs = 0.0001 + 0.5 * 0.0099 * (
                1 + np.cos(np.arange(0, 10) * np.pi / 10.))
            for n, _ in enumerate(iter(self.dataloader_train)):
                self.assertAlmostEqual(lr_scheduler.get_last_lr()[0], expected_lrs[n])
                lr_scheduler.step()

        def test_012_warmup_cosine_scheduler(self):
            model_dict = self._get_model()
            scheduler_func = medaka.torch_ext.linear_warmup_cosine_decay(
                end_ratio=0.01,
                warmup_steps=3)
            lr_scheduler = scheduler_func(
                model_dict["optimizer"], self.dataloader_train, 1, 0)
            expected_lrs = [0.001, 0.004, 0.007,] + \
                list(0.0001 + 0.5 * 0.0099 * (1 + np.cos(np.arange(0, 7) * np.pi / 7.)))
            for n, _ in enumerate(iter(self.dataloader_train)):
                self.assertAlmostEqual(lr_scheduler.get_last_lr()[0], expected_lrs[n])
                lr_scheduler.step()


class TestCallbacksTorch(TestNamespace.ModelAndData, unittest.TestCase):

    CheckpointType = medaka.torch_ext.ModelMetaCheckpoint
    ModelType = medaka.datastore.ModelStoreTGZ
    val_model_fname = tempfile.mkdtemp()
    val_model_load_fname = val_model_fname + '/model-0.tar.gz'
    train_model_fname = tempfile.mkdtemp()
    train_model_load_fname = train_model_fname + '/model-0.tar.gz'


training_features = os.path.join(os.path.dirname(__file__), 'data', 'training_features.hdf5')


class TestSequence(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.batch_size = 5

        self.tb = medaka.training.TrainBatcher(
            [training_features], batch_size=self.batch_size)

        self.n_training_samples = len(self.tb.train_samples)
        self.n_valid_samples = len(self.tb.valid_samples)

        self.seq_train = medaka.torch_ext.Sequence(
            self.tb.train_samples, self.tb.sample_to_x_y, "train", mini_epochs=1)
        self.seq_valid = medaka.torch_ext.Sequence(
            self.tb.valid_samples, self.tb.sample_to_x_y, "validation")

        self.mini_epochs = 3
        self.seq_mini = medaka.torch_ext.Sequence(
            self.tb.train_samples, self.tb.sample_to_x_y, "train", mini_epochs=self.mini_epochs)

    def test_000_correct_num_train_samples(self):
        self.assertEqual(len(self.seq_train), self.n_training_samples)
        self.assertEqual(self.seq_train.samples_per_epoch, self.n_training_samples)

    def test_001_correct_num_validation_batches(self):
        self.assertEqual(len(self.seq_valid), self.n_valid_samples)

    def test_002_shuffle_occurs_on_epoch_end(self):
        pre_shuffle = copy.deepcopy(self.seq_train.data)
        self.seq_train.on_epoch_end()
        post_shuffle = self.seq_train.data
        self.assertNotEqual(pre_shuffle, post_shuffle)
        self.assertEqual(set(pre_shuffle), set(post_shuffle))

    def test_003_correct_num_samples_with_miniepochs(self):
        self.assertEqual(len(self.seq_mini),
            self.n_training_samples // self.mini_epochs)
        self.assertEqual(self.seq_mini.samples_per_epoch,
            self.n_training_samples // self.mini_epochs)

    def test_004_shuffle_only_on_epoch_end(self):
        pre_shuffle = copy.deepcopy(self.seq_mini.data)
        for _ in range(self.mini_epochs - 1):
            self.seq_mini.on_epoch_end()
            post_shuffle = self.seq_train.data
            self.assertEqual(pre_shuffle, post_shuffle)
        self.seq_mini.on_epoch_end()
        post_shuffle = self.seq_train.data
        self.assertNotEqual(pre_shuffle, post_shuffle)

    def test_005_epoch_start_coord_set(self):
        expected_starts = [
                *[i * self.n_training_samples // self.mini_epochs for i in range(0, self.mini_epochs)],
                0
        ]
        starts = []
        for _ in range(self.mini_epochs + 1):
            self.seq_mini.on_epoch_start()
            starts.append(self.seq_mini.epoch_start)
            self.seq_mini.on_epoch_end()
        self.assertEqual(expected_starts, starts)

    def test_006_validation_init_with_miniepochs_gt_1_fails(self):
        with self.assertRaises(ValueError) as context:
            seq = medaka.torch_ext.Sequence(
                self.tb.valid_samples, self.tb.sample_to_x_y, dataset='validation', mini_epochs=3)

    def test_007_seq_init_with_incorrect_dataset_name_fails(self):
        with self.assertRaises(ValueError) as context:
            seq = medaka.torch_ext.Sequence(
                self.tb.train_samples, self.tb.sample_to_x_y, dataset='sasquatch') # only 'train' or 'validation' permitted

    def test_008_seq_with_random_seed(self):
        seq = medaka.torch_ext.Sequence(
            self.tb.train_samples, self.tb.sample_to_x_y, seed=2)


class TestSequenceBatcher(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.batch_size = 5

        self.tb = medaka.training.TrainBatcher(
            [training_features], batch_size=self.batch_size)

        self.n_training_samples = len(self.tb.train_samples)
        self.n_valid_samples = len(self.tb.valid_samples)

        self.sb_train = self.tb.train_loader()
        self.sb_valid = self.tb.valid_loader()

        self.mini_epochs = 3
        self.sb_mini = self.tb.train_loader(mini_epochs=self.mini_epochs)
 
    def test_000_correct_batch_size(self):
        first_batch = next(iter(self.sb_train))
        self.assertEqual(len(first_batch[0]), self.batch_size) # features
        self.assertEqual(len(first_batch[1]), self.batch_size) # labels

    def test_001_correct_num_train_batches(self):
        self.assertEqual(len(self.sb_train),
            self.n_training_samples // self.batch_size)

    def test_002_correct_num_validation_batches(self):
        self.assertEqual(len(self.sb_valid),
            self.n_valid_samples // self.batch_size)

    def test_003_correct_num_samples_with_miniepochs(self):
        self.assertEqual(len(self.sb_mini),
            self.n_training_samples // self.batch_size // self.mini_epochs)


class TestMetrics(unittest.TestCase):

    @classmethod
    def setUpClass(self) -> None:
        self.y_true = torch.Tensor([[0, 1, 1, 0, 0]])
        self.y_pred = torch.Tensor([
            [0.4, 0.6], # incorrect
            [0.1, 0.9], # correct
            [0.2, 0.8], # correct
            [0.7, 0.3], # correct
            [0.1, 0.9], # incorrect
        ])

    def test_000_sparse_cat_acc(self):
        # Test sparse categorical accuracy using class labels
        expected = [0., 1., 1., 1., 0.]
        metric = medaka.torch_ext.SparseCategoricalAccuracy()
        self.assertSequenceEqual(
            list(metric(self.y_true, self.y_pred).squeeze()),
            expected
        )

    def test_001_cat_acc(self):
        # Test categorical accuracy using one-hot classes
        expected = [0., 1., 1., 1., 0.]
        y_true = torch.nn.functional.one_hot(torch.as_tensor(self.y_true, dtype=torch.long))
        metric = medaka.torch_ext.CategoricalAccuracy()
        self.assertSequenceEqual(
            list(metric(y_true, self.y_pred).squeeze()),
            expected
        )

    def test_010_qscore(self):
        expected = -10.0 * np.log10(0.4)
        metric = medaka.torch_ext.QScore(accuracy=medaka.torch_ext.SparseCategoricalAccuracy)
        self.assertAlmostEqual(
            metric(self.y_true, self.y_pred).item(),
            expected,
            places=6
        )

    def test_011_qscore_clipped(self):
        eps = 1e-6
        expected = 60.0
        metric = medaka.torch_ext.QScore(eps=eps, accuracy=medaka.torch_ext.SparseCategoricalAccuracy)
        self.assertAlmostEqual(
            metric(torch.Tensor([1, 1, 1, 0, 1]), self.y_pred).item(),
            expected,
            places=6
        )
