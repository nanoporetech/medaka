import copy
import functools
import os
import tempfile
import unittest

import numpy as np
import tarfile
import torch

import medaka.common
import medaka.datastore
import medaka.training
import medaka.torch_ext
import medaka.models

from medaka.test.test_sample import get_test_samples

class TestNamespace:

    class ModelAndData:
        CheckpointType = medaka.torch_ext.ModelMetaCheckpoint
        ModelType = medaka.datastore.ModelStoreTGZ

        val_model_fname = tempfile.mkdtemp()
        val_model_load_fname = os.path.join(val_model_fname,'model-0.tar.gz')
        train_model_fname = tempfile.mkdtemp()
        train_model_load_fname = os.path.join(train_model_fname, 'model-0.tar.gz')
        input_dim = 2
        num_hidden = 4
        num_classes = 2
        batch_size = 2
        num_train = 20
        num_test = 20


        @classmethod
        def setUpClass(self):
            data_train, data_test = self._get_counts_matrix_data(self)

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
            meta_data = {
                'model_function': model_dict.pop('model_function'),
                'feature_encoder': self.FeatureEncoder,
                'label_scheme': self.LabelScheme,
            }
            
            _ = medaka.torch_ext.run_epoch(
                dataloader=self.dataloader_train,
                is_training_epoch=False,
                **model_dict
            )
            train_checkpoint = self.CheckpointType(meta_data, self.train_model_fname)
            train_checkpoint.on_epoch_end(0, model_dict['model'])

            _ = medaka.torch_ext.run_epoch(
                dataloader=self.dataloader_test,
                is_training_epoch=False,
                **model_dict
            )
            valid_checkpoint = self.CheckpointType(meta_data, self.val_model_fname)
            valid_checkpoint.on_epoch_end(0, model_dict['model'])

        def _get_model(self):
            model_function = functools.partial(medaka.models.build_default_architecture, (1000, 10))
            model = model_function()
            loss = torch.nn.CrossEntropyLoss()
            optimizer = torch.optim.RMSprop(model.parameters())

            return {
                'model_function': model_function,
                'model': model,
                'loss_fn': loss,
                'optimizer': optimizer,
                'scaler': None,
                'clip_grad': None,
            }

        def _get_counts_matrix_data(self):
            return get_test_samples(
                num_train=self.num_train,
                num_test=self.num_test,
                read_level_features=False)
        
        def _get_read_level_features_data(self):
            return get_test_samples(
                num_train=self.num_train,
                num_test=self.num_test,
                read_level_features=True)

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
            start_lr = 0.01
            end_ratio = 0.01
            scheduler_func = medaka.torch_ext.linear_warmup_cosine_decay(
                end_ratio=end_ratio,
                warmup_steps=0)
            lr_scheduler = scheduler_func(
                model_dict["optimizer"], self.dataloader_train, 1, 0)
            expected_lrs = end_ratio*start_lr + 0.5 * (1-end_ratio)* start_lr * (
                1 + np.cos(np.linspace(0,1,len(self.dataloader_train)+1) * np.pi))
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

        def test_013_clip_grad(self):
            clip_grad = medaka.torch_ext.ClipGrad(quantile=0.5, factor=2)
            buffer_vals = torch.ones(len(clip_grad.buffer))

            for i in range(len(clip_grad.buffer)):
                clip_grad.append(buffer_vals[i])
            
            random_parameters = torch.nn.Parameter(1000*torch.ones(10))
            y = random_parameters[None] @ torch.ones(10)
            y.backward()
            clip_grad(random_parameters)
            clipped_grad_norm = random_parameters.grad.norm()
            
            # check gradient of clipped parameters is good
            self.assertTrue(clipped_grad_norm <= 2*torch.quantile(buffer_vals, 0.5))

class TestCallbacksTorch(TestNamespace.ModelAndData, unittest.TestCase):

    CheckpointType = medaka.torch_ext.ModelMetaCheckpoint
    LabelScheme = medaka.labels.HaploidLabelScheme
    FeatureEncoder = medaka.features.CountsFeatureEncoder
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
            self.tb.train_samples, self.tb.load_sample_worker, "train")
        self.seq_valid = medaka.torch_ext.Sequence(
            self.tb.valid_samples, self.tb.load_sample_worker, "validation")

    def test_000_correct_num_train_samples(self):
        self.assertEqual(len(self.seq_train), self.n_training_samples)

    def test_001_correct_num_validation_batches(self):
        self.assertEqual(len(self.seq_valid), self.n_valid_samples)

    def test_007_seq_init_with_incorrect_dataset_name_fails(self):
        with self.assertRaises(ValueError) as context:
            seq = medaka.torch_ext.Sequence(
                self.tb.train_samples, self.tb.load_sample_worker, dataset='sasquatch') # only 'train' or 'validation' permitted

    def test_008_seq_with_random_seed(self):
        seq = medaka.torch_ext.Sequence(
            self.tb.train_samples, self.tb.load_sample_worker, seed=2)


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
 
    def test_000_correct_batch_size(self):
        first_batch = next(iter(self.sb_train))
        self.assertEqual(len(first_batch.counts_matrix), self.batch_size) # features
        self.assertEqual(len(first_batch.labels), self.batch_size) # labels

    def test_001_correct_num_train_batches(self):
        self.assertEqual(len(self.sb_train),
            self.n_training_samples // self.batch_size)

    def test_002_correct_num_validation_batches(self):
        self.assertEqual(len(self.sb_valid),
            self.n_valid_samples // self.batch_size)

class TestModelExport(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.root_dir = os.path.abspath(os.path.dirname(__file__))
        self.toml_config = os.path.join(self.root_dir, 'data', 'test_export_config.toml')
        self.model_dict = {
            'type': "GRUModel",
            'kwargs': {
                'num_features': 10,
                'num_classes': 5,
                'gru_size': 128}}
        self.label_scheme = medaka.labels.HaploidLabelScheme
        self.feature_encoder = medaka.features.CountsFeatureEncoder

    def test_export_model(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            # build model
            meta_data = {
                'model_function': functools.partial(medaka.models.model_from_dict, self.model_dict),
                'feature_encoder': self.feature_encoder(),
                'label_scheme': self.label_scheme(),
            }
            model = meta_data['model_function']()
            # save model to tarball
            checkpoint = medaka.torch_ext.ModelMetaCheckpoint(
                meta_data, os.path.join(tmpdir, 'model'))
            checkpoint.on_epoch_end(0, model)

            model_saved_fp = os.path.join(tmpdir, 'model', 'model-0.tar.gz')
            export_name = model_saved_fp.replace('-0.tar.gz','_export')

            class DummyArgs:
                model = model_saved_fp
                output = export_name 
                force = True
                script = True
                supported_basecallers = ["dna_r10.4.1_e8.2_400bps_hac@v5.0.0",]

            medaka.torch_ext.export_model(DummyArgs())

            assert os.path.exists(export_name + '.tar.gz')
            # check that there is a model file in the export directory

            with tarfile.open(export_name + '.tar.gz', 'r:gz') as tar:
                files=tar.getnames()
                assert 'model/weights.pt' in files
                assert 'model/config.toml' in files
                assert 'model/model.pt' in files


