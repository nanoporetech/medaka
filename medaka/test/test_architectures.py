import os
import functools
import unittest

import toml
import torch

import medaka.architectures
import medaka.architectures.base_classes
import medaka.features
import medaka.models
from medaka.test.test_sample import get_test_samples
import medaka.torch_ext


__model_config__ = os.path.join(os.path.dirname(__file__), 'data', 'test_model_spec.toml')

class TestArchitectures(unittest.TestCase):
    batch_size = 10
    @classmethod
    def setUpClass(cls):
        counts_matrix_samples, _ = get_test_samples(num_train=cls.batch_size, num_test=1, read_level_features=False)
        read_level_samples, _ = get_test_samples(num_train=cls.batch_size, num_test=1, read_level_features=True)

        cls.counts_matrix_dataloader = medaka.torch_ext.SequenceBatcher(
            medaka.torch_ext.Sequence(
                counts_matrix_samples,
                sampler_func=lambda x: x,
                dataset="train",
            ),
            batch_size=cls.batch_size,
            collate_fn=functools.partial(medaka.torch_ext.Batch.collate,
                                         counts_matrix=True),
        )
        cls.read_level_dataloader = medaka.torch_ext.SequenceBatcher(
            medaka.torch_ext.Sequence(
                read_level_samples,
                sampler_func=lambda x: x,
                dataset="train",
            ),
            batch_size=cls.batch_size,
            collate_fn=functools.partial(
                medaka.torch_ext.Batch.collate,
                counts_matrix=True),
        )

        class DummyInvalidFeatureEncoder:
            pass


        cls.feature_encoders = {
            'counts_matrix': medaka.features.CountsFeatureEncoder(),
            'rl': medaka.features.ReadAlignmentFeatureEncoder(include_dwells=False),
            'rl_dwells': medaka.features.ReadAlignmentFeatureEncoder(include_dwells=True),
            'other': DummyInvalidFeatureEncoder(),
        }

    def _run_iteration(self, dataloader, model):
        """Check that the model can run both inference and training passes on batch"""
        batch = next(iter(dataloader))
        loss_fn = torch.nn.CrossEntropyLoss()
        model.process_batch(batch, loss_fn)
        probs = model.predict_on_batch(batch)
        assert probs.shape == (batch.labels.shape[0], batch.labels.shape[1], model.num_classes)

    def _test_model_passes(self, model):
        """Test that the model can run on both read level and counts matrix data (if applicable)"""
        # run iteration on counts matrix data
        if isinstance(model, medaka.architectures.base_classes.CountsMatrixModel):
            self._run_iteration(self.counts_matrix_dataloader, model)
            # check that the counts matrix model can also be run on read level data
            self._run_iteration(self.read_level_dataloader, model)
        elif isinstance(model, medaka.architectures.base_classes.ReadLevelFeaturesModel):
            self._run_iteration(self.read_level_dataloader, model)

    def _test_fenc_compat(self, model, valid_fencs, invalid_fencs):
        for vfenc in valid_fencs:
            model.check_feature_encoder_compatibility(self.feature_encoders[vfenc])
        for ivfenc in invalid_fencs:
            with self.assertRaises(ValueError):
                model.check_feature_encoder_compatibility(self.feature_encoders[ivfenc])

    def test_model_from_toml(self):
        model_dict = toml.load(__model_config__)
        model = medaka.models.model_from_dict(model_dict)
        self._test_model_passes(model)

    def test_majority_model(self):
        model = medaka.architectures.MajorityVoteModel()
        self._test_fenc_compat(model, ['counts_matrix','rl'], ['other'])
        self._test_model_passes(model)

    def test_gru_model(self):
        model = medaka.architectures.GRUModel()
        self._test_fenc_compat(model, ['counts_matrix','rl'], ['other'])
        self._test_model_passes(model)

    def test_latent_gru_model(self):
        model = medaka.architectures.LatentSpaceLSTM()
        self._test_fenc_compat(model, ['rl'], ['counts_matrix','other'])
        self._test_model_passes(model)
    
    def test_latent_gru_model_w_dwells(self):
        model = medaka.architectures.LatentSpaceLSTM(use_dwells=True)
        self._test_fenc_compat(model, ['rl_dwells'], ['counts_matrix','other','rl'])
        self._test_model_passes(model)
