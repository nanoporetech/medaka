import os
import tempfile
import unittest
import zipfile

import medaka.export
import medaka.features
import medaka.labels
import medaka.models


class TestDictNesting(unittest.TestCase):
    def test_001_expand_keys(self):
        flat_dict = {"a.b.C": 0, "a.b.D": ["1a", "1b"], "B": 2}
        expected = {'a': {'b': {'C': 0, 'D': ['1a', '1b']}}, 'B': 2}
        nested_dict = medaka.export.parse_params_to_nested_dict(flat_dict)
        self.assertDictEqual(nested_dict, expected)


class TestSpecLevelValidation(unittest.TestCase):
    def test_001_passes(self):
        config = {
            "config_version": 4,
            "supported_basecallers": ["hac", "sup"],
            "chunk_size": 100,
            "chunk_overlap": 10,
            "candidate_filtering": False,
            "model": {},
            "feature_encoder": {},
            "label_scheme": {},
        }
        medaka.export.validate_config_subsection(config, medaka.export.top_level_spec, 4)

    def test_002_fail_missing_required_key(self):
        # should pass for version <4 but fail for version 4
        config = {
            "config_version": 4,
            "supported_basecallers": ["hac", "sup"],
            # "chunk_size": 100,
            "chunk_overlap": 10,
            "candidate_filtering": False,
            "model": {},
            "feature_encoder": {},
            "label_scheme": {},
        }
        with self.assertRaises(medaka.export.ModelValidationError):
            medaka.export.validate_config_subsection(config, medaka.export.top_level_spec, 4)

        config["config_version"] = 3
        medaka.export.validate_config_subsection(config, medaka.export.top_level_spec, 3)
   
    def test_003_fail_extra_key(self):
        config = {
            "config_version": 4,
            "supported_basecallers": ["hac", "sup"],
            "chunk_size": 100,
            "chunk_overlap": 10,
            "candidate_filtering": False,
            "name": "spam",
            "model": {},
            "feature_encoder": {},
            "label_scheme": {},
        }
        with self.assertRaises(medaka.export.ModelValidationError):
            medaka.export.validate_config_subsection(config, medaka.export.top_level_spec, 4)
    
    def test_004_fail_type_mismatch(self):
        config = {
            "config_version": 4,
            "supported_basecallers": ["hac", "sup"],
            "chunk_size": "100",
            "chunk_overlap": 10,
            "candidate_filtering": False,
            "model": {},
            "feature_encoder": {},
            "label_scheme": {},
        }
        with self.assertRaises(medaka.export.ModelValidationError):
            medaka.export.validate_config_subsection(config, medaka.export.top_level_spec, 4)


class TestModelExport(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.root_dir = os.path.abspath(os.path.dirname(__file__))
        self.toml_config = os.path.join(self.root_dir, 'data', 'test_export_config.toml')
        self.dwells_model = medaka.models.resolve_model('r1041_e82_400bps_hac_v6.0.0_rl_lstm384_dwells')
        self.counts_model = medaka.models.resolve_model('r1041_e82_400bps_hac_v6.0.0')
        self.label_scheme = medaka.labels.HaploidLabelScheme
        self.feature_encoder = medaka.features.CountsFeatureEncoder

    def test_01_export_model_v3(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            export_name = os.path.join(tmpdir, 'model_export')

            class DummyArgs:
                model = self.counts_model
                output = export_name
                force = True
                script = True
                compress = True
                params = {
                    "supported_basecallers": [
                        "dna_r10.4.1_e8.2_400bps_hac@v5.0.0",
                        "dna_r10.4.1_e8.2_400bps_sup@v5.0.0"]
                }
                config_version = 3

            medaka.export.export_model(DummyArgs())

            archive_name = export_name + '.zip'
            assert os.path.exists(archive_name)
            assert zipfile.is_zipfile(archive_name)
            # check that there is a model file in the export directory

            with zipfile.ZipFile(archive_name, 'r') as zip:
                assert zipfile.Path(zip, 'model_export/weights.pt').is_file()
                assert zipfile.Path(zip, 'model_export/config.toml').is_file()
                assert zipfile.Path(zip, 'model_export/model.pt').is_file()

    def test_02_export_model_v4(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            export_name = os.path.join(tmpdir, 'model_export')

            class DummyArgs:
                model = self.counts_model
                output = export_name
                force = True
                script = False
                compress = True
                params = {
                    "supported_basecallers": ["dna_r10.4.1_e8.2_400bps_hac@v5.0.0"],
                    "chunk_size": 100,
                    "chunk_overlap": 10,
                    "candidate_filtering": False
                }
                config_version = 4

            medaka.export.export_model(DummyArgs())

            archive_name = export_name + '.zip'
            assert os.path.exists(archive_name)
            # check that there is a model file in the export directory

            with zipfile.ZipFile(archive_name, 'r') as zip:
                assert zipfile.Path(zip, 'model_export/weights.pt').is_file()
                assert zipfile.Path(zip, 'model_export/config.toml').is_file()
                assert not zipfile.Path(zip, 'model_export/model.pt').exists()  # script=False

    def test_03_export_model_v4_invalid_args(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            export_name = os.path.join(tmpdir, 'model_export')

            class DummyArgs:
                model = self.dwells_model
                output = export_name
                force = True
                script = False
                compress = True
                params = {
                    "supported_basecallers": ["dna_r10.4.1_e8.2_400bps_hac@v5.0.0"],
                    "chunk_size": 100,
                    # "chunk_overlap": 10,
                    "candidate_filtering": False
                }
                config_version = 4

            with self.assertRaises(medaka.export.ModelValidationError):
                medaka.export.export_model(DummyArgs())

    def test_04_export_model_v3_uncompressed(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            export_name = os.path.join(tmpdir, 'model_export')

            class DummyArgs:
                model = self.dwells_model
                output = export_name
                force = True
                script = False
                compress = False
                params = {"supported_basecallers": ["dna_r10.4.1_e8.2_400bps_hac@v5.0.0"]}
                config_version = 3

            medaka.export.export_model(DummyArgs())

            self.assertTrue(os.path.isdir(export_name))
            self.assertTrue(os.path.isfile(os.path.join(export_name, "weights.pt")))
            self.assertTrue(os.path.isfile(os.path.join(export_name, "config.toml")))
            self.assertFalse(os.path.exists(os.path.join(export_name, "model.pt")))

