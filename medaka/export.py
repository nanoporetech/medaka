"""Export models with a config.toml that can be loaded by dorado."""
import collections
import os
import shutil
import zipfile

import toml
import torch

import medaka.common
import medaka.models


# Specifications of dorado model configuration files
# see dorado/dorado/secondary/architectures/model_config_spec.md

class ModelValidationError(Exception):
    """Exception class for model config validation."""

    pass


ConfigEntry = collections.namedtuple(
    'ConfigEntry', ['name', 'min_version', 'type'])


LATEST_CONFIG = 4


top_level_spec = [
    ConfigEntry(name="config_version", min_version=1, type=int),
    ConfigEntry(name="supported_basecallers", min_version=2, type=list),
    ConfigEntry(name="chunk_size", min_version=4, type=int),
    ConfigEntry(name="chunk_overlap", min_version=4, type=int),
    ConfigEntry(name="candidate_filtering", min_version=4, type=bool),
    ConfigEntry(name="model", min_version=1, type=dict),
    ConfigEntry(name="feature_encoder", min_version=1, type=dict),
    ConfigEntry(name="label_scheme", min_version=1, type=dict),
]


subsection_specs = {
    "model": [
        ConfigEntry(name="type", min_version=1, type=str),
        ConfigEntry(name="kwargs", min_version=1, type=dict),
    ],
    "feature_encoder": [
        ConfigEntry(name="type", min_version=1, type=str),
        ConfigEntry(name="kwargs", min_version=1, type=dict),
    ],
    "label_scheme": [
        ConfigEntry(name="type", min_version=1, type=str),
        ConfigEntry(name="kwargs", min_version=None, type=dict),
    ]
}


feature_encoder_specs = {
    "CountsFeatureEncoder": [
        ConfigEntry(name="normalise", min_version=1, type=str),
        ConfigEntry(name="tag_keep_missing", min_version=1, type=bool),
        ConfigEntry(name="min_mapq", min_version=1, type=int),
        ConfigEntry(name="sym_indels", min_version=1, type=bool),
        ConfigEntry(name="dtypes", min_version=1,
                    type=collections.abc.Collection)
    ],
    "ReadAlignmentFeatureEncoder": [
        ConfigEntry(name="dtypes", min_version=1,
                    type=collections.abc.Collection),
        ConfigEntry(name="tag_keep_missing", min_version=1, type=bool),
        ConfigEntry(name="min_mapq", min_version=1, type=int),
        ConfigEntry(name="max_reads", min_version=1, type=int),
        ConfigEntry(name="row_per_read", min_version=1, type=bool),
        ConfigEntry(name="include_dwells", min_version=1, type=bool),
        ConfigEntry(name="include_haplotype", min_version=1, type=bool),
        ConfigEntry(name="include_snp_qv", min_version=4, type=bool),
        ConfigEntry(name="right_align_insertions", min_version=4, type=bool),
        ConfigEntry(name="region_split", min_version=None, type=int),
    ]
}


label_scheme_specs = {
    # no keys are required
    "HaploidLabelScheme": [
        ConfigEntry(name="ploidy", min_version=None, type=int),
        ConfigEntry(name="ordered", min_version=None, type=bool),
        ConfigEntry(name="right_align_insertions", min_version=None,
                    type=bool),
    ],
    "DiploidLabelScheme": [
        ConfigEntry(name="ploidy", min_version=None, type=int),
        ConfigEntry(name="ordered", min_version=None, type=bool),
        ConfigEntry(name="right_align_insertions", min_version=None,
                    type=bool),
    ],
}


model_specs = {
    "GRUModel": [
        ConfigEntry(name="num_features", min_version=1, type=int),
        ConfigEntry(name="num_classes", min_version=1, type=int),
        ConfigEntry(name="gru_size", min_version=1, type=int),
        ConfigEntry(name="n_layers", min_version=1, type=int),
        ConfigEntry(name="bidirectional", min_version=1, type=bool),
    ],
    "LatentSpaceLSTM": [
        ConfigEntry(name="num_classes", min_version=1, type=int),
        ConfigEntry(name="lstm_size", min_version=1, type=int),
        ConfigEntry(name="cnn_size", min_version=1, type=int),
        ConfigEntry(name="pooler_type", min_version=1, type=str),
        ConfigEntry(name="bases_alphabet_size", min_version=1, type=int),
        ConfigEntry(name="bases_embedding_size", min_version=1, type=int),
        ConfigEntry(name="kernel_sizes", min_version=1,
                    type=collections.abc.Collection),
        ConfigEntry(name="use_dwells", min_version=1, type=bool),
        ConfigEntry(name="bidirectional", min_version=1, type=bool),
        ConfigEntry(name="pooler_args", min_version=None, type=dict),
    ]
}


def validate_config(test_dict):
    """Validate a configuration.

    :param test_dict: nested dictionary for a model export.
    """
    if "config_version" not in test_dict.keys():
        raise ModelValidationError("Configuration is missing a version number")

    config_version = test_dict["config_version"]
    # validate top-level keys
    validate_config_subsection(test_dict, top_level_spec, config_version)

    # validate subsections
    for key, spec_dict in [
            ("model", model_specs), ("feature_encoder", feature_encoder_specs),
            ("label_scheme", label_scheme_specs)]:
        validate_config_subsection(
            test_dict[key], subsection_specs[key], config_version)
        if (
            "type" in test_dict and
            test_dict[key]["type"] not in spec_dict.keys()
        ):
            msg = "Unknown type {} for {}, allowed values are {}"
            raise ModelValidationError(
                msg.format(test_dict[key]["type"], key, spec_dict.keys()))
        if "kwargs" in test_dict[key]:
            validate_config_subsection(
                test_dict[key]["kwargs"], spec_dict[test_dict[key]["type"]],
                config_version)


def validate_config_subsection(test_dict, spec, config_version):
    """Validate a config dictionary according to a specification.

    This function performs the following checks:
    - all keys in test_dict are names of config entries in spec,
    - all required entries for the given config_version in spec are in the
    keys of test_dict,
    - the types of all values in test_dict match with the type of the
    corresponding config entry in spec.

    :param test_dict: dictionary to test
    :param spec: list of ConfigEntrys
    :param config_version: int, config version for spec to validate against.
    """
    spec_keys = set([c.name for c in spec])
    required_keys = set([c.name for c in spec if (
        c.min_version is not None and c.min_version <= config_version)])
    provided_keys = set(test_dict.keys())

    # check that all the provided keys are present in the spec
    if not provided_keys.issubset(spec_keys):
        msg = "Unknown keys present: {}"
        raise ModelValidationError(msg.format(provided_keys - spec_keys))
    # check that all required keys for the given config version are present
    if not required_keys.issubset(provided_keys):
        msg = "Required keys are missing: {}"
        raise ModelValidationError(msg.format(required_keys - provided_keys))
    # check that types match for all present keys
    for config_entry in spec:
        name = config_entry.name
        if name in provided_keys and not isinstance(
                test_dict[name], config_entry.type):
            msg = "Type mismatch for key '{}': provided '{}',  requires '{}'"
            raise ModelValidationError(
                msg.format(type(test_dict[name]), name, config_entry.type))


def parse_params_to_nested_dict(flat_dict):
    """Expand dotted keys in `flat_dict` into a nested dictionary."""
    nested_dict = dict()

    def _assign_leaf(target, key, val):
        if key in target:
            if isinstance(target[key], dict):
                raise ValueError
            else:
                raise IndexError
        else:
            target[key] = val

    for key, value in flat_dict.items():
        if not isinstance(key, str):
            raise ValueError("Key must be a string, received {}".format(key))

        parts = key.split(".")
        current = nested_dict

        for n, part in enumerate(parts[:-1]):
            if part not in current:
                current[part] = {}
            elif not isinstance(current[part], dict):
                msg = (
                    "Cannot create level {}, already has a non-dictionary "
                    "value.")
                raise ValueError(msg.format(".".join(parts[:n+1])))
            current = current[part]

        try:
            _assign_leaf(current, parts[-1], value)
        except IndexError:
            msg = "Cannot assign to existing key: {}"
            raise IndexError(msg.format(key))
        except ValueError:
            msg = "Cannot assign to key {}, already exists as a nested level."
            raise ValueError(msg.format(key))

    return nested_dict


def export_model(args):
    """Export a model to a torchscript model and config file."""
    logger = medaka.common.get_named_logger('ModelExport')
    model_fp = args.model
    output_fp = args.output
    force = args.force
    script = args.script

    if not os.path.exists(model_fp):
        raise FileNotFoundError(f"Model file not found: {model_fp}")

    params_dict = parse_params_to_nested_dict(args.params)

    model_store = medaka.models.open_model(model_fp)
    model = model_store.load_model()
    label_scheme = model_store.get_meta("label_scheme")
    feature_encoder = model_store.get_meta("feature_encoder")

    config = {
        'config_version': args.config_version,
        'model': model.to_dict(),
        'feature_encoder': feature_encoder.to_dict(),
        'label_scheme': label_scheme.to_dict()}
    config = config | params_dict

    validate_config(config)

    if output_fp is None:
        output_fp = os.path.basename(model_fp).replace(".tar.gz", "_export")

    logger.info(f"Exporting model from {model_fp} to {output_fp}.tar.gz")
    if force:
        os.makedirs(output_fp, exist_ok=True)
    else:
        os.mkdir(output_fp)

    msg = "Output file path cannot be the same as the model file path"
    assert model_fp != output_fp, msg

    with open(os.path.join(output_fp, 'config.toml'), 'w') as f:
        toml.dump(config, f)

    if script:
        for name, module in model.named_modules():
            if (
                isinstance(module, torch.nn.LSTM) or
                isinstance(module, torch.nn.GRU)
            ):
                module.flatten_parameters()
        scripted_model = torch.jit.script(model)
        scripted_model.save(os.path.join(output_fp, 'model.pt'))

    # save state dict
    w = {k: v.cpu() for k, v in model.state_dict().items()}
    torch.save(w, os.path.join(output_fp, 'weights.pt'))

    # now add to tarfile
    if args.compress:
        output_zipfile = output_fp + ".zip"
        with zipfile.ZipFile(output_zipfile, "w") as zip:
            for file in os.listdir(output_fp):
                zip.write(
                    os.path.join(output_fp, file),
                    arcname=os.path.join(os.path.basename(output_fp), file))
        shutil.rmtree(output_fp)
