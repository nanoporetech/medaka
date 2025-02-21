"""Creation and loading of models."""

import abc
import importlib
import inspect
import itertools
import os
import pathlib
import tempfile

import pysam
import requests
import torch

import medaka.common
import medaka.datastore
import medaka.options


logger = medaka.common.get_named_logger('ModelLoad')


class DownloadError(ValueError):
    """Raised when model is unsuccessfully downloaded."""


model_suffixes = ["_model_pt.tar.gz", ]

DEFAULT_MODEL_DICT = {
    "type": "GRUModel",
    "kwargs": {
        "num_features": 10,
        "num_classes": 5,
        "gru_size": 256,
    },
}


def resolve_model(model):
    """Resolve a model filepath, downloading known models if necessary.

    :param model_name: str, model filepath or model ID

    :returns: str: filepath to hdf/tar.gz model file or model directory.
    """
    if os.path.exists(model):  # model is path to model file
        return model
    elif model in medaka.options.deprecated_models:
        raise medaka.options.DeprecationError(model)
    elif model not in medaka.options.known_models:
        # try to resolve as a basecaller model <model>:<variety>
        err_msg = f"Failed to interpret '{model}' as a basecaller model."
        try:
            bc_model, variety = model.split(":")
            if variety in {'consensus', 'variant'} and \
                    bc_model in medaka.options.basecaller_models.keys():
                consensus, var = medaka.options.basecaller_models[bc_model]
                mod = var if variety == 'variant' else consensus
                return resolve_model(mod)
            else:
                raise ValueError(err_msg)
        except medaka.options.DeprecationError as ex:
            raise ex
        except Exception:
            logger.warning(err_msg)

        extra_msg = ""
        if "fast" in model:
            extra_msg = (
                " 'fast' basecalling models are not supported by medaka, "
                "if you have used a 'fast' model to perform basecalling you "
                "will need to first perform basecalling with a "
                "'high accuracy' or 'super accuracy' model before using "
                "medaka.")
        raise ValueError(
            f"The model '{model}' is not a recognised basecaller model or "
            "existant file. This could indicate a malformed input file (for "
            "which medaka was unable to identify correctly the basecaller "
            "meta-information) or simply be that the model is not supported "
            f"by medaka.{extra_msg}")

    else:
        # check for model in model stores
        for suffix in model_suffixes:
            fname = '{}{}'.format(model, suffix)
            for ms in medaka.options.model_stores:
                fp = os.path.join(ms, fname)
                if os.path.exists(fp):
                    return fp

        # try to download model
        download_errors = 0
        data = None
        for suffix in model_suffixes:
            fname = '{}{}'.format(model, suffix)
            url = medaka.options.model_url_template.format(
                pkg=__package__, subdir=medaka.options.model_subdir,
                fname=fname)
            try:
                data = requests.get(url).content
                with tempfile.TemporaryDirectory() as tmpdir:
                    # write the data and check it looks like a model
                    tmp_file = os.path.join(tmpdir, "tmp{}".format(suffix))
                    with open(tmp_file, 'wb') as tmp_model:
                        tmp_model.write(data)
                    # this will call ourself recursively, but we give a
                    # filepath now so will return immediately.
                    with open_model(tmp_file) as model_store:
                        model_store.get_meta('model_function')
            except Exception:
                download_errors += 1
            else:
                break
        if download_errors == len(model_suffixes):
            raise DownloadError(
                "The model file for {} is not already installed and "
                "could not be downloaded. Check you are connected to "
                "the internet and try again.".format(model))
        else:
            # save the model, try all locations
            for ms in medaka.options.model_stores:
                fp = os.path.join(ms, fname)
                try:
                    d = os.path.dirname(fp)
                    pathlib.Path(d).mkdir(parents=True, exist_ok=True)
                    with open(fp, 'wb') as fh:
                        fh.write(data)
                    return fp
                except Exception:  # we might not have write access
                    pass
            msg = (
                "The model file for {} is not installed and could not be "
                "installed to any of {}. If you cannot gain write "
                "permissions, download the model file manually from {} and "
                "use the downloaded model as the --model option.")
            raise RuntimeError(
                msg.format(
                    model, ' or '.join(medaka.options.model_stores), url))
    raise RuntimeError("Model resolution failed")


def model_from_basecaller(fname, variant=False, bacteria=False):
    """Determine correct medaka model from basecaller output file.

    :param fname: a basecaller output (.sam/.bam/.cram/.fastq).
    :param variant: whether to return variant model (otherwise consensus).
    :param bacteria: whether to override basecaller consensus model with
        bacterial methylation model (for 5khz hac and sup models only).

    There are slight differences is the search strategy for .bam and .fastq
    files due to differences in what information is available in each.

    For .bam (and related) files the DS subfield of the read group header
    is examined to find the "basecall_model=" key=value entry. The found
    model is returned vebatim without further checks.

    For .fastq (and related) a RG:Z key=value in record comments is searched
    in the first 100 records. Due to ambiguities in the representation the
    search looks explicitely for known models.
    """
    logger = medaka.common.get_named_logger("MdlInspect")
    logger.info("Trying to find model")
    try:
        models = _model_from_bam(fname)
    except Exception:
        try:
            models = _model_from_fastq(fname)
        except Exception:
            raise IOError(
                "Failed to parse basecaller models from input file.")

    if len(models) != 1:
        # TODO: this potentially conflicts with medaka's ability to use
        #       multiple data types. In that case a user can just
        #       explicitely provide the correct model.
        raise ValueError(
            "Input file did not contain precisely 1 basecaller "
            "model reference.")

    basecaller = list(models)[0]
    if basecaller not in medaka.options.basecaller_models.keys():
        raise KeyError(
            "Unknown basecaller model. Please provide a medaka model "
            "explicitely using --model.")

    consensus, var = medaka.options.basecaller_models[basecaller]
    model = var if variant else consensus
    if model is None:
        txt = "variant" if variant else "consensus"
        raise ValueError(
            f"No {txt} model available for basecaller {basecaller}.")

    # replace with bacterial consensus model if needed
    if bacteria and not variant:
        if model in medaka.options.bact_methyl_compatible_models:
            model = medaka.options.bact_methyl_model
        else:
            logger.warning(
                "WARNING: --bacteria specified but input data was not "
                "compatible. Using default model {}.".format(model))
    return model


def _model_from_bam(fname):
    """Search for basecaller models listed in a .bam."""
    models = set()
    with pysam.AlignmentFile(fname, check_sq=False) as bam:
        callers = [rg['DS'] for rg in bam.header['RG']]
        logger.info(f"Found basecall models: {callers}")
        for caller in callers:
            models.add(caller.split("basecall_model=")[1].split()[0])
    return models


def _model_from_fastq(fname):
    """Search for model files listed in a .fastq."""
    known_models = list(medaka.options.basecaller_models.keys())
    models = set()
    with pysam.FastxFile(fname, 'r') as fastq:
        for rec in itertools.islice(fastq, 100):
            try:
                # dorado SAM converted to FASTQ with e.g. samtools fastq
                # model is embedded in RG:Z: tag of comment as
                # <run_id>_<model>_<barcode>, but model has _
                # characters in also so search for known models
                read_group = rec.comment.split("RG:Z:")[1].split()[0]
                for model in known_models:
                    if model in read_group:
                        models.add(model)
            except Exception:
                # minknow/guppy
                # basecall_model_version_id=<model>
                try:
                    model = rec.comment.split(
                        "basecall_model_version_id=")[1].split()[0]
                    models.add(model)
                except Exception:
                    try:
                        # rogue non-spec conforming minknow versions
                        model = rec.comment.split(
                            "model_version_id=")[1].split()[0]
                        models.add(model)
                    except Exception:
                        pass
    if len(models) > 1:
        # filter out any models without an `@`. These are likely FPs of
        # the search above (there are unversioned models whose name
        # is a substring of the versioned models).
        unversioned = {m for m in models if '@' not in m}
        versioned = {m for m in models if '@' in m}
        remove = set()
        for unver, ver in itertools.product(unversioned, versioned):
            if unver in ver:
                remove.add(unver)
        models = models - remove
    return models


def open_model(fname):
    """Determine model type from model name.

    :param fname: model filepath

    : returns: model store object
    """
    fname = resolve_model(fname)
    ext = os.path.splitext(fname)[-1].lower()
    if ext == ".gz":
        return medaka.datastore.ModelStoreTGZ(fname)
    else:
        raise ValueError(
            f"Model {fname} has an unknown extension."
            "Only .gz model extensions are supported in medaka>=2.0.0"
        )


class TorchModel(torch.nn.Module):
    """Base class for pytorch models."""

    def __init__(self):
        """Initialise underlying Module."""
        super().__init__()
        self.half_precision = False
        self.logger = medaka.common.get_named_logger('TorchModel')

    @abc.abstractmethod
    def forward(self, *args, **kwargs):
        """Model forward pass, to be overwritted by the specific model."""
        raise NotImplementedError

    def device(self):
        """Device where model has been loaded."""
        try:
            return next(self.parameters()).device
        except StopIteration:
            return torch.device("cpu")

    def half(self):
        """Set model to half precision."""
        super().half()
        self.half_precision = True

    def predict_on_batch(self, batch):
        """Run inference on a feature batch.

        Handles moving features to the correct device and type conversion
        if required. Returns a cpu tensor.
        """
        x = self.get_model_input_features(batch).to(self.device())
        with torch.inference_mode():
            with torch.amp.autocast("cuda", enabled=self.half_precision):
                x = self.forward(x).detach().cpu()
        return x

    def process_batch(self, batch, loss_fn,):
        """Process the training batch, returning loss and metrics.

        :param batch (medaka.torch_ext.Batch): A batch of samples
        :param loss_fn (torch.nn.Module): callable loss function to use,
            with signature loss_fn(logits, labels) -> loss
        """
        labels = batch.labels.to(self.device())
        feature_matrix = self.get_model_input_features(batch).to(self.device())

        with torch.amp.autocast("cuda", enabled=self.half_precision):
            logits = self.forward(feature_matrix)
            labels = labels.to(logits.device)

            # calculate loss
            loss = loss_fn(logits.flatten(0, 1), labels.flatten())

        model_preds = logits.detach().argmax(dim=-1)

        # now calculate metrics of interest
        metrics = {}
        metrics["n_model_correct"] = (model_preds == labels).sum().item()

        if batch.majority_vote_probs is not None:
            maj_vote_probs = batch.majority_vote_probs.to(self.device())
            maj_vote_preds = maj_vote_probs.argmax(dim=-1)
            maj_vote_correct = (maj_vote_preds == labels)
            metrics["n_argmax_correct"] = maj_vote_correct.sum().item()
        metrics["n_positions"] = labels.numel()

        return loss, metrics

    def to_dict(self):
        """Return a dict of the model name and args.

        Used for serialisation to recreate the model.
        """
        kwargs = inspect.signature(self.__class__.__init__).parameters
        out_kwargs = {}
        for k, v in kwargs.items():
            if k == "self":
                continue
            elif hasattr(self, k):
                out_kwargs[k] = getattr(self, k)
            elif v.default != inspect.Parameter.empty:
                out_kwargs[k] = v.default
            else:
                raise ValueError(
                    f"Model parameter {k} not set, Cannot serialise model."
                )
        return {"type": self.__class__.__name__, "kwargs": out_kwargs}

    @abc.abstractmethod
    def get_model_input_features(self, batch):
        """Return the input features for the model."""
        pass

    def count_parameters(self):
        """Return number of trainable parameters."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)

    def check_feature_encoder_compatibility(self, fenc):
        """Check feature encoder is valid for this model."""
        clsname = type(self).__name__
        msg = (
            "Feature encoder compatibility check not currently implemented for"
            f"{clsname}, be careful!"
        )
        self.logger.warn(msg)


def build_default_architecture(feature_shape):
    """Build the default model architecture."""
    model = model_from_dict(DEFAULT_MODEL_DICT)
    return model


def model_from_dict(dict, time_steps=None):
    """Create a model from a config dict."""
    name = dict["type"]
    kwargs = dict["kwargs"]
    if 'read_majority_threshold' in kwargs:
        kwargs.pop('read_majority_threshold')
    symbol = importlib.import_module("medaka.architectures")
    model = getattr(symbol, name)(**kwargs)
    return model


def build_model_torch(
        feature_len, num_classes, gru_size=128,
        classify_activation='softmax', time_steps=None):
    """Build a bidirectional GRU model.

    This function allows models saved from previous generations of medaka
    to be loaded and used with the current version of medaka.

    :param feature_len: int, number of features for each pileup column.
    :param num_classes: int, number of output class labels.
    :param gru_size: int, size of each GRU layer.
    :param classify_activation: str, currently ignored (always softmax).
    :param time_steps: int, number of pileup columns in a sample (ignored).

    :returns: `torch.nn.Module` object.

    """
    # create old style medaka model
    model_dict = {
        "type": "GRUModel",
        "kwargs": {
            "num_features": feature_len,
            "num_classes": num_classes,
            "gru_size": gru_size,
        },
    }

    model = model_from_dict(model_dict)
    return model


def build_model(*args, **kwargs):
    """Catch old format TF model builder commands."""
    raise ValueError("Model format is not supported by medaka v2.x.")
