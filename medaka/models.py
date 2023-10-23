"""Creation and loading of models."""

import itertools
import os
import pathlib
import tempfile

import pysam
import requests

import medaka.common
import medaka.datastore
import medaka.options

logger = medaka.common.get_named_logger('ModelLoad')


class DownloadError(ValueError):
    """Raised when model is unsuccessfully downloaded."""


def resolve_model(model):
    """Resolve a model filepath, downloading known models if necessary.

    :param model_name: str, model filepath or model ID

    :returns: str: filepath to hdf model file or TF model directory.
    """
    suffixes = ("_model.tar.gz", "_model.hdf5")
    if os.path.exists(model):  # model is path to model file
        return model
    elif model not in medaka.options.allowed_models:
        raise ValueError(
            "Model {} is not a known model or existant file.".format(model))
    else:
        # check for model in model stores
        for suffix in suffixes:
            fname = '{}{}'.format(model, suffix)
            for ms in medaka.options.model_stores:
                fp = os.path.join(ms, fname)
                if os.path.exists(fp):
                    return fp

        # try to download model
        download_errors = 0
        data = None
        for suffix in suffixes:
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
        if download_errors == len(suffixes):
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


def model_from_basecaller(fname, variant=False):
    """Determine correct medaka model from basecaller output file.

    :param fname: a basecaller output (.sam/.bam/.cram/.fastq).
    :param variant: whether to return variant model (otherwise consensus).

    There are slight differences is the search strategy for .bam and .fastq
    files due to differences in what information is available in each.

    For .bam (and related) files the DS subfield of the read group header
    is examined to find the "basecall_model=" key=value entry. The found
    model is returned vebatim without fursther checks.

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
            # model is embedded in RG:Z: tag of comment as
            # <run_id>_<model>_<barcode>, but model has _
            # characters in also so search for known models
            try:
                read_group = rec.comment.split("RG:Z:")[1].split()[0]
                for model in known_models:
                    if model in read_group:
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
    if ext == ".hdf5":
        return medaka.datastore.ModelStore(fname)
    elif ext == ".gz":
        return medaka.datastore.ModelStoreTF(fname)
    else:
        raise ValueError(
            "Model {} does not have .hdf5 or .gz extension.".format(fname))


def build_model(
        feature_len, num_classes, gru_size=128,
        classify_activation='softmax', time_steps=None):
    """Build a bidirectional GRU model with CuDNNGRU support.

    CuDNNGRU implementation is claimed to give speed-up on GPU of 7x.
    The function will build a model capable of running on GPU with
    CuDNNGRU provided a) a GPU is present, b) the arguments to the
    keras layer meet the CuDNN kernal requirements for cudnn = True;
    otherwise a compatible (but not CuDNNGRU accelerated model) is built.

    :param feature_len: int, number of features for each pileup column.
    :param num_classes: int, number of output class labels.
    :param gru_size: int, size of each GRU layer.
    :param classify_activation: str, activation to use in classification layer.
    :param time_steps: int, number of pileup columns in a sample.

    :returns: `keras.models.Sequential` object.

    """
    import tensorflow as tf
    from tensorflow.keras.layers import \
        Bidirectional, Dense, GRU
    from tensorflow.keras.models import Sequential
    #  Tensorflow2 uses a fast cuDNN implementation if a GPU is available
    #  and the arguments to the layer meet the CuDNN kernal requirements
    if tf.config.list_physical_devices('GPU'):
        logger.info("GPU available: building model with cudnn optimization")

    model = Sequential()
    input_shape = (time_steps, feature_len)
    for i in [1, 2]:
        name = 'gru{}'.format(i)
        gru = GRU(
            gru_size, reset_after=True, recurrent_activation='sigmoid',
            return_sequences=True, name=name)
        model.add(Bidirectional(gru, input_shape=input_shape))

    # see keras #10417 for why we specify input shape
    model.add(Dense(
        num_classes, activation=classify_activation, name='classify',
        input_shape=(time_steps, 2 * gru_size)
    ))

    return model


def build_majority(
        feature_len, num_classes, gru_size=128,
        classify_activation='softmax', time_steps=None):
    """Build a mock model that simply sums counts.

    :param feature_len: int, number of features for each pileup column.
    :param num_classes: int, number of output class labels.
    :param gru_size: int, size of each GRU layer.
    :param classify_activation: str, activation to use in classification layer.
    :param time_steps: int, number of pileup columns in a sample.

    :returns: `keras.models.Sequential` object.

    """
    import tensorflow as tf
    from tensorflow.keras.layers import \
        Activation, Lambda
    from tensorflow.keras.models import Sequential

    def sum_counts(f):
        """Sum forward and reverse counts."""
        # TODO write to handle multiple dtypes
        # acgtACGTdD
        # sum base counts
        b = f[:, :, 0:4] + f[:, :, 4:8]
        # sum deletion counts (indexing in this way retains correct shape)
        d = f[:, :, 8:9] + f[:, :, 9:10]
        return tf.concat([d, b], axis=-1)

    model = Sequential()
    model.add(Lambda(sum_counts, output_shape=(time_steps, num_classes)))
    model.add(Activation('softmax'))
    return model


default_model = 'two_layer_bidirectional_CuDNNGRU'
model_builders = {
    'two_layer_bidirectional_CuDNNGRU': build_model,
    'majority_vote': build_majority,
}
