"""Creation and loading of models."""

import os
import pathlib
import tempfile

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
    suffixes = ("_model.hdf5", "_model.tar.gz")
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


def build_model(feature_len, num_classes, gru_size=128,
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
    from tensorflow.keras.models import Sequential
    from tensorflow.keras.layers import Dense, GRU, Bidirectional

    #  Tensorflow2 uses a fast cuDNN implementation if a GPU is available
    #  and the arguments to the layer meet the CuDNN kernal requirements
    if len(tf.test.gpu_device_name()) > 0:
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


def build_majority(feature_len, num_classes, gru_size=128,
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
    from tensorflow.keras.models import Sequential
    from tensorflow.keras.layers import Lambda, Activation

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
