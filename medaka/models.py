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

    :returns: str, filepath to model file.
    """
    if os.path.exists(model):  # model is path to model file
        return model
    elif model not in medaka.options.allowed_models:
        raise ValueError(
            "Model {} is not a known model or existant file.".format(model))
    else:
        # check for model in model stores
        fname = '{}_model.hdf5'.format(model)
        fps = [
            os.path.join(ms, fname)
            for ms in medaka.options.model_stores]
        for fp in fps:
            if os.path.exists(fp):
                return fp

        # try to download model
        url = medaka.options.model_url_template.format(
            pkg=__package__, subdir=medaka.options.model_subdir, fname=fname)
        try:
            data = requests.get(url).content
            # check data is a model
            with tempfile.TemporaryDirectory() as tmpdir:
                tmp_file = os.path.join(tmpdir, "tmp_model.hdf5")
                with open(tmp_file, 'wb') as tmp_model:
                    tmp_model.write(data)
                with medaka.datastore.DataStore(tmp_file) as ds:
                    ds.get_meta('model_function')
        except Exception:
            raise DownloadError(
                "The model file for {} is not already installed and "
                "could not be downloaded. Check you are connected to"
                " the internet and try again.".format(model))
        else:
            # save the model
            for fp in fps:  # try saving the model
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
            raise RuntimeError(msg.format(model, ' or '.join(fps), url))
    raise RuntimeError("Model resolution failed")


def load_model(fname, time_steps=None, allow_cudnn=True):
    """Load a model from an .hdf file.

    :param fname: .hdf file containing model (or model name).
    :param time_steps: number of time points in RNN, `None` for dynamic.
    :param allow_cudnn: allow use of CuDNN optimizations.

    ..note:: keras' `load_model` cannot handle CuDNNGRU layers, hence this
        function builds the model then loads the weights.

    """
    fname = resolve_model(fname)
    with medaka.datastore.DataStore(fname) as ds:
        model_partial_function = ds.get_meta('model_function')
        model = model_partial_function(
            time_steps=time_steps, allow_cudnn=allow_cudnn)
        try:
            model.load_weights(fname)
        except ValueError():
            pass
        finally:
            return model


def build_model(feature_len, num_classes, gru_size=128,
                classify_activation='softmax', time_steps=None,
                allow_cudnn=True):
    """Build a bidirectional GRU model with CuDNNGRU support.

    CuDNNGRU implementation is claimed to give speed-up on GPU of 7x.
    The function will build a model capable of running on GPU with
    CuDNNGRU provided a) a GPU is present, b) the option has been
    allowed by the `allow_cudnn` argument; otherwise a compatible
    (but not CuDNNGRU accelerated model) is built.

    :param feature_len: int, number of features for each pileup column.
    :param num_classes: int, number of output class labels.
    :param gru_size: int, size of each GRU layer.
    :param classify_activation: str, activation to use in classification layer.
    :param time_steps: int, number of pileup columns in a sample.
    :param allow_cudnn: bool, opt-in to cudnn when using a GPU.

    :returns: `keras.models.Sequential` object.

    """
    import tensorflow as tf
    from tensorflow.keras.models import Sequential
    from tensorflow.keras.layers import Dense, GRU, CuDNNGRU, Bidirectional

    # Determine whether to use CuDNNGRU or not
    cudnn = False
    if tf.test.is_gpu_available(cuda_only=True) and allow_cudnn:
        cudnn = True
    logger.info("Building model with cudnn optimization: {}".format(cudnn))

    model = Sequential()
    input_shape = (time_steps, feature_len)
    for i in [1, 2]:
        name = 'gru{}'.format(i)
        # Options here are to be mutually compatible: train with CuDNNGRU
        # but allow inference with GRU (on cpu).
        # https://gist.github.com/bzamecnik/bd3786a074f8cb891bc2a397343070f1
        if cudnn:
            gru = CuDNNGRU(gru_size, return_sequences=True, name=name)
        else:
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
                   classify_activation='softmax', time_steps=None,
                   allow_cudnn=True):
    """Build a mock model that simply sums counts.

    :param feature_len: int, number of features for each pileup column.
    :param num_classes: int, number of output class labels.
    :param gru_size: int, size of each GRU layer.
    :param classify_activation: str, activation to use in classification layer.
    :param time_steps: int, number of pileup columns in a sample.
    :param allow_cudnn: bool, opt-in to cudnn when using a GPU.

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
