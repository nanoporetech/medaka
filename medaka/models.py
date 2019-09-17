"""Creation and loading of models."""
import medaka.common
import medaka.datastore

logger = medaka.common.get_named_logger('ModelLoad')


def load_model(fname, time_steps=None, allow_cudnn=True):
    """Load a model from an .hdf file.

    :param fname: .hdf file containing model.
    :param time_steps: number of time points in RNN, `None` for dynamic.
    :param allow_cudnn: allow use of CuDNN optimizations.

    ..note:: keras' `load_model` cannot handle CuDNNGRU layers, hence this
        function builds the model then loads the weights.

    """
    with medaka.datastore.DataStore(fname) as ds:
        model_partial_function = ds.get_meta('model_function')
        model = model_partial_function(time_steps=time_steps,
                                       allow_cudnn=allow_cudnn)
        model.load_weights(fname)
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


default_model = 'two_layer_bidirectional_CuDNNGRU'
model_builders = {
    'two_layer_bidirectional_CuDNNGRU': build_model,
}
