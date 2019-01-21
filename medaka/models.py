def build_legacy_model(chunk_size, feature_len, num_classes, gru_size=128):
    """Builds a bidirectional GRU model
    :param chunk_size: int, number of pileup columns in a sample.
    :param feature_len: int, number of features for each pileup column.
    :param num_classes: int, number of output class labels.
    :param gru_size: int, size of each GRU layer.
    :returns: `keras.models.Sequential` object.
    """

    from keras.models import Sequential
    from keras.layers import Dense, GRU
    from keras.layers.wrappers import Bidirectional

    model = Sequential()
    input_shape=(chunk_size, feature_len)
    for i in [1, 2]:
        name = 'gru{}'.format(i)
        gru = GRU(gru_size, activation='tanh', return_sequences=True, name=name)
        model.add(Bidirectional(gru, input_shape=input_shape))

    # see keras #10417 for why we specify input shape
    model.add(Dense(
        num_classes, activation='softmax', name='classify',
        input_shape=(chunk_size, 2 * feature_len)
    ))
    return model


def build_model(chunk_size, feature_len, num_classes, gru_size=128):
    """Builds a bidirectional GRU model. Uses CuDNNGRU for additional
    speed-up on GPU (claimed 7x).

    :param chunk_size: int, number of pileup columns in a sample.
    :param feature_len: int, number of features for each pileup column.
    :param num_classes: int, number of output class labels.
    :param gru_size: int, size of each GRU layer.

    :returns: `keras.models.Sequential` object.
    """

    from keras import backend as K
    from keras.models import Sequential
    from keras.layers import Dense, GRU, CuDNNGRU, Bidirectional

    # if we can see a gpu, use CuDNNGRU for speed
    cudnn = False
    if len(K.tensorflow_backend._get_available_gpus()) > 1:
        cudnn = True

    model = Sequential()
    input_shape=(chunk_size, feature_len)
    for i in [1, 2]:
        name = 'gru{}'.format(i)
        # Options here are to be mutually compatible: train with CuDNNGRU
        # but allow inference with GRU (on cpu).
        # https://gist.github.com/bzamecnik/bd3786a074f8cb891bc2a397343070f1
        if cudnn:
            gru = CuDNNGRU(gru_size, return_sequences=True, name=name)
        else:
            gru = GRU(gru_size, reset_after=True, recurrent_activation='sigmoid', return_sequences=True, name=name)
        model.add(Bidirectional(gru, input_shape=input_shape))

    # see keras #10417 for why we specify input shape
    model.add(Dense(
        num_classes, activation='softmax', name='classify',
        input_shape=(chunk_size, 2 * feature_len)
    ))

    return model

model_builders = {
    'two_layer_bidirectional_CuDNNGRU': build_model,
    'two_layer_bidirectional_GRU': build_legacy_model,
}

default_model = 'two_layer_bidirectional_CuDNNGRU'
