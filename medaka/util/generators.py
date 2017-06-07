import itertools
import numpy as np


def serve_sample(data, label, window_size):
    """Generate windowed (data, label) tuples infinitely

    :param data: 2D feature array (each row is a feature)
    :param label: 1D label array
    :param window_size: int (odd) window of positions flanking predicted label
    :returns: tuple (data, label)

        - data: 2D array (shape=(window_size, feature_length))
        - label: 1D array (shape=window_size,)

    """
    assert window_size % 2 != 0
    i = 0
    while True:
        yield data[i:i + window_size], label[i + window_size // 2]
        i += 1
        if i == data.shape[0] - window_size + 1:
            i = 0


def serve_sample_batch(data, label, batch_size, window_size):
    """Generate batches of windowed (data, label) samples

    :param data: full feature array (each row is a feature)
    :param label: full label array
    :param batch_size: int batch size
    :param window_size: int (odd) window of positions flanking predicted label
    :returns: tuple (data batch, label batch)

        - data batch: 3D array (shape=(batch_size, window_size, feature_length))
        - label batch: 2D array (shape=batch_size, window_size,)

    """
    dataiter = serve_sample(data, label, window_size)
    while True:
        nb1, nb2 = itertools.tee(itertools.islice(dataiter, 0, batch_size))
        data_out = np.stack(d[0] for d in nb1)
        label_out = np.stack(d[1] for d in nb2)
        yield data_out, label_out


def serve_data(data, window_size):
    """Generate windowed data items infinitely

    :param data: full feature array (each row is a feature)
    :param window_size: int (odd) window of positions flanking predicted label
    :returns: 2D data array (shape=(window_size, feature_length))

    """
    assert window_size % 2 != 0
    i = 0
    while True:
        yield data[i:i + window_size]
        i += 1
        if i == data.shape[0] - window_size + 1:
            i = 0


def serve_data_batch(data, batch_size, window_size):
    """Generate batches of windowed data items

    :param data: full feature array (each row is a feature)
    :param label: full label array
    :param window_size: int (odd) window of positions flanking predicted label
    :returns: 3D data batch array (shape=(batch_size, window_size, feature_length))

    """
    dataiter = serve_data(data, window_size)
    while True:
        sl = itertools.islice(dataiter, 0, batch_size)
        data_out = np.stack(sl)
        yield data_out
