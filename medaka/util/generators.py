import itertools
import functools
import threading
import numpy as np


class ThreadsafeIter:
    """Takes an iterator and makes it thread-safe by
    serializing call to `next` method.

    `<http://anandology.com/blog/using-iterators-and-generators/>`_
    """
    def __init__(self, it):
        self.it = it
        self.lock = threading.Lock()

    def __iter__(self):
        return self

    def __next__(self):
        with self.lock:
            return next(self.it)


def threadsafe_generator(f):
    """A decorator that takes a generator function and makes it
    thread-safe.

    `<http://anandology.com/blog/using-iterators-and-generators/>`_
    """
    @functools.wraps(f)
    def g(*args, **kwargs):
        return ThreadsafeIter(f(*args, **kwargs))
    return g


def stack_slice(it, n):
    """stack first n array-like elements from iterable
    along first dimension (equal size elements required)

    :param it: iterable
    :param n: upper limit of slice into iterable
    :returns: numpy ndarray
    """
    while True:
        yield np.stack(itertools.islice(it, n))


def _serve_windows_once(arr, window_size):
    """generator yielding subarrays constituting
    a sliding window along the first dimension of an array

    :param arr: numpy ndarray
    :param window_size: int (odd) sliding window size
    :returns: numpy ndarray

    """
    assert window_size % 2 != 0
    for i in range(arr.shape[0] - window_size + 1):
        yield arr[i: i + window_size]


@threadsafe_generator
def serve_windows(arr, window_size):
    """returns generator infinitely yielding subarrays constituting
    a sliding window along the first dimension of the array

    :param arr: numpy ndarray
    :param window_size: int (odd) sliding window size
    """
    return itertools.cycle(_serve_windows_once(arr, window_size))


def serve_centrals(label, window_size):
    """returns generator infinitely yielding the central value
    from a sliding window along the first dimension of an array

    :param arr: numpy ndarray
    :param window_size: int (odd) sliding window size
    """
    assert window_size % 2 != 0
    centrals_iter = iter(label[window_size // 2: -window_size // 2 + 1])
    return itertools.cycle(centrals_iter)


@threadsafe_generator
def serve_sample(data, label, window_size):
    """return generator infinitely yielding windowed (data, label) tuples

    :param data: 2D numpy ndarray feature array (each row is a feature)
    :param label: 2D numpy ndarray label array (each row is a label)
    :param window_size: int (odd) sliding window size
    :returns: tuple (data, label)

        - data: numpy ndarray (shape=(window_size, feature dim))
        - label: numpy ndarray (shape=(window_size, label dim))
    """
    return zip(serve_windows(data, window_size),
               serve_centrals(label, window_size))


@threadsafe_generator
def serve_data_batch(data, batch_size, window_size):
    """return generator infinitely yielding batches of windowed data

    :param arr: 2D numpy ndarray feature array (each row is a feature)
    :param batch_size: int batch size
    :param window_size: int (odd) sliding window size
    :returns: numpy ndarray (shape=(batch_size, window_size, feature dim))
    """
    return stack_slice(serve_windows(data, window_size), batch_size)


def serve_label_batch(label, batch_size, window_size):
    """return generator infinitely yielding batches of central
    values from a sliding window along first dimension of array

    :param arr: 2D numpy ndarray label array (each row is a value)
    :param batch_size: int batch size
    :param window_size: int (odd) sliding window size
    :returns: numpy ndarray (shape=(batch_size, value dim))
    """
    return stack_slice(serve_centrals(label, window_size), batch_size)


@threadsafe_generator
def serve_sample_batch(data, label, batch_size, window_size):
    """return generator infinitely yielding tuples of
    (window data batch, central label value batch)

    :param data: 2D numpy ndarray feature array (each row is a feature)
    :param label: 2D numpy ndarray label array (each row is a label)
    :param batch_size: int batch size
    :param window_size: int (odd) sliding window size
    :returns: tuple (data, label)

        - data: numpy ndarray (shape=(batch_size, window_size, feature dim))
        - label: numpy ndarray (shape=(window_size, label dim))

    """
    return zip(serve_data_batch(data, batch_size, window_size),
               serve_label_batch(label, batch_size, window_size))
