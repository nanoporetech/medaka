import errno
import itertools
import os
import numpy as np

import logging
logger = logging.getLogger(__name__)


def segment_limits(start, end, segment_len=20000, overlap_len=1000):
    """Generate segments of a range [0, end_point].

    :param start: startpoint of range.
    :param end: endpoint of range.
    :param segment_len: length of resultant segments.
    :param overlap_len: length of overlap between segments.
    :yields: tuples (start, end) indicating overlapping segments of
        the input range.
    """
    for n in range(start, end, segment_len):
        yield max(start, n - overlap_len), min(n + segment_len, end - 1)


def rle(array):
    """Calculate a run length encoding (rle), of an input vector.

    returns: structured array with fields `start`, `length`, and `value`.
    """
    if len(array.shape) != 1:
        raise TypeError("Input array must be one dimensional.")

    def _gen():
        start = 0
        for key, group in itertools.groupby(array):
            length = sum(1 for x in group)
            yield length, start, key
            start += length

    dtype = [('length', int), ('start', int), ('value', array.dtype)]
    return np.fromiter(_gen(), dtype=dtype)


def sliding_window(a, window=3, step=1, axis=0):
    """Sliding window across an array.

    :param a: input array.
    :param window: window length.
    :param step: step length between consecutive windows.
    :param axis: axis over which to apply window.

    :yields: windows of the input array.
    """
    slicee = [slice(None)] * a.ndim
    for i in range(0, a.shape[axis] - window + 1, step):
        slicee[axis] = slice(i, i + window)
        yield a[slicee]


def get_common_index(arrays):
    """Create a common index for a list of arrays.

    :param arrays: input arrays of an indentical type.

    :returns: array of index.
    """
    dtype = arrays[0].dtype
    for arr in arrays:
        if arr.dtype != dtype:
            raise TypeError("Arrays do not have a matching type.")
    return np.sort(np.unique(np.concatenate(arrays)))


def mkdir_p(path, info=None):
    """Make a directory if it doesn't exist."""
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            if info is not None:
                info = " {}".format(info)
            logging.warn("The path {} exists.{}".format(path, info))
            pass
        else:
            raise
