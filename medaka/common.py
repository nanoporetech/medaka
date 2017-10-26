import errno
import itertools
import os
import numpy as np
from collections import OrderedDict, namedtuple

import logging
logger = logging.getLogger(__name__)

# Codec for converting tview output to ints.
_gap_ = '*'
_ref_gap_ = '#'
_read_sep_ = ' '
decoding = _gap_ + 'acgtACGT' + _read_sep_
# store encoding in ordered dict as the order will always be the same
# (we need order to be the same for training and inference)
encoding = OrderedDict(((a, i) for i, a in enumerate(decoding)))

# Encoded (and possibly reindexed) pileup
Pileup = namedtuple('Pileup', ['bam', 'ref_name', 'reads', 'positions'])
LabelledPileup = namedtuple('LabelledPileup', ['pileups', 'labels', 'ref_seq'])
Sample = namedtuple('Sample', ['ref_name', 'features', 'labels', 'ref_seq', 'positions', 'label_probs'])

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


def rle(array, low_mem=False):
    """Calculate a run length encoding (rle), of an input vector.

    :param array: 1D input array.
    :param low_mem: use a lower memory implementation

    returns: structured array with fields `start`, `length`, and `value`.
    """
    if len(array.shape) != 1:
        raise TypeError("Input array must be one dimensional.")
    dtype = [('length', int), ('start', int), ('value', array.dtype)]

    if not low_mem:
        pos = np.where(np.diff(array) != 0)[0]
        pos = np.concatenate(([0], pos+1, [len(array)]))
        return np.fromiter(
            ((length, start, array[start]) for (length, start) in zip(pos[1:], pos[:-1])),
            dtype, count=len(pos) - 1,
        )
    else:
        def _gen():
            start = 0
            for key, group in itertools.groupby(array):
                length = sum(1 for x in group)
                yield length, start, key
                start += length
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
    end = 0
    for start in range(0, a.shape[axis] - window + 1, step):
        end = start + window
        slicee[axis] = slice(start, end)
        yield a[slicee]
    # yield the remainder with the same window size
    if a.shape[axis] > end:
        start = a.shape[axis] - window
        slicee[axis] = slice(start, a.shape[axis])
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

