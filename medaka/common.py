import errno
import h5py
import itertools
import os
import numpy as np
import re
from Bio import SeqIO
from collections import OrderedDict, namedtuple, defaultdict

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

# We might consider creating a Sample class with methods to encode, decode sample name
# and perhaps to write/read Samples to/from hdf etc
Sample = namedtuple('Sample', ['ref_name', 'features', 'labels', 'ref_seq', 'positions', 'label_probs'])
_sample_name_decoder_ = re.compile(r"(?P<ref_name>\w+):(?P<start>\d+\.\d+)-(?P<end>\d+\.\d+)")


def encode_sample_name(sample):
    """Encode a `Sample` object into a str key.

    :param sample: `Sample` object.
    :returns: str
    """
    p = sample.positions
    key = '{}:{}.{}-{}.{}'.format(sample.ref_name,
                                  p['major'][0] + 1, p['minor'][0],
                                  p['major'][-1] + 1, p['minor'][-1])
    return key


def decode_sample_name(key):
    """Decode a the result of `encode_sample_name` into a dict.

    :param key: `Sample` object.
    :returns: dict
    """
    d = None
    m = re.match(_sample_name_decoder_, key)
    if m is not None:
        d = m.groupdict()
    return d


def write_sample_to_hdf(s, hdf_fh):
    grp = encode_sample_name(s)
    for field in s._fields:
        if getattr(s, field) is not None:
            data = getattr(s, field)
            if isinstance(data, np.ndarray) and isinstance(data[0], np.unicode):
                data = np.char.encode(data)
            hdf_fh['{}/{}'.format(grp, field)] = data


def load_sample_from_hdf(key, hdf_fh):
    s = {}
    for field in Sample._fields:
        pth = '{}/{}'.format(key, field)
        if pth in hdf_fh:
            s[field] = hdf_fh[pth][()]
            if isinstance(s[field], np.ndarray) and isinstance(s[field][0], type(b'')):
                s[field] = np.char.decode(s[field])
        else:
            s[field] = None
    return Sample(**s)


def get_sample_overlap(s1, s2):
    ovl_start_ind1 = np.searchsorted(s1.positions, s2.positions[0])
    if ovl_start_ind1 == len(s1.positions):
        # they don't overlap
        print('{} and {} do not overlap'.format(encode_sample_name(s1), encode_sample_name(s2)))
        return None, None

    ovl_end_ind2 = np.searchsorted(s2.positions, s1.positions[-1], side='right')
    pos1_ovl = s1.positions[ovl_start_ind1:]
    pos2_ovl = s2.positions[0:ovl_end_ind2]
    assert len(pos1_ovl) == len(pos2_ovl)
    overlap_len = len(pos1_ovl)
    pad_1 = overlap_len // 2
    pad_2 = overlap_len - pad_1
    end_1_ind = ovl_start_ind1 + pad_1
    start_2_ind = ovl_end_ind2 - pad_2

    contr_1 = s1.positions[ovl_start_ind1:end_1_ind]
    contr_2 = s2.positions[start_2_ind:ovl_end_ind2]
    assert len(contr_1) + len(contr_2) == overlap_len

    return end_1_ind, start_2_ind


def get_sample_index_from_files(fnames, filetype='hdf'):
    ref_names = defaultdict(list)
    if filetype == 'hdf':
        fhs = (h5py.File(f) for f in fnames)
    elif filetype == 'fasta':
        fhs = (SeqIO.index(f, 'fasta') for f in fnames)
    for fname, fh in zip(fnames, fhs):
        for key in fh:
            d = decode_sample_name(key)
            if d is not None:
                d['key'] = key
                d['filename'] = fname
                ref_names[d['ref_name']].append(d)
    for fh in fhs:
        fh.close()

    # sort dicts so that refs are in order and within a ref, chunks are in order
    ref_names_ordered = OrderedDict()
    for ref_name in sorted(ref_names.keys()):
        sorter = lambda x: float(x['start'])
        ref_names[ref_name].sort(key=sorter)
        ref_names_ordered[ref_name] = ref_names[ref_name]

    return ref_names_ordered


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

