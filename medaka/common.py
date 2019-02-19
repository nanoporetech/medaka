from collections import OrderedDict, namedtuple
import errno
import functools
import itertools
import logging
import os
from pkg_resources import resource_filename
import queue
import re
import threading

import numpy as np
import pysam

# don't import keras here, as it means a slow import of tensorflow
#from keras.utils.np_utils import to_categorical

# Codec for converting tview output to ints.
#TODO: this can likely be renomved
_gap_ = '*'
_ref_gap_ = '#'
_read_sep_ = ' '
_alphabet_ = 'ACGT'
_extra_bases_ = 'NMRWSKYWS' #TODO: put this in htslib order
#TODO: change name of these
decoding = _gap_ + _alphabet_.lower() + _alphabet_.upper() + _read_sep_ + _extra_bases_
# store encoding in ordered dict as the order will always be the same
# (we need order to be the same for training and inference)
encoding = OrderedDict(((a, i) for i, a in enumerate(decoding)))

AlignPos = namedtuple('AlignPos', ('qpos', 'qbase', 'rpos', 'rbase'))
ComprAlignPos = namedtuple('AlignPos', ('qpos', 'qbase', 'qlen', 'rpos', 'rbase', 'rlen'))

#TODO: refactor all this..
# We might consider creating a Sample class with methods to encode, decode sample name
# and perhaps to write/read Samples to/from hdf etc
_Sample = namedtuple('Sample', ['ref_name', 'features', 'labels', 'ref_seq', 'positions', 'label_probs'])
class Sample(_Sample):

    def _get_pos(self, index):
        p = self.positions
        return p['major'][index], p['minor'][index]

    @property
    def first_pos(self):
        """Zero-based first reference co-ordinate."""
        return self._get_pos(0)

    @property
    def last_pos(self):
        """Zero-based (end inclusive) last reference co-ordinate."""
        return self._get_pos(-1)

    @property
    def span(self):
        """Size of sample in terms of reference positions."""
        return self.last_pos[0] - self.first_pos[0]

    @property
    def is_empty(self):
        return self.size == 0

    @property
    def size(self):
        return len(self.positions)

    @property
    def name(self):
        """Create zero-based (end inclusive) samtools-style region string."""
        fmaj, fmin = self.first_pos
        lmaj, lmin = self.last_pos
        return '{}:{}.{}-{}.{}'.format(
            self.ref_name, fmaj, fmin, lmaj + 1, lmin)


    @staticmethod
    def decode_sample_name(name):
        """Decode a the result of Sample.name into a dict.

        :param key: `Sample` object.
        :returns: dict.
        """
        d = None
        name_decoder = re.compile(r"(?P<ref_name>.+):(?P<start>\d+\.\d+)-(?P<end>\d+\.\d+)")
        m = re.match(name_decoder, name)
        if m is not None:
            d = m.groupdict()
        return d

    def chunks(self, chunk_len=1000, overlap=200):
        """Create overlapping chunks of self.

        :param chunk_len: chunk length (number of columns)
        :param overlap: overlap length.

        :yields: chunked :py:class:`Sample` instances.
        """
        chunker = functools.partial(sliding_window,
            window=chunk_len, step=chunk_len - overlap, axis=0)
        sample = self._asdict()
        chunkers = {
            field: chunker(sample[field])
            if sample[field] is not None else itertools.repeat(None)
            for field in sample.keys()
        }

        for i, pos in enumerate(chunkers['positions']):
            fields = set(sample.keys()) - set(['positions', 'ref_name'])
            new_sample = {
                'positions':pos, 'ref_name':sample['ref_name']
            }
            for field in fields:
                new_sample[field] = next(chunkers[field])
            yield Sample(**new_sample)


#TODO: refactor this
_Region = namedtuple('Region', 'ref_name start end')
class Region(_Region):

    @property
    def name(self):
        """Samtools-style region string, zero-base end exclusive."""
        return self.__str__()


    def __str__(self):
        # This will be zero-based, end exclusive
        return '{}:{}-{}'.format(self.ref_name, self.start, self.end)


    @property
    def size(self):
        return self.end - self.start


    @classmethod
    def from_string(cls, region):
        """Parse region strings into `Region` objects.

        :param regions: iterable of str

        >>> parse_regions(['Ecoli'])[0]
        Region(ref_name='Ecoli', start=None, end=None)
        >>> parse_regions(['Ecoli:1000-2000'])[0]
        Region(ref_name='Ecoli', start=1000, end=2000)
        >>> parse_regions(['Ecoli:1000'])[0]
        Region(ref_name='Ecoli', start=0, end=1000)
        >>> parse_regions(['Ecoli:500-'])[0]
        Region(ref_name='Ecoli', start=500, end=None)

        """
        if ':' not in region:
            ref_name, start, end = region, None, None
        else:
            start, end = None, None
            ref_name, bounds = region.split(':')
            if bounds[0] == '-':
                start = 0
                end = int(bounds[1:])
            elif bounds[-1] == '-':
                start = int(bounds[:-1])
                end = None
            else:
                start, end = [int(b) for b in bounds.split('-')]
        return cls(ref_name, start, end)


    def split(region, size, overlap=0):
        """Split region into sub-regions.

        :param size: size of sub-regions.
        :param overlap: overlap between ends of sub-regions.

        :returns: a list of sub-regions.

        """
        regions = [
            Region(region.ref_name, start, stop) for (start, stop) in
            segment_limits(region.start, region.end, segment_len=size, overlap_len=overlap)
        ]
        # correct end co-ordinate of the last
        last = regions[-1]
        regions[-1] = Region(last.ref_name, last.start, last.end + 1)
        return regions


def get_regions(bam, region_strs=None):
    """Create `Region` objects from a bam and region strings.

    :param bam: `.bam` file.
    :param region_strs: iterable of str in zero-based (samtools-like)
        region format e.g. ref:start-end or filepath containing a
        region string per line.

    :returns: list of `Region` objects.
    """
    with pysam.AlignmentFile(bam) as bam_fh:
        ref_lengths = dict(zip(bam_fh.references, bam_fh.lengths))
    if region_strs is not None:
        if os.path.isfile(region_strs[0]):
            with open(region_strs[0]) as fh:
                region_strs = [l.strip() for l in fh.readlines()]

        regions = []
        for r in (Region.from_string(x) for x in region_strs):
            start = r.start if r.start is not None else 0
            end = r.end if r.end is not None else ref_lengths[r.ref_name]
            regions.append(Region(r.ref_name, start, end))
    else:
        regions = [Region(ref_name, 0, end) for ref_name, end in ref_lengths.items()]

    return regions


def lengths_to_rle(lengths):
    runs = np.empty(len(lengths), dtype=[('start', int), ('length', int)])
    runs['length'][:] = lengths
    runs['start'][0] = 0
    runs['start'][1:] = np.cumsum(lengths)[:-1]
    return runs


def get_pairs(aln):
    """Return generator of pairs.

    :param aln: `pysam.AlignedSegment` object.
    :returns: generator of `ComprAlignPos` objects.
    """
    seq = aln.query_sequence
    pairs = (ComprAlignPos(qp, seq[qp], 1, rp, rb, 1) if qp is not None
            else ComprAlignPos(qp, None, 1, rp, rb, 1)
            for qp, rp, rb in aln.get_aligned_pairs(with_seq=True))

    return pairs


def seq_to_hp_lens(seq):
    """Return array of HP lengths for every position in the sequence

    :param seq: str, sequence
    :returns: `np.ndarray` containing HP lengths (the same length as seq).
    """
    q_rle = rle(np.fromiter(seq, dtype='U1', count=len(seq)), low_mem=True)
    # get rle encoding for every position (so we can slice with qb instead of
    # searching)
    qlens = np.repeat(q_rle['length'], q_rle['length'])
    return qlens


def get_pairs_with_hp_len(aln, ref_seq):
    """Return generator of pairs in which the qbase is not just the base,
    but is decorated with extra information. E.g. an A which is part of a basecalled
    6mer HP would be AAAAAA.

    :param aln: `pysam.AlignedSegment` object.
    :param ref_seq: str containing reference sequence or result of
        `seq_to_hp_lens` (ref_seq).

    :yields: `ComprAlignPos` objects.
    """
    seq = aln.query_sequence
    qlens = seq_to_hp_lens(seq)

    rlens = seq_to_hp_lens(ref_seq) if isinstance(ref_seq, str) else ref_seq

    h = lambda ar, i, alt: ar[i] if i is not None else alt
    for qp, rp, rb in aln.get_aligned_pairs(with_seq=True):
        a = ComprAlignPos(qpos=qp, qbase=h(seq, qp, None), qlen=h(qlens, qp, 1),
                          rpos=rp, rbase=rb, rlen=h(rlens, rp, 1))
        yield a


def yield_compressed_pairs(aln, ref_rle):
    """Yield ComprAlignPos objects for aligned pairs of an AlignedSegment.

    :param aln: pysam AlignedSegment.
    :param ref_rle: structured array with columns 'start' and 'length' encoding ref.
    :yields: `ComprAlignPos` objects.
    """
    seq = aln.query_sequence
    qlens = aln.query_qualities
    if qlens is None:
        raise ValueError('Found an alignment without qscores, try filter to only primary alignments.')
    qrle = lengths_to_rle(qlens)
    h = lambda ar, i, alt: ar[i] if i is not None else alt
    for qp, rp, rb in aln.get_aligned_pairs(with_seq=True):
        a = ComprAlignPos(qpos=qp, qbase=h(seq, qp, None), qlen=h(qlens, qp, 1),
                          rpos=rp, rbase=rb, rlen=h(ref_rle['length'], rp, 1))
        yield a


def get_sample_overlap(s1, s2):
    """Calculate overlap (in pileup columns) between two `Sample` objects.

    :param s1: First `Sample` object.
    :param s2: Second `Sample` object.
    :returns: (int end1, int start2) such that s1.positions[:end1] could be
        concatenated with s2.positions[start2:] to join them without gaps or overlap.
    """
    ovl_start_ind1 = np.searchsorted(s1.positions, s2.positions[0])
    if ovl_start_ind1 == len(s1.positions):
        # they don't overlap
        print('{} and {} do not overlap'.format(s1.name, s2.name))
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
        yield a[tuple(slicee)]
    # yield the remainder with the same window size
    if a.shape[axis] > end:
        start = a.shape[axis] - window
        slicee[axis] = slice(start, a.shape[axis])
        yield a[tuple(slicee)]


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


def background_generator(generator, max_size, daemon=True):
    """Run a generator in background thread.
    
    :param max_size: maximum number of items from the generator to cache.
    :param daemon: run generator in daemon thread

    """
    results = queue.Queue(maxsize=max_size)
    have_data = threading.Event()
    have_data.set()
    def gen_to_queue():
        for item in generator:
            results.put(item)
        have_data.clear()

    thread = threading.Thread(target=gen_to_queue)
    thread.daemon = daemon
    thread.start()
    # yield results concurrently whilst filling queue
    while have_data.is_set():
        yield results.get()

    # empty remaining items in queue when input has finished
    while not results.empty():
        yield results.get()

    thread.join(5)


def grouper(gen, batch_size=4):
    """Group together elements of an iterable, yielding remainder without padding."""
    while True:
        batch = []
        for i in range(batch_size):
            try:
                batch.append(next(gen))
            except StopIteration:
                if len(batch) > 0:
                    yield batch
                raise StopIteration
        yield batch


def print_data_path():
    """Print data directory containing models"""
    print(resource_filename(__package__, 'data'))


def get_named_logger(name):
    logger = logging.getLogger('{}.{}'.format(__package__, name))
    logger.name = name
    return logger
