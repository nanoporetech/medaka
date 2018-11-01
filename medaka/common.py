from collections import OrderedDict, namedtuple, defaultdict, Counter
import errno
import functools
import itertools
import os
from pkg_resources import resource_filename
import re
import threading
from timeit import default_timer as now
import yaml

from Bio import SeqIO
import h5py
from keras.utils.np_utils import to_categorical
import numpy as np
import pysam


import logging
logger = logging.getLogger(__name__)

# Codec for converting tview output to ints.
#TODO: this can likely be renomved
_gap_ = '*'
_ref_gap_ = '#'
_read_sep_ = ' '
_alphabet_ = 'ACGT'
_extra_bases_ = 'N'
decoding = _gap_ + _alphabet_.lower() + _alphabet_.upper() + _read_sep_ + _extra_bases_
# store encoding in ordered dict as the order will always be the same
# (we need order to be the same for training and inference)
encoding = OrderedDict(((a, i) for i, a in enumerate(decoding)))

AlignPos = namedtuple('AlignPos', ('qpos', 'qbase', 'rpos', 'rbase'))
ComprAlignPos = namedtuple('AlignPos', ('qpos', 'qbase', 'qlen', 'rpos', 'rbase', 'rlen'))

_feature_opt_path_ = 'medaka_features_kwargs'
_model_opt_path_ = 'medaka_model_kwargs'
_sample_path_ = 'samples'
_label_decod_path_ = 'medaka_label_decoding'
_feature_decoding_path_ = 'medaka_feature_decoding'
_label_counts_path_ = 'medaka_label_counts'
_feature_batches_path_ = 'medaka_feature_batches'
_label_batches_path_ = 'medaka_label_batches'


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
        """Create overlapping chunks of
        
        :param chunk_len: chunk length (number of columns)
        :param overlap: overlap length.

        :yields: chunked `Sample`s.
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


#TODO refactor the below two function into single interface
def write_sample_to_hdf(s, hdf_fh):
    """Write a sample to HDF5.

    :param s: `Sample` object.
    :param hdf_fh: `h5py.File` object.
    """
    grp = s.name
    for field in s._fields:
        if getattr(s, field) is not None:
            data = getattr(s, field)
            if isinstance(data, np.ndarray) and isinstance(data[0], np.unicode):
                data = np.char.encode(data)
            hdf_fh['{}/{}/{}'.format(_sample_path_, grp, field)] = data


def write_samples_to_hdf(fname, samples):
    """Write samples to hdf, ensuring a sample is not written twice and maintaining
    a count of labels seen.

    :param fname: str, output filepath.
    :param samples: iterable of `Sample` objects.
    """
    logging.info("Writing samples to {}".format(fname))
    samples_seen = set()
    labels_counter = Counter()
    with h5py.File(fname, 'w') as hdf:
        for s in samples:
            s_name = s.name
            logging.debug("Written sample {}".format(s_name))
            if s_name not in samples_seen:
                write_sample_to_hdf(s, hdf)
            else:
                logging.debug('Not writing {} as it is present already'.format(s_name))
            samples_seen.add(s_name)
            if s.labels is not None:
                if len(s.labels.dtype) == 2:
                    labels_counter.update([tuple(l) for l in s.labels])
                else:
                    labels_counter.update(s.labels)
        hdf[_label_counts_path_] = yaml.dump(labels_counter)

    h = lambda l: (decoding[l[0]], l[1]) if type(l) == tuple else l
    logging.info("Label counts:\n{}".format('\n'.join(
        ['{}: {}'.format(h(label), count) for label, count in labels_counter.items()]
    )))


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


def load_yaml_data(fname, group):
    """Load a yml str either from hdf or .yml file"""
    data = None
    if os.path.splitext(fname)[-1] == '.yml':
        with open(fname) as fh:
            yml_str = fh.read()
        d = yaml.load(yml_str)
        if group in d:
            data = d[group]
    else:
        with h5py.File(fname) as hdf:
            if group in hdf:
                yml_str = hdf[group][()]
                data = yaml.load(yml_str)
    return data


def write_yaml_data(fname, data):
    """Save a data structure to a yml str either within a hdf or .yml file"""
    if os.path.splitext(fname)[-1] == '.yml':
        with open(fname, 'w') as fh:
            fh.write(yaml.dump(data))
    else:
        with h5py.File(fname) as hdf:
            for group, d in data.items():
                if group in hdf:
                    del hdf[group]
                hdf[group] = yaml.dump(d)


def load_sample_from_hdf(key, hdf_fh):
    """Load `Sample` object from HDF5

    :param key: str, sample name.
    :param hdf_fh: `h5.File` object.
    :returns: `Sample` object.
    """
    s = {}
    for field in Sample._fields:
        pth = '{}/{}/{}'.format(_sample_path_, key, field)
        if pth in hdf_fh:
            s[field] = hdf_fh[pth][()]
            if isinstance(s[field], np.ndarray) and isinstance(s[field][0], type(b'')):
                s[field] = np.char.decode(s[field])
        else:
            s[field] = None
    return Sample(**s)


def yield_from_feature_files(fnames, ref_names=None, index=None, samples=None):
    """Yield `Sample` objects from one or more feature files.

    :param fnames: iterable of str of filepaths.
    :ref_names: iterable of str, only process these references.
    :index: result of `get_sample_index_from_files`
    :samples: iterable of sample names to yield (in order in which they are supplied).
    :yields: `Sample` objects.
    """
    handles = { fname: h5py.File(fname, 'r') for fname in fnames}

    if samples is not None:
        # yield samples in the order they are asked for
        for sample, fname in samples:
            yield load_sample_from_hdf(sample, handles[fname])
    else:
        # yield samples sorted by ref_name and start
        if index is None:
            index = get_sample_index_from_files(fnames, 'hdf')
        if ref_names is None:
            ref_names = sorted(index.keys())
        for ref_name in ref_names:
            for d in index[ref_name]:
                yield load_sample_from_hdf(d['key'], handles[d['filename']])

    for fh in handles.values():
        fh.close()


def load_from_feature_files(fnames, key_fnames):
    """Yield `Sample` objects in order from the result of `save_feature_file`.

    :param fnames: iterable of str of filepaths.
    :key_fnames: iterable of (str sample key, fname) tuples.
    :returns: list of `Sample` objects.
    """
    t0 = now()
    handles = { fname: h5py.File(fname, 'r') for fname in fnames}
    t1 = now()
    logging.info('Creating h5py.File objs took {}s'.format(t1 - t0))

    t0 = now()
    samples = [ load_sample_from_hdf(key, handles[fname]) for (key, fname) in key_fnames]
    t1 = now()
    logging.info('Loading samples took {}s'.format(t1 - t0))
    for fh in handles.values():
        fh.close()
    return samples


def load_feature_file(fname, ref_names=None):
    """Load `Sample` objects from HDF5.

    :param fname: str, HDF5 filepath.
    :param ref_names: iterable of str, only fetch samples for this reference.
    :returns: tuple of `Sample` objects.
    """
    return tuple([s for s in yield_from_feature_files([fname], ref_names)])


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


def get_sample_index_from_files(fnames, filetype='hdf', max_samples=np.inf):
    """Create index of samples from one or more HDF5 or fasta files.

    :param fnames: iterable of str of filepaths.
    :param filetype: str, 'hdf' or 'fasta'.
    :param max_samples: int, stop after max_samples.
    :returns: dict of lists, indexed by reference.
        each list contains a list of samples dicts, and is sorted by sample start.
    """

    count = 0
    ref_names = defaultdict(list)
    if filetype == 'hdf':
        fhs = (h5py.File(f) for f in fnames)
    elif filetype == 'fasta':
        fhs = (SeqIO.index(f, 'fasta') for f in fnames)
    for fname, fh in zip(fnames, fhs):
        if count > max_samples:
            break
        keys = (k for k in fh[_sample_path_]) if filetype == 'hdf' else (k for k in fh)
        for key in keys:
            d = Sample.decode_sample_name(key)
            if d is not None:
                d['key'] = key
                d['filename'] = fname
                ref_names[d['ref_name']].append(d)
                count += 1
            if count > max_samples:
                break
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


def chunk_samples(samples, chunk_len=1000, overlap=200):
    """Return generator of chunked `Sample` objs"""
    return (c for s in samples for c in chunk_sample(s, chunk_len=chunk_len, overlap=overlap))


def chunk_sample(sample, chunk_len=1000, overlap=200):
    """Chunk `Sample` objects into smaller overlapping chunks.

    :param sample: `Sample` obj
    :param chunk_len: chunk length (number of columns)
    :param overlap: overlap length.

    :yields: `Sample` objs.
    """
    yield from sample.chunk(chunk_len, overlap)


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


def sample_to_x_y(s, encoding, max_label_len=np.inf):
    """Convert a `Sample` object into an x,y tuple for training.

    :param s: `Sample` object.
    :param encoding: dict of label encodings.
    :max_label_len: int, maximum label length, longer labels will be truncated.
    :returns: (np.ndarray of inputs, np.ndarray of labels)
    """
    if s.labels is None:
        raise ValueError("Cannot train without labels.")
    x = s.features
    # labels can either be unicode strings or (base, length) integer tuples
    if isinstance(s.labels[0], np.unicode):
        y = np.fromiter((encoding[l[:min(max_label_len, len(l))]]
                     for l in s.labels), dtype=int, count=len(s.labels))
    else:
        y = np.fromiter((encoding[tuple((l['base'],
                                     min(max_label_len, l['run_length'])))]
                     for l in s.labels), dtype=int, count=len(s.labels))
    y = y.reshape(y.shape + (1,))
    return x, y


@threadsafe_generator
def chain_thread_safe(*gens):
    """Threadsafe version of itertools.chain"""
    gen = itertools.chain(*gens)
    while True:
        yield next(gen)


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


@threadsafe_generator
def gen_train_batch(xy_gen, batch_size, name=''):
    """Yield training batches.
    """
    count=0
    while True:
        batch = [next(xy_gen) for i in range(batch_size)]
        xs, ys = zip(*batch)
        logging.debug("Yielding {} batch {}".format(name, count))
        count += 1
        yield np.stack(xs), np.stack(ys)


@threadsafe_generator
def yield_batches_from_hdf(h5, keys):
    """Yield training batches from HDF5.

    :param h5: `h5.File` object.
    :param keys: iterable of keys to include in the batch.
    :yields: (np.ndarray of inputs, np.ndarray of labels).
    """
    for i in keys:
        yield (h5['{}/{}'.format(_feature_batches_path_, i)][()],
               h5['{}/{}'.format(_label_batches_path_, i)][()])


@threadsafe_generator
def yield_batches_from_hdfs(handles, batches=None, sparse_labels=True, n_classes=None):
    """Yield batches of training features and labels indefinitely, shuffling between epochs.

    :param handles: {str filename: `h5.File`}.
    :param batches: iterable of (str filename, str batchname).
    :param sparse_labels: bool, if False, labels will be one-hot encoded.
    :param n_classes: int, number of classes for one-hot encoding.
    :yields: (np.ndarray of inputs, np.ndarray of labels).
    """
    if batches is None:
        batches = [(fname, k) for (fname, fh) in handles.items() for k in fh[_feature_batches_path_]]

    np.random.shuffle(batches)
    epoch = 0
    while True:
        batch = 0
        for (fname, i) in batches:
            h5 = handles[fname]
            t0 = now()
            xs = h5['{}/{}'.format(_feature_batches_path_, i)][()]
            ys = h5['{}/{}'.format(_label_batches_path_, i)][()]
            if not sparse_labels:
                ys = to_categorical(ys, num_classes=n_classes)
            t1 = now()
            logging.debug("Took {:5.3}s to load batch {} (epoch {})".format(t1-t0, batch, epoch))
            yield xs, ys
            batch += 1
        epoch += 1
        # shuffle between epochs
        np.random.shuffle(batches)


def print_data_path():
    """Print data directory containing models"""
    print(resource_filename(__package__, 'data'))
