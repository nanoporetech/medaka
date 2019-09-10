from collections import OrderedDict, namedtuple
from distutils.version import LooseVersion
import enum
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


ComprAlignPos = namedtuple('ComprAlignPos', ('qpos', 'qbase', 'qlen', 'rpos', 'rbase', 'rlen'))


class OverlapException(Exception):
    pass


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
            self.ref_name, fmaj, fmin, lmaj, lmin)

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


    @staticmethod
    def from_samples(samples):
        """Create a sample by concatenating an iterable of `Sample` objects.

        :param samples: iterable of `Sample` objects.
        :returns: `Sample` obj
        """
        samples = list(samples)
        for s1, s2 in zip(samples[0:-1], samples[1:]):
            rel = Sample.relative_position(s1, s2)
            if rel is not Sample.Relationship.forward_abutted:
                msg = 'Refusing to concatenate unordered/non-abutting samples {} and {} with relationship {}.'
                raise ValueError(msg.format(s1.name, s2.name, repr(rel)))

        # Sample.Relationship.forward_abutted guarantees all samples have the
        # same ref_name
        non_concat_fields = {'ref_name'}
        def concat_attr(attr):
            vals = [getattr(s, attr) for s in samples]
            if attr not in non_concat_fields:
                all_none = all([v is None for v in vals])
                c = np.concatenate(vals) if not all_none else None
            else:
                assert len(set(vals)) == 1
                c = vals[0]
            return c

        return Sample(**{attr:  concat_attr(attr) for attr in Sample._fields})


    class Relationship(enum.Enum):
            different_ref_name = 'Samples come from different reference contigs.'
            forward_overlap = 'The end of s1 overlaps the start of s2.'
            reverse_overlap = 'The end of s2 overlaps the start of s1.'
            forward_abutted = 'The end of s1 abuts the start of s2.'
            reverse_abutted = 'The end of s2 abuts the start of s1.'
            forward_gapped = 's2 follows s1 with a gab inbetween.'
            reverse_gapped = 's1 follows s2 with a gab inbetween.'
            s2_within_s1 = 's2 is fully contained within s1.'
            s1_within_s2 = 's1 is fully contained within s2.'


    @staticmethod
    def relative_position(s1, s2):
        """Check the relative position of two samples in genomic coordinates.

        :param s1, s2: `medaka.common.Sample` objs.

        :returns: `Relationship` enum member.
        """

        def ordered_abuts(s1, s2):
            """Check if end of s1 abuts the start of s2 (i.e. is adjacent but not overlapping)
            """
            s1_end_maj, s1_end_min = s1.last_pos
            s2_start_maj, s2_start_min = s2.first_pos
            if s2_start_maj == s1_end_maj + 1 and s2_start_min == 0:
                abuts = True
            elif s2_start_maj == s1_end_maj and s2_start_min == s1_end_min + 1:
                abuts = True
            else:
                abuts = False
            return abuts

        def ordered_contained(s1, s2):
            """Check if s2 is contained within s1.
            """
            return s2.first_pos >= s1.first_pos and s2.last_pos <= s1.last_pos

        def ordered_overlaps(s1, s2):
            """Check if end of s1 overlaps start of s2.
            """
            s1_end_maj, s1_end_min = s1.last_pos
            s2_start_maj, s2_start_min = s2.first_pos
            if s2_start_maj < s1_end_maj:  # we have overlap of major coordinates
                overlaps = True
            elif s2_start_maj == s1_end_maj and s2_start_min < s1_end_min + 1:
                # we have overlap of minor coordinates
                overlaps = True
            else:
                overlaps = False
            return overlaps

        def ordered_gapped(s1, s2):
            """Check if there is a gap between the end of s1 and the start of s2.
            """
            s1_end_maj, s1_end_min = s1.last_pos
            s2_start_maj, s2_start_min = s2.first_pos
            gapped = False
            if s2_start_maj > s1_end_maj + 1:  # gap in major
                gapped = True
            elif s2_start_maj > s1_end_maj and s2_start_min > 0:  # missing minors
                gapped = True
            elif s2_start_maj == s1_end_maj and s2_start_min > s1_end_min + 1:  # missing minors
                gapped = True
            return gapped

        if s1.ref_name != s2.ref_name:  # different ref_names
            return Sample.Relationship.different_ref_name

        s1_ord, s2_ord = sorted((s1, s2), key=lambda x: (x.first_pos, -x.size))
        is_ordered = s1_ord.name == s1.name

        # is one sample within the other?
        if ordered_contained(s1_ord, s2_ord):
            return Sample.Relationship.s2_within_s1 if is_ordered else Sample.Relationship.s1_within_s2

        # do samples abut?
        elif ordered_abuts(s1_ord, s2_ord):
            return Sample.Relationship.forward_abutted if is_ordered else Sample.Relationship.reverse_abutted

        # do samples overlap?
        elif ordered_overlaps(s1_ord, s2_ord):
            return Sample.Relationship.forward_overlap if is_ordered else Sample.Relationship.reverse_overlap

        # if we got this far there should be a gap between s1_ord and s2_ord
        elif ordered_gapped(s1_ord, s2_ord):
            return Sample.Relationship.forward_gapped if is_ordered else Sample.Relationship.reverse_gapped

        else:
            raise RuntimeError('Something went wrong checking the relative position of {} and {}'.format(s1.name, s2.name))


    @staticmethod
    def overlap_indices(s1, s2):
        """Find indices to trim end of s1 which overlaps start of s2 to allow concatenation without gaps or overlap,::

                       end1
                       |
            s1 ............
                       |
            s2      ...............
                       |
                       start2

        such that
        .. code: python

            Sample.from_samples([s1.slice(slice(0, end1)), s2.slice(slice(start2, None))])

        would join them without gaps or overlap.

        :param s1: First `Sample` object.
        :param s2: Second `Sample` object.

        :returns: (int end1: start of overlap in s1, int start2: end of overlap in s2,
            such that s1.positions[:end1] could be concatenated with
            s2.positions[start2:] to join them without gaps or overlap.
        :raises: `OverlapException` if samples do not overlap nor abut.

        """

        rel = Sample.relative_position(s1, s2)

        if rel is Sample.Relationship.forward_abutted:  # they can be concatenated witout slicing
            return None, None

        if rel is not Sample.Relationship.forward_overlap:
            msg = 'Cannot overlap samples {} and {} with relationhip {}'
            raise OverlapException(msg.format(s1.name, s2.name, repr(rel)))

        # find where the overlap starts in s1 indices
        ovl_start_ind1 = np.searchsorted(s1.positions, s2.positions[0])

        # find where the overlap ends in s2 indices
        ovl_end_ind2 = np.searchsorted(s2.positions, s1.positions[-1], side='right')
        pos1_ovl = s1.positions[ovl_start_ind1:]
        pos2_ovl = s2.positions[0:ovl_end_ind2]
        assert len(pos1_ovl) == len(pos2_ovl)  # check overlap length is the same in both
        overlap_len = len(pos1_ovl)
        # find mid-point (handling case where we have an odd number of overlap columns)
        pad_1 = overlap_len // 2
        pad_2 = overlap_len - pad_1
        end_1_ind = ovl_start_ind1 + pad_1
        start_2_ind = ovl_end_ind2 - pad_2

        # check the length of the overlap obtained using our indices is the same
        # in each sample
        contr_1 = s1.positions[ovl_start_ind1:end_1_ind]
        contr_2 = s2.positions[start_2_ind:ovl_end_ind2]
        assert len(contr_1) + len(contr_2) == overlap_len

        return end_1_ind, start_2_ind


    def chunks(self, chunk_len=1000, overlap=200):
        """Create overlapping chunks of self.

        :param chunk_len: chunk length (number of columns)
        :param overlap: overlap length.

        :yields: chunked `Sample` instances.
        """
        # TODO - could refactor this to use Sample.slice?
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


    def slice(self, key):
        """Slice fields along the genomic axis

        :param key: slice or index (on array indices)
        :returns: `Sample` obj with views of slices of the original `Sample`.

        >>> pos = np.array([(0, 0), (1, 0), (1, 1), (2, 0)], dtype=[('major', int), ('minor', int)])
        >>> feat = np.arange(len(pos))
        >>> s = Sample('contig1', feat , None, None, pos, None)
        >>> s.slice(2)
        Sample(ref_name='contig1', features=2, labels=None, ref_seq=None, positions=(1, 1), label_probs=None)
        >>> s.slice(slice(1,3))
        Sample(ref_name='contig1', features=array([1, 2]), labels=None, ref_seq=None, positions=array([(1, 0), (1, 1)], dtype=[('major', '<i8'), ('minor', '<i8')]), label_probs=None)
        """
        non_slice_fields = {'ref_name'}
        def slice_attr(attr):
            a = getattr(self, attr)
            if attr not in non_slice_fields:
                a = a[key] if a is not None else None
            return a
        return Sample(**{attr: slice_attr(attr) for attr in self._fields})


    def __eq__(self, other):
        for field in self._fields:
            s = getattr(self, field)
            o = getattr(other, field)
            if type(s) != type(o):
                return False
            elif isinstance(s, np.ndarray):
                if (s.shape != o.shape or np.any(s != o)):
                    return False
            elif s != o:
                return False
        return True


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
        """Parse region string into `Region` objects.

        :param region: region str

        >>> Region.from_string('Ecoli') == Region(ref_name='Ecoli', start=None, end=None)
        True
        >>> Region.from_string('Ecoli:1000-2000') == Region(ref_name='Ecoli', start=1000, end=2000)
        True
        >>> Region.from_string('Ecoli:1000') == Region(ref_name='Ecoli', start=0, end=1000)
        True
        >>> Region.from_string('Ecoli:-1000') == Region(ref_name='Ecoli', start=0, end=1000)
        True
        >>> Region.from_string('Ecoli:500-') == Region(ref_name='Ecoli', start=500, end=None)
        True
        """
        if ':' not in region:
            ref_name, start, end = region, None, None
        else:
            start, end = None, None
            ref_name, bounds = region.split(':')
            if bounds[0] == '-' or '-' not in bounds:
                start = 0
                end = int(bounds.replace('-', ''))
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


def ref_name_from_region_str(region_str):
    """Parse region strings, returning a list of reference names.

    :param regions: iterable of region strings.

    :returns: tuple of reference name str.
    """
    ref_names = [Region.from_string(r).ref_name for r in region_str]
    return tuple(set(ref_names))


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


def loose_version_sort(it, key=None):
    """Try to sort iterable with distutils.version.LooseVersion, falling back on regular sort.

    >>> loose_version_sort(['chr10', 'chr2', 'chr1'])
    ['chr1', 'chr2', 'chr10']
    >>> sorted(['chr10', 'chr2', 'chr1'])
    ['chr1', 'chr10', 'chr2']
    >>> loose_version_sort(['chr{}c{}'.format(i,j) for i,j in itertools.product([1, 10, 2] , [1,10,2])])
    ['chr1c1', 'chr1c2', 'chr1c10', 'chr2c1', 'chr2c2', 'chr2c10', 'chr10c1', 'chr10c2', 'chr10c10']
    """
    def version_sorter(x):
        return LooseVersion(x) if key is None else LooseVersion(key(x))
    it = list(it)
    try:
        result = sorted(it, key=version_sorter)
    except:
        logger = get_named_logger("VariantSort")
        logger.debug("Could not sort with LooseVersion")
        result = sorted(it, key=key)
    return result


comp = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'X': 'X', 'N': 'N',
    'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'x': 'x', 'n': 'n',
    #'-': '-'
}
comp_trans = str.maketrans(''.join(comp.keys()), ''.join(comp.values()))


def reverse_complement(seq):
    """Reverse complement sequence.

    :param: input sequence string.

    :returns: reverse-complemented string.

    """
    return seq.translate(comp_trans)[::-1]


if __name__ == "__main__":
    import doctest
    doctest.testmod()
