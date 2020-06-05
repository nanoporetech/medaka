"""Commonly used data structures and functions."""
import collections
from distutils.version import LooseVersion
import enum
import errno
import functools
import itertools
import logging
import os
import re

import numpy as np
from pkg_resources import resource_filename
import pysam

import libmedaka

ComprAlignPos = collections.namedtuple(
    'ComprAlignPos',
    ('qpos', 'qbase', 'qlen', 'rpos', 'rbase', 'rlen'))


class OverlapException(Exception):
    """Exception class used when examining range overlaps."""

    pass


class Relationship(enum.Enum):
    """Enumeration of types of overlap."""

    different_ref_name = 'Samples come from different reference contigs.'
    forward_overlap = 'The end of s1 overlaps the start of s2.'
    reverse_overlap = 'The end of s2 overlaps the start of s1.'
    forward_abutted = 'The end of s1 abuts the start of s2.'
    reverse_abutted = 'The end of s2 abuts the start of s1.'
    forward_gapped = 's2 follows s1 with a gab inbetween.'
    reverse_gapped = 's1 follows s2 with a gab inbetween.'
    s2_within_s1 = 's2 is fully contained within s1.'
    s1_within_s2 = 's1 is fully contained within s2.'


# provide read only access to key sample attrs
_Sample = collections.namedtuple(
    'Sample',
    ['ref_name', 'features', 'labels', 'ref_seq', 'positions', 'label_probs'])


class Sample(_Sample):
    """Represents a pileup range."""

    def _asdict(self):
        # https://bugs.python.org/issue24931
        return collections.OrderedDict(zip(self._fields, self))

    def amend(self, **kwargs):
        """Create new `Sample` with some attributes changed."""
        d = self._asdict()
        for k, v in kwargs.items():
            if k not in self._fields:
                raise KeyError('Invalid key for Sample: {}'.format(k))
            d[k] = v
        return Sample(**d)

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
        """Is pileup empty, synomymous to `sample.size == 0`."""
        return self.size == 0

    @property
    def size(self):
        """Return number of columns of pileup."""
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
        name_decoder = re.compile(
            r"(?P<ref_name>.+):(?P<start>\d+\.\d+)-(?P<end>\d+\.\d+)")
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
            if rel is not Relationship.forward_abutted:
                msg = (
                    'Refusing to concatenate unordered/non-abutting '
                    'samples {} and {} with relationship {}.')
                raise ValueError(msg.format(s1.name, s2.name, repr(rel)))

        # Relationship.forward_abutted guarantees all samples have the
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

    @staticmethod
    def relative_position(s1, s2):
        """Check the relative position of two samples in genomic coordinates.

        :param s1, s2: `medaka.common.Sample` objs.

        :returns: `Relationship` enum member.
        """
        def ordered_abuts(s1, s2):
            """Check if end of s1 abuts the start of s2.

            i.e. is adjacent but not overlapping
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
            """Check if s2 is contained within s1."""
            return s2.first_pos >= s1.first_pos and s2.last_pos <= s1.last_pos

        def ordered_overlaps(s1, s2):
            """Check if end of s1 overlaps start of s2."""
            s1_end_maj, s1_end_min = s1.last_pos
            s2_start_maj, s2_start_min = s2.first_pos
            if s2_start_maj < s1_end_maj:  # overlap of major coordinates
                overlaps = True
            elif s2_start_maj == s1_end_maj and s2_start_min < s1_end_min + 1:
                # we have overlap of minor coordinates
                overlaps = True
            else:
                overlaps = False
            return overlaps

        def ordered_gapped(s1, s2):
            """Check for grap between end of s1 and start of s2."""
            s1_end_maj, s1_end_min = s1.last_pos
            s2_start_maj, s2_start_min = s2.first_pos
            gapped = False
            if s2_start_maj > s1_end_maj + 1:  # gap in major
                gapped = True
            elif (s2_start_maj > s1_end_maj and
                    s2_start_min > 0):  # missing minors
                gapped = True
            elif (s2_start_maj == s1_end_maj and
                    s2_start_min > s1_end_min + 1):  # missing minors
                gapped = True
            return gapped

        if s1.ref_name != s2.ref_name:  # different ref_names
            return Relationship.different_ref_name

        s1_ord, s2_ord = sorted((s1, s2), key=lambda x: (x.first_pos, -x.size))
        is_ordered = s1_ord.name == s1.name

        # is one sample within the other?
        if ordered_contained(s1_ord, s2_ord):
            if is_ordered:
                return Relationship.s2_within_s1
            else:
                return Relationship.s1_within_s2

        # do samples abut?
        elif ordered_abuts(s1_ord, s2_ord):
            if is_ordered:
                return Relationship.forward_abutted
            else:
                return Relationship.reverse_abutted

        # do samples overlap?
        elif ordered_overlaps(s1_ord, s2_ord):
            if is_ordered:
                return Relationship.forward_overlap
            else:
                return Relationship.reverse_overlap

        # if we got this far there should be a gap between s1_ord and s2_ord
        elif ordered_gapped(s1_ord, s2_ord):
            if is_ordered:
                return Relationship.forward_gapped
            else:
                return Relationship.reverse_gapped

        else:
            raise RuntimeError(
                'Could not calculate relative position of {} and {}'.format(
                    s1.name, s2.name))

    @staticmethod
    def overlap_indices(s1, s2):
        """Find indices by which to trim samples to allow concatenation.

        ::
            #          | end1
            s1 ............
            s2      ...............
            #          | start2


        For example:

        .. code: python

            Sample.from_samples([
                s1.slice(slice(0, end1)),
                s2.slice(slice(start2, None))])

        would join the samples without gaps or overlap.

        :param s1: First `Sample` object.
        :param s2: Second `Sample` object.

        :returns: (end1, start2
        :raises: `OverlapException` if samples do not overlap nor abut.

        """
        heuristic = False
        rel = Sample.relative_position(s1, s2)

        # trivial case
        if rel is Relationship.forward_abutted:
            return None, None, heuristic

        if rel is not Relationship.forward_overlap:
            msg = 'Cannot overlap samples {} and {} with relationhip {}'
            raise OverlapException(msg.format(s1.name, s2.name, repr(rel)))

        # find where the overlap starts (ends) in s1 (s2) indices
        ovl_start_ind1 = np.searchsorted(s1.positions, s2.positions[0])
        ovl_end_ind2 = np.searchsorted(
            s2.positions, s1.positions[-1], side='right')

        end_1_ind, start_2_ind = None, None
        pos1_ovl = s1.positions[ovl_start_ind1:]
        pos2_ovl = s2.positions[0:ovl_end_ind2]
        try:
            # the nice case where everything lines up
            if not np.array_equal(pos1_ovl['minor'], pos2_ovl['minor']):
                raise OverlapException("Overlaps are not equal in structure")
            overlap_len = len(pos1_ovl)
            # take mid point as break point
            pad_1 = overlap_len // 2
            pad_2 = overlap_len - pad_1
            end_1_ind = ovl_start_ind1 + pad_1
            start_2_ind = ovl_end_ind2 - pad_2

            contr_1 = s1.positions[ovl_start_ind1:end_1_ind]
            contr_2 = s2.positions[start_2_ind:ovl_end_ind2]
            if len(contr_1) + len(contr_2) != overlap_len:
                raise OverlapException(
                    "Resultant is not same length as overlap.")
        except OverlapException:
            heuristic = True
            # Some sample producing methods will not create 1-to-1 mappings
            # in their sets of columns, e.g. where chunking has affected the
            # reads used. Here we find a split point near the middle where
            # the two halfs have the same number of minor positions
            # (i.e. look similar).
            # Require seeing a number of major positions
            UNIQ_MAJ = 3
            end_1_ind, start_2_ind = None, None
            if (len(np.unique(pos1_ovl['major'])) > UNIQ_MAJ and
                    len(np.unique(pos2_ovl['major'])) > UNIQ_MAJ):

                start, end = pos1_ovl['major'][0], pos1_ovl['major'][-1]
                mid = start + (end - start) // 2
                offset = 1
                while end_1_ind is None:

                    if (mid + offset > max(s1.positions['major']) and
                            mid - offset < min(s2.positions['major'])):
                        # run off the edge
                        break
                    for test in (+offset, -offset):
                        left = np.where(
                            s1.positions['major'] == mid + test)[0]
                        right = np.where(
                            s2.positions['major'] == mid + test)[0]
                        if len(left) == len(right):
                            # found a nice junction
                            end_1_ind = left[0]
                            start_2_ind = right[0]
                            break
                    offset += 1
            if end_1_ind is None or start_2_ind is None:
                raise OverlapException(
                    "Could not find viable junction for {} and {}".format(
                        s1.name, s2.name))

        return end_1_ind, start_2_ind, heuristic

    def chunks(self, chunk_len=1000, overlap=200):
        """Create overlapping chunks of self.

        :param chunk_len: chunk length (number of columns)
        :param overlap: overlap length.

        :yields: chunked `Sample` instances.
        """
        # TODO - could refactor this to use Sample.slice?
        chunker = functools.partial(
            sliding_window,
            window=chunk_len, step=chunk_len - overlap, axis=0)
        chunkers = {
            field: chunker(getattr(self, field))
            if getattr(self, field) is not None else itertools.repeat(None)
            for field in self._fields
        }

        for pos in chunkers['positions']:
            fields = set(self._fields) - set(['positions', 'ref_name'])
            new_sample = {
                'positions': pos, 'ref_name': self.ref_name}
            for field in fields:
                new_sample[field] = next(chunkers[field])
            yield Sample(**new_sample)

    def slice(self, key):
        """Slice fields along the genomic axis.

        :param key: slice or index (on array indices)
        :returns: `Sample` obj with views of slices of the original `Sample`.

        >>> pos = np.array(
        ...     [(0, 0), (1, 0), (1, 1), (2, 0)],
        ...     dtype=[('major', int), ('minor', int)])
        >>> feat = np.arange(len(pos))
        >>> s = Sample('contig1', feat , None, None, pos, None)
        >>> s.slice(2)  #doctest: +ELLIPSIS
        Sample(...features=2, ..., positions=(1, 1), label_probs=None)
        >>> s.slice(slice(1,3)) #doctest: +ELLIPSIS
        Sample(..., features=array([1, 2]),...)
        """
        non_slice_fields = {'ref_name'}

        def slice_attr(attr):
            a = getattr(self, attr)
            if attr not in non_slice_fields:
                a = a[key] if a is not None else None
            return a
        return Sample(**{attr: slice_attr(attr) for attr in self._fields})

    def __eq__(self, other):
        """Test equality."""
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


# provide read only access to key region attrs
_Region = collections.namedtuple('Region', 'ref_name start end')


class Region(_Region):
    """Represents a genomic region."""

    @property
    def name(self):
        """Samtools-style region string, zero-base end exclusive."""
        return self.__str__()

    def __str__(self):
        """Return string representation of region."""
        # This will be zero-based, end exclusive
        start = 0 if self.start is None else self.start
        end = '' if self.end is None else self.end
        return '{}:{}-{}'.format(self.ref_name, start, end)

    @property
    def size(self):
        """Return size of region."""
        return self.end - self.start

    @classmethod
    def from_string(cls, region):
        """Parse region string into `Region` objects.

        :param region: region str

        >>> Region.from_string('Ecoli') == Region(
        ...     ref_name='Ecoli', start=None, end=None)
        True
        >>> Region.from_string('Ecoli:1000-2000') == Region(
        ...     ref_name='Ecoli', start=1000, end=2000)
        True
        >>> Region.from_string('Ecoli:1000') == Region(
        ...     ref_name='Ecoli', start=1000, end=None)
        True
        >>> Region.from_string('Ecoli:-1000') == Region(
        ...     ref_name='Ecoli', start=0, end=1000)
        True
        >>> Region.from_string('Ecoli:500-') == Region(
        ...     ref_name='Ecoli', start=500, end=None)
        True
        >>> Region.from_string('A:B:c:500-') == Region(
        ...     ref_name='A:B:c', start=500, end=None)
        True
        """
        if ':' not in region:
            ref_name, start, end = region, None, None
        else:
            start, end = None, None
            ref_name, bounds = region.rsplit(':', 1)
            if bounds[0] == '-':
                start = 0
                end = int(bounds.replace('-', ''))
            elif '-' not in bounds:
                start = int(bounds)
                end = None
            elif bounds[-1] == '-':
                start = int(bounds[:-1])
                end = None
            else:
                start, end = [int(b) for b in bounds.split('-')]
        return cls(ref_name, start, end)

    def split(region, size, overlap=0, fixed_size=True):
        """Split region into sub-regions of a given length.

        :param size: size of sub-regions.
        :param overlap: overlap between ends of sub-regions.
        :param fixed_size: ensure all sub-regions are equal in size. If `False`
            then the final chunk will be created as the smallest size to
            conform with `overlap`.

        :returns: a list of sub-regions.

        """
        regions = list()
        for start in range(region.start, region.end, size - overlap):
            end = min(start + size, region.end)
            regions.append(Region(region.ref_name, start, end))
        if len(regions) > 1:
            if fixed_size and regions[-1].size < size:
                del regions[-1]
                end = region.end
                start = end - size
                if start > regions[-1].start:
                    regions.append(Region(region.ref_name, start, end))
        return regions

    def overlaps(self, other):
        """Determine if a region overlaps another.

        :param other: a second Region to test overlap.

        :returns: True if regions overlap.

        """
        if self.ref_name != other.ref_name:
            return False

        def _limits(x):
            x0 = x.start if x.start is not None else -1
            x1 = x.end if x.end is not None else float('inf')
            return x0, x1

        a0, a1 = _limits(self)
        b0, b1 = _limits(other)
        return (
            (a0 < b1 and a1 > b0) or
            (b0 < a1 and b1 > a0))


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
                region_strs = [line.strip() for line in fh.readlines()]

        regions = []
        for r in (Region.from_string(x) for x in region_strs):
            start = r.start if r.start is not None else 0
            end = r.end if r.end is not None else ref_lengths[r.ref_name]
            regions.append(Region(r.ref_name, start, end))
    else:
        regions = [
            Region(ref_name, 0, end)
            for ref_name, end in ref_lengths.items()]

    return regions


def ref_name_from_region_str(region_str):
    """Parse region strings, returning a list of reference names.

    :param regions: iterable of region strings.

    :returns: tuple of reference name str.
    """
    ref_names = [Region.from_string(r).ref_name for r in region_str]
    return tuple(set(ref_names))


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
    """Group together elements of an iterable without padding remainder."""
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


def roundrobin(*iterables):
    """Take items from iterables in a round-robin."""
    pending = len(iterables)
    nexts = itertools.cycle(iter(it).__next__ for it in iterables)
    while pending:
        try:
            for next in nexts:
                yield next()
        except StopIteration:
            pending -= 1
            nexts = itertools.cycle(itertools.islice(nexts, pending))


def print_data_path():
    """Print data directory containing models."""
    print(resource_filename(__package__, 'data'))


def get_named_logger(name):
    """Create a logger with a name."""
    logger = logging.getLogger('{}.{}'.format(__package__, name))
    logger.name = name
    return logger


def loose_version_sort(it, key=None):
    """Sort an iterable.

    Items will be sorted with `distutils.version.LooseVersion`,
    falling back to regular sort when this fails.

    >>> loose_version_sort(['chr10', 'chr2', 'chr1'])
    ['chr1', 'chr2', 'chr10']
    >>> sorted(['chr10', 'chr2', 'chr1'])
    ['chr1', 'chr10', 'chr2']
    >>> loose_version_sort([
    ...     'chr{}c{}'.format(i,j)
    ...     for i,j in itertools.product(
    ...        [1, 10, 2] , [1,10,2])])  # doctest: +ELLIPSIS
    ['chr1c1', 'chr1c2', 'chr1c10', 'chr2c1', ..., 'chr10c2', 'chr10c10']
    """
    def version_sorter(x):
        return LooseVersion(x) if key is None else LooseVersion(key(x))
    it = list(it)
    try:
        result = sorted(it, key=version_sorter)
    except Exception:
        logger = get_named_logger("VariantSort")
        logger.debug("Could not sort with LooseVersion")
        result = sorted(it, key=key)
    return result


comp = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'X': 'X', 'N': 'N',
    'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'x': 'x', 'n': 'n',
    # '-': '-'
}
comp_trans = str.maketrans(''.join(comp.keys()), ''.join(comp.values()))


def reverse_complement(seq):
    """Reverse complement sequence.

    :param: input sequence string.

    :returns: reverse-complemented string.

    """
    return seq.translate(comp_trans)[::-1]


def read_key_value_tsv(fname):
    """Read a dictionary from a .tsv file.

    :param fname: a .tsv file with two columns: keys and values.

    :returns: a dictionary.

    """
    try:
        as_string = libmedaka.ffi.string
        text = libmedaka.lib.read_key_value(fname.encode())
        data = dict()
        for i in range(0, text.n, 2):
            key = as_string(text.strings[i]).decode()
            value = as_string(text.strings[i+1]).decode()
            data[key] = value
        libmedaka.lib.destroy_string_set(text)
    except Exception:
        raise RuntimeError(
            'Failed to parse {} as two-column .tsv file.'.format(fname))
    return data


def initialise_alignment(
        query_name, reference_id, reference_start,
        query_sequence, cigarstring, flag, mapping_quality=60,
        query_qualities=None, tags=None):
    """Create a `Pysam.AlignedSegment` object.

    :param query_name: name of the query sequence
    :param reference_id: index to the reference name
    :param reference_start: 0-based index of first leftmost reference
        coordinate
    :param query_sequence: read sequence bases, including those soft clipped
    :param cigarstring: cigar string representing the alignment of query
        and reference
    :param flag: bitwise flag representing some properties of the alignment
        (see SAM format)
    :param mapping_quality: optional quality of the mapping or query to
        reference
    :param query_qualities: optional base qualities of the query, including
        soft-clipped ones!

    :returns: `pysam.AlignedSegment` object
    """
    if tags is None:
        tags = dict()

    a = pysam.AlignedSegment()
    a.query_name = query_name
    a.reference_id = reference_id
    a.reference_start = reference_start
    a.query_sequence = query_sequence
    a.cigarstring = cigarstring
    a.flag = flag
    a.mapping_quality = mapping_quality
    if query_qualities is not None:
        a.query_qualities = query_qualities

    for tag_name, tag_value in tags.items():
        a.set_tag(tag_name, tag_value)

    return a


def yield_from_bed(bedfile):
    """Yield chrom, start, stop tuples from a bed file.

    :param bedfile: str, filepath.
    :yields: (str chrom, int start, int stop).

    """
    with open(bedfile) as fh:
        for line in fh:
            split_line = line.split()
            if split_line[0] in {'browser', 'track'} or len(split_line) < 3:
                continue
            chrom = split_line[0]
            start = int(split_line[1])
            stop = int(split_line[2])
            yield chrom, start, stop
