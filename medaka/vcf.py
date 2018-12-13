from copy import deepcopy
from collections import defaultdict
import itertools
from threading import Lock
import logging

from intervaltree import IntervalTree

logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
logger = logging.getLogger(__name__)


def self_return(x):
    return x


# Source: Table1 in 'The Variant Call Format Specification VCFv4.3', Table 1
# Tuples below are (number, type), where number can be:
#   A: The field has one value per alternate allele
#   R: The field has one value for each possible allele, including the reference
#   G: The field has one value for each possible genotype
#   .(dot): The number of possible values varies, is unknown or unbounded
reserved_info_fields = {'AA': (1, str), 'AC': ('A', int), 'AD': ('R', int), 'ADF': ('R', int),
                        'ADR': ('R', int), 'AF': ('A', float), 'AN': (1, int), 'BQ': (1, float),
                        'CIGAR': ('A', str), 'DB': (0, self_return), 'DP': (1, int), 'END': (1, int),
                        'H2': (0, self_return), 'H3': (0, self_return), 'MQ': (1, self_return),
                        'MQ0': (1, int), 'NS': (1, int), 'SB': ('.', self_return),
                        'SOMATIC': (0, self_return), 'VALIDATED': (0, self_return), '1000G': (0, self_return)
                        }
own_info_fields = {'SCORES': ('R', float)}
all_info_fields = reserved_info_fields.copy()
all_info_fields.update(own_info_fields)

def parse_tags_to_string(tags):
    str_tags = []
    for key, value in tags.items():
        # If key is of type 'Flag', print only key, else 'key=value'
        if value is True:
            str_tags.append(key)
        else:
            if isinstance(value, (tuple, list)):
                value = ','.join((str(x) for x in value))
            str_tags.append('{}={}'.format(key, value))
    return ';'.join(str_tags)


def parse_string_to_tags(string, splitter=','):
    tags = {}
    for field in string.split(';'):
        try:
            tag, value = field.split('=')
            if tag in all_info_fields.keys():
                _type = all_info_fields[tag][1]
                value = [_type(x) for x in value.split(splitter)]
                if len(value) == 1:
                    value = value[0]
        except ValueError:
            tag = field
            value = True

        tags[tag] = value
    return tags


class Variant(object):
    # TODO: ref/alt could be a symbolic allele "<ID>".
    # TODO: alt could contain breakends.
    # TODO: Handle genomic fields.

    def __init__(self, chrom, pos, ref, alt='.', id='.', qual='.', filter='.', info='.'):
        self.chrom = chrom
        self.pos = int(pos)
        self.ref = ref.upper()
        # self.alt should be a list/tuple of alternatives
        self.alt = alt.split(',') if isinstance(alt, str) else alt
        self.id = str(id)
        self.qual = float(qual) if qual != '.' else qual
        self.filter = filter.split(';') if ';' in filter else filter
        if isinstance(info, dict):
            self.info = info
        else:
            self.info = parse_string_to_tags(info)


    def __eq__(self, other):
        for field in ('chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info'):
            if getattr(self, field) != getattr(other, field):
                return False
        return True


    def __ne__(self, other):
        return not self.__eq__(other)


    @property
    def info_string(self):
        return parse_tags_to_string(self.info)


    @classmethod
    def from_text(cls, line):
        chrom, pos, id, ref, alt, qual, filter, info, *others = line.split('\t')
        pos = int(pos)
        pos -= 1 # VCF is 1-based, python 0-based
        return cls(chrom, pos, ref, alt=alt, id=id, qual=qual, filter=filter, info=info)


    def add_tag(self, tag, value=None):
        self.info[tag] = value

        # Remove default value if more than one exists
        if len(self.info.keys()) > 0:
            self.info.pop('.', None)


    def get_tag(self, tag):
        return self.info[tag]


    def __repr__(self):
        attributes = {}
        for field in ('chrom', 'pos', 'ref', 'alt', 'id', 'qual', 'filter', 'info_string'):
            attributes[field] = getattr(self, field)
        return ("Variant('{chrom}', {pos}, '{ref}', alt={alt}, id={id}, qual={qual},"
                " filter={filter}, info='{info_string}')".format(**attributes))


    def deep_copy(self):
        return deepcopy(self)


class VCFWriter(object):
    def __init__(self, filename, mode='w',
                 header=('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'),
                 meta_info=[]
                 ):

        self.filename = filename
        self.mode = mode
        self.header = header
        self.meta = ['fileformat=VCFv4.3'] + meta_info


    def __enter__(self):
        self.handle = open(self.filename, self.mode, encoding='utf-8')
        self.handle.write('\n'.join('##' + line for line in self.meta) + '\n')
        self.handle.write('#' + '\t'.join(self.header) + '\n')
        return self


    def __exit__(self, exc_type, exc_val, exc_tb):
        self.handle.close()


    def write_variant(self, variant):
        # Some fields can be multiple
        for attribute in ('alt', 'filter'):
            value = getattr(variant, attribute)
            if isinstance(value, (tuple, list)):
                setattr(variant, attribute, ','.join(str(x) for x in value))

        # Convert info dictionary to string
        variant.info = variant.info_string

        elements = [getattr(variant, field.lower()) for field in self.header]
        # VCF POS field is 1-based
        elements[self.header.index('POS')] += 1
        line = '\t'.join([str(x) for x in elements])
        self.handle.write('{}\n'.format(line))


class VCFReader(object):
    def __init__(self, filename, cache=True):
        """Basic VCF parser.

        :param filename: .vcf file.
        :param cache: if True, all parsed variants are stored in memory for
            faster subsequent access.

        """

        self.filename = filename
        self.cache = cache
        self._indexed = False
        self._tree = None
        self._parse_lock = Lock()

        # Read both metadata and header
        self.meta = []
        self.header = None
        with open(filename, encoding='utf-8') as handle:
            for line in handle:
                line = line.replace('\n', '')
                if line.startswith('##'):
                    self.meta.append(line[2:])
                elif line.startswith('#'):
                    line = line[1:]
                    self.header = line.split('\t')
                    break


    def _parse(self):
        # just parsing the file to yield records
        last_pos = [None, None]
        with open(self.filename, encoding='utf-8') as handle:
            for index, line in enumerate(handle):
                line = line.replace('\n', '')

                # Already read meta and header in self.__init__
                if line.startswith('#'):
                    continue

                try:
                    variant = Variant.from_text(line)
                except Exception as e:
                    raise IOError(
                        'Exception while reading variant #{}.\n'
                        'Line: {}'.format(index, line)) from e

                if variant.chrom != last_pos[0]:
                    last_pos = [variant.chrom, None]
                elif last_pos[1] is not None and last_pos[1] > variant.pos:
                    raise IOError('.vcf is unsorted at index #{}.'.format(index))
                yield variant
                last_pos[1] = variant.pos


    def index(self):
        """Index the input file for faster fetches."""

        # calling this method implies caching
        self.cache = True
        if self._indexed or not self.cache:
            return

        if self._parse_lock.acquire(blocking=False):
            try:
                # clear out an incomplete parse, actually this doesn't matter since
                #    the values in the tree are set-like.
                self._tree = defaultdict(IntervalTree)
                for variant in self._parse():
                    self._tree[variant.chrom][variant.pos:variant.pos + len(variant.ref)] = variant
            except Exception:
                raise
            else:
                # record we've done a complete parse
                self._indexed = True
            finally:
                self._parse_lock.release()

        else:
            # wait for lock to be released, then return
            self._parse_lock.acquire(blocking=True)
            if not self._indexed:
                raise IOError("Waited for parsing, but parsing did not occur.")


    def fetch(self, ref_name=None, start=None, end=None, strict=True):
        """Yield all variants spanned by a region.

        :param ref_name: reference name (CHROM field of .vcf).
        :param start: inclusive start co-ordinate (0-based).
        :param end: exclusive end co-ordinate (0-based).
        :param strict: if False variants overlapping the region, but not
            contained enitrely within the region are yielded also.

        :yields: :py:class:`Variant` instances.

        """
        if start is None:
            start = float('-inf')
        if end is None:
            end = float('inf')

        if not self.cache:
            # if not using a cache, just keep re-reading the file
            for variant in self._parse():
                if not all([ref_name is None or variant.chrom == ref_name,
                            start is None or variant.pos > start,
                            end is None or variant.pos + len(variant.ref) < end]):
                    continue
                yield variant
        else:
            self.index()
            if ref_name is not None:
                results = self._tree[ref_name].search(start, end, strict=strict)
            else:
                results = itertools.chain(*(
                     x.search(start, end, strict=True)
                     for x in self._tree.values()
                ))
            # spec says .vcf is sorted, lets follow
            for interval in sorted(results):
                yield interval.data


