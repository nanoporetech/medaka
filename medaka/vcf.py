from copy import deepcopy
from collections import defaultdict, OrderedDict
import itertools
from threading import Lock

from intervaltree import IntervalTree

from medaka.common import get_named_logger


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
    _format_fields_ = ('gt', 'gq',)

    def __init__(self, chrom, pos, ref, alt='.', id='.', qual='.', filter='.', info='.', sample_dict=None):
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
        self.sample_dict = sample_dict if sample_dict is not None else OrderedDict()


    def __eq__(self, other):
        for field in ('chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info'):
            if getattr(self, field) != getattr(other, field):
                return False
        return True


    def __ne__(self, other):
        return not self.__eq__(other)


    @property
    def format(self):
        return ':'.join((str(v) for v in self.sample_dict.keys()))


    @property
    def sample(self):
        return ':'.join((str(v) for v in self.sample_dict.values()))


    @property
    def info_string(self):
        return parse_tags_to_string(self.info)


    @classmethod
    def from_text(cls, line):
        chrom, pos, id, ref, alt, qual, filter, info, sample_fields, sample_data, *others = line.split('\t')
        pos = int(pos)
        pos -= 1 # VCF is 1-based, python 0-based
        sample_dict = OrderedDict(zip(sample_fields.split(':'), sample_data.split(':')))
        return cls(chrom, pos, ref, alt=alt, id=id, qual=qual, filter=filter, info=info, sample_dict=sample_dict)


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
        attributes['sample_repr'] = ';'.join('{}={}'.format(k,v) for k,v in self.sample_dict.items())
        return ("Variant('{chrom}', {pos}, '{ref}', alt={alt}, id={id}, qual={qual},"
                " filter={filter}, info='{info_string}', sample='{sample_repr}')".format(**attributes))


    def deep_copy(self):
        return deepcopy(self)


class VCFWriter(object):
    # some tools don't like VCFv4.3, preferring VCFv4.1 - so we should be able to
    # write VCFv4.1 files. VCFv4.3 has a few extra reserved fields ('AD', 'ADF', and
    # 'ADR') but there is no harm in including those files written as VCFv4.1 - they
    # just won't be recognised and used as reserved fields.
    version_options = {'4.3', '4.1'}
    def __init__(self, filename, mode='w',
                 header=('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'),
                 meta_info=[],
                 version='4.3'
                 ):

        self.filename = filename
        self.mode = mode
        self.header = header
        if version not in self.version_options:
            raise ValueError('version must be one of {}'.format(self.version_options))
        self.version = version
        self.meta = ['fileformat=VCFv{}'.format(self.version)] + meta_info
        self.logger = get_named_logger('VCFWriter')


    def __enter__(self):
        self.handle = open(self.filename, self.mode, encoding='utf-8')
        self.handle.write('\n'.join('##' + line for line in self.meta) + '\n')
        self.handle.write('#' + '\t'.join(self.header) + '\n')
        return self


    def __exit__(self, exc_type, exc_val, exc_tb):
        self.handle.close()


    def write_variant(self, variant):
        variant = variant.deep_copy()
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
        self.logger = get_named_logger('VCFReader')

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

        def _tree_search(tree, start, end, strict):
            return tree.overlap(start, end) if strict else tree.envelop(start, end)

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
                results = _tree_search(self._tree[ref_name], start, end, strict)
            else:
                results = itertools.chain(*(
                     _tree_search(x, start, end, strict=True)
                     for x in self._tree.values()
                ))
            # spec says .vcf is sorted, lets follow
            for interval in sorted(results):
                yield interval.data
