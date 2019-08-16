import collections
import contextlib
from copy import deepcopy
import itertools
import os
from threading import Lock

import intervaltree
import numpy as np
import pysam

import medaka.common


def self_return(x):
    return x


# Source: Table1 in 'The Variant Call Format Specification VCFv4.3', Table 1
# Tuples below are (number, type), where number can be:
#   A: The field has one value per alternate allele
#   R: The field has one value for each possible allele,
#      including the reference
#   G: The field has one value for each possible genotype
#   .(dot): The number of possible values varies, is unknown or unbounded
reserved_info_fields = {
    'AA': (1, str), 'AC': ('A', int), 'AD': ('R', int), 'ADF': ('R', int),
    'ADR': ('R', int), 'AF': ('A', float), 'AN': (1, int), 'BQ': (1, float),
    'CIGAR': ('A', str), 'DB': (0, self_return), 'DP': (1, int),
    'END': (1, int), 'H2': (0, self_return), 'H3': (0, self_return),
    'MQ': (1, self_return), 'MQ0': (1, int), 'NS': (1, int),
    'SB': ('.', self_return), 'SOMATIC': (0, self_return),
    'VALIDATED': (0, self_return), '1000G': (0, self_return)}
own_info_fields = {'SCORES': ('R', float)}
all_info_fields = reserved_info_fields.copy()
all_info_fields.update(own_info_fields)


def parse_tags_to_string(tags):
    str_tags = []
    for key, value in sorted(tags.items()):
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


class MetaInfo(object):

    __valid_groups__ = ('INFO', 'FILTER', 'FORMAT')
    __valid_group_sort__ = {v: k for k, v in enumerate(__valid_groups__)}
    __valid_non_int_nums__ = {'A', 'R', 'G', '.'}
    __valid_types__ = {'Integer', 'Float', 'Flag', 'Character', 'String'}

    def __init__(self, group, id, number, typ, descr):
        """Class for representing meta info for VCF header

        :param group: str, one of {'INFO', 'FILTER', 'FORMAT'}
        :param id: str, short name as it occurs in a VCF data line.
        :param number: int or one of {'A', 'R', 'G', '.'}.
        :param type: one of {'Integer', 'Float', 'Flag', 'Character', 'String'}
        :param descr: str, free form description.
        """
        if group not in self.__valid_groups__:
            raise ValueError(
                'Group {} is not one of {}'.format(
                    group, self.__valid_groups__))

        if not isinstance(number, int) \
                and not (isinstance(number, str) and number.isdigit()) \
                and number not in self.__valid_non_int_nums__:
            raise ValueError(
                'Number {} is not an int, digit str or one of {}'.format(
                    number, self.__valid_non_int_nums__))

        if typ not in self.__valid_types__:
            raise ValueError('typ {} is not one of {}'.format(
                typ, self.__valid_types__))

        self.group = group
        self.id = id
        self.number = number
        self.typ = typ
        self.descr = descr

    def __repr__(self):
        return '{}=<ID={},Number={},Type={},Description="{}">'.format(
            self.group, self.id, self.number, self.typ, self.descr)

    def __str__(self):
        return self.__repr__()


class Variant(object):
    # TODO: ref/alt could be a symbolic allele "<ID>".
    # TODO: alt could contain breakends.
    # TODO: Handle genomic fields.

    def __init__(
            self, chrom, pos, ref,
            alt='.', id='.', qual='.', filter='.', info='.',
            sample_dict=None):
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
        if sample_dict is not None:
            self.sample_dict = sample_dict
        else:
            self.sample_dict = collections.OrderedDict()

    def __eq__(self, other):
        for field in (
                'chrom', 'pos', 'id', 'ref', 'alt',
                'qual', 'filter', 'info', 'sample_dict'):
            if getattr(self, field) != getattr(other, field):
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    @property
    def _sorted_format_keys(self):
        start_with = ['GT', 'GQ']
        sorted_keys = [k for k in start_with if k in self.sample_dict.keys()]
        sorted_keys.extend([
            k for k in sorted(self.sample_dict.keys())
            if k not in sorted_keys])
        return sorted_keys

    @property
    def format(self):
        return ':'.join((str(v) for v in self._sorted_format_keys))

    @property
    def sample(self):
        return ':'.join((
            str(self.sample_dict[k]) for k in self._sorted_format_keys))

    @property
    def info_string(self):
        return parse_tags_to_string(self.info)

    @property
    def gt(self):
        if 'GT' in self.sample_dict:
            gts = self.sample_dict['GT'].replace('|', '/').split('/')
            return tuple(int(x) for x in gts)
        else:
            return None

    @property
    def alleles(self):
        all_alleles = [self.ref] + self.alt
        if self.gt is None:
            return None
        else:
            return tuple([all_alleles[i] for i in self.gt])

    @classmethod
    def from_text(cls, line):
        chrom, pos, ident, ref, alt, qual, filt, info, \
            sample_fields, sample_data, *others = line.split('\t')
        pos = int(pos)
        pos -= 1  # VCF is 1-based, python 0-based
        sample_dict = collections.OrderedDict(
            tuple(zip(sample_fields.split(':'), sample_data.split(':'))))
        instance = cls(
            chrom, pos, ref,
            alt=alt, id=ident, qual=qual, filter=filt, info=info,
            sample_dict=sample_dict)
        return instance

    def add_tag(self, tag, value=None):
        self.info[tag] = value

        # Remove default value if more than one exists
        if len(self.info.keys()) > 0:
            self.info.pop('.', None)

    def get_tag(self, tag):
        return self.info[tag]

    def __repr__(self):
        attributes = {}
        for field in (
                'chrom', 'pos', 'ref', 'alt', 'id',
                'qual', 'filter', 'info_string'):
            attributes[field] = getattr(self, field)
        attributes['sample_repr'] = ';'.join(
            '{}={}'.format(k, self.sample_dict[k])
            for k in self._sorted_format_keys)
        return (
            "Variant('{chrom}', {pos}, '{ref}', alt={alt}, id={id}, "
            "qual={qual}, filter={filter}, info='{info_string}', "
            "sample='{sample_repr}')".format(**attributes))

    def deep_copy(self):
        return deepcopy(self)

    def to_dict(self):
        """Return a dictionary representation."""
        d = dict(alt=','.join(self.alt))
        for attr in ['chrom', 'pos', 'qual', 'id', 'filter', 'ref']:
            d[attr] = getattr(self, attr)
        d.update(self.info)
        d.update(self.sample_dict)
        return d

    def trim(self):
        """Return new trimmed Variant with minimal ref and alt sequence."""
        def get_trimmed_start_ref_alt(seqs):
            def trim_start(seqs):
                min_len = min([len(s) for s in seqs])
                trim_start = 0
                for bases in zip(*seqs):
                    bases = list(bases)
                    bases_same = len(set(bases)) == 1
                    if not bases_same or trim_start == min_len - 1:
                        break
                    if bases_same:
                        trim_start += 1
                return trim_start, [s[trim_start:] for s in seqs]

            # trim ends
            rev_seqs = [s[::-1] for s in seqs]
            _, trimmed_rev_seqs = trim_start(rev_seqs)
            seqs = [s[::-1] for s in trimmed_rev_seqs]

            trim_start, seqs = trim_start(seqs)
            return trim_start, seqs

        trimmed = self.deep_copy()
        seqs = [trimmed.ref] + trimmed.alt
        trim_start, (ref, *alt) = get_trimmed_start_ref_alt(seqs)
        trimmed.pos += trim_start
        trimmed.ref = ref
        trimmed.alt = alt
        return trimmed

    def split_haplotypes(self):
        """Split multiploid variants into list of non-ref haploid variants.

        :returns: (haplotype number, `vcf.Variant` or None for ref allele)
        """
        if 'GT' not in self.sample_dict:
            return tuple()

        vs = []
        sample_dict = self.sample_dict.copy()
        sample_dict['GT'] = '1/1'
        for hap_n, n in enumerate(self.gt, 1):
            if n == 0:
                v = None
            else:
                v = Variant(
                    self.chrom, self.pos, self.ref, self.alt[n - 1],
                    qual=self.qual, info=self.info.copy(),
                    sample_dict=sample_dict)
            vs.append((hap_n, v))
        return tuple(vs)


class VCFWriter(object):
    # some tools don't like VCFv4.3, preferring VCFv4.1 - so we should be able
    # to write VCFv4.1 files. VCFv4.3 has a few extra reserved fields ('AD',
    # 'ADF', and 'ADR') but there is no harm in including those files written
    # as VCFv4.1 - they just won't be recognised and used as reserved fields.
    version_options = {'4.3', '4.1'}

    def __init__(self, filename, mode='w',
                 header=(
                     'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                     'FILTER', 'INFO', 'FORMAT', 'SAMPLE'),
                 contigs=None,
                 meta_info=None,
                 version='4.1'
                 ):

        self.filename = filename
        self.mode = mode
        self.header = header
        if version not in self.version_options:
            raise ValueError('version must be one of {}'.format(
                self.version_options))
        self.version = version
        self.meta = ['fileformat=VCFv{}'.format(self.version)]

        if contigs is not None:
            self.meta.extend(['contig=<ID={}>'.format(c) for c in contigs])

        if meta_info is not None:
            # try to sort so we get INFO, FILTER, FORMAT in that order
            try:
                meta_info.sort(
                    key=lambda x: MetaInfo.__valid_group_sort__[x.group])
            except Exception:
                # we probably have a pre-formed meta str we assume are in order
                pass
            meta_info = [str(m) for m in meta_info]
            # remove version if this is present in meta_info
            meta_info = [m for m in meta_info if 'fileformat=VCFv' not in m]
            self.meta.extend(meta_info)

        self.logger = medaka.common.get_named_logger('VCFWriter')

    def __enter__(self):
        self.handle = open(self.filename, self.mode, encoding='utf-8')
        self.handle.write('\n'.join('##' + line for line in self.meta) + '\n')
        self.handle.write('#' + '\t'.join(self.header) + '\n')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.handle.close()

    def write_variants(self, variants, sort=True):
        """Write variants to file, optionally sorting before writing."""
        if sort:
            variants = medaka.common.loose_version_sort(
                variants, key=lambda v: '{}-{}'.format(v.chrom, v.pos))

        for variant in variants:
            self.write_variant(variant)

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
        self.chroms = list()  # keep record of chroms in order they were read
        self._indexed = False
        self._tree = None
        self._parse_lock = Lock()
        self.logger = medaka.common.get_named_logger('VCFReader')

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
                    raise IOError(
                        '.vcf is unsorted at index #{}.'.format(index))
                if variant.chrom not in self.chroms:
                    self.chroms.append(variant.chrom)
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
                # clear out an incomplete parse, actually this doesn't matter
                # since the values in the tree are set-like.
                self._tree = collections.defaultdict(intervaltree.IntervalTree)
                for variant in self._parse():
                    self._tree[variant.chrom][
                        variant.pos:variant.pos + len(variant.ref)] = variant
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

        :yields: `Variant` instances.

        """
        if start is None:
            start = float('-inf')
        if end is None:
            end = float('inf')

        def _tree_search(tree, start, end, strict):
            search = tree.overlap if strict else tree.envelop
            return search(start, end)

        if not self.cache:
            # if not using a cache, just keep re-reading the file
            for variant in self._parse():
                if not all([
                        ref_name is None or variant.chrom == ref_name,
                        start is None or variant.pos > start,
                        end is None or variant.pos + len(variant.ref) < end]):
                    continue
                yield variant
        else:
            self.index()
            # spec says .vcf is sorted, lets follow. Keep ordering of
            # chromosomes as they are in the file, and sort positions.
            if ref_name is not None:
                results = sorted(_tree_search(
                    self._tree[ref_name], start, end, strict))
            else:
                results = itertools.chain(*(
                     sorted(_tree_search(
                         self._tree[chrom], start, end, strict=True))
                     for chrom in self.chroms
                ))
            for interval in results:
                yield interval.data


def _get_hap(v, trees):
    for hap, tree in enumerate(trees, 1):
        at_pos = tree.at(v.pos)
        for vp in at_pos:
            if vp.data is v:
                return hap


def _merge_variants(
        interval, trees, ref_seq, detailed_info=False, discard_phase=True):
    """Merge variants in an interval into a `Variant` obj

    .. note::

        It is assumed that variants in each haplotype have a single alt (an
        exception will be raised if this is not the case) and that that if two
        overlapping variants have the same alt, the GT it 1/1, else if the alts
        are different, the GT is 1/2 (or the phased equivalents if
        discard_phase is False)

    :param interval: `intervaltree.Interval` with .data containing list of
        `medaka.vcf.Variant` objs to be merged
    :param trees: iterable of `intervaltree.IntervalTree` objs containing the
        `medaka.vcf.Variant` objs of each haplotype (used to determine which
        haplotype variants in `interval` belong to).
    :param ref_seq: str, reference sequence
    :param detailed_info: bool, whether to add more detail to Variant info.
    :param discard_phase: bool, if False, preserve phase, else return unphased
        variants.

    :returns: `medaka.vcf.Variant` obj
    """
    if interval.end > len(ref_seq):
        raise ValueError(
            'A variant occurs after the end of the reference sequence.')
    ref = ref_seq[interval.begin: interval.end]
    alts_dict = collections.OrderedDict()
    info = {}
    mixed_vars = collections.defaultdict(list)
    for v in interval.data:
        mixed_vars[str(_get_hap(v, trees))].append(v)
    qual = 0.0
    for hap, hap_vars in sorted(mixed_vars.items()):
        alt = list(ref)
        for v in hap_vars:
            if len(v.alt) > 1:
                raise ValueError(
                    'Only single-allele variants from two vcfs can be merged')
            start_i = v.pos - interval.begin
            end_i = start_i + len(v.ref)
            if v.ref != ref[start_i:end_i]:
                msg = 'Variant ref {} does not match ref {} at {}:{}'
                raise ValueError(
                    msg.format(v.ref, ref[start_i:end_i], v.chrom, v.pos))
            # also check ref is correct within unsliced ref seq
            assert ref_seq[v.pos:v.pos + len(v.ref)] == v.ref
            alt[start_i:end_i] = [''] * len(v.ref)
            alt[start_i] = v.alt[0]
        # calculate mean GQ for each haplotype, and take mean of these for
        # overall qual.
        # Use mean otherwise we might need very different thresholds for
        # short vs long variants and homozygous vs heterozygous variants.
        info['q{}'.format(hap)] = sum((
            float(v.qual) for v in hap_vars)) / len(hap_vars)
        info['pos{}'.format(hap)] = ','.join(str(v.pos + 1) for v in hap_vars)
        if detailed_info:
            # + 1 as VCF is 1-based, v.pos is 0 based
            info['ref{}'.format(hap)] = ','.join((v.ref for v in hap_vars))
            info['alt{}'.format(hap)] = ','.join((v.alt[0] for v in hap_vars))
        qual += info['q{}'.format(hap)] / len(mixed_vars)
        alts_dict[hap] = ''.join(alt)

    haps = list(alts_dict.keys())
    alts = list(alts_dict.values())
    gt_sep = '/' if discard_phase else '|'
    if len(alts) == 2:
        if alts[0] == alts[1]:  # homozygous 1/1
            gt = gt_sep.join(len(haps) * '1')
            alts = alts[:1]
        else:  # heterozygous 1/2
            gt = gt_sep.join(map(str, haps))
    else:  # heterozygous 0/1
        assert len(haps) == 1
        gts = [0, 1]  # appropriate if hap with variant is 2
        if not discard_phase:
            if int(haps[0]) == 1:
                gts = [1, 0]
        gt = gt_sep.join(map(str, gts))

    sample = {'GT': gt, 'GQ': qual}
    return medaka.vcf.Variant(
        v.chrom, interval.begin, ref, alt=alts,
        filter='PASS', info=info, qual=qual,
        sample_dict=sample).trim()


class Haploid2DiploidConverter(object):

    def __init__(
            self, vcf1, vcf2, ref_fasta,
            only_overlapping=True, discard_phase=True, detailed_info=False):
        """Merge variants from two haploid VCFs into a diploid vcf. Variants in
        one file which overlap with variants in the other will have their alts
        padded.

        .. note::

            Variants in a single vcf file should not overlap with each other.

        :param vcf1, vcf2: paths to haploid vcf files.
        :param ref_fasta: path to reference.fasta file.
        :param only_overlapping: bool, merge only overlapping variants (not
            adjacent ones).
        :param discard_phase: bool, if False, preserve phase, else output
            unphased variants.
        """

        self.only_overlapping = only_overlapping
        self.discard_phase = discard_phase
        self.detailed_info = detailed_info

        self.logger = medaka.common.get_named_logger('VCFMERGE')

        self.vcfs = [medaka.vcf.VCFReader(vcf) for vcf in (vcf1, vcf2)]
        for vcf in self.vcfs:
            vcf.index()  # create tree
        self.fasta = pysam.FastaFile(ref_fasta)
        self.chroms = set(itertools.chain(*[v.chroms for v in self.vcfs]))

    def variants(self):
        """Yield diploid variants.

        :yields `medaka.vcf.Variant` objs
        """

        for chrom in medaka.common.loose_version_sort(self.chroms):
            self.logger.info('Merging variants in chrom {}'.format(chrom))
            merged = []
            trees = [vcf._tree[chrom] for vcf in self.vcfs]
            # assign haplotype so that otherwise identical variants in both
            # trees are not treated as identical (we need to be able to
            # distinguish between 0/1 and 1/1)
            for h, tree in enumerate(trees):
                for i in tree.all_intervals:
                    i.data.info['mhap'] = h
            comb = intervaltree.IntervalTree(
                trees[0].all_intervals.union(trees[1].all_intervals))
            # if strict, merge only overlapping intervals (not adjacent ones)
            comb.merge_overlaps(
                strict=self.only_overlapping,
                data_initializer=list(),
                data_reducer=lambda x, y: x + [y])
            ref_seq = self.fasta.fetch(chrom)
            for interval in comb.all_intervals:
                merged.append(_merge_variants(
                    interval, trees, ref_seq,
                    detailed_info=self.detailed_info,
                    discard_phase=self.discard_phase))
            yield from sorted(merged, key=lambda x: x.pos)

    @property
    def meta_info(self):
        m = []
        for h in 1, 2:
            m.append(MetaInfo(
                'INFO', 'pos{}'.format(h), '.', 'Integer',
                'POS of incorporated variants from haplotype {}'.format(h)))
            m.append(MetaInfo(
                'INFO', 'q{}'.format(h), 1, 'Float',
                'Combined qual score for haplotype {}'.format(h)))
        if self.detailed_info:
            for h in 1, 2:
                m.append(MetaInfo(
                    'INFO', 'ref{}'.format(h), '.', 'String',
                    'ref alleles of incorporated variants '
                    'from haplotype {}'.format(m)))
                m.append(MetaInfo(
                    'INFO', 'alt{}'.format(h), '.', 'String',
                    'alt alleles of incorporated variants '
                    'from haplotype {}'.format(m)))

        m.append(MetaInfo(
            'FORMAT', 'GT', 1, 'String', 'Genotype'))
        m.append(MetaInfo(
            'FORMAT', 'GQ', 1, 'Float', 'Genotype quality score'))
        return m


def haploid2diploid(args):
    """Entry point for merging two haploid vcfs into a diploid vcf"""

    convertor = Haploid2DiploidConverter(args.vcf1, args.vcf2, args.ref_fasta,
                                         only_overlapping=not args.adjacent)

    with VCFWriter(
            args.vcfout, 'w', version='4.1', contigs=convertor.chroms,
            meta_info=convertor.meta_info) as vcf_writer:
        for v in convertor.variants():
            vcf_writer.write_variant(v)


def split_variants(vcf_fp, trim=True):
    """Split diploid vcf into two haploid vcfs.

    :param vcf: str, path to input vcf.
    :param trim: bool, trim variants to minimal alt and ref and update pos.

    :returns: tuple of output files written
    """

    vcf = VCFReader(vcf_fp, cache=False)
    q = collections.defaultdict(list)
    for v in vcf.fetch():
        for k, v in v.split_haplotypes():
            if v is not None:
                v = v.trim() if trim else v
                q[k].append(v)

    basename, ext = os.path.splitext(vcf_fp)
    output_files = []
    for k, variants in q.items():
        output_files.append('{}_hap{}{}'.format(basename, k, ext))
        with medaka.vcf.VCFWriter(
                output_files[-1], meta_info=vcf.meta) as vcf_writer:
            # output variants in the same order they occured in the input file
            vcf_writer.write_variants(variants, sort=False)
    return tuple(output_files)


def diploid2haploid(args):
    """Entry point to split a diploid VCF into two haploid VCFs."""
    output_files = split_variants(args.vcf, trim=not args.notrim)
    logger = medaka.common.get_named_logger('Diploid2Haploid')
    logger.info('VCFs files written to: {}'.format('\n'.join(output_files)))


def classify_variant(v):
    """Classify a variant.

    :param v: `medaka.vcf.Variant` obj

    :returns: typ: str, variant classification.
    """

    def is_start_same(x):
        return all([a[0] == x.ref[0] for a in x.alt])

    def is_end_same(x):
        return all([a[-1] == x.ref[-1] for a in x.alt])

    def is_len_same(x):
        return all(len(a) == len(x.ref) for a in x.alt)

    def is_longer(x):
        return all([len(a) > len(v.ref) for a in x.alt])

    def is_shorter(x):
        return all([len(a) < len(v.ref) for a in x.alt])

    typ = 'other'
    if is_len_same(v):  # sub
        typ = 'snp' if len(v.ref) == 1 else 'mnp'
    elif is_start_same(v) or is_end_same(v):  # i.e. not a sub and indel
        if is_longer(v):  # ins
            if all(len(a) == len(v.ref) + 1 for a in v.alt):
                typ = 'sni'
            else:
                typ = 'mni'
        elif is_shorter(v):  # del
            if all(len(a) == len(v.ref) - 1 for a in v.alt):
                typ = 'snd'
            else:
                typ = 'mnd'
        elif all(len(a) != len(v.ref) for a in v.alt):
            typ = 'indel'  # could be mix of ins and del
    return typ


def classify_variants(args):
    """Entry point to classify variants and write out a new VCF for each class
    and class-groups.

    Example classifications for subsitutions are 'snp' and 'mnp', which belong
    to the same class-group 'sub'.
    """
    path, ext = os.path.splitext(args.vcf)

    base_typ = ['snp', 'mnp', 'sni', 'mni', 'snd', 'mnd', 'other', 'all']
    typs = {x: {x} for x in base_typ}
    typs['sub'] = {'snp', 'mnp'}
    typs['del'] = {'snd', 'mnd'}
    typs['ins'] = {'sni', 'mni'}
    typs['snid'] = {'sni', 'snd'}
    typs['mnid'] = {'mni', 'mnd'}
    typs['indel'] = set.union({'indel'}, typs['del'], typs['ins'])
    typs['all'] = set.union(*[s for s in typs.values()])

    vcfs = {s: '{}.{}{}'.format(path, s, ext) for s in typs.keys()}

    vcf = VCFReader(args.vcf, cache=True)

    type_meta = str(MetaInfo('INFO', 'type', 1, 'String', 'Type of variant'))
    if args.replace_info:
        meta_info = [type_meta]
    else:
        meta_info = vcf.meta + [type_meta]

    with contextlib.ExitStack() as stack:
        writers = {
            s: stack.enter_context(VCFWriter(
                vcfs[s], header=vcf.header, meta_info=meta_info))
            for s in typs.keys()}
        for v in vcf.fetch():
            v_typ = classify_variant(v)
            d = {'type': v_typ}
            if args.replace_info:
                v.info = d
            else:
                v.info.update(d)
            for typ, typ_set in typs.items():
                if v_typ in typ_set:
                    writers[typ].write_variant(v)


def vcf2tsv(args):
    """Entry point to convert vcf to tsv, unpacking info and sample fields"""
    try:
        import pandas as pd
    except ImportError:
        print("Converting vcf to tsv requires pandas to be installed.")
        print("Try running 'pip install pandas'")

    vcf = VCFReader(args.vcf, cache=False)
    df = pd.DataFrame((v.to_dict() for v in vcf.fetch()))
    df.to_csv(args.vcf + '.tsv', sep='\t', index=False)


def get_homozygous_regions(args):
    """Entry point to find homozygous regions from an input diploid VCF.
    """

    vcf = VCFReader(args.vcf, cache=False)
    reg = medaka.common.Region.from_string(args.region)
    if reg.start is None or reg.end is None:
        raise ValueError('Region start and end must be specified')

    def get_hetero_pos(v):
        gt = v.to_dict()['GT']
        pos = list()
        if gt[0] != gt[-1]:  # heterozygous
            pos = [p for p in range(v.pos, v.pos + len(v.ref))]
        return pos

    def get_homo_regions(ref_name, pos, min_len=1000):
        regions = []
        hetero_gaps = np.ediff1d(pos)
        sort_inds = np.argsort(hetero_gaps)[::-1]
        homo_len = 0
        for i in sort_inds:
            if hetero_gaps[i] < min_len:
                break
            homo_len += hetero_gaps[i]
            start = pos[i]
            end = start + hetero_gaps[i]
            regions.append(medaka.common.Region(ref_name, start, end))
        return regions, homo_len

    pos = [reg.start]
    for v in vcf.fetch(ref_name=reg.ref_name, start=reg.start, end=reg.end):
        pos.extend(get_hetero_pos(v))
    pos.append(reg.end)

    homo_regions, homo_len = get_homo_regions(
        reg.ref_name, pos, min_len=args.min_len)

    # sort regions by start
    homo_regions.sort(key=lambda r: r.start)

    # invert homozygous regions to get regions containing heterozygous calls
    hetero_regions = []
    start = reg.start
    hetero_len = 0
    full_reg = [medaka.common.Region(reg.ref_name, start=reg.end, end=None)]
    for homo_reg in homo_regions + full_reg:
        end = homo_reg.start
        if end - start > args.min_len:
            hetero_regions.append(
                medaka.common.Region(reg.ref_name, start, end))
            hetero_len += end - start
        start = homo_reg.end

    homo_fp = 'homozygous_' + args.suffix
    with open(homo_fp, 'w') as fh:
        fh.write('\n'.join((r.name for r in homo_regions)))
    hetero_fp = 'heterozygous_' + args.suffix
    with open(hetero_fp, 'w') as fh:
        fh.write('\n'.join((r.name for r in hetero_regions)))

    reg_len = reg.end - reg.start
    logger = medaka.common.get_named_logger('REGIONS')
    msg = 'Found {} {} regions > {} within {}, spanning {}/{} ({:.2%})'
    logger.info(msg.format(
        len(homo_regions), 'homozygous', args.min_len,
        reg.name, reg_len, homo_len, homo_len / reg_len))
    logger.info(msg.format(
        len(hetero_regions), 'heterozygous', args.min_len,
        reg.name, reg_len, hetero_len, hetero_len / reg_len))
