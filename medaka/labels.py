"""Schemes to encode and decode truth labels and network outputs."""
import abc
import collections
from copy import copy
import functools
import itertools
from operator import attrgetter

import intervaltree
import numpy as np
import pysam

import medaka.common
import medaka.rle
import medaka.vcf


class TruthAlignment(object):
    """Process truth alignments."""

    def __init__(self, alignment):
        """Create a `TruthAlignment` list from an `AlignedSegment`.

        :param alignment: `pysam.AlignedSegment`.
        """
        self.aln = alignment  # so we can get positions and labels later
        # initialise start and end (which might be moved)
        self.start = self.aln.reference_start  # zero-based
        self.end = self.aln.reference_end
        self.is_kept = True
        self.logger = medaka.common.get_named_logger('TruthAlign')

    def _get_overlap_with(self, other):
        first, second = sorted((self, other),
                               key=attrgetter('aln.reference_start'))
        overlap = first.aln.reference_end > second.aln.reference_start
        if overlap:
            # return positions of start and end of overlapping region
            return second.aln.reference_start, first.aln.reference_end
        else:
            return None

    @staticmethod
    def _filter_alignments(alignments, region, min_length=1000,
                           length_ratio=2.0, overlap_fraction=0.5):
        """Filter alignments to yield only segments suitable for training.

        :param alignments: iterable of `TruthAlignment` objects.
        :param region: `medaka.common.Region` obj. all alignments will be
            trimmed to this region.
        :param min_length: int, minimum length of alignment segment.
        :param length_ratio: float, in cases of overlap, ratio of
            longer segment / shorter segment
            above which the longer segment is assumed to be more reliable.
        :param overlap_fraction: float, length_ratio: float, in cases of
            overlap, fraction of shorter segment overlapping with the longer
            segment above which the segments are considered highly overlapping.
        :returns: list of `TruthAlignment` objects

        """
        # don't want to modify original alignments
        filtered_alignments = [copy(a) for a in alignments]

        # git rid of any alignments with ambiguity bases in reference or query
        def only_valid_symbols(al):
            """`TruthAlignment` is free of ambiguous bases in ref and query."""
            symbols = set(list('ACGT'))
            ref = al.aln.get_reference_sequence().upper()
            query = al.aln.query_sequence.upper()
            result = all([set(ref).issubset(symbols),
                          set(query).issubset(symbols)])
            return result

        filtered_alignments = [al for al in filtered_alignments
                               if only_valid_symbols(al)]

        for al_i, al_j in itertools.combinations(filtered_alignments, 2):
            first, second = sorted((al_i, al_j),
                                   key=attrgetter('aln.reference_start'))
            overlap = first._get_overlap_with(second)
            if overlap is None:
                continue
            ovlp_start, ovlp_end = overlap
            shorter, longer = sorted((al_i, al_j),
                                     key=attrgetter('aln.reference_length'))
            length_ratio_ij = (longer.aln.reference_length /
                               shorter.aln.reference_length)
            overlap_fraction_ij = ((ovlp_end - ovlp_start) /
                                   shorter.aln.reference_length)
            # 4 cases
            # we don't trust one more than the other
            if length_ratio_ij < length_ratio:
                if overlap_fraction_ij >= overlap_fraction:
                    # 1) large overlap; significant ambiguity, discard both
                    shorter.is_kept = False
                    longer.is_kept = False
                else:
                    # 2) small overlap; trim overlapping portions of alignments
                    first.end = ovlp_start
                    second.start = ovlp_end
            else:  # we trust the longer one more than the shorter one
                if overlap_fraction_ij >= overlap_fraction:
                    # 3) large overlap; discard shorter alignment
                    shorter.is_kept = False
                else:
                    # 4) small overlap; trim overlapping portion of
                    # shorter alignment
                    second.start = ovlp_end

        # trim starts and ends if needed
        if region.start > 0 or region.end is not None:
            for al in filtered_alignments:
                if region.start > 0:
                    al.start = max(region.start, al.start)
                if region.end is not None:
                    al.end = min(region.end, al.end)
        # do filtering
        filtered_alignments = [al for al in filtered_alignments
                               if (al.is_kept and
                                   al.end - al.start >= min_length)]
        filtered_alignments.sort(key=attrgetter('start'))
        return filtered_alignments

    @staticmethod
    def bam_to_alignments(truth_bam, region, haplotag=None, min_length=1000):
        """Get processed truth alignments.

        :param truth_bam: (sorted indexed) bam with true sequence
            aligned to reference
        :param region: `medaka.common.Region` obj,
            (all alignments with any overlap with the interval
            start:end will be retrieved)
        :param haplotag: bam tag specifying which haplotype the alignment
            belongs to (for polyploid labels)
        :param min_length: minimum length for valid alignments.

        :returns: list of tuples where each tuple contains `TruthAlignment`
            for each haplotype trimmed to common genomic window.

        """
        algns = TruthAlignment._load_alignments(truth_bam, region, haplotag)
        # filter truth alignments to restrict ourselves to regions of the
        # ref where the truth is unambiguous in each haplotype
        algns = {
            h: TruthAlignment._filter_alignments(
                h_algns, region=region, min_length=min_length)
            for h, h_algns in algns.items()}
        # Group truth alignments from the haplotypes by genomic window and trim
        # to common genomic range
        if len(algns) == 0:
            return []
        else:
            grouped = TruthAlignment._group_and_trim_by_haplotype(algns)
            return grouped

    @staticmethod
    def _group_and_trim_by_haplotype(alignments):
        """Group alignments by haplotype tag and trim to common genomic window.

        :param alignments: {haplotype: [`TruthAlignment`]}
        :returns: list of tuples where each tuple contains `TruthAlignment`
            for each haplotype trimmed to common genomic window.

        .. note:: We should avoid the situation of staggered alignments
             which could occur by independently chunking each haplotype
             by chunking the draft and aligning to both haplotypes, then
             chunking both haplotypes according to draft-chunks, then realining
             haplotype chunks to back to the draft - this should minimize
             staggering of truth alignments and hence the number of labels
             discarded.
        """
        logger = medaka.common.get_named_logger("Group_and_trim")
        haplotypes = sorted(list(alignments.keys()))
        if len(haplotypes) == 1:  # haploid
            grouped = [(a,) for a in alignments[haplotypes[0]]]
        else:
            # create interval trees for other haplotypes
            trees = {}
            for h in haplotypes[1:]:
                trees[h] = intervaltree.IntervalTree(
                    [intervaltree.Interval(a.start, a.end, a)
                        for a in alignments[h]])
            # loop over alignments in first haplotype and find overlapping
            # alignments in other haplotypes. If there are multiple overlapping
            # alignments, take the one with the longest overlap.
            grouped = []
            for a in alignments[haplotypes[0]]:
                group = [a]
                common_start, common_end = a.start, a.end
                for h, tree in trees.items():
                    h_algns = list(tree.overlap(common_start, common_end))
                    # no alignments on this haplotype, skip this alignment
                    if len(h_algns) == 0:
                        break
                    # find most overlapping alignment
                    elif len(h_algns) > 1:
                        ovlps = [min(common_end, o.end)
                                 - max(common_start, o.begin)
                                 for o in h_algns]
                        h_algn = h_algns[np.argmax(ovlps)].data
                    # keep this alignment
                    elif len(h_algns) == 1:
                        h_algn = h_algns[0].data

                    common_start = max(common_start, h_algn.start)
                    common_end = min(common_end, h_algn.end)
                    group.append(h_algn)
                # skip this group
                if len(group) != len(haplotypes):
                    msg = 'Skipping {}:{}-{}; missing alignment for haplotype'
                    logger.info(msg.format(a.aln.reference_name,
                                           a.start, a.end))
                    continue
                else:
                    # trim all alignments to common start/end
                    for i in group:
                        i.start = common_start
                        i.end = common_end
                    grouped.append(tuple(group))
        return grouped

    @staticmethod
    def _load_alignments(truth_bam, region, haplotag=None):
        """Create list of `TruthAlignment` s from a truth bam.

        :param truth_bam: (sorted indexed) bam with true sequence(s) aligned
            to reference
        :param region: `medaka.common.Region`
        :param haplotag: bam tag specifying which haplotype the alignment
            belongs to (for polyploid labels)

        :returns: {haplotype: [`TruthAlignment`]}

        """
        alignments = collections.defaultdict(list)
        with pysam.AlignmentFile(truth_bam, 'rb') as bamfile:
            aln_reads = bamfile.fetch(
                reference=region.ref_name,
                start=region.start, end=region.end)
            for r in aln_reads:
                if (r.is_unmapped or r.is_secondary):
                    continue
                else:
                    hap = r.get_tag(haplotag) if haplotag is not None else None
                    alignments[hap].append(TruthAlignment(r))

        logger = medaka.common.get_named_logger("TruthAlign")
        for hap in alignments.keys():
            alignments[hap].sort(key=attrgetter('start'))
            logger.info("Retrieved {} alignments for haplotype {}.".format(
                len(alignments[hap]), hap))
        return alignments


label_schemes = dict()


class LabelSchemeRegistrar(type):
    """Class for registering label schemes."""

    def __new__(cls, clsname, bases, attrs):
        """Register class to `label_schemes` dict upon instantiation."""
        newclass = super(LabelSchemeRegistrar, cls).__new__(
            cls, clsname, bases, attrs)
        cls.register_label_scheme(clsname, newclass)
        return newclass

    def register_label_scheme(clsname, cls):
        """Add `LabelScheme` to `label_schemes` dict."""
        # do not display base class as command line option
        if clsname != 'BaseLabelScheme':
            label_schemes[clsname] = cls


class LabelSchemeMeta(abc.ABC, LabelSchemeRegistrar):
    """Metaclass facilitating registration of `LabelScheme` s."""

    pass


class BaseLabelScheme(metaclass=LabelSchemeMeta):
    """Base class for labelling schemes.

    Labels are 'truth' or 'y' vectors used to train models.
    Models aim to predict a label (y) from a feature (X).

    A LabelScheme contains all logic for creating label vectors
    from a source of truth (typically an alignment of true reference
    sequence(s) to a draft sequence).

    A LabelScheme class also contains all logic for converting the
    model output into useful output (usually  a consensus
    sequence or a set of variant calls).

    Encoding (conversion of truth alignments to training vectors):

    The public method `self.encode` converts alignments
    to arrays of integer-encoded intermediate representations that
    are saved to file during feature creation.
    These private method must be implemented.

        - `_alignment_to_pairs` (called by _alignment_to_labels)
        - `_labels_to_encoded_labels`

    At training time, it is necessary to call the following public
    method for the conversion of intermediate forms saved to file to
    training vectors input to the model training process.

        - `encoded_labels_to_training_vectors`

    Decoding (interpretation of network outputs):

    For SNP calling, the public method is `self.decode_snps`
    This private method must be implemented.

        - `_prob_to_snp` (for SNP decoding)

    Public methods may also be defined for other types of decoding.

    For example, the concrete HaploidLabelScheme defines:

        - `decode_consensus`
        - `decode_variants`

    """

    # default set of symbols used throughout all LabelSchemes
    symbols = '*ACGT'

    # flag to turn on decoding of variants with extra information
    verbose = True

    @property
    @abc.abstractmethod
    def n_elements(self):
        """Return number of elements provided by truth alignment.

        (Mostly) synonymous with ploidy. i.e. n_elements = 2 where
        two symbols are provided for the two haplotypes of a truth
        alignment for diploid training.
        """
        return 1

    @property
    @abc.abstractmethod
    def num_classes(self):
        """Return number of elements in output layer of neural network."""

    @property
    @abc.abstractmethod
    def padding_vector(self):
        """Return the integer encoding used to denote a gap.

        Used to inform calling programs that need to expand
        label arrays to align with feature arrays,
        where reads introduce minor positions.
        """

    @staticmethod
    def _singleton(it):
        """Test whether iterable contains one unique element.

        :param it: iterable

        :returns: bool, iterable contains one unique element
        """
        return len(set(it)) == 1

    @staticmethod
    def _phred(err, cap=70.0):
        """Calculate phred quality score.

        :param err: float, error probability
        :param cap: float, maximum q score (prevent inf)

        :returns: float, phred quality score
        """
        # add smallest positive usable number to err to avoid
        # RuntimeWarning: divide by zero encountered in log10
        err += np.finfo(float).tiny
        q = -10 * np.log10(err)
        return np.minimum(q, cap)

    @staticmethod
    def _pfmt(p, dp=3):
        """Cast float to string rounding to dp decimal places.

        Used to format probabilities and quality scores for vcf output.
        """
        return '{:.{dp}f}'.format(round(p, dp), dp=dp)

    @abc.abstractmethod
    def _alignment_to_pairs(self, aln):
        """Convert `pysam.AlignedSegment` to aligned pairs."""

    def _alignments_to_labels(self, truth_alns):
        """Convert `TruthAlignment` s to nd.array of positions and labels.

        :param truth_alns: tuple of `TruthAlignment` s for each haplotype
            spanning the same genomic range

        :returns: tuple(positions, label_array)
            - positions: numpy structured array with 'major'
            (reference position index) and 'minor'
            (trailing insertion index) fields.

            - label_array: numpy array of self.n_elements tuples

        """
        n_haps = len(truth_alns)
        if n_haps != self.n_elements:
            raise ValueError(
                '{} alignments were passed to {}, requires {}'.format(
                    n_haps, type(self), self.n_elements))

        # ensure all alignments have same start and end
        if not (self._singleton(a.start for a in truth_alns) and
                self._singleton(a.end for a in truth_alns)):

            raise ValueError('Alignments must have identical \
                genomic start and end.')

        # list of dicts mapping position to symbol for haplotypes
        pos_maps = list()

        for aln in truth_alns:
            # default to gap on lookup
            pos_to_symbol = collections.defaultdict(lambda: '*')
            ins_count = 0

            labels = self._alignment_to_pairs(aln.aln)
            trimmed_labels = itertools.dropwhile(
                lambda x: (x[0] is None) or (x[0] < aln.start), labels)
            for rpos, label in trimmed_labels:
                if rpos is not None and rpos >= aln.end:
                    break
                if rpos is None:
                    ins_count += 1
                else:
                    ins_count = 0
                    current_pos = rpos

                pos = (current_pos, ins_count)
                pos_to_symbol[pos] = label

            pos_maps.append(pos_to_symbol)

        # get ordered positions present in any haplotype mapping
        positions = sorted(set(itertools.chain.from_iterable(
            h.keys() for h in pos_maps)))

        # n_element tuples of symbols from haplotypes
        labels = [tuple(h[pos] for h in pos_maps)
                  for pos in positions]
        positions = np.array(positions,
                             dtype=[('major', int), ('minor', int)])

        return (positions, labels)

    @abc.abstractmethod
    def _labels_to_encoded_labels(self, labels):
        """Convert list of label tuples to array of integer encoded labels.

        The logic of many to one mappings, where multiple labels
        map to a common integer encoding is specified here.

        :param labels: list of label tuples

        :returns: np.ndarray of integer (tuples)
        """

    @abc.abstractmethod
    def encoded_labels_to_training_vectors(self, enc_labels):
        """Convert integer (tuple) encoded labels to truth vectors.

        (e.g. one-hot encoded truth vectors) that represent the truth
        for comparison with network output logits in
        e.g. metric(truth, pred) or loss(truth, pred)) functions.

        :param enc_labels: np.ndarray of integer (tuples)

        :returns: np.ndarray of training vectors
        """

    @property
    @abc.abstractmethod
    def _encoding(self):
        """Return a dictionary mapping from label tuple to integer (tuple).

        This property is accessed by methods that are responsible
        for encoding. Also, it can be written to file in order to
        log the encoding scheme used.
        """

    @property
    @functools.lru_cache(1)
    def _decoding(self):
        """Return a dictionary mapping from integer (tuple) to label tuple.

        Inverse of encoding.
        """
        return {v: k for k, v in self._encoding.items()}

    @property
    def _unitary_encoding(self):
        """Return a dictionary mapping from all symbol 1-tuples to integers."""
        return {v: k for k, v in enumerate(self._unitary_labels())}

    @property
    def _unitary_decoding(self):
        """Return a dictionary mapping from integers to all symbol 1-tuples."""
        return {v: k for k, v in self._unitary_encoding.items()}

    def encode(self, truth_alns):
        """Convert truth alignment(s) to array of intermediate representation.

        In most cases the intermediate representation consists of integers.

        :param truth_alns: tuple of `pysam.AlignedSegment` s for each haplotype
            spanning the same genomic range.

        :returns: tuple(positions, training_vectors)

            - positions: numpy structured array with 'major'
              (reference position index) and 'minor'
              (trailing insertion index) fields.

            - encoded: nd.array of encoded labels

        .. note ::
            It is generally the case that the returned encoded labels must be
            padded with encoded gap labels when aligned to corresponding
            training feature data.

        """
        # Labels is a list of tuples with alleles ('A', ), ('A', 'C'), ('C', 3)
        positions, labels = self._alignments_to_labels(truth_alns)

        # Encoded is an array of integers
        encoded = self._labels_to_encoded_labels(labels)

        return positions, encoded

    def _unitary_labels(self):
        """Return all symbol 1-tuples."""
        return tuple((s,) for s in self.symbols)

    def _unordered_label_combinations(self):
        """Generate all combinations of n_elements tuples; order is ignored.

        ('A','T') == ('T','A').
        """
        return tuple(itertools.combinations_with_replacement(
            self.symbols, self.n_elements))

    def _decode_snps(self, sample):
        """Convert network output in sample to a set of medaka.vcf.Variants.

        Recording SNPs but NOT indels.

        :param sample: medaka.common.Sample

        :returns: list of medaka.vcf.Variant objects for SNPs
        """
        ref_name = sample.ref_name
        pos = sample.positions
        probs = sample.label_probs

        if self.ref_vcf is not None:
            # process (all) loci in ref_vcf
            loci = {v.pos for v in self.ref_vcf.fetch(
                ref_name=ref_name, start=sample.first_pos,
                end=sample.end_pos)}

            indices = list()
            ref_symbols = list()
            for i in range(len(probs)):
                reference_symbol = self.ref_seq[pos['major'][i]]
                # insertions are ignored
                # ref_symbol must be in standard set
                if not all((pos['major'][i] in loci,
                            pos['minor'][i] == 0,
                            reference_symbol in self.symbols)):
                    continue
                indices.append(i)
                ref_symbols.append(reference_symbol)

            snps = self._prob_to_snp(
                probs[indices], pos['major'][indices], ref_name,
                ref_symbols, return_all=True)

        else:
            # process all loci but only retain variant loci
            indices = list()
            ref_symbols = list()
            for i in range(len(probs)):
                reference_symbol = self.ref_seq[pos['major'][i]]
                if not all((pos['minor'][i] == 0,
                            reference_symbol in self.symbols)):
                    continue
                indices.append(i)
                ref_symbols.append(reference_symbol)

            snps = self._prob_to_snp(
                probs[indices], pos['major'][indices], ref_name,
                ref_symbols, return_all=False)

        return snps

    def decode_snps(self, sample, ref_seq, ref_vcf=None, threshold=0.04):
        """Decode network outputs to medaka.vcf.Variant objects recording SNPs.

        Optionally, SNPs are returned at ALL loci specified in suplied ref_vcf.

        :param sample: `medaka.common.Sample`
        :param ref_seq: reference sequence, should be upper case.
        :param ref_vcf: `.vcf` file
        :param threshold: threshold for acceptance of secondary call.

        :returns: list of `medaka.vcf.Variant` objects

        """
        self.ref_seq = ref_seq
        self.secondary_threshold = threshold
        self.ref_vcf = medaka.vcf.VCFReader(ref_vcf) \
            if ref_vcf else None

        return self._decode_snps(sample)

    @abc.abstractmethod
    def _prob_to_snp(
            self, outputs, positions, ref_name, ref_symbols,
            return_all=False):
        """Convert network output to `medaka.common.Variant` s.

        This method contains all logic for converting network
        output into SNP Variants.

        :param outputs: np.ndarray of network outputs
            for multiple loci
        :param pos: array of genomic indices
        :param ref_name: reference name
        :param ref_symbols: the symbols at ref_seq[pos]
        :param return_all: return a `medaka.vcf.Variant` even if
            there is no SNP at a locus.

        :returns: a list of `medaka.vcf.Variant` s
        """

    @property
    def snp_metainfo(self):
        """Return meta data for use in `.vcf` header."""
        MI = medaka.vcf.MetaInfo
        m = [MI('FORMAT', 'GT', 1, 'String', 'Medaka genotype'),
             MI('FORMAT', 'GQ', 1, 'Integer', 'Medaka genotype quality score')]
        if self.verbose:
            m.extend([MI('INFO', 'ref_prob', 1, 'Float',
                         'Medaka probability for reference allele'),
                      MI('INFO', 'primary_prob', 1, 'Float',
                         'Medaka probability of primary call'),
                      MI('INFO', 'primary_call', 1, 'String',
                         'Medaka primary call'),
                      MI('INFO', 'secondary_prob', 1, 'Float',
                         'Medaka probability of secondary call'),
                      MI('INFO', 'secondary_call', 1, 'String',
                         'Medaka secondary call'), ])

        return m


class HaploidLabelScheme(BaseLabelScheme):
    """A single-element label per genomic position.

    Unitary symbols map to integers, which map to one-hot encoded
    training vectors.

    Encoding ('A',) -> 1 -> [0, 1, 0, 0, 0]

    Consensus decoding is simple argmax.
    [[0.02, 0.9, 0.02, 0.01, 0.05], [0.01, 0.02, 0.01, 0.9, 0.05]] -> "AG"

    SNP decoding utilises a hard threshold to define secondary calls.
    """

    @property
    def n_elements(self):
        """Return number of elements provided by truth alignment.

        Synonymous with ploidy.
        """
        return 1

    @property
    def num_classes(self):
        """Return number of elements in output layer of neural network."""
        return len(self._decoding)

    @property
    def padding_vector(self):
        """Return the integer encoding used to denote a gap."""
        gap = [('*',)]
        return self._labels_to_encoded_labels(gap)[0]

    @property
    @functools.lru_cache(1)
    def _encoding(self):
        """Return a dictionary mapping from label tuple to integer (tuple).

        This property is accessed by methods that are responsible
        for encoding. Also, it can be written to file in order to
        log the encoding scheme used.
        """
        return self._unitary_encoding

    def _alignment_to_pairs(self, aln):
        """Convert `pysam.AlignedSegment` to aligned pairs."""
        seq = aln.query_sequence
        for qpos, rpos in aln.get_aligned_pairs():
            yield rpos, seq[qpos].upper() if qpos is not None else '*'

    def _labels_to_encoded_labels(self, labels):
        """Convert list of label tuples to array of integer encoded labels."""
        return np.fromiter(
            (self._encoding[x] for x in labels), dtype=int)

    def encoded_labels_to_training_vectors(self, enc_labels):
        """Convert integer (tuple) encoded labels to truth vectors.

        (e.g. sparse one-hot encoded truth vectors) that represent the truth
        for comparison with network output logits in
        e.g. metric(truth, pred) or loss(truth, pred)) functions.
        """
        # TODO remove once legacy files extinct.
        # legacy features had (base, runlength) with the encoding:
        # gap, lowercase bases, uppercase bases, other stuff
        if len(enc_labels.dtype) == 2:
            enc_labels = np.array(
                [max(0, x[0] - 4) for x in enc_labels], dtype='int64')
        return np.expand_dims(enc_labels, axis=1)  # sparse 1-hot

    def _prob_to_snp(
            self, outputs, positions, ref_name, ref_symbols,
            return_all=False):
        """Convert networkout output to `medaka.common.Variant` s.

        A threshold is used to define significant secondary calls
        allowing the prediciton of homozygous and heterozygous
        diploid calls.

        Where a significant secondary call is a deletion, the
        SNP is considered homozygous in the primary call.
        """
        results = list()
        for network_output, pos, ref_symbol in zip(
                outputs, positions, ref_symbols):
            # TODO: some optimisation here?
            secondary_call, primary_call = (
                self._decoding[p][0] for p in np.argsort(network_output)[-2:])

            secondary_prob, primary_prob = np.sort(network_output)[-2:]
            ref_prob = network_output[self._encoding[(ref_symbol,)]]

            if self.verbose:
                info = {
                    'ref_prob': self._pfmt(ref_prob),
                    'primary_prob': self._pfmt(primary_prob),
                    'primary_call': primary_call,
                    'secondary_prob': self._pfmt(secondary_prob),
                    'secondary_call': secondary_call}
            else:
                info = {}

            # logical tests
            primary_is_reference = primary_call == ref_symbol
            primary_is_deletion = primary_call == '*'
            secondary_is_deletion = secondary_call == '*'
            secondary_exceeds_threshold = \
                secondary_prob >= self.secondary_threshold

            # homozygous snp
            if all((not primary_is_reference,
                    not primary_is_deletion,
                    not secondary_exceeds_threshold)):

                alt = primary_call
                qual = self._phred(1 - primary_prob)
                genotype = {'GT': '1/1', 'GQ': self._pfmt(qual, 0)}
                results.append(medaka.vcf.Variant(
                    ref_name, pos, ref_symbol, alt, filt='PASS', info=info,
                    qual=self._pfmt(qual), genotype_data=genotype))

            # heterozygous, no deletions
            elif all((not primary_is_deletion,
                      not secondary_is_deletion,
                      secondary_exceeds_threshold)):

                err = 1 - (primary_prob + secondary_prob)
                qual = self._phred(err)
                # filtering by list comp maintains order
                alt = [c for c in [primary_call, secondary_call]
                       if not c == ref_symbol]
                gt = '0/1' if len(alt) == 1 else '1/2'
                genotype = {'GT': gt, 'GQ': self._pfmt(qual, 0)}

                results.append(medaka.vcf.Variant(
                    ref_name, pos, ref_symbol, alt, filt='PASS',
                    info=info, qual=self._pfmt(qual), genotype_data=genotype))

            # heterozygous, secondary deletion
            elif all((not primary_is_reference,
                      not primary_is_deletion,
                      secondary_is_deletion,
                      secondary_exceeds_threshold)):

                alt = primary_call
                qual = self._phred(1 - primary_prob)
                genotype = {'GT': '1/1', 'GQ': self._pfmt(qual, 0)}
                results.append(medaka.vcf.Variant(
                    ref_name, pos, ref_symbol, alt, filt='PASS',
                    info=info, qual=self._pfmt(qual), genotype_data=genotype))

            # no snp at this location
            else:
                if return_all:
                    # return variant even though it is not a snp
                    qual = self._phred(1 - primary_prob)
                    genotype = {'GT': '0/0', 'GQ': self._pfmt(qual, 0)}
                    results.append(medaka.vcf.Variant(
                        ref_name, pos, ref_symbol, alt='.', filt='PASS',
                        info=info, qual=self._pfmt(qual),
                        genotype_data=genotype))
        return results

    def decode_variants(self, sample, ref_seq):
        """Convert network output in sample to a set of `medaka.vcf.Variant` s.

        A consensus sequence is decoded and compared with a reference sequence.
        Both substitution and indel variants that may be multi-base will be
        reported in the output `medaka.vcf.Variant` s.

        :param sample: `medaka.common.Sample`.
        :param ref_seq: reference sequence, should be upper case.

        :returns: list of `medaka.vcf.Variant` objects.
        """
        pos = sample.positions
        probs = sample.label_probs

        assert sample.positions['minor'][0] == 0

        # array of symbols retaining gaps
        predicted = np.array(list(
            self.decode_consensus(sample, with_gaps=True)))

        # get reference sequence with insertions marked as '*'

        def get_symbol(p):
            return ref_seq[p['major']] if p['minor'] == 0 else '*'

        reference = np.fromiter((get_symbol(p) for p in pos),
                                dtype='|U1', count=len(pos))

        # find variants by looking for runs of labels which differ.
        # If both labels are gap, we don't want to consider this a
        # match so as to avoid splitting up long insertions if
        # there happens to be a gap label in ref and pred.
        mismatch = predicted != reference
        both_gap = np.logical_and(predicted == '*',
                                  reference == '*')
        is_variant = np.logical_or(mismatch, both_gap)

        # medaka.common.rle requires numeric input
        runs = medaka.rle.rle(is_variant)
        variant_runs = runs[np.where(runs['value'])]
        variants = []
        encoding = self._encoding
        for run in variant_runs:
            start = run['start']
            end = start + run['length']
            # if call or ref starts with gap, pad left
            # (or right if we are at start of genome).
            pad_right = False
            # if chunks start with a deletion, need to pad (see comment below).
            pad_del_at_start_of_chunk = False
            while ((not pad_right and
                    (reference[start] == '*' or
                     predicted[start] == '*'))
                   or
                   (pad_right and
                    (reference[end] == '*' or
                     predicted[end] == '*'))):

                if tuple(pos[start]) == (0, 0) or pad_right:
                    end += 1
                    # avoid getting stuck in loop if there is a run
                    # of dels at start of ref
                    pad_right = True
                elif (start == 0 and pos[start]['minor'] == 0
                      and predicted[start] == '*'):
                    # If variant calling is being performed on a region (rather
                    # than an entire chr), it is possible that a chunk will
                    # start with a deletion. In which case, the VCF spec says
                    # one should pad left.  This means creating a variant which
                    # asserts that the base at the previous position was not a
                    # variant - something we have no evidence of.
                    pad_del_at_start_of_chunk = True
                    break
                else:
                    assert start != 0
                    start -= 1

            var_ref_with_gaps = ''.join(s for s in reference[start:end])
            var_ref = var_ref_with_gaps.replace('*', '')

            var_pred_with_gaps = ''.join(s for s in predicted[start:end])
            var_pred = var_pred_with_gaps.replace('*', '')

            if var_ref == var_pred:
                # not a variant
                continue
            elif not set(var_ref).issubset(set(self.symbols)):
                # don't call where reference is ambiguous
                continue

            var_probs = probs[start:end]

            var_ref_encoded = (encoding[(s,)] for s in var_ref_with_gaps)
            var_pred_encoded = (encoding[(s,)] for s in var_pred_with_gaps)

            ref_probs = np.array([var_probs[i, j]
                                 for i, j in enumerate(var_ref_encoded)])
            pred_probs = np.array([var_probs[i, j]
                                  for i, j in enumerate(var_pred_encoded)])

            ref_quals = [self._phred(1 - p) for p in ref_probs]
            pred_quals = [self._phred(1 - p) for p in pred_probs]

            if self.verbose:
                info = {
                    'ref_seq': var_ref_with_gaps,
                    'pred_seq': var_pred_with_gaps,
                    'ref_qs': ','.join((self._pfmt(q) for q in ref_quals)),
                    'pred_qs': ','.join((self._pfmt(q) for q in pred_quals)),
                    'ref_q': self._pfmt(sum(ref_quals)),
                    'pred_q': self._pfmt(sum(pred_quals)),
                    'n_cols': len(pred_quals)
                }
            else:
                info = {}

            # log likelihood ratio
            qual = sum(pred_quals) - sum(ref_quals)
            genotype = {'GT': '1', 'GQ': self._pfmt(qual, 0)}

            var_pos = pos['major'][start]

            if pad_del_at_start_of_chunk:
                var_pos -= 1
                var_ref = ref_seq[var_pos] + var_ref
                var_pred = ref_seq[var_pos] + var_pred

            variant = medaka.vcf.Variant(sample.ref_name, var_pos, var_ref,
                                         alt=var_pred, filt='PASS', info=info,
                                         qual=self._pfmt(qual),
                                         genotype_data=genotype)
            variant = variant.trim()
            variants.append(variant)

        return variants

    @property
    def variant_metainfo(self):
        """Return meta data for use in `.vcf` header."""
        MI = medaka.vcf.MetaInfo

        m = [MI('FORMAT', 'GT', 1, 'String',
                'Medaka genotype.'),
             MI('FORMAT', 'GQ', 1, 'Integer',
                'Medaka genotype quality score'), ]
        if self.verbose:
            m.extend([MI('INFO', 'ref_seq', 1, 'String',
                      'Medaka reference sequence'),
                      MI('INFO', 'pred_seq', 1, 'String',
                      'Medaka predicted sequence'),
                      MI('INFO', 'ref_qs', '.', 'Float',
                      'Medaka quality score for reference'),
                      MI('INFO', 'pred_qs', '.', 'Float',
                      'Medaka quality score for prediction'),
                      MI('INFO', 'ref_q', 1, 'Float',
                      'Medaka per position quality score for reference'),
                      MI('INFO', 'pred_q', 1, 'Float',
                      'Medaka per position quality score for prediction'),
                      MI('INFO', 'n_cols', 1, 'Integer',
                      'Number of medaka pileup columns in variant call')])
        return m

    def decode_consensus(self, sample, with_gaps=False):
        """Convert network output to consensus sequence by argmax decoding.

        :param sample: medaka.common.Sample
        :param with_gaps: include gap ("*") characters in output.

        :returns: str, consensus sequence
        """
        # property access is slow
        decode = self._decoding
        # most probable class
        mp = np.argmax(sample.label_probs, -1)
        seq = ''.join((decode[x][0] for x in mp))
        # delete gap symbol from sequence
        if not with_gaps:
            seq = seq.replace('*', '')
        return seq


class DiploidLabelScheme(BaseLabelScheme):
    """LabelScheme defines a two-element label per genomic position.

    Each label maps to an integer, and we make a direct diploid call.
    A categorical crossentropy loss is appropriate.

    Encoding: ('A', 'C') -> 6 -> [6] (a sparse one-hot encoded training vector)

    SNP decoding: [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0] -> ('A', 'C')
    (network output is not constrained to be 0 or 1, but for simplicity
    we show the argmax probability as 1 and all others as 0)

    Depending on the reference, this may result in a homozygous or
    heterozygous SNP call.

    Where reference is ('A',), we output ref='A', alt=['C'], GT='0/1'
    Where reference is ('T',), we output ref='T' alt=['A', 'C'], GT='1/2'

    When predicion is homozygous e.g. ('A', 'A'):
    Where reference is ('T',), we output ref='T', alt=['A'], GT='1/1'

    When the prediction contains a gap symbol e.g. ('A', '*'):
    Where reference is ('T',), we output ref='T', alt=['A'], GT='1/1'
    """

    @property
    def n_elements(self):
        """Return number of elements provided by truth alignment.

        Synonymous with ploidy.
        """
        return 2

    @property
    def num_classes(self):
        """Return number of elements in output layer of neural network."""
        return len(self._decoding)

    @property
    def padding_vector(self):
        """Return the integer encoding used to denote a gap."""
        gap = [('*', '*')]
        return self._labels_to_encoded_labels(gap)[0]

    @property
    @functools.lru_cache(1)
    def _encoding(self):
        return {v: k for k, v in enumerate(
            self._unordered_label_combinations())}

    def _alignment_to_pairs(self, aln):
        """Convert `pysam.AlignedSegment` to aligned pairs."""
        seq = aln.query_sequence
        for qpos, rpos in aln.get_aligned_pairs():
            yield rpos, seq[qpos].upper() if qpos is not None else '*'

    def _labels_to_encoded_labels(self, labels):
        """Convert a list of labels to array of integer encoded labels."""
        return np.fromiter(
            (self._encoding[tuple(sorted(x))] for x in labels), dtype=int)

    def encoded_labels_to_training_vectors(self, enc_labels):
        """Convert integer (tuple) encoded labels to truth vectors.

        (e.g. sparse one-hot encoded truth vectors) that represent the truth
        for comparison with network output logits in
        e.g. metric(truth, pred) or loss(truth, pred)) functions.
        """
        return np.expand_dims(enc_labels, axis=1)  # sparse 1-hot

    def _prob_to_snp(
            self, outputs, positions, ref_name,
            ref_symbols, return_all=False):
        """Convert network output single locus to medaka.common.Variant."""
        argmax = outputs.argmax(axis=1)
        probs = outputs[np.arange(outputs.shape[0]), argmax]
        quals = self._phred(1 - probs)

        def _make_info(rs, p, c):
            if not self.verbose:
                return {}
            # helper for variant info formating
            rp = network_output[self._encoding[(rs, rs)]]
            info = {
                'ref_prob': self._pfmt(rp), 'prob': self._pfmt(p), 'call': c}
            return info

        results = list()
        data = zip(outputs, argmax, probs, quals, positions, ref_symbols)
        for network_output, amax, prob, qual, pos, ref_symbol in data:
            call = self._decoding[amax]

            # notes: we output only variants with one or more substitutions,
            #        and deletions are masked from the records. TODO: is this
            #        really the behaviour we want?
            if call == (ref_symbol, ref_symbol):  # is reference
                if return_all:
                    # return variant even though it is not a snp
                    q = self._pfmt(qual)
                    genotype = {'GT': '0/0', 'GQ': self._pfmt(qual, 0)}
                    info = _make_info(ref_symbol, prob, call)
                    results.append(medaka.vcf.Variant(
                        ref_name, pos, ref_symbol, alt='.', filt='PASS',
                        info=info, qual=q, genotype_data=genotype))
            else:  # is variant
                contains_deletion = '*' in call
                if not self._singleton(call):  # heterozygous
                    if not contains_deletion:
                        alt = [s for s in call if s != ref_symbol]
                        gt = '0/1' if len(alt) == 1 else '1/2'
                        q = self._pfmt(qual)
                        genotype = {'GT': gt, 'GQ': self._pfmt(qual, 0)}
                        info = _make_info(ref_symbol, prob, call)
                        results.append(medaka.vcf.Variant(
                            ref_name, pos, ref_symbol, alt, filt='PASS',
                            info=info, qual=q, genotype_data=genotype))
                    else:  # with deletion
                        # disallow (ref, *), and record (alt, *) as (alt, alt)
                        contains_nonref_nondel = len(
                            [s for s in call
                                if s != ref_symbol and s != '*']) > 0
                        if contains_nonref_nondel:
                            alt = [s for s in call if s != '*']
                            gt = '1/1'
                            q = self._pfmt(qual)
                            genotype = {'GT': gt, 'GQ': self._pfmt(qual, 0)}
                            info = _make_info(ref_symbol, prob, call)
                            results.append(medaka.vcf.Variant(
                                ref_name, pos, ref_symbol, alt, filt='PASS',
                                info=info, qual=q, genotype_data=genotype))
                elif not contains_deletion:  # homozygous (alt, alt)
                    alt = call[0]
                    q = self._pfmt(qual)
                    genotype = {'GT': '1/1', 'GQ': self._pfmt(qual, 0)}
                    info = _make_info(ref_symbol, prob, call)
                    results.append(medaka.vcf.Variant(
                        ref_name, pos, ref_symbol, alt, filt='PASS',
                        info=info, qual=q, genotype_data=genotype))
        return results

    @property
    def snp_metainfo(self):
        """Return meta data for use in `.vcf` header."""
        MI = medaka.vcf.MetaInfo
        m = [MI('FORMAT', 'GT', 1, 'String', 'Medaka genotype'),
             MI('FORMAT', 'GQ', 1, 'Float', 'Medaka genotype quality score')]
        if self.verbose:
            m.extend([
                MI('INFO', 'ref_prob', 1, 'Float',
                   'Medaka probability of reference'),
                MI('INFO', 'prob', 1, 'Float',
                   'Medaka probability of variant'),
                MI('INFO', 'call', 1, 'String', 'Medaka variant call'),
                ])

        return m


class DiploidZygosityLabelScheme(DiploidLabelScheme):
    """Diploid label scheme with an explicit heterozygosity element.

    The intermediate integer encoding is identical to the
    DiploidLabelScheme.

    The training vector encoding is multi-hot. For example, to encode
    ('A', C') at a locus, three elements of the training vectors representing
    'A', 'C' and 'is_heterogygous' are set to 1.

    ('A', 'C') -> 6 -> [0, 1, 1, 0, 0, 1]

    During inference, we predict the probability of each base existing
    independently, and also (independently) the probability of the locus
    being heterozygous. A binary crossentropy loss is appropriate.
    The output vector will not sum to 1.

    Importantly, during inference, there may be an is_heterozygous
    prediction > 0.5 which conflicts with independent predictions of which
    (how many) symbols are present (which symbols have independent
    probabilities > 0.5. In this case, the is_heterozygous call dominates
    and we return the top two most probable bases.

    SNP decoding:

    [0, 0.6, 0.7, 0, 0, 0.8] -> ('C', 'A') (in order of most probable)

    Depending on the reference, this may result in a homozygous or
    heterozygous SNP call.

    Where reference is ('A',), we output ref='A', alt=['C'], GT='0/1'
    Where reference is ('T',), we output ref='T' alt=['A', 'C'], GT='1/2'

    When predicion is homozygous e.g. ('A', 'A'):
    Where reference is ('T',), we output ref='T', alt=['A'], GT='1/1'

    When the prediction contains a gap symbol e.g. ('A', '*'):
    Where reference is ('T',), we output ref='T', alt=['A'], GT='1/1'

    If the is_heterozygous had been < 0.5, this would alter the call.
    [0, 0,6, 0.7, 0, 0, 0.4] -> ('C', 'C')
    """

    @property
    def n_elements(self):
        """Return number of elements provided by truth alignment.

        Synonymous with ploidy.
        """
        return 2

    @property
    def num_classes(self):
        """Return number of elements in output layer of neural network."""
        return len(self.symbols) + 1

    @property
    def padding_vector(self):
        """Return the integer encoding used to denote a gap."""
        gap = [('*', '*')]
        return self._labels_to_encoded_labels(gap)[0]

    def _is_het(self, x):
        return 1 if not self._singleton(x) else 0

    def encoded_labels_to_training_vectors(self, enc_labels):
        """Convert integer (tuple) encoded labels to truth vectors.

        (e.g. one-hot encoded truth vectors) that represent the truth
        for comparison with network output logits in
        e.g. metric(truth, pred) or loss(truth, pred)) functions.
        """
        # we hardcode the conversion between intermediate integer
        # encodings and the desired multi-hot vectors
        to_training_vector = {0:  np.array([1, 0, 0, 0, 0, 0]),  # ('*', '*')
                              1:  np.array([1, 1, 0, 0, 0, 1]),  # ('*', 'A')
                              2:  np.array([1, 0, 1, 0, 0, 1]),  # ('*', 'C')
                              3:  np.array([1, 0, 0, 1, 0, 1]),  # ('*', 'G')
                              4:  np.array([1, 0, 0, 0, 1, 1]),  # ('*', 'T')
                              5:  np.array([0, 1, 0, 0, 0, 0]),  # ('A', 'A')
                              6:  np.array([0, 1, 1, 0, 0, 1]),  # ('A', 'C')
                              7:  np.array([0, 1, 0, 1, 0, 1]),  # ('A', 'G')
                              8:  np.array([0, 1, 0, 0, 1, 1]),  # ('A', 'T')
                              9:  np.array([0, 0, 1, 0, 0, 0]),  # ('C', 'C')
                              10: np.array([0, 0, 1, 1, 0, 1]),  # ('C', 'G')
                              11: np.array([0, 0, 1, 0, 1, 1]),  # ('C', 'T')
                              12: np.array([0, 0, 0, 1, 0, 0]),  # ('G', 'G')
                              13: np.array([0, 0, 0, 1, 1, 1]),  # ('G', 'T')
                              14: np.array([0, 0, 0, 0, 1, 0])}  # ('T', 'T')

        vectors = np.stack([to_training_vector[x] for x in enc_labels])

        return vectors

    def _prob_to_snp(
            self, outputs, positions, ref_name, ref_symbols,
            return_all=False):
        """Convert network output to `medaka.common.Variant` s."""
        # where probability of heterozygosity > 0.5,
        # we make a heterozygous call, taking the two most
        # probable symbols.
        # TODO: optimise this!

        results = list()
        for network_output, pos, ref_symbol in zip(
                outputs, positions, ref_symbols):
            het_prob = network_output[-1]
            symbol_probs = network_output[:-1]

            is_het = het_prob > 0.5

            secondary_call, primary_call = (
                self._unitary_decoding[i][0]
                for i in np.argsort(symbol_probs)[-2:])
            secondary_prob, primary_prob = np.sort(symbol_probs)[-2:]

            ref_prob = symbol_probs[self._unitary_encoding[(ref_symbol,)]]

            if not is_het:
                call = tuple((primary_call, primary_call))
            else:
                call = tuple((primary_call, secondary_call))

            if self.verbose:
                info = {
                    'ref_prob': self._pfmt(ref_prob),
                    'primary_prob': self._pfmt(primary_prob),
                    'primary_call': primary_call,
                    'secondary_prob': self._pfmt(secondary_prob),
                    'secondary_call': secondary_call}
            else:
                info = {}

            # logical tests
            is_reference = call == (ref_symbol, ref_symbol)
            contains_deletion = '*' in call
            contains_nonref = len(
                [s for s in call if s != ref_symbol and s != '*']) > 0

            # homozygous
            if all((not is_reference,
                    not is_het,
                    not contains_deletion)):

                qual = self._phred(1 - primary_prob)
                alt = call[0]
                genotype = {'GT': '1/1', 'GQ': self._pfmt(qual, 0)}
                results.append(medaka.vcf.Variant(
                    ref_name, pos, ref_symbol, alt, filt='PASS',
                    info=info, qual=self._pfmt(qual), genotype_data=genotype))

            # heterozygous, no deletions
            elif all((is_het,
                      not contains_deletion)):

                err = 1 - 0.5 * (primary_prob + secondary_prob)
                qual = self._phred(err)
                alt = [s for s in call if s != ref_symbol]
                gt = '0/1' if len(alt) == 1 else '1/2'
                genotype = {'GT': gt, 'GQ': self._pfmt(qual, 0)}

                results.append(medaka.vcf.Variant(
                    ref_name, pos, ref_symbol, alt, filt='PASS',
                    info=info, qual=self._pfmt(qual), genotype_data=genotype))

            # heterozygous, one deletion
            elif all((is_het,
                      contains_nonref,
                      contains_deletion)):

                qual = self._phred(1 - primary_prob)
                alt = [s for s in call if s != '*']
                gt = '1/1'
                genotype = {'GT': gt, 'GQ': self._pfmt(qual, 0)}

                results.append(medaka.vcf.Variant(
                    ref_name, pos, ref_symbol, alt, filt='PASS',
                    info=info, qual=self._pfmt(qual), genotype_data=genotype))

            # no snp at this location
            else:
                if return_all:
                    qual = self._phred(ref_prob)
                    # return variant even though it is not a snp
                    genotype = {'GT': '0/0', 'GQ': self._pfmt(qual, 0)}
                    results.append(medaka.vcf.Variant(
                        ref_name, pos, ref_symbol, alt='.', filt='PASS',
                        info=info, qual=self._pfmt(qual),
                        genotype_data=genotype))
        return results


class RLELabelScheme(HaploidLabelScheme):
    """Class for RLE labelling schemes.

    The true length of the runs is encoded in the query scores.
    """

    def __init__(self, max_run=12):
        """Initialise class.

        :param max_run: Maximum run length (inclusive) to be
            considered. This will determine, amongst other things, the size
            of the labels created.

        """
        self.max_run = max_run

    @property
    def padding_vector(self):
        """Return the integer encoding used to denote a gap."""
        gap = [(('*', 1), )]
        return self._labels_to_encoded_labels(gap)[0]

    @property
    @functools.lru_cache(1)
    def _encoding(self):
        """Create a dictionary mapping from label tuple to integer."""
        encoding = dict()
        encoding[(('*', 1), )] = 0
        bases = [s for s in self.symbols if s != '*']
        lengths = range(1, self.max_run + 1)

        for i, (b, l) in enumerate(itertools.product(bases, lengths), 1):
            encoding[((b, l), )] = i

        return encoding

    def _alignment_to_pairs(self, aln):
        """Convert `pysam.AlignedSegment` to aligned pairs."""
        seq = aln.query_sequence
        # pysam gives back 0-based, but we skip the 0 for sanity
        run_lengths = aln.query_qualities
        for qpos, rpos in aln.get_aligned_pairs():
            qbase = seq[qpos] if qpos is not None else '*'
            # A deletion will have length 1
            qlen = run_lengths[qpos] if qpos is not None else 1
            # A larger run length that our maximum run length will be clipped
            qlen = min(qlen, self.max_run)
            yield rpos, (qbase, qlen)

    def _labels_to_encoded_labels(self, labels):
        """Convert a list of tuple labels to array of int encoded labels."""
        return np.fromiter(
            (self._encoding[x] for x in labels), dtype=int)

    def decode_consensus(self, sample):
        """Convert network output to consensus sequence by argmax decoding.

        :param sample: medaka.common.Sample

        :returns: str, consensus sequence
        """
        # property access is slow
        decode = self._decoding
        # most probable class
        mp = np.argmax(sample.label_probs, -1)

        def _get_substrings():
            for x in mp:
                ((base, run), ) = decode[x]
                if base != '*':
                    yield base * run
        seq = ''.join(_get_substrings())

        return seq

    def _prob_to_snp(self):
        """Convert network output to medaka.common.Variant."""
        raise NotImplementedError
