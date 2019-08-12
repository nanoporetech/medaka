import abc
import collections
import inspect
import itertools
import operator

import numpy as np
import pysam

import medaka.common
import medaka.datastore
import medaka.rle
import medaka.vcf

variant_decoders = {}

# when calculating qscores, we can get inf, but cannot write that to vcf.
__max_qscore__ = 9999


def register_variant_decoder(clsname, cls):
    """Register a variant decoder."""
    if clsname != 'BaseDecoder':
        variant_decoders[clsname] = cls


class DecoderRegistrar(type):
    """metaclass for registering variant decoders."""

    def __new__(cls, clsname, bases, attrs):
        newclass = super(DecoderRegistrar, cls).__new__(cls, clsname, bases, attrs)
        register_variant_decoder(clsname, newclass)  # register variant decoders
        return newclass


def yield_trimmed_consensus_chunks(sample_gen):
    """Yield adjacent consensus chunks trimmed to remove overlap.

    :param sample_gen: `generator of medaka.common.Sample` obj

    :yields: (`Sample`, bool is_last_in_contig)
    """

    logger = medaka.common.get_named_logger('TrimOverlap')
    s1 = next(sample_gen)
    start_1_ind = None  # don't trim beginning of s1
    start_2_ind = None  # initialise in case we just have 1 sample
    for s2 in itertools.chain(sample_gen, (None,)):

        s1_name = 'Unknown' if s1 is None else s1.name
        s2_name = 'Unknown' if s2 is None else s2.name

        is_last_in_contig = False
        if s2 is None:  # s1 is last chunk
            end_1_ind = None  # go to the end of s1
            is_last_in_contig = True
        else:
            rel = medaka.common.Sample.relative_position(s1, s2)
            # skip s2 if it is contained within s1 (i.e. s2 is redundant)
            if rel is medaka.common.Sample.Relationship.s2_within_s1:
                logger.info('{} is contained within {}, skipping.'.format(
                    s2_name, s1_name
                ))
                continue
            elif rel is medaka.common.Sample.Relationship.forward_overlap:
                end_1_ind, start_2_ind = medaka.common.Sample.overlap_indices(s1, s2)
            elif rel is medaka.common.Sample.Relationship.forward_gapped:
                is_last_in_contig = True
                end_1_ind, start_2_ind = (None, None)
                msg = '{} and {} cannot be concatenated as there is no overlap and they do not abut.'
                logger.info(msg.format(s1_name, s2_name))
            else:
                raise RuntimeError('Unexpected sample relationship {} between {} and {}'.format(repr(rel), s1.name, s2.name))

        yield s1.slice(slice(start_1_ind, end_1_ind)), is_last_in_contig

        s1 = s2
        start_1_ind = start_2_ind


def join_chunked_variants(gen, ref_seq_encoded, gap_encoded):
    """Process stream of trimmed consensus chunks ensuring a variant is
    not split accross multiple chunks.

    :param gen: stream of (`medaka.common.Sample` obj, bool is_last_in_contig)
    :param ref_seq_encoded: `np.ndarray` of the ref sequence (i.e. only major positions) encoded as medaka labels.
    :param gap_encoded: int, encoding of gaps.

    :yields: `medaka.common.Sample`
    """
    queue = []
    for s, is_last_in_contig in gen:
        if is_last_in_contig:
            # there are no further samples in this contig, so all remaining variants must
            # be in this contig or in the queue, so process and empty queue
            queue.append(s)
            yield medaka.common.Sample.from_samples(queue)
            queue = []
            continue

        # find last non-variant call, i.e. a major position which has no minor
        # positions and is the same as the ref

        call = np.argmax(s.label_probs, axis=1)

        # get encoded ref seq containing insertions as gaps
        get_encoded = lambda p: ref_seq_encoded[p['major']] if p['minor'] == 0 else gap_encoded
        s_ref_encoded = np.fromiter((get_encoded(p) for p in s.positions), int, count=len(s.positions))
        is_diff = call != s_ref_encoded
        # don't count a minor position called as gap as a match
        both_gap = np.logical_and(call == gap_encoded, s_ref_encoded == gap_encoded)
        is_var = np.logical_or(is_diff, both_gap)

        # if the entire chunk is a variant, queue it up
        if np.all(is_var):
            queue.append(s)
            continue

        major_inds = np.where(s.positions['minor'] == 0)
        major_pos = s.positions['major'][major_inds]
        major_call = call[major_inds]
        is_diff = major_call != ref_seq_encoded[(major_pos,)]
        for offset, is_diff_i in enumerate(is_diff[::-1]):
            # Even if the last position does not look like a variant, we can't know
            # if the first pos in next chunk is not an insertion => we need
            # last non-variant position before last position.
            if not is_diff_i:
                break
        last_non_var_pos = major_pos[-1 - offset]
        last_non_var_start = np.searchsorted(s.positions['major'], last_non_var_pos, side='left')
        to_yield = s.slice(slice(None, last_non_var_start))
        to_queue = s.slice(slice(last_non_var_start, None))
        yield medaka.common.Sample.from_samples(queue + [to_yield])
        queue = [to_queue]
    if len(queue) > 0:
        raise ValueError('Reached end of generator at {} without is_last_in_contig being True'.format(s.name))


class Meta(abc.ABC, DecoderRegistrar):
    pass


class BaseDecoder(metaclass=Meta):
    """Abstract base class to decode diploid SNPs"""

    def __init__(self, meta, ref_vcf=None):
        """
        :param meta: dict containing 'medaka_label_decoding', 'medaka_feature_decoding' and 'medaka_multi_label'
        :param ref_vcf: input vcf to force evaluation only at these loci, even if we would
                not otherwise call a SNP there.
        """
        self.logger = medaka.common.get_named_logger(self.__class__.__name__)

        self.label_decoding = meta['medaka_label_decoding']
        self.augmented_label_encoding = self._label_encoding_with_extra_bases()

        # to be able to add feature date to the VCF info fields, we need a label
        # for each feature row, which we can get from
        # meta['medaka_feature_decoding']
        # which is a tuple (dtype, is_rev, base, base_len)
        # if we have a del, base is None.
        # Wrap in a try except as it might fail with any non-default features
        try:
            fmt_feat = lambda x: '{}{}{}'.format(x[0], 'rev' if x[1] else 'fwd', x[3] * (x[2] if x[2] is not None else '-'))
            self.feature_row_names = [fmt_feat(x) for x in meta['medaka_feature_decoding']]
        except:
            self.logger.debug('Could not create feature info for VCF')
            self.feature_row_names = None

        self.logger.debug("Label decoding is:\n{}".format('\n'.join(
            '{}: {}'.format(i, x) for i, x in enumerate(self.label_decoding)
        )))

        self.ref_vcf = medaka.vcf.VCFReader(ref_vcf) if ref_vcf is not None else None
        if self.ref_vcf is not None:
            self.ref_loci = collections.defaultdict(set)
            for v in self.ref_vcf.fetch():
                self.ref_loci[v.chrom].add(v.pos)

        # if the model was trained with the multi_label option
        # label probabilities are independent and do not sum up to 1, so we need
        # to calculate qualties differently
        if 'medaka_multi_label' in meta:
            self.multi_label = meta['medaka_multi_label']
        else:
            self.multi_label = False
        self.logger.info('Decoding with multi_label set to {}'.format(self.multi_label))


    def _label_encoding_with_extra_bases(self):
        """Get map from label decodings (eg. 'A') to int labels with extra bases so we can encode e.g. reference N's.
        """
        # include extra bases so we can encode e.g. N's in reference and create
        # our own label encoding rather than use that in the meta, as
        # meta['medaka_label_encoding'] could map from pairs of raw labels to
        # one or more encoded labels while we want a simple map from
        # label_decoding (eg. 'A') to int label classification.

        # This is used to encode the reference with the (augmented) medaka
        # labels once per reference contig, and this encoded ref is compared to
        # the consensus calls in each sample chunk to check for variants. This
        # is more efficient than having to convert every sample chunk to str
        # arrays to compare with a reference string array
        return {label: ind for ind, label in enumerate(
                itertools.chain(self.label_decoding, medaka.common._extra_bases_))}


    @abc.abstractmethod
    def decode_variants(self, samples, ref_seq):
        """Abstract method to decode variants.

        :param samples: stream of `medaka.common.Sample` objects for a given contig.
        :param ref_seq: reference sequence.

        :returns: sorted list of `medaka.vcf.Variant` objects.
        """


class SNPDecoder(BaseDecoder):
    """Class to decode SNPs using a threshold"""


    def __init__(self, meta, threshold=0.04, ref_vcf=None):
        """
        :param meta: dict containing 'medaka_label_decoding', 'medaka_feature_decoding' and 'medaka_multi_label'
        :param threshold: threshold below which a secondary call (which would make
                for a heterozygous call) is deemed insignificant.
        :param ref_vcf: input vcf to force evaluation only at these loci, even if we would
                not otherwise call a SNP there.
        """
        super().__init__(meta, ref_vcf=ref_vcf)
        self.threshold = threshold


    def _get_ref_variant(self, ref_name, pos):
        """Check if there is a variant at the given ref_name and pos in the input reference vcf.

        :param ref_name: str, genomic contig.
        :param pos: int, genomic position.

        :returns: dict containing 'ref_alt' and 'ref_gt' fields which are both 'na' if there is no variant at that loci.
        """
        ref_variants = list(self.ref_vcf.fetch(ref_name=ref_name, start=pos, end=pos+1))
        assert len(ref_variants) < 2
        if len(ref_variants) == 0:
            ref_info = {'ref_alt': 'na', 'ref_gt': 'na'}
        else:
            ref_info = {'ref_alt': ','.join(ref_variants[0].alt),
                        'ref_gt': ref_variants[0].sample_dict['GT']}
        return ref_info


    def _process_sample(self, s, ref_seq):

        #TODO one could refactor this to use specialised VariantSample classes that
        # inherit from medaka.common.Sample and implement methods like this

        data = self._process_sample_probs(s, ref_seq)
        homozygous_snp_inds, heterozygous_snp_inds = self._get_hom_het_snp_inds_(data)

        return homozygous_snp_inds, heterozygous_snp_inds, data


    def _preprocess_sample_probs(self, s):
        """Preprocess sample label probabilities
        :param s: `medaka.common.Sample` instance
        :returns: dict containing sample fields"""
        #TODO one could refactor this to use specialised VariantSample classes that
        # inherit from medaka.common.Sample and implement methods like this
        data = s._asdict()
        return data


    def _process_sample_probs(self, s, ref_seq):

        #TODO one could refactor this to use specialised VariantSample classes that
        # inherit from medaka.common.Sample and implement methods like this

        data = self._preprocess_sample_probs(s)

        # for snps, we discard minor positions (insertions)
        data = self._discard_minor(data)

        data['ref_seq_encoded'] = np.fromiter(
            (self.augmented_label_encoding[ref_seq[i]] for i in data['positions']['major']),
            int, count=len(data['positions'])
        )
        sorted_prob_inds = np.argsort(data['label_probs'], -1)
        sorted_probs = np.take_along_axis(data['label_probs'], sorted_prob_inds, axis=-1)
        data['primary_labels'] = sorted_prob_inds[:, -1]
        data['secondary_labels'] = sorted_prob_inds[:, -2]
        data['primary_probs'] = sorted_probs[:, -1]
        data['secondary_probs'] = sorted_probs[:, -2]

        # we'll want to skip positions where ref is not a label (*ATCG)
        data['is_ref_valid_label'] = np.isin(data['ref_seq_encoded'], np.arange(len(self.label_decoding)))

        data['is_primary_diff_to_ref'] = np.not_equal(data['primary_labels'], data['ref_seq_encoded'])
        data['is_secondary_diff_to_ref'] = np.not_equal(data['secondary_labels'], data['ref_seq_encoded'])
        data['is_primary_not_del'] = data['primary_labels'] != self.augmented_label_encoding[medaka.common._gap_]
        data['is_secondary_del'] = data['secondary_labels'] == self.augmented_label_encoding[medaka.common._gap_]

        return data


    def _discard_minor(self, data):
        """Discard minor positions (insertions) from processed sample data (modifies data in-place)

        :param data: dict of `medaka.common.Sample._fields` containing at least sample positions array.
        :returns data: dict in which all arrays have had minor positions removed. Non-array fields will be unchanged.
        """

        #TODO one could refactor this to use specialised VariantSample classes that
        # inherit from medaka.common.Sample and implement methods like this

        major_inds = np.where(data['positions']['minor'] == 0)

        for field in data.keys():
            if data[field] is not None and isinstance(data[field], np.ndarray):
                data[field] = data[field][major_inds]
        return data


    def _get_hom_het_snp_inds_(self, data):
        """Get indices of homozygous and heterozygous SNPs using a threshold"""

        #TODO one could refactor this to use specialised VariantSample classes that
        # inherit from medaka.common.Sample and implement methods like this

        # for homozygous SNP max_prob_label not in {ref, del} and
        # (2nd_max_prob < threshold or 2nd_max_prob_label is del)
        # for heterozygous SNP 2nd_max_prob > threshold and del not in {max_prob_label, 2nd_prob_label}
        # this catches both SNPs where the genotype contains the
        # reference, and where both copies are mutated.

        is_secondary_prob_lt_thresh = data['secondary_probs'] < self.threshold
        is_not_secondary_call = np.logical_or(data['is_secondary_del'], is_secondary_prob_lt_thresh)

        is_homozygous_snp = np.logical_and.reduce([
            data['is_primary_diff_to_ref'], data['is_primary_not_del'],
            is_not_secondary_call, data['is_ref_valid_label']
        ])
        homozygous_snp_inds = np.where(is_homozygous_snp)[0]  # to get indices in 0th axis

        # heterozygous SNPs
        is_secondary_prob_ge_thresh = np.logical_not(is_secondary_prob_lt_thresh)
        is_secondary_not_del = np.logical_not(data['is_secondary_del'])
        is_heterozygous_snp = np.logical_and.reduce([
            is_secondary_prob_ge_thresh, is_secondary_not_del,
            data['is_primary_not_del'], data['is_ref_valid_label']
        ])
        heterozygous_snp_inds = np.where(is_heterozygous_snp)[0]  # to get indices in 0th axis

        return homozygous_snp_inds, heterozygous_snp_inds


    def _get_variant_data(self, data, i):
        pos = data['positions']['major'][i]
        ref_base_encoded = data['ref_seq_encoded'][i]
        info = {'ref_prob': data['label_probs'][i][ref_base_encoded],
                'primary_prob': data['primary_probs'][i],
                'primary_label': self.label_decoding[data['primary_labels'][i]],
                'secondary_prob': data['secondary_probs'][i],
                'secondary_label': self.label_decoding[data['secondary_labels'][i]],
                }
        if self.ref_vcf is not None:
            ref_info = self._get_ref_variant(data['ref_name'], pos)
            info.update(ref_info)
        if data['features'] is not None and self.feature_row_names is not None:
            info.update(dict(zip(self.feature_row_names, data['features'][i])))

        ref = self.label_decoding[ref_base_encoded]

        return pos, ref, info


    def _yield_homozygous_snps(self, homozygous_snp_inds, data, called_loci):

        for i in homozygous_snp_inds:
            pos, ref, info = self._get_variant_data(data, i)
            alt = self.label_decoding[data['primary_labels'][i]]

            qual = -10 * np.log10(1 - data['primary_probs'][i])
            sample = {'GT': '1/1', 'GQ': qual,}
            called_loci.add(pos)
            yield medaka.vcf.Variant(data['ref_name'], pos, ref, alt,
                                     filter='PASS', info=info, qual=qual,
                                     sample_dict=sample)


    def _yield_heterozygous_snps(self, heterozygous_snp_inds, data, called_loci):

        for i in heterozygous_snp_inds:
            pos, ref, info = self._get_variant_data(data, i)
            if self.multi_label:
                # label probabilities are independent and not normalised
                # so calculate error for each haplotype as 1-p(call)
                # average the error so that a similar threshold can be used
                # for homozygous and heterozygous SNPs (were one to sum, if
                # one assumes that the network can be equally confident at
                # calling homozygous and heterozygous loci, the heterozygous
                # quals would always have double the uncertainty.
                err = 1 - 0.5 * (data['primary_probs'][i] + data['secondary_probs'][i])
            else:
                # ideally one will have almost 0.5 probability in each
                # label, with any residual error in the other labels
                # (as compared to homozygous loci where one would have
                # almost unit probability in a single label
                err = 1 - (data['primary_probs'][i] + data['secondary_probs'][i])

            qual = -10 * np.log10(err)
            ref_base_encoded = data['ref_seq_encoded'][i]
            alt = [self.label_decoding[l] for l in
                   (data['primary_labels'][i], data['secondary_labels'][i]) if l != ref_base_encoded]
            gt = '0/1' if len(alt) == 1 else '1/2'  # / => unphased
            sample = {'GT': gt, 'GQ': qual,}
            called_loci.add(pos)
            yield medaka.vcf.Variant(data['ref_name'], pos, ref, alt=alt,
                                     filter='PASS', info=info, qual=qual,
                                     sample_dict=sample)


    def _yield_reference_snps(self, missing_snp_inds, data, called_loci):

        for i in missing_snp_inds:
            pos = data['positions']['major'][i]
            ref_base_encoded = data['ref_seq_encoded'][i]
            info = {'ref_prob': data['label_probs'][i][ref_base_encoded],
                    'primary_prob': data['primary_probs'][i],
                    'primary_label': self.label_decoding[data['primary_labels'][i]],
                    'secondary_prob': data['secondary_probs'][i],
                    'secondary_label': self.label_decoding[data['secondary_labels'][i]],
                    }
            ref_info = self._get_ref_variant(ref_name, data['positions']['major'][i])
            info.update(ref_info)
            if data['features'] is not None:
                info.update(dict(zip(self.feature_row_names, data['features'][i])))
            qual = -10 * np.log10(1 - data['primary_probs'][i])
            sample = {'GT': 0, 'GQ': qual,}
            yield medaka.vcf.Variant(ref_name, pos,
                                        self.label_decoding[ref_base_encoded],
                                        alt='.', filter='PASS', info=info,
                                        qual=qual, sample_dict=sample)


    def decode_variants(self, samples, ref_seq):
        """Decode SNPs.

        :param samples: stream of `medaka.common.Sample` objects for a given contig.
        :param ref_seq: reference sequence.

        :yields: `medaka.vcf.Variant` objects.
        """

        # For SNPS, we assume we just want the label probabilities at major positions
        # We need to handle boundary effects simply to make sure we take the best
        # probability (which is furthest from the boundary).
        called_loci = set()
        ref_name = None
        for s, _ in yield_trimmed_consensus_chunks(samples):
            ref_name = s.ref_name if ref_name is None else ref_name
            assert ref_name == s.ref_name

            homozygous_snp_inds, heterozygous_snp_inds, data = self._process_sample(s, ref_seq)
            yield from self._yield_homozygous_snps(homozygous_snp_inds, data, called_loci)
            yield from self._yield_heterozygous_snps(heterozygous_snp_inds, data, called_loci)

            if self.ref_vcf is not None:
                # if we provided a vcf, check which positions are missing
                major_pos = data['positions']['major']
                missing_loci = self.ref_loci[ref_name] - called_loci
                missing_loci_in_chunk = missing_loci.intersection(major_pos)
                missing_loci_in_chunk = np.fromiter(missing_loci_in_chunk, int,
                                                    count=len(missing_loci_in_chunk))
                is_missing = np.isin(major_pos, missing_loci_in_chunk)
                missing_snp_inds = np.where(is_missing)[0]
                yield from self._yield_reference_snps(missing_snp_inds, data)


    @property
    def meta_info(self):
        m = [
            medaka.vcf.MetaInfo('INFO', 'ref_prob', 1, 'Float', 'Medaka label probability for the reference allele'),
            medaka.vcf.MetaInfo('INFO', 'primary_prob', 1, 'String', 'Medaka label probability of primary call'),
            medaka.vcf.MetaInfo('INFO', 'primary_label', 1, 'String', 'Medaka label of primary call'),
            medaka.vcf.MetaInfo('INFO', 'secondary_prob', 1, 'Float', 'Medaka label probability of secondary call'),
            medaka.vcf.MetaInfo('INFO', 'secondary_label', 1, 'String', 'Medaka label of secondary call'),
            medaka.vcf.MetaInfo('FORMAT', 'GT', 1, 'String', 'Medaka genotype'),
            medaka.vcf.MetaInfo('FORMAT', 'GQ', 1, 'Float', 'Medaka genotype quality score'),
        ]
        return m


class FactoredBaseZygositySNPDecoder(SNPDecoder):
    """Class to decode SNPs and deletions using factored base and zygosity labels."""
    _heterozygous_label_decoding_ = 'heterozygous'
    _homozygous_label_decoding_ = 'homozygous'

    def __init__(self, meta, ref_vcf=None):
        """
        :param meta: dict containing 'medaka_label_decoding', 'medaka_feature_decoding' and 'medaka_multi_label'
        :param ref_vcf: input vcf to force evaluation only at these loci, even if we would
                not otherwise call a SNP there.
        """
        super().__init__(meta, ref_vcf=ref_vcf)

        # tuple of labels e.g. ('*', 'A', 'AA', 'C', 'G', 'T', 'homozygous', 'heterozygous')
        decoding = meta['medaka_label_decoding']
        # indices of zygosity labels
        self._zygosity_labels_ = (decoding.index(self._homozygous_label_decoding_),
                                  decoding.index(self._heterozygous_label_decoding_))
        # indices of base labels
        self._base_labels_ = tuple(set(range(len(decoding))) - set(self._zygosity_labels_))

        # create label_decoding and label_encoding only for base labels
        self.label_decoding = operator.itemgetter(*self._base_labels_)(decoding)
        self.augmented_label_encoding = self._label_encoding_with_extra_bases()
        # create zygosity label decoding
        self.label_decoding_zygosity = operator.itemgetter(*self._zygosity_labels_)(decoding)


    def _preprocess_sample_probs(self, s):
        """Preprocess sample label probabilities to split out zygosity labels from base labels"""

        #TODO one could refactor this to use specialised VariantSample classes that
        # inherit from medaka.common.Sample and implement methods like this

        data = s._asdict()

        data['zygosity_probs'] = s.label_probs[:, self._zygosity_labels_]
        # replace label_probs with just base labels
        data['label_probs'] = s.label_probs[:, self._base_labels_]

        # add a is_heterozygous field
        heterozygous_label = self.label_decoding_zygosity.index(self._heterozygous_label_decoding_)
        predicted_zygosity_label = np.argmax(data['zygosity_probs'], axis=1)
        data['is_heterozygous'] = predicted_zygosity_label == heterozygous_label

        return data


    def _get_hom_het_snp_inds_(self, data):
        """Get indices of homozygous and heterozygous SNPs using zygosity labels"""

        #TODO one could refactor this to use specialised VariantSample classes that
        # inherit from medaka.common.Sample and implement methods like this

        is_heterozygous_snp = np.logical_and.reduce([
            data['is_ref_valid_label'],
            data['is_heterozygous'],
            np.logical_or(data['is_primary_diff_to_ref'], data['is_secondary_diff_to_ref']),
            np.logical_and(data['is_primary_not_del'], np.logical_not(data['is_secondary_del'])),
        ])

        is_homozygous_snp = np.logical_and.reduce([
            data['is_ref_valid_label'],
            data['is_primary_not_del'],
            np.logical_not(data['is_heterozygous']),
            data['is_primary_diff_to_ref'],
        ])

        heterozygous_snp_inds = np.where(is_heterozygous_snp)[0]
        homozygous_snp_inds = np.where(is_homozygous_snp)[0]

        return homozygous_snp_inds, heterozygous_snp_inds


def parse_regions(regions):
    """Parse region strings, returning a list of ref_names.

    :param regions: iterable of region strings.

    :returns: tuple of str of ref_names.
    """
    #TODO: respect entire region specification
    ref_names = list()
    for region in (medaka.common.Region.from_string(r) for r in regions):
        if region.start is not None or region.end is not None:
            logger.warning("Ignoring start:end for '{}'.".format(region))
        if region.ref_name not in ref_names:
            ref_names.append(region.ref_name)
    return tuple(ref_names)


class HaploidVariantDecoder(SNPDecoder):
    """Class to decode haploid variants (snps and indels)."""

    def decode_variants(self, samples, ref_seq):
        """Decode haploid variants.

        Multi-base or otherwise complex variants will be decoded in their entirety.

        :param samples: stream of `medaka.common.Sample` objects for a given contig.
        :param ref_seq: reference sequence.

        :returns: sorted list of `medaka.vcf.Variant` objects.
        """

        ref_seq_encoded = np.fromiter((self.augmented_label_encoding[b] for b in ref_seq), int, count=len(ref_seq))

        gap_encoded = self.augmented_label_encoding[medaka.common._gap_]

        for s in join_chunked_variants(yield_trimmed_consensus_chunks(samples), ref_seq_encoded, gap_encoded):

            #TODO parallelize? one could dispatch each sample s to a different process..

            assert s.positions['minor'][0] == 0
            pred_labels = np.argmax(s.label_probs, -1)

            # get encoded ref seq containing insertions as gaps
            get_label = lambda p: ref_seq_encoded[p['major']] if p['minor'] == 0 else gap_encoded
            ref_labels = np.fromiter((get_label(p) for p in s.positions), int, count=len(s.positions))

            # find variants by looking for runs of labels which differ. If both
            # labels are gap, we don't want to consider this a match so as to
            # avoid splitting up long insertions if there happens to be a gap
            # label in ref and pred.
            is_diff = pred_labels != ref_labels
            both_gap = np.logical_and(pred_labels == gap_encoded, ref_labels == gap_encoded)
            is_var = np.logical_or(is_diff, both_gap)
            runs = medaka.rle.rle(is_var)
            var_runs = runs[np.where(runs['value'])]
            # TODO: optionally load in var_runs from an input vcf -
            # we can then decode quals for what medaka had, the ref and input alt

            for var_run in var_runs:
                start_i = var_run['start']
                stop_i = start_i + var_run['length']
                # if call or ref starts with gap, pad left (or right if we are
                # at start of genome).
                pad_right = False
                while ((not pad_right and (ref_labels[start_i] == gap_encoded or pred_labels[start_i] == gap_encoded))
                        or
                       (pad_right and (ref_labels[stop_i] == gap_encoded or pred_labels[stop_i] == gap_encoded))
                      ):
                    if tuple(s.positions[start_i]) == (0, 0) or pad_right:
                        stop_i += 1
                        # avoid getting stuck in loop if there is a run of dels at start of ref
                        pad_right = True
                    else:
                        assert start_i != 0
                        start_i -= 1

                var_ref_labels = ref_labels[start_i:stop_i]
                if var_ref_labels.max() >= len(self.label_decoding):
                    # Ref seq contains non-label bases (e.g. N's)
                    continue
                    self.logger.debug('Ref seq contains non-label bases, skipping {}'.format())

                var_ref_decoded = ''.join([self.label_decoding[x] for x in var_ref_labels])
                var_ref_seq = var_ref_decoded.replace(medaka.common._gap_, '')
                var_pred_labels = pred_labels[start_i:stop_i]
                var_pred_decoded = ''.join([self.label_decoding[x] for x in var_pred_labels])
                var_pred_seq = var_pred_decoded.replace(medaka.common._gap_, '')

                if var_ref_seq == var_pred_seq:
                    # this is not a variant
                    continue
                probs = s.label_probs[start_i:stop_i]
                var_probs = np.fromiter((probs[i,j] for i, j in enumerate(var_pred_labels)),
                                        dtype=probs.dtype, count=len(var_pred_labels))
                ref_probs = np.fromiter((probs[i,j] for i, j in enumerate(var_ref_labels)),
                                        dtype=probs.dtype, count=len(var_ref_labels))
                var_quals = -10 * np.log10(1 - var_probs)
                ref_quals = -10 * np.log10(1 - ref_probs)
                qfmt = lambda q : '{:.3f}'.format(q if q < np.inf else __max_qscore__)
                info = {'ref_seq': var_ref_decoded,
                        'pred_seq': var_pred_decoded,
                        'pred_qs': ','.join((qfmt(q) for q in var_quals)),
                        'ref_qs': ','.join((qfmt(q) for q in ref_quals)),
                        'pred_q': qfmt(var_quals.sum()),
                        'ref_q': qfmt(ref_quals.sum()),
                        'n_cols': len(var_quals),
                        }

                qual = qfmt(var_quals.sum() - ref_quals.sum())  # log likelihood ratio
                # TODO: test if normalising by n_cols_diff makes it easier to
                # pick a single threshold that will work equally well for short
                # and long variants (at the moment, longer variants have a
                # higher qual)
                # n_cols_diff = sum(p != r for p, r in zip(info['pred_seq'], info['ref_seq']))
                # qual /= n_cols_diff
                sample = {'GT': '1/1', 'GQ': qual,}

                self.logger.debug('Yielding variant {}:{}'.format(s.ref_name, s.positions[start_i]))
                yield medaka.vcf.Variant(s.ref_name, s.positions['major'][start_i],
                                         var_ref_seq, alt=var_pred_seq,
                                         filter='PASS', info=info, qual=qual,
                                         sample_dict=sample).trim()


    @property
    def meta_info(self):
        m = [
            medaka.vcf.MetaInfo('INFO', 'ref_seq', 1, 'String', 'Medaka reference label for each pileup column'),
            medaka.vcf.MetaInfo('INFO', 'pred_seq', 1, 'String', 'Medaka predicted-label for each pileup column'),
            medaka.vcf.MetaInfo('INFO', 'pred_qs', '.', 'Float', 'Medaka predicted-label quality score for each pileup column'),
            medaka.vcf.MetaInfo('INFO', 'ref_qs', '.', 'Float', 'Medaka reference-label quality score for each pileup column'),
            medaka.vcf.MetaInfo('INFO', 'pred_q', 1, 'Float', 'Medaka total predicted-label quality score'),
            medaka.vcf.MetaInfo('INFO', 'ref_q', 1, 'Float', 'Medaka total reference-label quality score'),
            medaka.vcf.MetaInfo('INFO', 'n_cols', 1, 'Integer', 'Number of Medaka pileup columns in the variant call'),
            medaka.vcf.MetaInfo('FORMAT', 'GT', 1, 'String', 'Medaka genotype.'),
            medaka.vcf.MetaInfo('FORMAT', 'GQ', 1, 'Float', 'Medaka genotype quality score'),
        ]
        return m


def variants_from_hdf(args):
    """Entry point for variant calling from HDF5 consensus probability files."""
    index = medaka.datastore.DataIndex(args.inputs)
    if args.regions is not None:
        ref_names = parse_regions(args.regions)
    else:
        ref_names = medaka.common.loose_version_sort(index.index.keys())

    decoder_cls = variant_decoders[args.decoder]

    opts = inspect.signature(decoder_cls.__init__).parameters.keys()
    opts = {k:getattr(args, k) for k in opts if k not in  {'self', 'meta', 'ref_vcf'}}

    decoder = decoder_cls(index.meta, ref_vcf=args.ref_vcf, **opts)

    logger = medaka.common.get_named_logger('Variants')

    with medaka.vcf.VCFWriter(args.output, 'w', version='4.1', contigs=ref_names,
                              meta_info=decoder.meta_info) as vcf_writer:
        for ref_name in ref_names:
            logger.info("Processing {}.".format(ref_name))
            ref_seq = pysam.FastaFile(args.ref_fasta).fetch(reference=ref_name)
            samples = index.yield_from_feature_files([ref_name])
            variants = decoder.decode_variants(samples, ref_seq)
            vcf_writer.write_variants(variants, sort=True)
