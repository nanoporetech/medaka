from collections import defaultdict
import inspect
import itertools

import numpy as np
import pysam

import medaka.common
import medaka.datastore
import medaka.vcf

variant_decoders = {}

# when calculating qscores, we can get inf, but cannot write that to vcf.
__max_qscore__ = 9999


def register_variant_decoder(clsname, cls):
    """Register a variant decoder."""
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


class SNPDecoder(object, metaclass=DecoderRegistrar):
    """Class to decode diploid SNPs and deletions using a threshold"""

    def __init__(self, meta, threshold=0.04, ref_vcf=None):
        """Class to decode SNPs and deletions from hfds using a threshold

        :param meta: dict containing 'medaka_label_decoding' and 'medaka_feature_decoding'i
        :param threshold: threshold below which a secondary call (which would make
                for a heterozygous call) is deemed insignificant.
        :param ref_vcf: input vcf to force evaluation only at these loci, even if we would
                not otherwise call a SNP there.
        """
        self.logger = medaka.common.get_named_logger('SNPs')


        self.label_decoding = meta['medaka_label_decoding']
        # include extra bases so we can encode e.g. N's in reference
        self.label_encoding = {label: ind for ind, label in enumerate(self.label_decoding + list(medaka.common._extra_bases_))}

        # to be able to add feature date to the VCF info fields, we need a label
        # for each feature row, which we can get from
        # meta['medaka_feature_decoding']
        # which is a tuple (dtype, is_rev, base, base_len)
        # if we have a del, base is None.
        fmt_feat = lambda x: '{}{}{}'.format(x[0], 'rev' if x[1] else 'fwd', x[3] * (x[2] if x[2] is not None else '-'))
        self.feature_row_names = [fmt_feat(x) for x in meta['medaka_feature_decoding']]

        self.logger.debug("Label decoding is:\n{}".format('\n'.join(
            '{}: {}'.format(i, x) for i, x in enumerate(self.label_decoding)
        )))

        self.ref_vcf = medaka.vcf.VCFReader(ref_vcf) if ref_vcf is not None else None
        if self.ref_vcf is not None:
            self.ref_loci = defaultdict(set)
            for v in self.ref_vcf.fetch():
                self.ref_loci[v.chrom].add(v.pos)


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


    def decode_variants(self, samples, ref_seq, threshold=0.04):
        """Decode SNPs using a threshold.

        :param samples: stream of `medaka.common.Sample` objects for a given contig.
        :param ref_seq: reference sequence.
        :param threshold: threshold below which a secondary call (which would make
                for a heterozygous call) is deemed insignificant.

        :returns: sorted list of `medaka.vcf.Variant` objects.
        """
        # For SNPS, we assume we just want the label probabilities at major positions
        # We need to handle boundary effects simply to make sure we take the best
        # probability (which is furthest from the boundary).
        called_loci = set()
        ref_name = None
        for s, _ in yield_trimmed_consensus_chunks(samples):
            ref_name = s.ref_name if ref_name is None else ref_name
            assert ref_name == s.ref_name
            pos = s.positions
            probs = s.label_probs
            # discard minor positions (insertions)
            major_inds = np.where(pos['minor'] == 0)
            major_pos = pos[major_inds]['major']
            major_probs = probs[major_inds]
            major_feat = s.features[major_inds] if s.features is not None else None

            # for homozygous SNP max_prob_label not in {ref, del} and
            # (2nd_max_prob < threshold or 2nd_max_prob_label is del)
            # for heterozygous SNP 2nd_max_prob > threshold and del not in {max_prob_label, 2nd_prob_label}
            # this catches both SNPs where the genotype contains the
            # reference, and where both copies are mutated.

            ref_seq_encoded = np.fromiter((self.label_encoding[ref_seq[i]] for i in major_pos), int, count=len(major_pos))
            sorted_prob_inds = np.argsort(major_probs, -1)
            sorted_probs = np.take_along_axis(major_probs, sorted_prob_inds, axis=-1)
            primary_labels = sorted_prob_inds[:, -1]
            secondary_labels = sorted_prob_inds[:, -2]
            primary_probs = sorted_probs[:, -1]
            secondary_probs = sorted_probs[:, -2]
            # skip positions where ref is not a label (ATGC)
            is_ref_valid_label = np.isin(ref_seq_encoded, np.arange(len(self.label_decoding)))

            # homozygous SNPs
            is_primary_diff_to_ref = np.not_equal(primary_labels, ref_seq_encoded)
            is_primary_not_del = primary_labels != self.label_encoding[medaka.common._gap_]
            is_secondary_del = secondary_labels == self.label_encoding[medaka.common._gap_]
            is_secondary_prob_lt_thresh = secondary_probs < threshold
            is_not_secondary_call = np.logical_or(is_secondary_del, is_secondary_prob_lt_thresh)

            is_homozygous_snp = np.logical_and.reduce(
                [is_primary_diff_to_ref, is_primary_not_del, is_not_secondary_call, is_ref_valid_label]
            )
            homozygous_snp_inds = np.where(is_homozygous_snp)

            # heterozygous SNPs
            is_secondary_prob_ge_thresh = np.logical_not(is_secondary_prob_lt_thresh)
            is_secondary_not_del = secondary_labels != self.label_encoding[medaka.common._gap_]

            is_heterozygous_snp = np.logical_and.reduce(
                [is_secondary_prob_ge_thresh, is_secondary_not_del, is_primary_not_del, is_ref_valid_label]
            )
            heterozygous_snp_inds = np.where(is_heterozygous_snp)
            for i in homozygous_snp_inds[0]:
                ref_base_encoded = ref_seq_encoded[i]
                info = {'ref_prob': major_probs[i][ref_base_encoded],
                        'primary_prob': primary_probs[i],
                        'primary_label': self.label_decoding[primary_labels[i]],
                        'secondary_prob': secondary_probs[i],
                        'secondary_label': self.label_decoding[secondary_labels[i]],
                        }
                if self.ref_vcf is not None:
                    ref_info = self._get_ref_variant(ref_name, major_pos[i])
                    info.update(ref_info)
                if major_feat is not None:
                    info.update(dict(zip(self.feature_row_names, major_feat[i])))

                qual = -10 * np.log10(1 - primary_probs[i])
                sample = {'GT': '1/1', 'GQ': qual,}
                called_loci.add(major_pos[i])
                yield medaka.vcf.Variant(ref_name, major_pos[i],
                                         self.label_decoding[ref_base_encoded],
                                         alt=self.label_decoding[primary_labels[i]],
                                         filter='PASS', info=info, qual=qual,
                                         sample_dict=sample)

            for i in heterozygous_snp_inds[0]:
                ref_base_encoded = ref_seq_encoded[i]
                info = {'ref_prob': major_probs[i][ref_base_encoded],
                        'primary_prob': primary_probs[i],
                        'primary_label': self.label_decoding[primary_labels[i]],
                        'secondary_prob': secondary_probs[i],
                        'secondary_label': self.label_decoding[secondary_labels[i]]
                        }
                if self.ref_vcf is not None:
                    ref_info = self._get_ref_variant(ref_name, major_pos[i])
                    info.update(ref_info)
                if major_feat is not None:
                    info.update(dict(zip(self.feature_row_names, major_feat[i])))

                qual = -10 * np.log10(1 - primary_probs[i] - secondary_probs[i])
                alt = [self.label_decoding[l] for l in (primary_labels[i], secondary_labels[i]) if l != ref_base_encoded]
                gt = '0/1' if len(alt) == 1 else '1/2'  # / => unphased
                sample = {'GT': gt, 'GQ': qual,}
                called_loci.add(major_pos[i])
                yield medaka.vcf.Variant(ref_name, major_pos[i],
                                         self.label_decoding[ref_base_encoded],
                                         alt=alt, filter='PASS', info=info,
                                         qual=qual, sample_dict=sample)

            if self.ref_vcf is not None:
                # if we provided a vcf, check which positions are missing
                missing_loci = self.ref_loci[ref_name] - called_loci
                missing_loci_in_chunk = missing_loci.intersection(major_pos)
                missing_loci_in_chunk = np.fromiter(missing_loci_in_chunk, int,
                                                    count=len(missing_loci_in_chunk))
                is_missing = np.isin(major_pos, missing_loci_in_chunk)
                missing_snp_inds = np.where(is_missing)

                for i in missing_snp_inds[0]:
                    ref_base_encoded = ref_seq_encoded[i]
                    info = {'ref_prob': major_probs[i][ref_base_encoded],
                            'primary_prob': primary_probs[i],
                            'primary_label': self.label_decoding[primary_labels[i]],
                            'secondary_prob': secondary_probs[i],
                            'secondary_label': self.label_decoding[secondary_labels[i]],
                            }
                    ref_info = self._get_ref_variant(ref_name, major_pos[i])
                    info.update(ref_info)
                    if major_feat is not None:
                        info.update(dict(zip(self.feature_row_names, major_feat[i])))
                    qual = -10 * np.log10(1 - primary_probs[i])
                    sample = {'GT': 0, 'GQ': qual,}
                    yield medaka.vcf.Variant(ref_name, major_pos[i],
                                             self.label_decoding[ref_base_encoded],
                                             alt='.', filter='PASS', info=info,
                                             qual=qual, sample_dict=sample)

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


class HaploidVariantDecoder(SNPDecoder, metaclass=DecoderRegistrar):
    """Class to decode haploid variants (snps and indels)."""

    def decode_variants(self, samples, ref_seq):
        """Decode haploid variants.

        Multi-base or otherwise complex variants will be decoded in their entirety.

        :param samples: stream of `medaka.common.Sample` objects for a given contig.
        :param ref_seq: reference sequence.

        :returns: sorted list of `medaka.vcf.Variant` objects.
        """

        ref_seq_encoded = np.fromiter((self.label_encoding[b] for b in ref_seq), int, count=len(ref_seq))

        gap_encoded = self.label_encoding[medaka.common._gap_]

        for s in join_chunked_variants(yield_trimmed_consensus_chunks(samples), ref_seq_encoded, gap_encoded):

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
            runs = medaka.common.rle(is_var)
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
    """Entry point for variant calling from hdf consensus probability files."""
    index = medaka.datastore.DataIndex(args.inputs)
    if args.regions is not None:
        ref_names = parse_regions(args.regions)
    else:
        ref_names = medaka.common.loose_version_sort(index.index.keys())

    decoder_cls = variant_decoders[args.decoder]
    decoder = decoder_cls(index.meta, ref_vcf=args.ref_vcf)

    opts = inspect.signature(decoder.decode_variants).parameters.keys()
    opts = {k:getattr(args, k) for k in opts if k not in  {'self', 'samples', 'ref_seq'}}
    logger = medaka.common.get_named_logger('Variants')

    with medaka.vcf.VCFWriter(args.output, 'w', version='4.1', contigs=ref_names,
                              meta_info=decoder.meta_info) as vcf_writer:
        for ref_name in ref_names:
            logger.info("Processing {}.".format(ref_name))
            ref_seq = pysam.FastaFile(args.ref_fasta).fetch(reference=ref_name)
            samples = index.yield_from_feature_files([ref_name])
            variants = decoder.decode_variants(samples, ref_seq, **opts)
            vcf_writer.write_variants(variants, sort=True)
