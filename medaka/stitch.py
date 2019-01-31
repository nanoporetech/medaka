from collections import defaultdict
from itertools import chain

import numpy as np
import pysam

import medaka
from medaka.datastore import DataIndex
from medaka.vcf import Variant, VCFWriter


def write_fasta(filename, contigs):
    with open(filename, 'w') as fasta:
        for name, seq in contigs:
            fasta.write('>{}\n{}\n'.format(name, seq))


def stitch_from_probs(probs_hdfs, regions=None, model_yml=None):
    """Decode and join overlapping label probabilities from hdf to assemble complete sequence

    :param probs_hdfs: iterable of hdf filepaths.
    :ref_names: iterable of ref_names to limit stitching to.
    :returns: list of (contig_name, sequence)

    """
    logger = medaka.common.get_named_logger('Stitch')

    index = DataIndex(probs_hdfs)

    label_decoding = index.meta['medaka_label_decoding']

    logger.debug("Label decoding is:\n{}".format('\n'.join(
        '{}: {}'.format(i, x) for i, x in enumerate(label_decoding)
    )))
    if regions is None:
        ref_names = index.index.keys()
    else:
        #TODO: respect entire region specification
        ref_names = list()
        for region in (medaka.common.Region.from_string(r) for r in regions):
            if region.start is None or region.end is None:
                logger.warning("Ignoring start:end for '{}'.".format(region))
            ref_names.append(region.ref_name)

    get_pos = lambda s, i: '{}.{}'.format(s.positions[i]['major'] + 1, s.positions[i]['minor'])
    ref_assemblies = []
    for ref_name in ref_names:
        logger.info("Processing {}.".format(ref_name))
        data_gen = index.yield_from_feature_files(ref_names=(ref_name,))
        seq=''
        s1 = next(data_gen)
        start = get_pos(s1, 0)
        start_1_ind = None
        for s2 in chain(data_gen, (None,)):
            if s2 is None:  # s1 is last chunk
                end_1_ind = None
            else:
                if s2.last_pos <= s1.last_pos:
                    logger.info('{} ends before {}, skipping.'.format(
                        s2.name, s1.name
                    ))
                    continue
                else:
                    end_1_ind, start_2_ind = medaka.common.get_sample_overlap(s1, s2)

            best = np.argmax(s1.label_probs[start_1_ind:end_1_ind], -1)
            seq += ''.join([label_decoding[x] for x in best]).replace(medaka.common._gap_, '')
            if end_1_ind is None:
                key = '{}:{}-{}'.format(s1.ref_name, start, get_pos(s1, -1))
                ref_assemblies.append((key, seq))
                seq = ''
                if start_2_ind is None:
                    msg = 'There is no overlap betwen {} and {}'
                    logger.info(msg.format(s1.name, s2.name))
                    start = get_pos(s2, 0)

            s1 = s2
            start_1_ind = start_2_ind
    return ref_assemblies


def find_snps(probs_hdfs, ref_fasta, out_file, regions=None, threshold=0.1, ref_vcf=None):
    """Find potential homozygous and heterozygous snps based on a probability threshold.

    :param probs_hdfs: iterable of hdf filepaths.
    :param ref_fasta: reference fasta.
    :param out_file: vcf output file.
    :param threshold: threshold below which a secondary call (which would make
            for a heterozygous call) is deemed insignificant.
    :param ref_vcf: input vcf to force evaluation only at these loci, even if we would
            not otherwise call a SNP there.
    """
    logger = medaka.common.get_named_logger('SNPs')

    index = DataIndex(probs_hdfs)

    label_decoding = index.meta['medaka_label_decoding']
    # include extra bases so we can encode e.g. N's in reference
    label_encoding = {label: ind for ind, label in enumerate(label_decoding + list(medaka.common._extra_bases_))}

    fmt_feat = lambda x: '{}{}'.format('rev' if x[0] else 'fwd', x[2] * (x[1] if x[1] is not None else '-'))
    feature_row_names = [fmt_feat(x) for x in index.meta['medaka_feature_decoding']]

    logger.debug("Label decoding is:\n{}".format('\n'.join(
        '{}: {}'.format(i, x) for i, x in enumerate(label_decoding)
    )))
    if regions is None:
        ref_names = index.index.keys()
    else:
        #TODO: respect entire region specification
        ref_names = list()
        for region in (medaka.common.Region.from_string(r) for r in regions):
            if region.start is None or region.end is None:
                logger.warning("Ignoring start:end for '{}'.".format(region))
            ref_names.append(region.ref_name)


    def _get_ref_variant(ref_vcf, ref_name, pos):
        ref_variants = list(ref_vcf.fetch(ref_name=ref_name, start=pos, end=pos+1))
        assert len(ref_variants) < 2
        if len(ref_variants) == 0:
            ref_info = {'ref_alt': 'na', 'ref_gt': 'na'}
        else:
            ref_info = {'ref_alt': ','.join(ref_variants[0].alt),
                        'ref_gt': ref_variants[0].sample_dict['GT']}
        return ref_info

    if ref_vcf is not None:
        ref_vcf = medaka.vcf.VCFReader(ref_vcf)
        vcf_loci = defaultdict(set)
        for v in ref_vcf.fetch():
            vcf_loci[v.chrom].add(v.pos)

    # For SNPS, we assume we just want the label probabilities at major positions
    # We need to handle boundary effects simply to make sure we take the best
    # probability (which is furthest from the boundary).
    with VCFWriter(out_file, 'w', version='4.1') as vcf_writer:
        for ref_name in ref_names:
            called_loci = set()
            ref_seq = pysam.FastaFile(ref_fasta).fetch(reference=ref_name)
            logger.info("Processing {}.".format(ref_name))
            # TODO: refactor this and stitch to use common func/generator to get
            # chunks with overlapping bits trimmed off
            data_gen = index.yield_from_feature_files(ref_names=(ref_name,))
            s1 = next(data_gen)
            start_1_ind = None  # don't trim beginning of s1
            for s2 in chain(data_gen, (None,)):
                if s2 is None:  # s1 is last chunk
                    end_1_ind = None  # go to the end of s2
                else:
                    end_1_ind, start_2_ind = medaka.common.get_sample_overlap(s1, s2)

                pos = s1.positions[start_1_ind:end_1_ind]
                probs = s1.label_probs[start_1_ind:end_1_ind]
                # discard minor positions (insertions)
                major_inds = np.where(pos['minor'] == 0)
                major_pos = pos[major_inds]['major']
                major_probs = probs[major_inds]
                major_feat = s1.features[start_1_ind:end_1_ind][major_inds] if s1.features is not None else None


                # for homozygous SNP max_prob_label not in {ref, del} and 2nd_max_prob < threshold
                # for heterozygous SNP 2nd_max_prob > threshold and del not in {max_prob_label, 2nd_prob_label}
                # this catches both SNPs where the genotype contains the
                # reference, and where both copies are mutated.

                # Note: if 2nd_max_prob > thresh and is a deletion, we will
                # ignore the SNP (it could be a SNP in one haplotype and a
                # deletion in the other).
                ref_seq_encoded = np.fromiter((label_encoding[ref_seq[i]] for i in major_pos), int, count=len(major_pos))
                sorted_prob_inds = np.argsort(major_probs, -1)
                sorted_probs = np.take_along_axis(major_probs, sorted_prob_inds, axis=-1)
                primary_labels = sorted_prob_inds[:, -1]
                secondary_labels = sorted_prob_inds[:, -2]
                primary_probs = sorted_probs[:, -1]
                secondary_probs = sorted_probs[:, -2]
                # skip positions where ref is not a label (ATGC)
                is_ref_valid_label = np.isin(ref_seq_encoded, np.arange(len(label_decoding)))

                # homozygous SNPs
                is_primary_diff_to_ref = np.not_equal(primary_labels, ref_seq_encoded)
                is_primary_not_del = primary_labels != label_encoding[medaka.common._gap_]
                is_secondary_prob_lt_thresh = secondary_probs < threshold
                is_homozygous_snp = np.logical_and(is_primary_diff_to_ref, is_primary_not_del)
                is_homozygous_snp = np.logical_and(is_homozygous_snp, is_secondary_prob_lt_thresh)
                is_homozygous_snp = np.logical_and(is_homozygous_snp, is_ref_valid_label)
                homozygous_snp_inds = np.where(is_homozygous_snp)
                # heterozygous SNPs
                is_secondary_prob_ge_thresh = np.logical_not(is_secondary_prob_lt_thresh)
                is_secondary_not_del = secondary_labels != label_encoding[medaka.common._gap_]
                is_heterozygous_snp = np.logical_and(is_secondary_prob_ge_thresh, is_secondary_not_del)
                is_heterozygous_snp = np.logical_and(is_heterozygous_snp, is_primary_not_del)
                is_heterozygous_snp = np.logical_and(is_heterozygous_snp, is_ref_valid_label)
                heterozygous_snp_inds = np.where(is_heterozygous_snp)
                variants = []
                for i in homozygous_snp_inds[0]:
                    ref_base_encoded = ref_seq_encoded[i]
                    info = {'ref_prob': major_probs[i][ref_base_encoded],
                            'primary_prob': primary_probs[i],
                            'primary_label': label_decoding[primary_labels[i]],
                            'secondary_prob': secondary_probs[i],
                            'secondary_label': label_decoding[secondary_labels[i]],
                            }
                    if ref_vcf is not None:
                        ref_info = _get_ref_variant(ref_vcf, ref_name, major_pos[i])
                        info.update(ref_info)
                    if major_feat is not None:
                        info.update(dict(zip(feature_row_names, major_feat[i])))

                    qual = -10 * np.log10(1 - primary_probs[i])
                    sample = {'GT': 1, 'GQ': qual,}
                    variants.append(Variant(ref_name, major_pos[i], label_decoding[ref_base_encoded],
                                      alt=label_decoding[primary_labels[i]],
                                      filter='PASS', info=info, qual=qual, sample_dict=sample))

                for i in heterozygous_snp_inds[0]:
                    ref_base_encoded = ref_seq_encoded[i]
                    info = {'ref_prob': major_probs[i][ref_base_encoded],
                            'primary_prob': primary_probs[i],
                            'primary_label': label_decoding[primary_labels[i]],
                            'secondary_prob': secondary_probs[i],
                            'secondary_label': label_decoding[secondary_labels[i]]
                            }
                    if ref_vcf is not None:
                        ref_info = _get_ref_variant(ref_vcf, ref_name, major_pos[i])
                        info.update(ref_info)
                    if major_feat is not None:
                        info.update(dict(zip(feature_row_names, major_feat[i])))

                    qual = -10 * np.log10(1 - primary_probs[i] - secondary_probs[i])
                    alt = [label_decoding[l] for l in (primary_labels[i], secondary_labels[i]) if l != ref_base_encoded]
                    gt = '0/1' if len(alt) == 1 else '1/2'  # / => unphased
                    sample = {'GT': gt, 'GQ': qual,}
                    variants.append(Variant(ref_name, major_pos[i], label_decoding[ref_base_encoded],
                                      alt=alt, filter='PASS', info=info, qual=qual, sample_dict=sample))

                if ref_vcf is not None:
                    # if we provided a vcf, check which positions are missing
                    called_loci.update({v.pos for v in variants})
                    missing_loci = vcf_loci[ref_name] - called_loci
                    missing_loci_in_chunk = missing_loci.intersection(major_pos)
                    missing_loci_in_chunk = np.fromiter(missing_loci_in_chunk, int, count=len(missing_loci_in_chunk))
                    is_missing = np.isin(major_pos, missing_loci_in_chunk)
                    missing_snp_inds = np.where(is_missing)

                    for i in missing_snp_inds[0]:
                        ref_base_encoded = ref_seq_encoded[i]
                        info = {'ref_prob': major_probs[i][ref_base_encoded],
                                'primary_prob': primary_probs[i],
                                'primary_label': label_decoding[primary_labels[i]],
                                'secondary_prob': secondary_probs[i],
                                'secondary_label': label_decoding[secondary_labels[i]],
                                }
                        ref_info = _get_ref_variant(ref_vcf, ref_name, major_pos[i])
                        info.update(ref_info)
                        if major_feat is not None:
                            info.update(dict(zip(feature_row_names, major_feat[i])))
                        qual = -10 * np.log10(1 - primary_probs[i])
                        sample = {'GT': 0, 'GQ': qual,}
                        variants.append(Variant(ref_name, major_pos[i], label_decoding[ref_base_encoded],
                                          alt='.',
                                          filter='PASS', info=info, qual=qual, sample_dict=sample))

                sorter = lambda v: v.pos
                variants.sort(key=sorter)
                for variant in variants:
                    vcf_writer.write_variant(variant)

                if end_1_ind is None:
                    if start_2_ind is None:
                        msg = 'There is no overlap betwen {} and {}'
                        logger.info(msg.format(s1.name, s2.name))
                s1 = s2
                start_1_ind = start_2_ind


def stitch(args):
    joined = stitch_from_probs(args.inputs, regions=args.regions)
    write_fasta(args.output, joined)


def snps(args):
    find_snps(args.inputs, args.ref_fasta, args.output, regions=args.regions, threshold=args.threshold, ref_vcf=args.ref_vcf)
