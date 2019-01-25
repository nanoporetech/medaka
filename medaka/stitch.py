from itertools import chain

import numpy as np
import pysam

from medaka.common import _gap_, Region, get_sample_overlap, get_named_logger
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
    logger = get_named_logger('Stitch')

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
        for region in (Region.from_string(r) for r in regions):
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
                end_1_ind, start_2_ind = get_sample_overlap(s1, s2)

            best = np.argmax(s1.label_probs[start_1_ind:end_1_ind], -1)
            seq += ''.join([label_decoding[x] for x in best]).replace(_gap_, '')
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


#def probs_at_panel_loci(probs_hdfs, loci_file, ref_fasta, out_file, regions=None):
#    """Get label probabilities at specifed reference loci
#
#    :param probs_hdfs: iterable of hdf filepaths.
#    :param loci_file: vcf/bed type file with first column chrom, second col pos (zero-based).
#    :param ref_fasta: reference fasta
#    :param output_file: out_file with chrom pos ref p(A) p(T) p(C) p(G) p(del).
#    """
#    logger = get_named_logger('Panel')
#
#    loci = defaultdict(set)
#    with open(loci_file) as fh:
#        for line in fh:
#            if line.startswith('#'):
#                continue
#            split_line = line.split()
#            chrom = split_line[0]
#            pos = int(split_line[1])
#            loci[chrom].add(pos)
#
#    index = DataIndex(probs_hdfs)
#
#    label_decoding = index.meta['medaka_label_decoding']
#
#    logger.debug("Label decoding is:\n{}".format('\n'.join(
#        '{}: {}'.format(i, x) for i, x in enumerate(label_decoding)
#    )))
#    if regions is None:
#        ref_names = index.index.keys()
#    else:
#        #TODO: respect entire region specification
#        ref_names = list()
#        for region in (Region.from_string(r) for r in regions):
#            if region.start is None or region.end is None:
#                logger.warning("Ignoring start:end for '{}'.".format(region))
#            ref_names.append(region.ref_name)
#
#    # For a panel, we assume we just want the label probability at given loci
#    # We need to handle boundary effects simply to make sure we take the best
#    # probability (which is furthest from the boundary).
#    headers = ['chrom', 'pos', 'ref'] + label_decoding
#    with open(out_file, 'w') as out_fh:
#        out_fh.write('\t'.join(headers) + '\n')
#        for ref_name in ref_names:
#            ref_seq = pysam.FastaFile(ref_fasta).fetch(reference=ref_name)
#            logger.info("Processing {}.".format(ref_name))
#            # TODO: refactor this and stitch to use common func/generator to get
#            # chunks with overlapping bits trimmed off
#            data_gen = index.yield_from_feature_files(ref_names=(ref_name,))
#            s1 = next(data_gen)
#            start = get_pos(s1, 0)
#            start_1_ind = None  # don't trim beginning of s1
#            for s2 in chain(data_gen, (None,)):
#                if s2 is None:  # s1 is last chunk
#                    end_1_ind = None  # go to the end of s2
#                else:
#                    end_1_ind, start_2_ind = get_sample_overlap(s1, s2)
#
#                trimmed_pos = s1.positions[start_1_ind:end_1_ind]
#                trimmed_probs = s1.label_probs[start_1_ind:end_1_ind]
#                loci_inds = np.where(np.isin(trimmed_pos['major'], loci[ref_name]))
#                loci_pos = trimmed_pos[loci_inds]
#                loci_probs = trimmed_probs[loci_inds]
#                # discard minor positions (insertions)
#                loci_major_inds = np.where(loci_pos['minor'] == 0)
#                loci_major_pos = loci_pos[loci_major_inds]
#                loci_major_probs = loci_probs[loci_major_inds]
#
#                for pos, probs in zip(loci_major_pos, loci_major_probs):
#                    out.write('\t'.join(map(str, [ref_name, pos, ref_seq[pos]] + list(probs))) + '\n')
#
#                if end_1_ind is None:
#                    if start_2_ind is None:
#                        msg = 'There is no overlap betwen {} and {}'
#                        logger.info(msg.format(s1.name, s2.name))
#                s1 = s2
#                start_1_ind = start_2_ind


def find_snps(probs_hdfs, ref_fasta, out_file, regions=None, threshold=0.1):
    """Find potential homozygous and heterozygous snps based on a probability threshold.

    :param probs_hdfs: iterable of hdf filepaths.
    :param ref_fasta: reference fasta.
    :param threshold: threshold below which a secondary call (which would make
            for a heterozygous call) is deemed insignificant.
    :param output_file: out_file with chrom pos ref p(A) p(T) p(C) p(G) p(del).
    """
    logger = get_named_logger('Panel')

    index = DataIndex(probs_hdfs)

    label_decoding = index.meta['medaka_label_decoding']
    label_encoding = {label: ind for ind, label in enumerate(label_decoding)}

    logger.debug("Label decoding is:\n{}".format('\n'.join(
        '{}: {}'.format(i, x) for i, x in enumerate(label_decoding)
    )))
    if regions is None:
        ref_names = index.index.keys()
    else:
        #TODO: respect entire region specification
        ref_names = list()
        for region in (Region.from_string(r) for r in regions):
            if region.start is None or region.end is None:
                logger.warning("Ignoring start:end for '{}'.".format(region))
            ref_names.append(region.ref_name)

    # For SNPS, we assume we just want the label probabilities at major positions
    # We need to handle boundary effects simply to make sure we take the best
    # probability (which is furthest from the boundary).
    with VCFWriter(out_file, 'w') as vcf_writer:
        for ref_name in ref_names:
            ref_seq = pysam.FastaFile(ref_fasta).fetch(reference=ref_name)
            logger.info("Processing {}.".format(ref_name))
            # TODO: refactor this and stitch to use common func/generator to get
            # chunks with overlapping bits trimmed off
            data_gen = index.yield_from_feature_files(ref_names=(ref_name,))
            s1 = next(data_gen)
            start = get_pos(s1, 0)
            start_1_ind = None  # don't trim beginning of s1
            for s2 in chain(data_gen, (None,)):
                if s2 is None:  # s1 is last chunk
                    end_1_ind = None  # go to the end of s2
                else:
                    end_1_ind, start_2_ind = get_sample_overlap(s1, s2)

                pos = s1.positions[start_1_ind:end_1_ind]
                probs = s1.label_probs[start_1_ind:end_1_ind]
                # discard minor positions (insertions)
                major_inds = np.where(pos['minor'] == 0)
                major_pos = pos[major_inds]
                major_probs = probs[major_inds]

                # for homozygous SNP max_prob_label not in {ref, del} and 2nd_max_prob < threshold
                # for heterozygous SNP 2nd_max_prob > threshold and del not in {max_prob_label, 2nd_prob_label}
                # this catches both SNPs where the genotype contains the
                # reference, and where both copies are mutated.

                # Note: if 2nd_max_prob > thresh and is a deletion, we will
                # ignore the SNP (it could be a SNP in one haplotype and a
                # deletion in the other).

                ref_seq_encoded = np.fromiter((label_encoding[ref_seq[i]] for i in major_pos), int, count=len(label_decoding))
                sorted_prob_inds = np.argsort(major_probs, -1)
                sorted_probs = np.take_along_axis(major_probs, sorted_prob_inds, axis=-1)
                primary_labels = sorted_prob_inds[:, -1]
                secondary_labels = sorted_prob_inds[:, -2]
                secondary_probs = sorted_probs[:, -2]
                # homozygous SNPs
                is_primary_diff_to_ref = np.not_equal(primary_labels, ref_seq_encoded)
                is_primary_not_del = primary_labels != label_encoding[_gap_]
                is_secondary_prob_lt_thresh = secondary_probs < threshold
                is_homozygous_snp = np.logical_and(is_primary_diff_to_ref, is_primary_not_del)
                is_homozygous_snp = np.logical_and(is_homozygous_snp, is_secondary_prob_lt_thresh)
                homozygous_snp_inds = np.where(is_homozygous_snp)
                # heterozygous SNPs
                is_secondary_prob_ge_thresh = np.logical_not(is_secondary_prob_lt_thresh)
                is_secondary_not_del = secondary_labels != label_encoding[_gap_]
                is_heterozygous_snp = np.logical_and(is_secondary_prob_ge_thresh, is_secondary_not_del)
                is_heterozygous_snp = np.logical_and(is_heterozygous_snp, is_primary_not_del)
                heterozygous_snp_inds = np.where(is_heterozygous_snp)

                for i in homozygous_snp_inds:
                    ref_base = ref_seq[i]
                    info = {ref_prob: major_probs[i][label_encoding[ref_base],
                            primary_prob: primary_prob[i],
                            secondary_prob: secondary_prob[i]
                            secondary_label: label_decoding[secondary_labels[i]]
                            }
                    variant = Variant(ref_name, major_pos[i], ref_seq[i],
                                      alt=label_decoding[primary_labels[i]],
                                      gt=1, info=info)
                    vcf_writer.write_variant(variant)

                for i in heterozygous_snp_inds:
                    ref_base = ref_seq[i]
                    info = {ref_prob: major_probs[i][label_encoding[ref_base],
                            primary_prob: primary_prob[i],
                            secondary_prob: secondary_prob[i]
                            secondary_label: label_decoding[secondary_labels[i]]
                            }
                    variant = Variant(ref_name, major_pos[i], ref_seq[i],
                                      alt=(label_decoding[primary_labels[i]], label_decoding[secondary_labels[i]]),
                                      gt='1/2', info=info)  # / => unphased
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
