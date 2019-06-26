import collections
import itertools
import intervaltree
import numpy as np
import pysam
from copy import copy
from operator import attrgetter

import medaka.common


class TruthAlignment(object):
    def __init__(self, alignment):
        """Create a TruthAlignment oblist from a `pysam.libcalignedsegment.AlignedSegment` object.

        :param alignment: `pysam.libcalignedsegment.AlignedSegment` object.
        """
        self.aln = alignment  # so we can get positions and labels later
        # initialise start and end (which might be moved)
        self.start = self.aln.reference_start  # zero-based
        self.end = self.aln.reference_end
        self.is_kept = True
        self.logger = medaka.common.get_named_logger('TruthAlign')

    def _get_overlap_with(self, other):
        first, second = sorted((self, other), key=attrgetter('aln.reference_start'))
        overlap = first.aln.reference_end > second.aln.reference_start
        if overlap:
            # return positions of start and end of overlapping region
            return second.aln.reference_start, first.aln.reference_end
        else:
            return None

    @staticmethod
    def _filter_alignments(alignments, region, min_length=1000, length_ratio=2.0, overlap_fraction=0.5):
        """Filter alignments to yield only segments suitable for training.

        :param alignments: iterable of `TruthAlignment` objects.
        :param region: `medaka.common.Region` obj. all alignments will be trimmed to this region.
        :param min_length: int, minimum length of alignment segment.
        :param length_ratio: float, in cases of overlap, ratio of longer segment / shorter segment
            above which the longer segment is assumed to be more reliable.
        :param overlap_fraction: float, length_ratio: float, in cases of overlap, fraction of
            shorter segment overlapping with the longer segment above which the segments are
            considered highly overlapping.
        :returns: list of `TruthAlignment` objects

        """
        filtered_alignments = [copy(a) for a in alignments]  # don't want to modify original alignments
        for al_i, al_j in itertools.combinations(filtered_alignments, 2):
            first, second = sorted((al_i, al_j), key=attrgetter('aln.reference_start'))
            overlap = first._get_overlap_with(second)
            if overlap is None:
                continue
            ovlp_start, ovlp_end = overlap
            shorter, longer = sorted((al_i, al_j), key=attrgetter('aln.reference_length'))
            length_ratio_ij = longer.aln.reference_length / shorter.aln.reference_length
            overlap_fraction_ij = (ovlp_end - ovlp_start) / shorter.aln.reference_length
            # 4 cases
            if length_ratio_ij < length_ratio:  # we don't trust one more than the other
                if overlap_fraction_ij >= overlap_fraction:
                    # 1) they overlap a lot; we have significant ambiguity, discard both
                    shorter.is_kept = False
                    longer.is_kept = False
                else:
                    # 2) they overlap a little; just trim overlapping portions of both alignments
                    first.end = ovlp_start
                    second.start = ovlp_end
            else:  # we trust the longer one more than the shorter one
                if overlap_fraction_ij >= overlap_fraction:
                    # 3) they overlap a lot; discard shorter alignment
                    shorter.is_kept = False
                else:
                    # 4) they overlap a little; trim overlapping portion of shorter alignment
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
                               if (al.is_kept and al.end - al.start >= min_length)]
        filtered_alignments.sort(key=attrgetter('start'))
        return filtered_alignments


    @staticmethod
    def bam_to_alignments(truth_bam, region, haplotag=None):
        """Get processed truth alignments

        :param truth_bam: (sorted indexed) bam with true sequence aligned to reference
        :param region: `medaka.common.Region` obj,
            (all alignments with any overlap with the interval start:end will be retrieved)
        :param haplotag: bam tag specifying which haplotype the alignment belongs to (for multi-ploidy labels)
        :returns: list of tuples where each tuple contains `TruthAlignment` objs
            for each haplotype trimmed to common genomic window.
        """
        algns = TruthAlignment._load_alignments(truth_bam, region, haplotag)
        # filter truth alignments to restrict ourselves to regions of the ref where the truth
        # in unambiguous in each haplotype
        algns = {h: TruthAlignment._filter_alignments(h_algns, region=region) for h, h_algns in algns.items()}
        # Group truth alignments from the haplotypes by genomic window and trim
        # to common genomic range
        grouped = TruthAlignment._group_and_trim_by_haplotype(algns)
        return grouped


    @staticmethod
    def _group_and_trim_by_haplotype(alignments):
        """Group alignments by haplotype tag and trim to common genomic window.
        :param alignments: {haplotype: [`TruthAlignment` objs]}
        :returns: list of tuples where each tuple contains `TruthAlignment` objs
            for each haplotype trimmed to common genomic window.

        .. note:: We should avoid the situation of staggered alignments
             which could occur by independently chunking each haplotype
             by chunking the draft and aligning to both haplotypes, then
             chunking both haplotypes according to draft-chunks, then realining
             haplotype chunks to back to the draft - this should minimize
             staggering of truth alignments and hence the number of labels
             discarded.
             -------------  ----------
                     --------------
        """
        haplotypes = sorted(list(alignments.keys()))
        if len(haplotypes) == 1:  # haploid
            grouped = [(a,) for a in alignments[haplotypes[0]]]
        else:
            # create interval trees for other haplotypes
            trees = {}
            for h in haplotypes[1:]:
                trees[h] = intervaltree.IntervalTree([intervaltree.Interval(a.start, a.end, a) for a in alignments[h]])
            # loop over alignments in first haplotype and find overlapping
            # alignments in other haplotypes. If there are multiple overlapping
            # alignments, take the one with the longest overlap.
            grouped = []
            for a in alignments[haplotypes[0]]:
                group = [a]
                common_start, common_end = a.start, a.end
                for h, tree in trees.items():
                    h_algns = list(tree.overlap(common_start, common_end))
                    if len(h_algns) == 0:  # no alignments on this haplotype, skip this alignment
                        break
                    elif len(h_algns) > 1:  # find most overlapping alignment
                        ovlps = [min(common_end, o.end) - max(common_start, o.begin)]
                        h_algn = h_algns[np.argmax(ovlps)].data
                    elif len(h_algns) == 1:  # keep this alignment
                        h_algn = h_algns[0].data

                    common_start = max(common_start, h_algn.start)
                    common_end = min(common_end, h_algn.end)
                    group.append(h_algn)
                if len(group) != len(haplotypes):  # skip this group
                    msg = 'Skipping {}:{}-{} as some haplotypes have no alignments'
                    self.logger.info(msg.format(a.alignment.ref_name, a.start, a.end))
                    continue
                else:  # trim all alignments to common start/end
                    for i in group:
                        i.start = common_start
                        i.end = common_end
                    grouped.append(tuple(group))
        return grouped


    @staticmethod
    def _load_alignments(truth_bam, region, haplotag=None):
        """Create list of TruthAlignment objects from a bam of Truth aligned to ref.

        :param truth_bam: (sorted indexed) bam with true sequence aligned to reference
        :param region: `medaka.common.Region` obj
        :param haplotag: bam tag specifying which haplotype the alignment belongs to (for multi-ploidy labels)

        :returns: {haplotype: [`TruthAlignment` objs]}

        """
        alignments = collections.defaultdict(list)
        with pysam.AlignmentFile(truth_bam, 'rb') as bamfile:
            aln_reads = bamfile.fetch(reference=region.ref_name, start=region.start, end=region.end)
            for r in aln_reads:
                if (r.is_unmapped or r.is_secondary):
                    continue
                else:
                    hap = r.get_tag(haplotag) if haplotag is not None else None
                    alignments[hap].append(TruthAlignment(r))

        logger = medaka.common.get_named_logger("TruthAlign")
        for hap in alignments.keys():
            alignments[hap].sort(key=attrgetter('start'))
            logger.info("Retrieved {} alignments for haplotype {}.".format(len(alignments[hap]), hap))
        return alignments


    @staticmethod
    def get_positions_and_labels(h_algns, aln_to_pairs):
        """Create labels and positions array.

        :param h_algns: tuple of alignments for each haplotype spanning the same genomic range.
        :param aln_to_pairs: callable, that gets aligned pairs from an `AlignedSegment` obj.
        :returns: tuple(positions, encoded_label_array)

            - positions: numpy structured array with 'ref_major'
              (reference position index) and 'ref_minor'
              (trailing insertion index) fields.

            - encoded_label_array: numpy array, with shape (pos, ploidy) and
                           dtype [('base, int), ('run_length', int)]
        """
        # ensure all provided alignments have the same start end
        if max(len(set((a.end for a in h_algns))), len(set((a.start for a in h_algns)))) != 1:
            raise ValueError('Alignments must be trimmed to common genomic region')

        start = h_algns[0].start
        end = h_algns[0].end

        pos2label_haps = []
        pad = (medaka.common.encoding[medaka.common._gap_], 1)

        for aln in h_algns:
            pairs = aln_to_pairs(aln.aln)
            # pad labels with encoded gap
            padder = lambda: pad
            pos2label = collections.defaultdict(padder)
            ins_count = 0
            for pair in itertools.dropwhile(lambda x: (x.rpos is None)
                                            or (x.rpos < start), pairs):
                if pair.rpos is not None and pair.rpos >= end:
                    break
                if pair.rpos is None:
                    ins_count += 1
                else:
                    ins_count = 0
                    current_pos = pair.rpos
                pos = (current_pos, ins_count)
                label = pair.qbase.upper() if pair.qbase else medaka.common._gap_
                if label == 'N':
                    self.logger.info('Found {} at pos {}'.format(label, pos))
                label = (medaka.common.encoding[label], pair.qlen)
                pos2label[pos] = label

            pos2label_haps.append(pos2label)

        # get positions common to all haplotypes
        positions = sorted(set(itertools.chain(*[p.keys() for p in pos2label_haps])))
        # create label array
        dtype = [('base', int), ('run_length', int)]
        ploidy = len(pos2label_haps)
        label_array = np.empty(shape=(len(positions), ploidy), dtype=dtype)
        # fill with pad so that insertions present on only one haplotype have
        # correct gap-label on other haplotype(s)
        label_array.fill(pad)

        for h, pos2label in enumerate(pos2label_haps):
            label_array[:, h] = np.fromiter((pos2label[p] for p in positions), dtype=dtype, count=len(positions))

        positions = np.array(positions, dtype=[('major', int), ('minor', int)])

        return positions, label_array
