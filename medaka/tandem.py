"""Variant calling in tandem repeats."""
import collections
import concurrent.futures
from copy import copy, deepcopy
import dataclasses
import functools
import os
import re
import sys
from timeit import default_timer as now


import edlib
import numpy as np
import pysam

from medaka import abpoa, parasail
import medaka.align
import medaka.common
import medaka.features
import medaka.rle
import medaka.smolecule
import medaka.vcf


class InsufficientCoverage(Exception):
    """Exception for tracking cases of insufficient read coverage."""

    pass


class RecordName(object):
    """Simple class for encoding/decoding sequence metadata to/from strings."""

    def __init__(
        self, query_name, ref_name, ref_start, ref_end, hap=0,
        strand='fwd', ref_start_padded=None, ref_end_padded=None,
    ):
        """Initialize repeat read analysis.

        :param query_name: str
        :param ref_name: str
        :param ref_start: int
        :param ref_end: int
        :param hap: int
        :param strand: str, {fwd, rev}
        :param ref_start_pad: padded reference start, int
        :param ref_end_pad: padded reference end, int.
        """
        self.query_name = query_name
        self.ref_name = ref_name
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.hap = hap
        self.strand = strand
        self.ref_start_padded = (
            ref_start if ref_start_padded is None else ref_start_padded)
        self.ref_end_padded = (
            ref_end if ref_end_padded is None else ref_end_padded)

    def __str__(self):
        """Encode as a string."""
        return (
            f"{self.query_name}_{self.ref_name}_{self.ref_start}_"
            f"{self.ref_end}_pad_{self.ref_start_padded}_"
            f"{self.ref_end_padded}_{self.strand}_hap{self.hap}"
        )

    @classmethod
    def from_str(cls, name):
        """Decode from string."""
        pattern = (
            r"(?P<query_name>.+)_(?P<ref_name>.+)_(?P<ref_start>\d+)_"
            r"(?P<ref_end>\d+)_pad_(?P<ref_start_padded>\d+)_"
            r"(?P<ref_end_padded>\d+)_(?P<strand>fwd|rev)_hap(?P<hap>\d+)"
        )
        m = re.match(pattern, name)
        if m is None:
            raise ValueError(f"Could not parse {name} using {pattern}")
        d = m.groupdict()
        for field in (
            'ref_start', 'ref_end', 'hap', 'ref_start_padded', 'ref_end_padded'
        ):
            d[field] = int(d[field])
        return cls(**d)

    def sorter(self):
        """Get sorting keys."""
        return self.ref_name, self.ref_start

    def to_padded_region(self):
        """Convert to padded `medaka.common.Region`."""
        return medaka.common.Region(
            self.ref_name, self.ref_start_padded, self.ref_end_padded)

    def to_unpadded_region(self):
        """Convert to unpadded `medaka.common.Region`."""
        return medaka.common.Region(
            self.ref_name, self.ref_start, self.ref_end)

#########################
# read-fetching functions
#########################


def get_trimmed_reads(bam, region, partial, read_filters):
    """Get trimmed reads without reference.

    :param bam: bam containing reads
    :param region: `medaka.common.Region` obj
    :param partial: bool, whether to keep reads which don't fully span region.
    :param read_filters: read filtering kwargs,
        see `medaka.features.get_trimmed_read`.

    :returns: [(bool is_rev, str trimmed read sequence)]
    """
    # TODO move to featues / common and remove vcf.get_trimmed_reads
    try:
        region_got, reads = next(
            medaka.features.get_trimmed_reads(
                region, bam, partial=partial,
                region_split=2*region.size,
                **read_filters)

        )
    except StopIteration:
        return list()
    if not region == region_got:  # check region was not e.g. split
        msg = 'Expected region {}, got region {}'
        raise ValueError(msg.format(region, region_got))
    # first read is ref, remove it.
    ref_is_rev, ref_seq = reads.pop(0)
    return ref_seq, reads


def get_subreads(bam_fp, rec, hap_tag_vals=(0, 1, 2), min_mapq=5):
    """Get reads separated by haplotype.

    :param bam_fp: str, path to bam file.
    :param rec: RecordName, record
    :param hap_tag_vals: HP tags to fetch reads.

    :returns: list of `medaka.smolecule.Subread` instances.
    """
    reg_padded = rec.to_padded_region()
    subreads = []

    # TODO rewrite C-code to fetch reads with multiple tags and return tag vals
    # to avoid querying bam three times.
    for h in hap_tag_vals:
        # 0 for untagged - no reads will have this tag so keep missing
        # will retrieve untagged
        read_filters = dict(
            tag_name='HP', tag_value=h, keep_missing=h == 0, min_mapq=min_mapq)

        _, reads = get_trimmed_reads(
            bam_fp, reg_padded, partial=False, read_filters=read_filters)

        rn_kwargs = {
            k: v for k, v in rec.__dict__.items()
            if k not in {'query_name', 'strand', 'hap'}
        }

        try:
            subreads.extend([
                medaka.smolecule.Subread(
                    str(RecordName(
                        query_name=f"read_{i:03d}",
                        strand='rev' if is_rev else 'fwd',
                        hap=h, **rn_kwargs)),
                    medaka.common.reverse_complement(seq) if is_rev else seq))
                for i, (is_rev, seq) in enumerate(reads, 1)]
            )
        except ValueError:
            continue
    return subreads


#####################
# alignment functions
#####################


def align_chunk_to_ref(chunk, ref_fasta, aln_header=None):
    """Align consensus chunk to reference using global alignment.

    :param chunk: `pysam.FastxRecord` or similar with name and sequence attrs.
        name should be the string representation of a `RecordName`.
    :param ref_fasta: `pysam.FastaFile`.
    :param aln_header: `pysam.AlignmentHeader`.

    :returns: `pysam.AlignedSegment` instance.
    """
    logger = medaka.common.get_named_logger("TR_ALIGN")
    if aln_header is None:
        aln_header = pysam.AlignmentHeader.from_references(
            ref_fasta.references, ref_fasta.lengths)

    rn = RecordName.from_str(chunk.name)

    ref_seq = ref_fasta.fetch(
        rn.ref_name, rn.ref_start_padded, rn.ref_end_padded)
    if rn.strand == 'fwd':
        query_seq = chunk.sequence
        flag = 0
    else:
        query_seq = medaka.common.reverse_complement(chunk.sequence)
        flag = 16

    # use global alignment
    parasail_aligner = functools.partial(
        parasail.nw_trace_scan_32,
        open=8, extend=4, matrix=parasail.dnafull
    )
    result = parasail_aligner(query_seq, ref_seq)
    rstart, cigar = medaka.align.parasail_to_sam(result, query_seq)
    # rstart may not be zero if alignent starts with a SNP / indel.
    # TODO - do we want to automatically increase padding and rerun?
    if rstart > 0:
        logger.warn(
            f"rstart not 0 when using global alignment for {chunk.name}. "
            "Consider rerunning this region with more padding.")

    aln = medaka.align.initialise_alignment(
        query_name=chunk.name,
        reference_id=aln_header.get_tid(rn.ref_name),
        reference_start=rn.ref_start_padded + rstart,
        query_sequence=query_seq,
        cigarstring=cigar,
        flag=flag,
        header=aln_header,
        tags=dict(HP=rn.hap),
    )
    return aln


def align_consensus_fx_to_ref(consensus_fx, bam_fp, ref_fasta):
    """Align consensus fastx to ref producing bam.

    :param consensus_fx: filepath to consensus fastx file, str.
    :param bam_fp: output filepath of bam.
    :param ref_fasta: `pysam.FastaFile`.
    """
    aln_header = pysam.AlignmentHeader.from_references(
        ref_fasta.references, ref_fasta.lengths)
    alns = [
        align_chunk_to_ref(polished_chunk, ref_fasta, aln_header=aln_header)
        for polished_chunk in pysam.FastxFile(consensus_fx)
    ]
    with pysam.AlignmentFile(bam_fp, 'wb', header=aln_header) as bam_fh:
        for aln in sorted(
            alns, key=lambda x: (x.reference_id, x.reference_start)
        ):
            bam_fh.write(aln)
    pysam.index(bam_fp)


#########################
# abpoa-related functions
#########################


def abpoa_consensus(aligner, subreads, max_n_cons=2):
    """Run abpoa consensus."""
    records = [RecordName.from_str(s.name) for s in subreads]
    seqs = [
        s.seq if r.strand == 'fwd' else medaka.common.reverse_complement(s.seq)
        for s, r in zip(subreads, records)]

    result = aligner.msa(
        seqs, out_cons=True, out_msa=False, max_n_cons=max_n_cons,
        min_freq=0.3)
    return result


def get_abpoa_clusters(subreads, record, put_bam_hp_in_name=True):
    """Cluster / phase reads using abpoa, checking for order dependence."""
    def rle_seq(seq):
        return ''.join(medaka.common.rle(seq)['value'])
    # run-length compress sequences
    subreads_rle = []
    for s in subreads:
        rn = RecordName.from_str(s.name)
        subreads_rle.append(medaka.smolecule.Subread(str(rn), rle_seq(s.seq)))
    # sort reads by length, see https://github.com/yangao07/abPOA/issues/48
    subreads_rle_asc = sorted(subreads_rle, key=lambda s: len(s.seq))
    subreads_rle_dsc = sorted(subreads_rle, key=lambda s: -len(s.seq))
    aligner = abpoa.msa_aligner(aln_mode='g')
    result_rle_asc = abpoa_consensus(aligner, subreads_rle_asc, max_n_cons=2)
    result_rle_dsc = abpoa_consensus(aligner, subreads_rle_dsc, max_n_cons=2)

    d = check_cluster_read_partitioning(
        result_rle_asc, result_rle_dsc, subreads_rle_asc, subreads_rle_dsc)

    subreads = {s.name: s for s in subreads}  # orig non-rle compressed reads
    clustered_subreads = {}
    for h in range(3):
        d.update(summarize_reads(d[f'hap{h}_reads'], prefix=f'hap{h}_'))
        clustered_subreads[h] = []
        for name in d[f'hap{h}_reads']:
            s = subreads[name]
            # optionally move SNP-based phasing HP tags into query name
            # and replace with cluster id with
            rn = RecordName.from_str(name)
            if put_bam_hp_in_name:
                rn.query_name += f'_BHP{rn.hap}'
            rn.hap = h
            clustered_subreads[h].append(
                medaka.smolecule.Subread(str(rn), s.seq)
            )
        del d[f'hap{h}_reads']

    return d, clustered_subreads


def process_record_abpoa(bam_fp, rec, ref_fasta, min_depth, min_mapq):
    """Process a record using abpoa for phasing."""
    subreads = get_subreads(bam_fp, rec, min_mapq=min_mapq)
    if len(subreads) < min_depth:
        raise InsufficientCoverage(
            f"{rec}: Retrieved too few reads ({len(subreads)} < {min_depth})")

    d, clustered_reads = get_abpoa_clusters(subreads, rec)
    cons = []
    for h in [1] if d['is_homozygous'] else [1, 2]:
        cons_rec = copy(rec)
        cons_rec.hap = h
        cons_rec.query_name += '_HOM' if d['is_homozygous'] else '_HET'
        cons.append(
            consensus_pileup_from_reads(cons_rec,
                                        clustered_reads[h],
                                        ref_fasta,
                                        min_depth=min_depth)
        )
    d['phasing_method'] = 'abpoa'
    return d, cons, clustered_reads[0]


def process_record_prephased(bam_fp, rec, ref_fasta, min_depth, min_mapq):
    """Process a record using prephased reads (using HP tags in bam)."""
    subreads = get_subreads(bam_fp, rec, min_mapq=min_mapq)
    if len(subreads) < 2 * min_depth:  # min_depth per haplotype
        raise InsufficientCoverage(
            f"{rec}: Retrieved too few reads "
            f"({len(subreads)} < 2 * {min_depth})")
    subreads_by_hap = collections.defaultdict(list)
    for s in subreads:
        rn = RecordName.from_str(s.name)
        subreads_by_hap[rn.hap].append(s)
    cons = []
    for h in [1, 2]:
        cons_rec = copy(rec)
        cons_rec.hap = h
        cons.append(
            consensus_pileup_from_reads(cons_rec,
                                        subreads_by_hap[h],
                                        ref_fasta,
                                        min_depth=min_depth))
    d = dict()
    for h in range(3):
        d.update(summarize_reads(
            [s.name for s in subreads_by_hap[h]],
            prefix=f'hap{h}_', bhp_counts=False))
    d['phasing_method'] = 'prephased'
    return d, cons, subreads_by_hap[0]


def process_record_haploid(bam_fp, rec, ref_fasta, min_depth, min_mapq):
    """Process haploid regions e.g. sex chromosomes."""
    subreads = get_subreads(bam_fp, rec, min_mapq=min_mapq)
    if len(subreads) < min_depth:
        raise InsufficientCoverage(
            f"{rec}: Retrieved too few reads ({len(subreads)} < {min_depth})")
    cons_rec = copy(rec)
    h = 1
    cons_rec.hap = h
    cons_rec.query_name += '_HOM'
    c = consensus_pileup_from_reads(
        cons_rec, subreads, ref_fasta, min_depth=min_depth)
    d = summarize_reads(
        [s.name for s in subreads], prefix='', bhp_counts=False)
    d['phasing_method'] = 'haploid'
    ambig_reads = []
    return d, [c], ambig_reads


def process_record_hybrid(*args):
    """Try to using existing phasing, fall back to abpoa."""
    try:
        m, haps, ambig_reads = process_record_prephased(*args)
        for hap in haps:  # hap is ConsensusResult
            if isinstance(hap.exception, InsufficientCoverage):
                # insufficient depth on a single haplotype
                return process_record_abpoa(*args)
    except InsufficientCoverage:
        # less than 2 * min_depth in total
        # still try to phase with abpoa as if abpoa decides region is haploid,
        # can get away with min_depth not 2 * min_depth
        # TODO revisit this.
        return process_record_abpoa(*args)

    return m, haps, ambig_reads


def check_cluster_read_partitioning(res_a, res_b, subreads_a, subreads_b):
    """Assess dependence of read ordering on phasing and generate metrics."""
    # TODO check if clu_read_ids of each cluster are superset of subreads
    # have seen abpoa error messages suggesting some reads are not assigned a
    # cluster - these should be included in counts of ambigous reads but are
    # probably not currently being used in consensus calculations.
    assert res_a.n_seq == res_b.n_seq
    # subreads should be in a different order but check we have same reads
    assert set(s.name for s in subreads_a) == set(s.name for s in subreads_b)
    assert all([r.n_cons <= 2 for r in (res_a, res_b)])
    # check consensus seq edit distances to see if cluster ordering has flipped
    cluster_edits = np.zeros(shape=(2, 2), dtype=int)
    cluster_ovlps = [[set(), set()], [set(), set()]]

    def edit_dist(seq1, seq2):
        return edlib.align(seq1, seq2, mode='NW')['editDistance']

    for a, (a_inds, a_cons) in enumerate(
        zip(res_a.clu_read_ids, res_a.cons_seq)
    ):
        for b, (b_inds, b_cons) in enumerate(
            zip(res_b.clu_read_ids, res_b.cons_seq)
        ):
            a_reads = [subreads_a[i].name for i in a_inds]
            b_reads = [subreads_b[i].name for i in b_inds]
            cluster_ovlps[a][b].update(set(a_reads).intersection(b_reads))
            cluster_edits[a][b] = edit_dist(a_cons, b_cons)

    diag_edits = cluster_edits.trace()
    off_diag_edits = cluster_edits.sum() - diag_edits

    consensus_is_flipped = off_diag_edits < diag_edits
    if not consensus_is_flipped:
        cluster1_reads = cluster_ovlps[0][0]
        cluster2_reads = cluster_ovlps[1][1]
        ambig_reads = cluster_ovlps[0][1] | cluster_ovlps[1][0]
    else:
        cluster1_reads = cluster_ovlps[0][1]
        cluster2_reads = cluster_ovlps[1][0]
        ambig_reads = cluster_ovlps[0][0] | cluster_ovlps[1][1]
        diag_edits, off_diag_edits = off_diag_edits, diag_edits

    # get fraction of reads that have jumped between clusters upon alternative
    # ordering of reads
    fraction_ambig = len(ambig_reads) / res_a.n_seq
    # in an ideal case of no POA order dependence, diagonal edits (i.e. between
    # a single haplotype with alternative orderings) will be zero, so edits
    # ratio will be zero. Large values imply stronger order depedence of the
    # differences in consensus sequence between haplotypes.
    edits_ratio = diag_edits / off_diag_edits if diag_edits else 0

    is_homozygous = all([r.n_cons == 1 for r in (res_a, res_b)])

    # if we have two clusters, but all reads in second cluster are ambigous and
    # hence removed  such that second cluster has no reads, make this haploid
    empty_second_cluster = False
    if (
        max(res_a.n_cons, res_b.n_cons) > 1 and
        min(len(cluster1_reads), len(cluster2_reads)) == 0
    ):
        is_homozygous = True
        empty_second_cluster = True
        cluster1_reads = cluster1_reads | cluster2_reads | ambig_reads
        cluster2_reads = set()
        ambig_reads = set()

    n_same_tags, n_switched = None, None
    if not is_homozygous:
        # if we have SNP-based bam HP phasing, use to phase clusters.
        cluster_bhp_ovlps = np.zeros(shape=(2, 2), dtype=int)
        subreads_by_bhp = {1: set(), 2: set()}
        for name in (cluster1_reads | cluster2_reads):
            rn = RecordName.from_str(name)
            if rn.hap in subreads_by_bhp:  # don't want to keep hap0
                subreads_by_bhp[rn.hap].add(name)
        for cid, cluster_names in enumerate([cluster1_reads, cluster2_reads]):
            for bhp_id, bhp_names in subreads_by_bhp.items():
                cluster_bhp_ovlps[cid, bhp_id - 1] = len(cluster_names
                                                         & bhp_names)

        n_same_tags = cluster_bhp_ovlps.trace()
        n_switched = cluster_bhp_ovlps.sum() - n_same_tags

        if n_switched > n_same_tags:  # switch cluster1/2 to match HP tags
            cluster1_reads, cluster2_reads = cluster2_reads, cluster1_reads
            n_same_tags, n_switched = n_switched, n_same_tags

    return dict(
        hap1_reads=cluster1_reads,
        hap2_reads=cluster2_reads,
        hap0_reads=ambig_reads,
        is_homozygous=is_homozygous,
        empty_second_cluster=empty_second_cluster,
        fraction_ambig_reads=round(fraction_ambig, 3),
        # smaller edits_ratio implies consensus seqs robust to read ordering
        edits_ratio=round(edits_ratio, 3),
        diag_edits=diag_edits,
        nreads_cluster_phasing_matches_bhp=n_same_tags,
        nreads_cluster_phasing_switched_wrt_bhp=n_switched,
    )


####################
# pipeline functions
####################


@dataclasses.dataclass
class ConsensusResult:
    """Simple class for organizing consensus results and alignments."""

    rec: RecordName
    subreads: tuple
    consensus_seq: str = None
    consensus_alignments: tuple = None
    ref_seq: str = None
    ref_alignments: tuple = None
    exception: Exception = None


def consensus_pileup_from_reads(
    rec, subreads, ref_fasta, min_depth=3, align_to_ref=False
):
    """Run consensus and generate read-consensus alignments."""
    if isinstance(rec, str):
        rec = RecordName.from_str(rec)

    res = ConsensusResult(rec, subreads=tuple(subreads))

    if len(res.subreads) < min_depth:  # definately insufficient coverage
        # though if partial == True could still be insufficient depth so need
        # extra checks later
        res.exception = InsufficientCoverage(
            f"{rec}: Too few reads ({len(res.subreads)} < {min_depth})")
        return res

    if isinstance(ref_fasta, str):
        ref_fasta = pysam.FastaFile(ref_fasta)

    res.ref_seq = ref_fasta.fetch(
        rec.ref_name, rec.ref_start_padded, rec.ref_end_padded)

    consensus_read = medaka.smolecule.Read(
        name=str(rec), subreads=res.subreads)
    # avoid realignment to discover strandedness
    consensus_read._orient = [
        RecordName.from_str(s.name).strand == 'fwd' for s in res.subreads
    ]
    consensus_read._initialized = True

    res.consensus_seq = consensus_read.poa_consensus(method='abpoa')

    # estimate depth of coverage from number of bp - and round to be generous
    depth = sum(len(s) for s in subreads)
    if round(depth) < min_depth:
        res.exception = InsufficientCoverage(
            f"{rec}: Estimated coverage too low ({depth:.2f} < {min_depth})")
        return res

    # use global aligner
    consensus_read.parasail_aligner_name = 'nw_trace_scan_32'

    # align sub-reads (i.e. trimmed pileup) to consensus
    res.consensus_alignments = consensus_read.align_to_template(
        template=res.consensus_seq,
        template_name=consensus_read.name,
    )

    # align sub-reads (i.e. trimmed pileup) to reference (for debug purposes)
    res.ref_alignments = []
    if align_to_ref:
        for aln in consensus_read.align_to_template(
            res.ref_seq, rec.ref_name,
        ):
            aln_dict = aln._asdict()
            # shift rstart to entire chromosome coords
            aln_dict['rstart'] += rec.ref_start_padded
            res.ref_alignments.append(medaka.smolecule.Alignment(**aln_dict))

    return res


def summarize_reads(names, prefix='', bhp_counts=False):
    """Generate counts of reads by strand and optionally bam HP tags."""
    records = [RecordName.from_str(n) for n in names]
    counts = collections.Counter()

    if bhp_counts:
        # we want to have 0 counts so initialize
        for h in range(3):
            counts[f'{prefix}n_reads_BHP{h}'] = 0
        counts.update([f'{prefix}n_reads_BHP{r.hap}' for r in records])

    for strand in 'fwd', 'rev':
        counts[f'{prefix}n_reads_{strand}'] = 0
    counts.update([f'{prefix}n_reads_{r.strand}' for r in records])

    counts[f'{prefix}n_reads'] += len(names)

    return dict(counts)


def bam_to_vcfs(bam_fp, ref_fasta):
    """Decode variants from alignments."""
    if isinstance(ref_fasta, str):
        ref_fasta = pysam.FastaFile(ref_fasta)

    logger = medaka.common.get_named_logger('BAM2VCF')
    contigs = [
        '{},length={}'.format(c, l)
        for c, l in zip(ref_fasta.references, ref_fasta.lengths)]

    meta = [
        medaka.vcf.MetaInfo(
            'INFO', 'rec', 1, 'String',
            'Medaka name for haplotype-specific consensus record.'),
        medaka.vcf.MetaInfo('FORMAT', 'GT', 1, 'String', 'Medaka genotype.'),
        medaka.vcf.MetaInfo(
            'FORMAT', 'PS', 1, 'Integer', 'Phase set identifier.'),
    ]

    prefix, ext = os.path.splitext(bam_fp)

    vcf = f'{prefix}.TR.vcf'
    logger.info(f'Writing variants to {vcf}')

    with medaka.vcf.VCFWriter(vcf, contigs=contigs, meta_info=meta) as vcfout:
        with pysam.AlignmentFile(bam_fp) as bam:
            for chrom in medaka.common.loose_version_sort(bam.references):
                rseq = ref_fasta.fetch(chrom)
                variants_by_pos = collections.defaultdict(list)
                for aln in bam.fetch(chrom):
                    # reference_start is 0 based, reference_end points to one
                    # past the last aligned residue, i.e. same as bed file
                    rn = RecordName.from_str(aln.query_name)
                    # clustered consensus results have _HOM/_HET suffixes and
                    # put all regions in their own phase blocks
                    if rn.query_name.endswith('_HOM'):
                        gts = ['0|1', '1|0']
                    elif rn.hap == 1:
                        gts = ['1|0']
                    elif rn.hap == 2:
                        gts = ['0|1']
                    else:
                        raise ValueError('Could not derive GT from {rec}')

                    for v in medaka.variant.yield_variants_from_aln(
                        aln, rseq, rn.to_unpadded_region(),
                    ):
                        for gt in gts:
                            v_gt = deepcopy(v)
                            v_gt.genotype_data['GT'] = gt
                            v_gt.info['rec'] = aln.query_name
                            v_gt.genotype_data['PS'] = rn.ref_start
                            variants_by_pos[v_gt.pos].append(v_gt)
                for pos, variants in sorted(variants_by_pos.items()):
                    for v in sorted(
                        variants, key=lambda x: x.gt, reverse=True
                    ):
                        vcfout.write_variant(v)


phasing_options = {
    'prephased': process_record_prephased,
    'abpoa': process_record_abpoa,
    'hybrid': process_record_hybrid,
}


def record_to_process_func(
    record, phasing_method, sex, par_regions, sex_chroms, logger
):
    """Determime which function should be used to process a given region."""
    region = record.to_unpadded_region()
    chr_x_name, chr_y_name = sex_chroms
    if sex == 'female':
        # exlude Y regions, but otherwise process everything as diploid
        if record.ref_name == chr_y_name:
            logger.info(
                (f'Skipping {record.to_unpadded_region()} on '
                 f'{chr_y_name} as sample is female.'))
            return
        return phasing_options[phasing_method]
    else:
        # for male samples, process X and Y as haploid with exception of
        # diploid PAR X regions.
        if any((region.overlaps(par_reg) for par_reg in par_regions)):
            logger.info(f'{region} is PAR, treating as diploid')
            return phasing_options[phasing_method]
        elif region.ref_name in sex_chroms:
            logger.info(f'{region} is sex chrom, so haploid for male.')
            return process_record_haploid
        else:
            return phasing_options[phasing_method]


def main(args):
    """Entry point for targeted tandem repeat variant calling."""
    logger = medaka.common.get_named_logger("TR")
    out_dir = args.output  # args.output will be later changed
    logger.info(f"Running medaka tr with options: {' '.join(sys.argv)}")
    if args.phasing in {'abpoa', 'hybrid'} and abpoa is None:
        # need to install pyabpoa which is not a formal requirement due to
        # https://github.com/yangao07/abPOA/issues/41#issuecomment-1400837516
        raise RuntimeError(
            'abpoa is not installed. Please install by running:\n'
            'pip install cython && pip install pyabpoa')

    medaka.common.mkdir_p(out_dir, info='Results will be overwritten.')

    ref_fasta = pysam.FastaFile(args.ref_fasta)

    contig_lengths = dict(zip(ref_fasta.references, ref_fasta.lengths))

    records = [RecordName(
            query_name='tr',
            ref_name=r.ref_name,
            ref_start=r.start,
            ref_end=r.end,
            ref_start_padded=max(r.start - args.pad, 0),
            ref_end_padded=min(
                r.end + args.pad, contig_lengths[r.ref_name]),
            hap=0,
            ) for r in args.regions
        ]

    header = {'HD': {'VN': 1.0}, 'SQ': []}
    consensuses = []
    all_consensus_alignments = []
    all_ref_alignments = []
    all_trimmed_reads = []
    all_ref_seqs = []
    records_skipped = []

    metrics_fhs = {
        k: open(os.path.join(out_dir, f'{k}_region_metrics.txt'), 'w')
        for k in ('prephased', 'abpoa', 'haploid')
    }

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=args.poa_threads
    ) as ex:
        futures_to_recs = {}
        for rec in records:
            # get appropriate function to process each region based on phasing
            # method with special handling of sex chromosomes
            process_func = record_to_process_func(
                rec, args.phasing, args.sex, args.PAR_regions,
                args.sex_chroms, logger)
            if process_func is not None:
                fut = ex.submit(
                    process_func, args.bam, rec, args.ref_fasta, args.depth,
                    args.min_mapq
                )
                futures_to_recs[fut] = rec
        t0 = now()
        t1, tlast = t0, t0
        n_done = 1
        for fut in concurrent.futures.as_completed(futures_to_recs.keys()):
            n_done += 1
            t1 = now()
            if t1 - tlast > 10:
                tlast = t1
                logger.info(
                    f"{n_done/len(records):.1%} Done "
                    f"({n_done}/{len(records)} regions) in {t1-t0:.1f}s"
                )
            rec = futures_to_recs[fut]
            e = fut.exception()
            if e is None:
                metrics = rec.to_unpadded_region()._asdict()
                (m, haps, ambig_reads) = fut.result()
                method = m.pop('phasing_method')
                metrics.update(m)
                if not metrics_fhs[method].tell():  # file empty
                    metrics_fhs[method].write("\t".join(metrics.keys()) + "\n")
                metrics_fhs[method].write(
                    "\t".join((str(v) for v in metrics.values())) + "\n")
                for hap in haps:  # hap is ConsensusResult
                    logger.debug(f"Generated poa consensus for {hap.rec}")
                    if isinstance(hap.exception, InsufficientCoverage):
                        # insufficient depth on a single haplotype
                        logger.info(hap.exception.args[0])
                        records_skipped.append(hap.rec)
                        continue
                    header['SQ'].append({
                        'LN': len(hap.consensus_seq),
                        'SN': str(hap.rec)})
                    all_consensus_alignments.append(hap.consensus_alignments)
                    # TODO write as much as possible on the fly rather than
                    # appending to lists. only thing need to accumulate is
                    # consensus alignments as don't have header
                    consensuses.append([str(hap.rec), hap.consensus_seq])
                    all_ref_alignments.append(hap.ref_alignments)
                    all_trimmed_reads += list(hap.subreads)
                    all_ref_seqs.append((str(hap.rec), hap.ref_seq))
                all_trimmed_reads += list(ambig_reads)
            elif isinstance(e, InsufficientCoverage):
                # insufficient depth on both haplotypes
                logger.info(e.args[0])
                rec = RecordName.from_str(e.args[0].split(':')[0])
                records_skipped.append(rec)
            else:
                logger.info(
                    f'Encountered exception whilst processing {rec}: {e}')
                records_skipped.append(rec)

        logger.info(
            f"Created {len(consensuses)} consensus with "
            f"{len(all_consensus_alignments)} alignments."
        )
    for fh in metrics_fhs.values():
        fh.close()

    # write bed file of regions skipped due to low depth
    skipped_bed_fp = os.path.join(out_dir, 'skipped.bed')
    with open(skipped_bed_fp, 'w') as fh:
        for r in records_skipped:
            fh.write(f"{r.ref_name}\t{r.ref_start}\t{r.ref_end}\t{r}\n")

    # write bam of trimmed reads aligned to poa consensus
    consensus_bam_file = os.path.join(out_dir, 'trimmed_reads_to_poa.bam')

    logger.info(
        "Writing trimmed reads to poa draft medaka input bam for "
        f"{len(all_consensus_alignments)} to {consensus_bam_file}."
    )

    medaka.smolecule.write_bam(
        consensus_bam_file, all_consensus_alignments, header)

    # write bam of trimmed reads aligned to reference
    # ref_bam_file = os.path.join(out_dir, 'trimmed_reads_to_ref.bam')
    # if isinstance(ref_fasta, str):
    #     ref_fasta = pysam.FastaFile(ref_fasta)
    #
    # ref_header = pysam.AlignmentHeader.from_references(
    #     ref_fasta.references, ref_fasta.lengths)
    #
    # sorted_alns = sorted(
    #     itertools.chain.from_iterable(all_ref_alignments),
    #     key=lambda x: (ref_header.get_tid(x.rname), x.rstart, x.qname))
    # medaka.smolecule.write_bam(ref_bam_file, [sorted_alns], ref_header)

    poa_file = os.path.join(out_dir, 'poa.fasta')
    logger.info(f"Writing poa consensus sequences to {poa_file}.")
    with open(poa_file, 'w') as fh:
        for rname, cons in consensuses:
            fh.write('>{}\n{}\n'.format(rname, cons))

    trimmed_reads_file = os.path.join(out_dir, 'trimmed_reads.fasta')
    logger.info(f"Writing trimmed reads to {trimmed_reads_file}.")
    with open(trimmed_reads_file, 'w') as fh:
        for read in all_trimmed_reads:
            fh.write(f'>{read.name}\n{read.seq}\n')

    ref_chunks_file = os.path.join(out_dir, 'ref_chunks.fasta')
    logger.info(f"Writing reference chunks to {ref_chunks_file}.")
    with open(ref_chunks_file, 'w') as fh:
        for reg_name, ref_seq in all_ref_seqs:
            fh.write(f'>{reg_name}\n{ref_seq}\n')

    # align poa consensus to ref
    poa_bam = os.path.join(out_dir, 'poa_to_ref.bam')
    align_consensus_fx_to_ref(poa_file, poa_bam, ref_fasta)
    bam_to_vcfs(poa_bam, ref_fasta)

    if args.poa_only:
        return

    logger.info("Running medaka consensus.")
    args.bam = consensus_bam_file
    args.output = os.path.join(out_dir, 'consensus.hdf')
    args.regions = None
    args.tag_name = None
    args.tag_value = None
    args.tag_keep_missing = False
    args.RG = None  # TODO plug this in and use in sample name

    # we run this in a subprocess so GPU resources are all cleaned
    # up when things are finished
    with concurrent.futures.ProcessPoolExecutor() as ex:
        fut = ex.submit(medaka.prediction.predict, args)
        _ = fut.result()

    logger.info("Running medaka stitch.")
    args.draft = poa_file
    args.inputs = [os.path.join(out_dir, 'consensus.hdf')]
    args.qualities = False
    out_ext = 'fasta'
    args.output = os.path.join(out_dir, 'consensus.{}'.format(out_ext))
    args.regions = None  # medaka consensus changes args.regions
    args.fillgaps = True
    args.fill_char = None  # fill with poa sequence where medaka did not polish
    args.min_depth = 0
    medaka.stitch.stitch(args)
    logger.info(f"Medaka consensus sequences written to {args.output}")
    # align consensus back to ref
    medaka_bam = os.path.join(out_dir, 'medaka_to_ref.bam')
    align_consensus_fx_to_ref(args.output, medaka_bam, ref_fasta)
    bam_to_vcfs(medaka_bam, ref_fasta)
