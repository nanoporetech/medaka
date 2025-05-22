"""IO modules for medaka tandem."""

from itertools import groupby
import os
from typing import List

import numpy as np
import pysam
import pysam.bcftools

import medaka.common
import medaka.features
import medaka.smolecule
from medaka.tandem.record_name import RecordName
import medaka.variant
import medaka.vcf


#########################
# read-fetching functions
#########################
# TODO move to featues / common and remove vcf.get_trimmed_reads

class SpanningReadsExtractor:
    """A class to extract and process reads spanning specified regions."""

    def __init__(self, bam_path: str, read_filters: dict):
        """Initialize the SpanningReadsExtractor with the BAM and read filters.

        Args:
            bam_path (str): Path to the BAM file.
            read_filters (dict): Filters to apply to reads, e.g., minimum mapq.

        """
        self.bam_handle = medaka.features.BAMHandler(bam_path, size=1)
        self.read_filters = read_filters
        self.bam_path = bam_path

    def get_subreads(self, rec) -> List[medaka.smolecule.Subread]:
        """Extract a sub-part of the read spanning the read record.

        Args:
            rec (RecordName): record descriping the tandem repeat region.

        Returns:
            List[medala.smolecule.Subread]: subreads spanning the record.

        """
        reg_padded = rec.to_padded_region()

        _, reads = self.get_trimmed_reads(reg_padded)
        if len(reads) == 0:
            return []

        rn_kwargs = {
            k: v
            for k, v in vars(rec).items()
            if k not in {"query_name", "strand", "hap", "phased_set"}
        }
        try:
            subreads = [
                medaka.smolecule.Subread(
                    str(
                        RecordName(
                            query_name=f"{read_name}",
                            strand="rev" if is_rev else "fwd",
                            hap=hap,
                            phased_set=phased_set,
                            **rn_kwargs,
                        )
                    ),
                    medaka.common.reverse_complement(seq) if is_rev else seq,
                )
                for is_rev, read_name, seq, hap, phased_set in reads
            ]

        except ValueError:
            subreads = []

        return subreads

    def get_trimmed_reads(self, region: medaka.common.Region):
        """Retrieve and trims reads spanning the specified region.

        Args:
            region (str): The genomic region in "chr:start-end" format.

        Returns:
            List[str]: A list of trimmed reads that span the specified region.

        """
        region_got, reads = next(
            medaka.features.get_trimmed_reads(
                region,
                self.bam_handle,
                partial=False,
                region_split=2 * region.size,
                include_empty_reads=True,
                **self.read_filters,
            ),
            (region, []),
        )
        if len(reads) == 0:
            raise ValueError(
                "No reads found for {} nor even reference sequence. "
                "Check Bam file {}".format(region, self.bam_path)
            )
        if not region == region_got:  # check region was not e.g. split
            msg = "Expected region {}, got region {}"
            raise ValueError(msg.format(region, region_got))
        # first read is ref, remove it.
        _, _, ref_seq, _, _ = reads.pop(0)
        return ref_seq, reads


#########################
# VCF export functions
#########################


def get_alt_from_aln(aln, record):
    """Extract alternate sequence covering the tandem repeat region.

    Args:
        aln: pysam.AlignedSegment .
        record: RecordName for the region.

    Returns:
        string containing alternate allele sequences  .

    """
    consensus_boundaries = [
        cons
        for cons, ref in aln.get_aligned_pairs()
        if cons is not None
        and ref is not None
        and record.ref_start <= ref <= record.ref_end
    ]
    if len(consensus_boundaries) == 0:
        return "<DEL>"
    elif len(consensus_boundaries) == 1:
        return aln.query_sequence[consensus_boundaries[0]]
    else:
        return aln.query_sequence[
            consensus_boundaries[0]: consensus_boundaries[-1]
        ]


def determine_gt_and_alleles(alignments, ref_seq):
    """Determine the genotype and alternate alleles from alignments.

    Compares the different alleles including the reference allele
      to determine the right gt and allele seuqnce for vcf output

    Args:
        alignments: array of one or two pysam.AlignedSegment.
        ref_seq: string containing the ref allele.

    Returns:
        A tuple of alternate allele sequences and a VCF genotype string.

    """
    if len(alignments) > 2:
        raise ValueError("More than two consensus sequences found.")
    rn = RecordName.from_str(alignments[0].query_name)
    alts = [get_alt_from_aln(aln, rn) for aln in alignments]

    alleles = set(alts + [ref_seq])
    # clustered consensus results have _HOM/_HET suffixes and
    if rn.query_name.endswith("_HOM"):
        if alts[0] == ref_seq:
            return ".", "0|0"
        else:
            return alts[0], "1|1"
    elif len(alleles) == 1:
        if len(alts) == 2:
            return (".", "0|0")
        else:
            return (".", "0|." if rn.hap == 1 else ".|0")
    elif len(alleles) == 2:
        if len(alts) == 1:
            return (alts, "1|." if rn.hap == 1 else ".|1")
        # Handle genotype based on reference presence in alternates
        genotype = f"{int(ref_seq != alts[0])}|{int(ref_seq != alts[1])}"
        return alts[1] if ref_seq == alts[0] else alts[0], genotype

    elif len(alleles) == 3:
        return alts, "1|2"
    else:
        raise ValueError("Impossible")


class RefFastaHandle:
    """Class to handle reference FASTA file."""

    def __init__(self, fasta_path: str):
        """Initialize the ReferenceFastaHandle class.

        Args:
            fasta_path (str): The path to the reference FASTA file.

        """
        self.fasta_path = fasta_path
        self.fasta_handle = pysam.FastaFile(fasta_path)
        self.sequence_cache = {}

    def get_reference_names_and_lengths(self):
        """Get the pair of reference names and lengths.

        Args:
            ref_fasta_handle (RefFastaHandle): A handle to the ref FASTA file.

        Returns:
            List[Tuple[str, int]]: A tuple of reference names and lengths.

        """
        return [
            (name, length)
            for name, length in zip(
                self.fasta_handle.references,
                self.fasta_handle.lengths
            )
        ]

    def get_sequence(self, name: str) -> str:
        """Get the sequence by name.

        Args:
            name (str): The name of the sequence.

        Returns:
            str: The sequence corresponding to the given name.

        """
        if name in self.sequence_cache:
            return self.sequence_cache[name]
        else:
            sequence = self.fasta_handle.fetch(name)
            self.sequence_cache[name] = sequence
            return sequence

    def clear_cache(self):
        """Clear the sequence cache."""
        self.sequence_cache.clear()

    @property
    def fasta_file(self) -> pysam.FastaFile:
        """Get the reference FASTA file handle."""
        return self.fasta_handle


def create_vcf_header_meta():
    """Create the VCF header meta information.

    Returns:
        list: A list of VCF header meta information.

    """
    return [
        medaka.vcf.MetaInfo(
            "INFO",
            "rec",
            ".",
            "String",
            "Medaka name for haplotype-specific consensus record.",
        ),
        medaka.vcf.MetaInfo("FORMAT", "GT", 1, "String", "Medaka genotype."),
        medaka.vcf.MetaInfo(
            "FORMAT", "PS", 1, "Integer", "Phase set identifier."
        ),
        medaka.vcf.MetaInfo(
            "FORMAT",
            "SD",
            ".",
            "Integer",
            "Number of spanning reads supporting each allele, "
            "reported separately per haplotype when phased. "
            "A single count is provided when haplotypes cannot "
            "be resolved (e.g., due to dropout or homozygous variants)."
        ),
        medaka.vcf.MetaInfo(
            "FORMAT",
            "MAD",
            ".",
            "Float",
            "Median Absolute Deviation of read lengths "
            "reported separately per haplotype when phased. "
            "A single count is provided when haplotypes cannot "
            "be resolved (e.g., due to dropout or homozygous variants)."
        ),
        medaka.vcf.MetaInfo(
            'FORMAT',
            'ALLR',
            ".",
            'String',
            "Allele length range "
            "reported separately per haplotype when phased. "
            "A single count is provided when haplotypes cannot "
            "be resolved (e.g., due to dropout or homozygous variants)."
        ),
        medaka.vcf.MetaInfo(
            "INFO",
            "read_names_hap1",
            "1",
            "String",
            "Names of supporting reads for hap1.",
        ),
        medaka.vcf.MetaInfo(
            "INFO",
            "read_names_hap2",
            "1",
            "String",
            "Names of supporting reads for hap2.",
        ),
        medaka.vcf.MetaInfo(
            "INFO",
            "read_names_hap0",
            "1",
            "String",
            "Names of supporting reads for sex chromosome",
        ),
    ]


def convert_alignments_to_variants(
    alignments,
    reads_bam,
    ref_fasta_handle,
    add_read_names,
    is_replacement_style,
):
    """Convert alignments to variants.

    Args:
        alignments (List[pysam.AlignedSegment]): alignments to be converted.
        reads_bam (pysam.AlignmentFile): The BAM path containing the reads.
        ref_fasta_handle (RefFastaHandle): A handle to the ref FASTA file.
        add_read_names (bool): Flag for adding read names to the vcf.
        is_replacement_style (bool): Flag for using replacement style.

    Returns:
        list: A list of variants converted from the alignments.

    Raises:
        None

    """
    alignments = sorted(
        alignments, key=lambda a: RecordName.from_str(a.query_name).hap
    )
    if is_replacement_style:
        return [
            convert_alignments_to_variants_replacement_style(
                alignments,
                reads_bam,
                add_read_names,
                ref_fasta_handle.fasta_file
            )
        ]
    else:
        ref_seq = ref_fasta_handle.get_sequence(alignments[0].reference_name)
        # I am using ref_fasta fetch in replacement style because it is faster.
        return convert_alignments_to_variants_decomposition(
            alignments, reads_bam, add_read_names, ref_seq
        )


def convert_alignments_to_variants_decomposition(
    alignments, reads_bam, add_read_names, rseq
):
    """Convert list of pysam.AlignedSegment to decomposed variants.

    Variants on the same positions but from different haplotypes
    will be written to the vcf on different records.

    Args:
        alignments (List[pysam.AlignedSegment]): A list of alignments.
        reads_bam (pysam.AlignmentFile): The BAM file containing the reads.
        add_read_names (bool): Flag to add read names to vcf.
        rseq (str): The reference sequence.

    Returns:
        List[Variant]: A list of variants.

    """
    results: List[medaka.vcf.Variant] = []
    for aln in alignments:
        rn = RecordName.from_str(aln.query_name)
        reads = list(reads_bam.fetch(aln.query_name))
        depth = len(reads)
        for v in medaka.variant.yield_variants_from_aln(
            aln,
            rseq,
            rn.to_unpadded_region(),
        ):
            v.genotype_data["SD"] = depth
            v.ident = (
                f"{rn.ref_name}_{rn.ref_start}_{rn.ref_end}_{v.pos}_"
                f"hap{rn.hap}"
            )
            if add_read_names:
                v.info[f"read_names_hap{rn.hap}"] = [
                    RecordName.from_str(trimmed_read.query_name).query_name
                    for trimmed_read in reads
                ]

            # Clusters variant by position,
            # using '.' for missing variants in haplotypes.
            # Assumes hap1 is processed first,
            # so initializes with '.' for hap2
            # to update if present in alignments.
            if rn.query_name.endswith("_HOM"):
                v.genotype_data["GT"] = "1|1"
            elif rn.hap == 1:
                v.genotype_data["GT"] = "1|0"
            elif rn.hap == 2:
                v.genotype_data["GT"] = "0|1"
            results.append(v)
    return results


def convert_alignments_to_variants_replacement_style(
    alignments, reads_bam, add_read_names, ref_fasta
):
    """Convert alignments to replacement style variants.

    Convert list of pysam.AlignedSegment to replacement style. Variants
    on the same positions but from different haplotypes
    will be printed on the same line.

    Args:
        alignments (List[pysam.AlignedSegment]): A list of alignments.
        reads_bam (pysam.AlignmentFile): The BAM file containing the reads.
        add_read_names (bool): Flag indicating whether to add read names
        to the variant information.
        ref_fasta (str): The reference sequence.

    Returns:
        List[Variant]: A list of variants.

    """
    format = {"SD": [], 'ALLR': []}
    info = {}
    depths = []
    reads_lengths = []
    mads = []
    chrom = alignments[0].reference_name
    for aln in alignments:
        hap = RecordName.from_str(aln.query_name).hap
        reads = list(reads_bam.fetch(aln.query_name))
        if add_read_names:
            info[f"read_names_hap{hap}"] = [
                RecordName.from_str(trimmed_read.query_name).query_name
                for trimmed_read in reads
            ]
        curr_reads_len = np.array([
            trimmed_read.infer_query_length()
            for trimmed_read in reads
        ])
        reads_lengths.append(
            f'{np.min(curr_reads_len)}'
            f'-{np.max(curr_reads_len)}'
        )
        med = np.median(curr_reads_len)
        mad = np.median(np.absolute(curr_reads_len - med))
        mads.append(f'{mad:.2f}')
        depths.append(str(len(reads)))

    format["SD"] = ",".join(depths)
    format['ALLR'] = ",".join(reads_lengths)
    format["MAD"] = ",".join(mads)

    rn = [RecordName.from_str(aln.query_name) for aln in alignments]
    ref = ref_fasta.fetch(chrom, start=rn[0].ref_start, end=rn[0].ref_end)
    alts, gt = determine_gt_and_alleles(alignments, ref)
    info["rec"] = [aln.query_name for aln in alignments]

    phase_sets = list({r.phased_set for r in rn})
    # check if both alleles have the same phaseset
    is_phased = len(phase_sets) == 1 and phase_sets[0] != 0
    is_phased &= not rn[0].query_name.endswith("_HOM")
    is_phased &= not rn[0].query_name.endswith("_HET")
    if is_phased:
        format["PS"] = phase_sets[0]
        format["GT"] = gt
    else:
        gt = gt.split("|")
        format["GT"] = "/".join(gt)

    variantID = f"{rn[0].ref_name}_{rn[0].ref_start}_{rn[0].ref_end}"

    return medaka.vcf.Variant(
        chrom=chrom,
        pos=rn[0].ref_start,
        ref=ref,
        alt=alts,
        ident=variantID,
        genotype_data=format,
        info=info,
    )


def bam_to_vcfs(
    bam_fp,
    ref_fasta,
    trimmed_reads_to_poa,
    *,
    replacement_style=False,
    add_read_names=False,
    sample_name="SAMPLE",
):
    """Decode variants from alignments."""
    ref_fasta_handle = RefFastaHandle(ref_fasta)
    logger = medaka.common.get_named_logger("BAM2VCF")

    contigs = [
        "{},length={}".format(chr, length)
        for chr, length in ref_fasta_handle.get_reference_names_and_lengths()
    ]
    meta = create_vcf_header_meta()
    prefix, ext = os.path.splitext(bam_fp)
    out_dir = os.path.dirname(bam_fp)
    vcf_temp = f"{prefix}.temp.vcf"
    vcf_final = f"{prefix}.TR.vcf"
    logger.info(f"Writing variants to {vcf_final}")
    header = (
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        sample_name,
    )
    with medaka.vcf.VCFWriter(
        vcf_temp, contigs=contigs, meta_info=meta, header=header
    ) as vcfout:
        with pysam.AlignmentFile(bam_fp) as bam:
            with pysam.AlignmentFile(trimmed_reads_to_poa) as reads_bam:
                for chrom in medaka.common.loose_version_sort(bam.references):
                    # Group alignments by tandem repeat regions
                    # so that different alleles will processed together
                    ref_fasta_handle.clear_cache()
                    for _, alignments in groupby(
                        bam.fetch(chrom),
                        key=lambda x:
                        (RecordName.from_str(x.query_name).ref_start,
                            RecordName.from_str(x.query_name).ref_end)
                    ):
                        variants = convert_alignments_to_variants(
                            alignments=alignments,
                            reads_bam=reads_bam,
                            ref_fasta_handle=ref_fasta_handle,
                            add_read_names=add_read_names,
                            is_replacement_style=replacement_style,
                        )
                        for v in variants:
                            vcfout.write_variant(v)

    pysam.bcftools.sort(
        "-o", vcf_final, "--temp-dir", out_dir, vcf_temp, catch_stdout=False
    )
    os.remove(vcf_temp)
