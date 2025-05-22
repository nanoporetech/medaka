"""alignments methods."""

import functools

import parasail
import pysam

import medaka.common
from medaka.tandem.record_name import RecordName


def align_chunk_to_ref(
    chunk: pysam.FastxRecord, ref_fasta: pysam.FastaFile, aln_header=None
) -> pysam.AlignedSegment:
    """Align consensus chunk to reference using global alignment.

    Args:
    chunk(pysam.FastxRecord): sequence to be aligned.
    name should be the string representation of a `RecordName`.
    ref_fasta(pysam.FastaFile): reference fasta file.
    aln_header(Optional[AligmentHeader]): optional alignment header.
    Otherwise, one will be created from the reference fasta.

    Returns :
    alignment : an instance of pysam.AlignedSegment.

    """
    logger = medaka.common.get_named_logger("TR_ALIGN")
    if aln_header is None:
        aln_header = pysam.AlignmentHeader.from_references(
            ref_fasta.references, ref_fasta.lengths
        )

    rn = RecordName.from_str(chunk.name)

    ref_seq = ref_fasta.fetch(
        rn.ref_name, rn.ref_start_padded, rn.ref_end_padded
    )
    if rn.strand == "fwd":
        query_seq = chunk.sequence
        flag = 0
    else:
        query_seq = medaka.common.reverse_complement(chunk.sequence)
        flag = 16

    # Case when the whole allele is deleted
    # and we can't find any sequence to align
    if query_seq == "":
        ref_size = rn.ref_end_padded - rn.ref_start_padded
        return medaka.align.initialise_alignment(
            query_name=chunk.name,
            reference_id=aln_header.get_tid(rn.ref_name),
            reference_start=rn.ref_start_padded,
            query_sequence=query_seq,
            cigarstring=f"{ref_size}D",
            flag=flag,
            header=aln_header,
            tags={"HP": rn.hap},
        )
    # use global alignment
    parasail_aligner = functools.partial(
        parasail.nw_trace_scan_32, open=8, extend=4, matrix=parasail.dnafull
    )
    result = parasail_aligner(query_seq, ref_seq)
    rstart, cigar = medaka.align.parasail_to_sam(result, query_seq)
    # rstart may not be zero if alignent starts with a SNP / indel.
    # TODO - do we want to automatically increase padding and rerun?
    if rstart > 0:
        logger.warning(
            f"rstart not 0 when using global alignment for {chunk.name}. "
            "Consider rerunning this region with more padding."
        )

    aln = medaka.align.initialise_alignment(
        query_name=chunk.name,
        reference_id=aln_header.get_tid(rn.ref_name),
        reference_start=rn.ref_start_padded + rstart,
        query_sequence=query_seq,
        cigarstring=cigar,
        flag=flag,
        header=aln_header,
        tags={"HP": rn.hap},
    )
    return aln


def align_consensus_fx_to_ref(
    consensus_fx: str, bam_fp: str, ref_fasta: pysam.FastaFile
) -> None:
    """Align consensus fastx to ref producing bam.

    Args:
    consensus_fx(string): filepath to consensus fastx file.
    bam_fp(string): output filepath of bam.
    ref_fasta(`pysam.FastaFile`): reference fasta file.

    Returns:
    None

    """
    ref_fasta = pysam.FastaFile(ref_fasta)
    aln_header = pysam.AlignmentHeader.from_references(
        ref_fasta.references, ref_fasta.lengths
    )
    temp_bam = bam_fp + "_temp.bam"
    with pysam.AlignmentFile(temp_bam, "wb", header=aln_header) as bam_fh:
        for polished_chunk in pysam.FastxFile(consensus_fx):
            aln = align_chunk_to_ref(
                polished_chunk, ref_fasta, aln_header=aln_header
            )
            bam_fh.write(aln)

    pysam.sort("-o", bam_fp, temp_bam)
    pysam.index(bam_fp)
