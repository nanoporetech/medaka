"""Variant calling in tandem repeats."""

import os
import re
import sys
from typing import List

from packaging.version import Version
import pkg_resources
import pysam

from medaka import abpoa
import medaka.common
import medaka.medaka
from medaka.tandem.consensus_generator import (
    ConsensusGenerator,
    ParallelConsensusGenerator,
)
from medaka.tandem.io import bam_to_vcfs
from medaka.tandem.record_name import RecordName
from medaka.tandem.spanning_read_clusterer import SpanningReadClusterFactory


def determine_ploidy(
    record, phasing, sex, sex_chromosomes, par_regions
) -> int:
    """Determine the ploidy of a genomic record.

    Args:
        record: The genomic record.
        phasing: The phasing information.
        sex: The sex of the individual.
        sex_chromosomes: The sex chromosomes.
        par_regions: The PAR (pseudoautosomal region) of the genome.

    Returns:
        int: The ploidy of the genomic record.

    """
    if phasing == "unphased":
        return 1
    if record.ref_name not in sex_chromosomes:
        return 2
    if sex == "female":
        _, chr_y_name = sex_chromosomes
        if record.ref_name == chr_y_name:
            raise Exception(
                f"Can't determine ploidy for {chr_y_name} for female samples"
            )
        return 2
    if sex == "male":
        if any(record.overlaps(par_reg) for par_reg in par_regions):
            logger = medaka.common.get_named_logger("TR")
            logger.debug(f"{record} is PAR, treating as diploid")
            return 2
        else:
            return 1


def check_abpoa_version():
    """
    Check if pyabpoa is installed.

    This function verifies that the pyabpoa library is installed and that its
    version is greater than or equal to the specified minimum version. If the
    library is not installed or the version is insufficient, a RuntimeError
    is raised with instructions for installation or updating.

    Raises:
        RuntimeError: If pyabpoa is not installed or its version is less than
                      the required minimum version.
    """
    abpoa_required_version = '1.5.1'
    if abpoa is None:
        raise RuntimeError(
            f'pyabpoa == {abpoa_required_version} is not installed. Refer to '
            'abpoa documentation at https://github.com/yangao07/abPOA '
            'for installation instructions.')
    else:
        version = pkg_resources.get_distribution("pyabpoa").version
        if Version(version) != Version(abpoa_required_version):
            raise RuntimeError(
                f"pyabpoa == {abpoa_required_version} required, got {version}")


def check_read_level_model(model: str) -> bool:
    """Check if the provided model is a read-level model.

    Args:
        model (str): The path to the model file.

    Returns:
        bool: True if the model is a read-level model, False otherwise.
    """
    read_level_regex = (
        r"r\d{4}_e\d{2}_\d{3}bps_[a-z]+(_v\d+\.\d+\.\d+)?_rl_lstm384_"
        r"(dwells|no_dwells)_model_pt\.tar\.gz"  # noqa: E501
    )
    return bool(re.search(read_level_regex, str(model)))


def main(args):
    """Entry point for targeted tandem repeat variant calling."""
    logger = medaka.common.get_named_logger("TR")
    out_dir = args.output  # args.output will be later changed
    logger.info(f"Running medaka tr with options: {' '.join(sys.argv)}")
    check_abpoa_version()
    if check_read_level_model(args.model):
        logger.error(
            "The provided model is incompatible "
            "with tandem repeat genotyping. "
            "Please use a compatible model."
        )
        logger.error(f"Model: {args.model}")
        return

    medaka.common.mkdir_p(out_dir, info="Results will be overwritten.")

    consensus_bam_file = os.path.join(out_dir, "trimmed_reads_to_poa.bam")
    poa_file = os.path.join(out_dir, "poa.fasta")
    output_fasta = os.path.join(out_dir, "consensus.{}".format("fasta"))

    with pysam.FastaFile(args.ref_fasta) as ref_fasta:
        contig_lengths = dict(zip(ref_fasta.references, ref_fasta.lengths))

    spanning_read_clusterer = SpanningReadClusterFactory.create_clusterer(
        args.phasing,
        min_depth=args.min_depth,
        remove_outliers=not args.disable_outlier_filter
    )

    if args.sex == "female":
        _, chr_y_name = args.sex_chrs
        args.regions = [
            region for region in args.regions if region.ref_name != chr_y_name
        ]
    # remove duplicates
    args.regions = set(str(r) for r in args.regions)
    args.regions = [medaka.common.Region.from_string(r) for r in args.regions]
    par_regions = [
        medaka.common.Region.from_string(r) for r in args.par_regions
    ]
    records: List[RecordName] = [
        RecordName(
            query_name="tr",
            ref_name=r.ref_name,
            ref_start=r.start,
            ref_end=r.end,
            ref_start_padded=max(r.start - args.padding, 0),
            ref_end_padded=min(
                r.end + args.padding, contig_lengths[r.ref_name]
            ),
            hap=0,
            ploidy=determine_ploidy(
                r, args.phasing, args.sex, args.sex_chrs, par_regions
            ),
        )
        for r in args.regions
    ]
    contig_lengths = None  # free memory
    args.regions = []
    read_filters = {"min_mapq": args.min_mapq}

    generator_args = {
        "regions": records,
        "bam": args.bam,
        "ref": args.ref_fasta,
        "reads_clusterer": spanning_read_clusterer,
        "min_depth": args.min_depth,
        "reads_filter": read_filters,
        "process_large_regions": args.process_large_regions,
        "output_prefix": out_dir,
        "model": args.model,
    }
    if args.workers == 1:
        generator = ConsensusGenerator(**generator_args)
        no_processed_regions = generator.process()
        success = no_processed_regions == len(records)
    else:
        generator = ParallelConsensusGenerator(
            **generator_args,
            num_processes=args.workers,
        )
        del records  # remove records to free memory
        success = generator.process()

    if not success:
        return

    if medaka.common.is_file_empty(poa_file) or \
            medaka.common.is_file_empty(output_fasta):
        logger.error(
            "Medaka failed to generate a consensus sequence "
            "for the input regions."
        )
        return

    medaka_bam = os.path.join(out_dir, "medaka_to_ref.bam")
    bam_to_vcfs(
        medaka_bam,
        args.ref_fasta,
        trimmed_reads_to_poa=consensus_bam_file,
        replacement_style=not args.decompose,
        add_read_names=args.add_read_names,
        sample_name=args.sample_name,
    )
    return
