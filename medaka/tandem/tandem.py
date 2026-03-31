"""Variant calling in tandem repeats."""

import importlib.metadata
import logging
import os
import re
import sys

from packaging.version import Version
import pysam

from medaka import abpoa
import medaka.common
import medaka.medaka
import medaka.models
from medaka.tandem.consensus_generator import (
    ConsensusGenerator,
    ParallelConsensusGenerator,
)
from medaka.tandem.io import bam_to_vcfs
from medaka.tandem.record_name import RecordName
from medaka.tandem.spanning_read_clusterer import SpanningReadClusterFactory

PYABPOA_REQUIRED_VERSION = "1.5.6"


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
    if abpoa is None:
        raise RuntimeError(
            "pyabpoa=={} is required for medaka tandem, but pyabpoa is not "
            "installed. See medaka/tandem/README.md "
            "(Requirements section).".format(PYABPOA_REQUIRED_VERSION)
        )

    version = importlib.metadata.distribution("pyabpoa").version
    if Version(version) != Version(PYABPOA_REQUIRED_VERSION):
        raise RuntimeError(
            "pyabpoa=={} is required for medaka tandem, but found "
            "pyabpoa=={}. See medaka/tandem/README.md "
            "(Requirements section)."
            .format(
                PYABPOA_REQUIRED_VERSION,
                version,
            )
        )


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


def normalize_model_name(model: str) -> str:
    """Normalize a model path/name for policy checks."""
    model_name = os.path.basename(str(model or "")).lower()
    for suffix in medaka.models.model_suffixes:
        if model_name.endswith(suffix):
            return model_name.removesuffix(suffix)
    return model_name


def _has_model_token(normalized_model: str, token: str) -> bool:
    """Check for model token using underscore-aware boundaries."""
    return bool(re.search(rf"(^|_){re.escape(token)}(_|$)", normalized_model))


def _is_read_level_model(model: str) -> bool:
    """Check if model belongs to the read-level family."""
    return check_read_level_model(model)


def _is_variant_model(normalized_model: str) -> bool:
    """Check if model belongs to the variant family."""
    return _has_model_token(normalized_model, "variant")


def _is_snp_model(normalized_model: str) -> bool:
    """Check if model belongs to the SNP family."""
    return _has_model_token(normalized_model, "snp")


def _is_bacterial_model(normalized_model: str) -> bool:
    """Check if model belongs to the bacterial family."""
    return _has_model_token(normalized_model, "bacterial")


def validate_tandem_model(model: str):
    """Validate that model is supported for tandem polishing."""
    normalized_model = normalize_model_name(model)
    if not any((
        _is_read_level_model(model),
        _is_variant_model(normalized_model),
        _is_snp_model(normalized_model),
        _is_bacterial_model(normalized_model),
    )):
        return
    raise RuntimeError(
        "The selected model is incompatible with Medaka Tandem. "
        "Medaka Tandem supports consensus models only; "
        f"detected unsupported model '{model}'. "
        "Please choose a consensus model name or a basecaller model "
        "suffixed with ':consensus'."
    )


def validate_tandem_inputs(
    bam_path: str, ref_fasta_path: str
) -> None:
    """Validate tandem BAM/CRAM and reference input paths."""
    if not os.path.exists(ref_fasta_path):
        raise RuntimeError(
            f"Reference FASTA for tandem does not exist: {ref_fasta_path}"
        )
    try:
        with pysam.FastaFile(ref_fasta_path):
            pass
    except Exception as exc:
        raise RuntimeError(
            "Reference FASTA for tandem could not be opened/read: {}. "
            "Check that the FASTA is valid and indexed (.fai; for "
            "bgzipped FASTA, also .gzi).".format(ref_fasta_path)
        ) from exc

    if not os.path.exists(bam_path):
        raise RuntimeError(
            f"Alignment file for tandem does not exist: {bam_path}"
        )

    bam_open_kwargs = {}
    if str(bam_path).lower().endswith(".cram"):
        bam_open_kwargs["reference_filename"] = ref_fasta_path
    try:
        with pysam.AlignmentFile(bam_path, **bam_open_kwargs) as bam:
            # Touch references to ensure file/header are readable.
            _ = bam.references
            try:
                has_index = bool(bam.check_index())
            except Exception:
                has_index = False
    except Exception as exc:
        raise RuntimeError(
            f"Alignment file for tandem could not be read: {bam_path}"
        ) from exc

    if not has_index:
        raise RuntimeError(
            "Alignment index for 'bam' is missing or unreadable: {}. "
            "Provide an indexed BAM (.bai/.csi) or CRAM (.crai).".format(
                bam_path
            )
        )


def cli_error(message: str, code: int = 2) -> None:
    """Print a concise CLI error and exit."""
    print(f"medaka tandem: error: {message}", file=sys.stderr)
    raise SystemExit(code)


def get_effective_regions(
    regions: list[medaka.common.Region], sex: str, sex_chromosomes: list
) -> list[medaka.common.Region]:
    """Apply tandem-specific region filtering and deduplicate."""
    if sex == "female":
        _, chr_y_name = sex_chromosomes
        regions = [
            region for region in regions if region.ref_name != chr_y_name
        ]

    if not regions:
        raise RuntimeError(
            "No tandem regions remain after input parsing/filtering; "
            "please provide a non-empty BED file."
        )
    # keep input order while removing duplicate intervals
    return list(dict.fromkeys(regions))


def main(args):
    """Entry point for targeted tandem repeat variant calling."""
    tandem_formatter = logging.Formatter(
        '[%(asctime)s - %(levelname)s - %(name)s] %(message)s',
        datefmt='%H:%M:%S')
    for handler in logging.getLogger().handlers:
        handler.setFormatter(tandem_formatter)

    logger = medaka.common.get_named_logger("TR")
    # Suppress medaka inference progress logs in tandem runs.
    medaka.common.get_named_logger("Predict").setLevel(logging.WARNING)
    out_dir = args.output  # args.output will be later changed
    logger.info(f"Running medaka tr with options: {' '.join(sys.argv)}")
    try:
        validate_tandem_inputs(args.bam, args.ref_fasta)
    except RuntimeError as exc:
        cli_error(str(exc))
    check_abpoa_version()
    validate_tandem_model(args.model)
    logger.info("Using polishing model: %s", args.model)
    with pysam.FastaFile(args.ref_fasta) as ref_fasta:
        contig_lengths = dict(zip(ref_fasta.references, ref_fasta.lengths))

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    elif not os.path.isdir(out_dir):
        raise NotADirectoryError(
            f"Output path exists and is not a directory: {out_dir}"
        )
    else:
        with os.scandir(out_dir) as entries:
            if next(entries, None) is not None:
                logger.warning(
                    "The path %s exists and is not empty. "
                    "Results may be overwritten.",
                    out_dir,
                )

    consensus_bam_file = os.path.join(out_dir, "trimmed_reads_to_poa.bam")
    poa_file = os.path.join(out_dir, "poa.fasta")
    output_fasta = os.path.join(out_dir, "consensus.{}".format("fasta"))

    spanning_read_clusterer = SpanningReadClusterFactory.create_clusterer(
        args.phasing,
        min_depth=args.min_depth,
        remove_outliers=not args.disable_outlier_filter
    )

    args.regions = get_effective_regions(args.regions, args.sex, args.sex_chrs)
    par_regions = [
        medaka.common.Region.from_string(r) for r in args.par_regions
    ]
    records: list[RecordName] = [
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
        raise RuntimeError(
            "Tandem processing failed. Please check the logs for details."
        )

    if medaka.common.is_file_empty(poa_file) or \
            medaka.common.is_file_empty(output_fasta):
        msg = (
            "Medaka failed to generate a consensus sequence "
            "for the input regions."
        )
        logger.error(msg)
        raise RuntimeError(msg)

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
