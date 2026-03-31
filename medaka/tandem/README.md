![Oxford Nanopore Technologies logo](https://github.com/nanoporetech/medaka/raw/master/images/ONT_logo_590x106.png)


Medaka Tandem
======

[![](https://img.shields.io/pypi/v/medaka.svg)](https://pypi.org/project/medaka/)
[![](https://img.shields.io/pypi/wheel/medaka.svg)](https://pypi.org/project/medaka/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/medaka)
[![](https://img.shields.io/conda/pn/bioconda/medaka.svg)](https://anaconda.org/bioconda/medaka)


Medaka Tandem is a tool for haplotype-aware genotyping of short tandem repeats (STRs) from Oxford Nanopore reads. It takes an aligned (and preferably haplotagged) BAM file and a BED file with tandem repeat regions as input and outputs VCF with STR alleles and auxiliary information (e.g., number of reads supporting each allele). It also produces additional outputs, e.g. FASTA files with predicted allele sequences and BAM files with reads aligned to respective predicted alleles for visual inspection.

Medaka Tandem carries out the following steps for each specified tandem repeat region:

1. **Haplotype separation** Reads spanning the region are identified and divided into haplotypes of origin (see the `--phasing` option description for more information).

2. **Draft consensus reconstruction** Then draft consensus sequence for each STR allele is reconstructed using abPOA library.

3. **Read re-alignment and consensus polishing.** Spanning reads are realigned to the corresponding draft allele sequence, followed by Medaka neural network based polishing.

4. **Allele reporting** Variety of output files are produced, including genotype information in VCF format. By default, the entire polished STR sequences are reported as alternative alleles. Optionally these sequences can be aligned back to the reference regions to derive a list of smaller variant calls (see the `--decompose` option).


Requirements
-----

In addition to usual Medaka dependencies, Medaka Tandem also requires [pyabpoa](https://github.com/yangao07/abPOA) v1.5.6 to be available:
- If medaka is being installed via pip, on Linux it can be installed with `pip install 'medaka[abpoa]'` but development libraries are required to build the `pyabpoa` wheel.
- If medaka was already installed via pip, `pyabpoa` can be added with `pip install pyabpoa==1.5.6`.
- If medaka was installed from conda, pyabpoa should have been installed as a dependency.
- If medaka was built from source with `make install`, pyabpoa was built as part of this installation and should already be present.


Usage
-----

`medaka tandem` command requires five positional arguments:
Positional arguments: bam, ref_fasta, regions, sex, output.
1. Read alignments (preferably haplotagged) in BAM/CRAM format; the alignment file must be indexed.
2. Reference genome in FASTA format.
3. Path to a BED or BED.gz file with genomic regions to be analysed.
4. Sample sex (`male|female`)
5. Output directory.

Adotto repeat catalogues can be downloaded from [here](https://github.com/ACEnglish/adotto/blob/main/regions/DataDescription.md).

Command-line example

```bash
medaka tandem \
  haplotagged.bam \
  ref.fna \
  adotto_TRregions_v1.2.bed \
  male \
  output_folder
```

## Additional options

### Performance-Related Options

1. `--workers` Number of parallel worker processes to use (default: 1).
2. `--process_large_regions` Process TRs with estimated length (of one or both alleles) exceeding 10kbp (default: False). Processing large regions can substantially increase RAM usage. With the default setting, the expected peak RAM consumption on Adotto repeat catalogue is approximately 14GB when using 8 workers, and 23GB when using 16 workers. Skipped regions will be output to `skipped_large.bed`.

### Behaviour-Related Options
1. `--phasing` Specify the strategy for dividing the reads between haplotypes. (default: `hybrid`)
    - `prephased` Rely on haplotype (HP) BAM tags for phasing.
    - `abpoa` Use abPOA clustering feature to identify haplotypes based on STR sequences in the reads.
    - `hybrid` Use haplotag assignments if both haplotypes have at least `min_depth` spanning reads assigned, otherwise fallback to the abPOA clustering.
    - `unphased` Assume the sample is haploid.
2. `--min_depth` Minimum number of spanning reads required for allele consensus reconstruction. (default: 3)
3. `--min_mapq` Minimum mapping quality (MAPQ) for alignments filtering. (default: 5)
4. `--disable_outlier_filter` Disable exclusion of reads with significantly divergent spanning region lengths. (default: False).
5. `--padding` Number of bases to pad spanning read regions and reference sequence. (default: 10)
6. `--sex_chrs` Names of X and Y chromosomes in reference FASTA as two arguments (e.g. `--sex_chrs chrX chrY`). (default: ['chrX', 'chrY'])
7. `--par_regions` Coordinates of pseudoautosomal regions (PARs) on the X chromosome. Will be treated as diploid in male samples. The analysis assumes that the corresponding PARs on chromosome Y have been hard-masked (i.e. replaced with Ns) in the reference to avoid ambiguous read alignments. (default: chrX:10000-2781479,chrX:155701382-156030895 assuming use of the GRCh38 analysis set, e.g. `GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta`)
8. `--model` Model to be used for polishing. Use a Medaka consensus model name (for example, `r1041_e82_400bps_hac_v5.2.0`), or a basecaller model name suffixed with `:consensus` (for example, `dna_r10.4.1_e8.2_400bps_hac@v4.1.0:consensus`).
9. `--auto_model TYPE INPUT` Automatically choose a model by inspecting `INPUT`. For Medaka Tandem, `TYPE` must be `consensus`. `INPUT` should be a basecaller output file containing basecaller-model metadata (BAM/CRAM/SAM with `@RG DS` header entries containing `basecall_model=...`, or FASTQ with basecaller model metadata in read comments).

Example (using the same BAM for both auto-model detection and tandem input):

```bash
BAM=haplotagged.bam
medaka tandem \
  --auto_model consensus "${BAM}" \
  "${BAM}" \
  ref.fna \
  adotto_TRregions_v1.2.bed \
  male \
  output_folder
```

`--model` and `--auto_model` are mutually exclusive.
Please review the [Medaka README](../../README.md) for more information about model selection.
Note: Medaka Tandem shares Medaka CLI help text for model arguments. Hence, while `--help` lists other model types (suffixed with `:variant` or `--auto_model variant`), Medaka Tandem can only be used with consensus models.

### Logging Options
1. `--quiet` Reduce logging verbosity to warnings and errors only.
2. `--debug` Enable verbose debug logging.
`--quiet` and `--debug` are mutually exclusive. By default, Medaka Tandem logs at INFO level.


### Output-Related Options

1. `--decompose` Align polished sequences back to reference region and extract a list of (left-aligned) variants. By default, Medaka Tandem reports entire haplotype-specific tandem repeats as alternative alleles.
2. `--add_read_names` Report names of spanning reads in the output VCF file.
3. `--sample_name` Sample name used in the output VCF file. (default: `SAMPLE`)

Output Folder Content
-----

1. `medaka_to_ref.TR.vcf`: VCF file with resulting tandem repeat genotypes. By default, Medaka Tandem reports entire haplotype-specific tandem repeats as alternative alleles. Use `--decompose` option to report smaller variants relative to the reference sequence.
2. `skipped.bed`: BED file containing the regions that skipped due to insufficient number of spanning reads.
3. `skipped_large.bed`: BED file containing skipped large regions (>10kbp). To analyse these regions, one can re-run Medaka with `--process_large_regions` option using this BED file to specify regions of interest.
4. `consensus.fasta`: FASTA file containing the polished consensus sequence covering the tandem repeat regions plus padding.
5. `medaka_to_ref.bam`: BAM file containing the mapping of `consensus.fasta` to the reference genome.
6. `{prephased,abpoa,unphased}_region_metrics.txt`: Text files providing the number of reads supporting each allele.
7. `trimmed_reads.fasta`: read segments spanning tandem repeat regions.
8. `poa.fasta`: Draft consensus allele sequences before polishing.
9. `trimmed_reads_to_poa.bam`: alignments of sequences from `trimmed_reads.fasta` to `poa.fasta`, used during polishing.

VCF Preparation for TDB
-----

A helper script is provided to preprocess tandem output VCF for TDB import:

```bash
./medaka/tandem/scripts/medaka_tandem_to_tdb medaka_to_ref.TR.vcf prepared_for_tdb.vcf.gz
```

It performs the preprocessing steps for TDB import and writes a `.vcf.gz` with `.tbi` index:
sort, add `AM` FORMAT field placeholder values, split multi-allelic records, and final sort.

Then create a TDB explicitly:

```bash
tdb create -o output.tdb prepared_for_tdb.vcf.gz
```

Help
----

**Licence and Copyright**

© 2018- Oxford Nanopore Technologies Ltd.

`medaka` is distributed under the terms of the Oxford Nanopore Technologies PLC. Public License Version 1.0

**Research Release**

Research releases are provided as technology demonstrators to provide early
access to features or stimulate Community development of tools. Support for
this software will be minimal and is only provided directly by the developers.
Feature requests, improvements, and discussions are welcome and can be
implemented by forking and pull requests. However much as we would
like to rectify every issue and piece of feedback users may have, the
developers may have limited resource for support of this software. Research
releases may be unstable and subject to rapid iteration by Oxford Nanopore
Technologies.
