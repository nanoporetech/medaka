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

In addition to usual Medaka dependencies, Medaka Tandem also requires [pyabpoa](https://github.com/yangao07/abPOA) v1.5.1 to be available.
On Linux, it can be installed as follows:
- if medaka was installed from conda, `conda install bioconda::pyabpoa==1.5.1`.
- if medaka was installed via pip, `pip install pyabpoa==1.5.1` but development libraries are required to build the wheel.
- if medaka was built from source with `make install`, pyabpoa was built as part of this installation and should already be present.


Usage
-----

`medaka tandem` command requires five positional arguments:
1. Reference genome in FASTA format.
2. Read alignments (preferably haplotagged) in BAM format.
3. Genomic regions to be analysed in BED format
4. Sample sex (‘male|female`)
5. Output directory
Addotto repeat catalogues can be downloaded from [here](https://github.com/ACEnglish/adotto/blob/main/regions/DataDescription.md).

Command line example

    medaka tandem \
            happlotagged.bam \
            ref.fna \
            adotto_TRregions_v1.2.bed \
            male \
            output_folder


## Additional options

### Performance-Related Options

1. `--workers` Number of parallel worker processes to use (default: 1).
2. `--process_large_regions` Process TRs with estimated length (of one or both alleles) exceeding 10kbp (default: False). Processing large regions can substantially increase RAM usage. With the default setting, the expected peak RAM consumption on Addotto repeat catalogue is approximately 14GB when using 8 workers, and 23GB when using 16 workers. Skipped regions will be output to `skipped_large.bed`.

### Behaviour-Related Options
1. `--phasing` Specify the strategy for dividing the reads between haplotypes. (default: `hybrid`)
    - `prephased` Rely on haplotype (HP) BAM tags for phasing.
    - `abpoa` Use abPOA clustering feature to identify haplotypes based on STR sequences in the reads.
    - `hybrid` Use haplotag assignments if both haplotypes have at least `min_depth` spanning reads assigned, otherwise fallback to the abPOA clustering.
    - `unphased` Assume the sample is haploid.
2. `--min_depth` Minimum number of spanning reads required for allele consensus reconstruction. (default: 3)
3. `--min_mapq` Minimum mapping quality (MAPQ) for alignments filtering. (default: 5)
4. `-- disable_outlier_filter` Disable exclusion of reads with significantly divergent spanning region lengths. (default: False).
5.  `--padding` Number of bases to pad spanning read regions and reference sequence. (default: 10)
6. `--sex_chrs` Comma separated names of X and Y chromosomes in reference FASTA. (default: ['chrX', 'chrY'])
7. `--par_regions` Coordinates of pseudoautosomal regions (PARs) on the X chromosome. Will be treated as diploid in male samples. The analysis assumes that the corresponding PARs on chromosome Y have been hard-masked (i.e. replaced with Ns) in the reference to avoid ambiguous read alignments. (default: chrX:10000-2781479,chrX:155701382-156030895 assuming use of the GRCh38 analysis set, e.g. `GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta`)
8. `--model` Model to be used for polishing. Can be a medaka model name or a basecaller model name suffixed with ':consensus'`
Please review the [Medaka Tandem README](../../README.md) for more information about the model choice.


### Output-Related Options

1. `--decompose` Align polished sequences back to reference region and extract a list of (left-aligned) variants. By default, Medaka Tandem reports entire haplotype-specific tandem repeats as alternative alleles.
2. `--add_read_names` Report names of spanning reads in the output VCF files.

3. `--sample_name` Sample name used in the output VCF file. (default: SAMPLE)

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



Licence and Copyright
======

© 2018- Oxford Nanopore Technologies Ltd.

`medaka` is distributed under the terms of the Mozilla Public License 2.0.

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

