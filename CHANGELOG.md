# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]
### Fixes
- Update documentation with `inference` and `sequence` command renaming.
### Added
- Models `dna_r10.4.1_e8.2_5khz_400bps_sup` and `dna_r10.4.1_e8.2_5khz_400bps_hac` added
  as aliases to those without `_5kz_` in their names. 
- Consensus models for v5.2.0 dorado models.
- Added `-B` option to `medaka_consensus` to allow passing a bed file or region to polish
  via `medaka inference --regions`.

## [v2.0.1]
### Fixed
- `medaka smolecule` was broken by change from `medaka consensus` to `medaka inference`.
### Changed
- Improved error message when model is not found.

## [v2.0.0]
Switched from tensorflow to pytorch.

Existing models for recent basecallers have been converted to the new format.
Pytorch format models contain a ``_pt`` suffix in the filename.

### Changed
- Inference is now performed using PyTorch instead of TensorFlow.
- The `medaka consensus` command has been renamed to `medaka inference` to reflect
  its function in running an arbitrary model and avoid confusion with `medaka_consensus`.
- The `medaka stitch` command has been renamed to `medaka sequence` to reflect its
  function in creating a consensus sequence.
- The `medaka variant` command has been renamed to `medaka vcf` to reflect its function
  in consolidating variants and avoid confusion with `medaka_variant`.
- Order of arguments to `medaka vcf` has been changed to be more consistent
  with `medaka sequence`.
- The helper script `medaka_haploid_variant` has been renamed `medaka_variant` to
  save typing.
- Make `--ignore_read_groups` option available to more medaka subcommands including `inference`.
### Removed
- The `medaka snp` command has been removed. This was long defunct as diploid SNP calling
  had been deprecated, and `medaka variant` is used to create VCFs for current models.
- Loading models in hdf format has been deprecated.
- Deleted minimap2 and racon wrappers in `medaka/wrapper.py`.
### Added
- Release conda packages for Linux (x86 and aarch64) and macOS (arm64).
- Option `--lr_schedule` allows using cosine learning rate schedule in training.
- Option `--max_valid_samples` to set number of samples in a training validation batch.
### Fixed
- Training models with DiploidLabelScheme uses categorical cross-entropy loss
  instead of binary cross-entropy.

## [v2.0.0a2]
### Changed
- Minor edits to README around model selection and package installation.
### Added
- Release conda packages for Linux (x86 and aarch64) and macOS (arm64).

## [v2.0.0a1]
Switched from tensorflow to pytorch.

Existing models for recent basecallers have been converted to the new format.
Pytorch format models contain a ``_pt`` suffix in the filename.
### Changed
- Inference is now performed using PyTorch instead of TensorFlow.
- The `medaka consensus` command has been renamed to `medaka inference` to reflect
  its function in running an arbitrary model and avoid confusion with `medaka_consensus`.
- The `medaka stitch` command has been renamed to `medaka sequence` to reflect its
  function in creating a consensus sequence.
- The `medaka variant` command has been renamed to `medaka vcf` to reflect its function
  in consolidating variants and avoid confusion with `medaka_variant`.
- Order of arguments to `medaka vcf` has been changed to be more consistent
  with `medaka sequence`.
- The helper script `medaka_haploid_variant` has been renamed `medaka_variant` to
  save typing.
### Removed
- The `medaka snp` command has been removed. This was long defunct as diploid SNP calling
  had been deprecated, and `medaka variant` is used to create VCFs for current models.
- Loading models in hdf format has been deprecated.
- Deleted minimap2 and racon wrappers in `medaka/wrapper.py`.
### Added
- Option `--lr_schedule` allows using cosine learning rate schedule in training.
- Option `--max_valid_samples` to set number of samples in a training validation batch.
### Fixed
- Training models with DiploidLabelScheme uses categorical cross-entropy loss
  instead of binary cross-entropy.

## [v1.12.1]
(Probably) final version of medaka using tensorflow. Future versions will use
pytorch instead.
### Fixed
- medaka_consensus: only keep bam tags if input file matches joint polishing pipeline.
- Pin numpy to <2.0.0.
### Added
- Consensus and variant models lookup for v3.5.1 Dorado models.

## [v1.12.0]
### Fixed
- tandem: Use haplotag 0 in unphased mode.
- tandem: Don't run consensus if regions set is empty.
### Added
- Models for version 5 basecaller models.
- Expose `sym_indels` option for training.
- Expose `--min_mapq` minimum  mapping quality alignment fitering option for medaka consensus. 
- tandem: Option `--ignore_read_groups` to ignore read groups present in input file.
- Wrapper script `medaka_consensus_joint` and convenience tools (`prepare_tagged_bam`,
  `get_model_dtypes`) to facilitate joint polishing with multiple datatypes.

## [v1.11.3]
### Added
- Consensus and variant models for v4.3.0 dorado models.

## [v1.11.2]
### Added
- Parsing model information from fastq headers output by Guppy and MinKNOW.
### Changed
- Additional explanatory information in VCF INFO fields concerning depth calculations.

## [v1.11.1]
### Fixed
- Do not exit if model cannot be interpreted, use the default instead.
- An issue with co-ordinate handling in computing variants from alignments.
### Added
- Ability to use basecaller model name as --model argument.
- Better handling or errors when running abpoa.

## [v1.11.0]
### Fixed
- Correct suffix of consensus file when `medaka_consensus` outputs a fastq.
### Added
- Choice of model file can be introspected from input files. For BAM files the
  read group (RG) headers are searched according to the dorado
  [specification](https://github.com/nanoporetech/dorado/blob/master/documentation/SAM.md),
  whilst for .fastq files the comment section of a number of reads are checked
  for corresponding read group information. In the latter case see README for
  information on correctly converting basecaller output to .fastq whilst
  maintaining the relevant meta information.
- `medaka tools resolve_model` can display the model that would automatically
  be used for a given input file.
### Changed
- If no model is provided on command-line interface (medaka consensus,
  medaka_consensus, and medaka_haploid_variant) automatic attempts will be made
  to choose the appropriate model.

## [v1.10.0]
### Changed
- Tensorflow logging level no longer set from Python.
- spoa and parasail are now strict requirements.
### Fixed
- Sort VCF before annotating in `medaka_haploid_variant`.
- Ignore errors when deleting temporary files.
- The output of the first POA run not being used in the second iteration in smolecule command.
### Added
- Support for Python 3.11.
- `--spoa_min_coverage` option to smolecule command.
### Removed
- Support for Python 3.7.

## [v1.9.1]
### Fixed
- A long-standing bug in pileup_counts that manifests for single-position pileups on ARM64.

## [v1.9.0]
### Added
- Added `medaka tandem` targeted tandem repeat variant calling.

## [v1.8.2]
### Added
* Updated features related to fetching of trimmed reads.
### Changed
* Refactored smolecule module.
* Faster inference and stitching of many short contigs.
* Tensorflow version 2.10 (allows for aarch64 wheels).

## [v1.8.1]
### Added
- Expose qualities parameter in medaka_consensus script with `-q` parameter.

## [v1.8.0]
### Added
- Consensus and variant models for v4.1 and v4.2 basecallers.
### Changed
- Changed default models to be r1041_e8.2_400bps_v4.2 models
- Clip probabilities in `_phred()` rather than adding smallest float.

## [v1.7.3]
### Added
- Consensus polishing models for Version 4 basecallers.
- Wheel builds for newer Python versions.
### Fixed
- Deprecated numpy.unicode use.
### Changed
- Set minimum Python version to 3.7.
- Updated tensorflow requirement to 2.8.
- Put lower bound on numpy requirement.
### Removed
- Dropped support for Python 3.6. Security support for Python 3.6 was ended on 23 Dec 2021;
  as such we have removed support for Python 3.6 and suggest users update their Python version.

## [v1.7.2]
### Added
- New models for R10.4.1 E8.2 260bps based sequencing chemistries.
### Changed
- Updated Hac and Fast models for R10.4.1 E8.2 400bps based sequencing chemistries.
- Removed models for Fast basecallers from pypi package

## [v1.7.1]
### Fixed
- medaka variant IndexError on long insertion

## [v1.7.0]
### Added
- capability to fill gaps in consensus sequence with a designated character 
  (e.g. 'N') instead of content from a reference sequence.
- option `-r` in `medaka_consensus` to set the designated fill character.
- option `--fill_char` in `medaka stitch` to set the designated fill character.
### Fixed
- CUDA initialization errors during `medaka smolecule`s stitch phase.

## [v1.6.1]
### Added
- New models for R10.4.1 E8.2 400bps based sequencing chemistries.

### Fixed
- DiploidZygosityLabelScheme renaming.

## [v1.6.0]
### Changed
- Updated to tensorflow~=2.7.0.
- Do not always force recreation of minimap2 index in helper scripts.
- PyPI wheel releases now built with libdeflate for faster BAM reading.
### Fixed
- Inclusion of inserted bases immediately after deletion in pileup counts.
### Added
- Makefile can now build environment for macOS M1.
- Publish ARMv8 wheels compatible with NVIDIA's [Jetpack 4.6.1 binary](https://developer.download.nvidia.com/compute/redist/jp/v461/tensorflow).
- `--qualities` option for `smolecule` and `stitch` to output consensus fastq.

## [v1.5.0]
### Changed
- Updated tensorflow requirement to ~=2.5.2.
- Spruced-up documentation.
### Added
- Light testing of Docker build.
### Removed
- Remove `medaka_variant` in deference of clair3.

## [v1.5.0.rc1]
### Changed
- Updated tensorflow requirement to ~=2.4.4.
### Added
- Light testing of Docker build.

## [v1.4.4]
### Changed
- tensorflow requirement to ~=2.2.2
### Added
- R10.4 E8.1 consensus models for Guppy version 5.0.15.

## [v1.4.3]
### Fixed
- `medaka tools` now displays its help rather than an error.
### Added
- `medaka tools download_models can download specific models.

## [v1.4.2]
### Fixed
- Missing sites in gVCF output.
### Changed
- Rewrittern algorithm for determining VCF records from RNN outputs for clarity
  and speed.

## [v1.4.1]
### Fixed
- Inclusion of select models in distributions.

## [v1.4.0]
### Added
- Models for Guppy version 5.0.7.

## [v1.3.4]
### Fixed
- Issue whereby tensorflow would spawn many threads that do not exit.

## [v1.3.3]
### Fixed
- Added missing default option to arparse instance in smolecule command.
- Copy across contigs with no aligned reads during `medaka stitch`.
- Quote strings in bash scripts to allow filenames with spaces.

## [v1.3.2]
### Fixed
- Typo in `medaka_consensus` causing a syntax error.

## [v1.3.1]
### Added
 - Option to output VCF record for all reference positions from
   `medaka variant`.

## [v1.3.0]
### Changed
 - Haploid variant calling reverted to old-style methodology.
### Fixed
 - Early exits on error in `medaka_consensus` and `medaka_variant`.
 - `INFO` field of VCFs is now correctly `.` when empty.
 

## [v1.2.6]
### Changed
 - Rewrote inference data loading code for clarity.
 - Removed pinned BioPython pin.
 - Formally update htslib program requirements to 1.11.
### Removed
 - Support for Python 3.5.
### Fixed
 - Corner case in consensus stitching.


## [v1.2.5]
### Fixed
 - Variant annotation when more than one CHROM record in VCF.


## [v1.2.4]
### Fixed
 - Variant annotation when counts matrix does not span variants.
### Changed
 - Updated Tensorflow requirement to 2.2.2


## [v1.2.3]
### Fixed
 - Off-by-one error during stitching of consensus chunks.


## [v1.2.2]

Minor release

### Fixed
 - Fixed incorrect read depth annotations in VCFs.
 - Fixed missing files in PyPI source distribution.
 - Fix `StopIteration` issues in newer Pythons.
### Added
 - Added `-n` option to `medaka_variant` to add a sample field to outputs.
 - Set `HDF5_USE_FILE_LOCKING=FALSE`, which some users report as useful.
 - Set `OMP_NUM_THREADS=1` required to make Tensorflow anaconda use CPU resource sensibly.


## [v1.2.1]

Minor release

### Fixed
 - Fix issue whereby variant ALTs were created equal to REF.
### Added
 - Build a medaka-cpu package depending on tensorflow-cpu.


## [v1.2.0]

Performance release

### Fixed
 - Fix long-standing issue where genome regions could be unprocessed.
### Added
 - Improve inference performance by 30%.
 - Add efficient multiprocessing to `medaka stitch`.


## [v1.1.3]

Bug fix and feature release

### Fixed
 - Fix iteration error in retrieving trimmed reads.
 - Work around tensorflow threading issue.
### Added
 - Add ability to `haploid2diploid` tool on VCFs generated by `medaka_haploid_variant`


## [v1.1.2]

Bug fix and feature release

### Fixed
### Added
 - Fix issues in command-line argument parsing.
 - Add true ploidy-1 variant caller.
 - Do not break contigs at unpolished regions (fill with input instead).
 - Add multi-nucleotide variant decomposition to be compatible with DeepVariant.


## [v1.1.1]

Bux fix release.

### Fixed
 - Remove python version check preventing Python >3.6 builds from running.


## [v1.1.0]

Update with new models and features.

### Fixed
 - Fix a few bugs in variant annotation program.
### Added
 - Add ARM builds to PyPI release.
 - Add Python 3.7 and 3.8 builds for x86-64.
 - Add PromethION model for Guppy 4.0.11.
### Changed
 - Upgrade to Tensorflow 2.2.
 - Option to split MNPs to independent SNPs (for compatibility with DeepVariant).
 - Single molecule consensus program now uses `pyspoa`.
### Removed
 - Remove methylation aggregation functionality.


## [v1.0.3]

Minor fixes release.

### Fixed
 - Fix occasional mangled sam output in guppy2sam.
### Changed
 - Update htslib ecosystem to 1.10 to fix conda installation issue.

## [v1.0.2]

Minor fixes and models release.

### Fixed
 - VCF GQ is now an integer in line with VCF spec.
 - Fixed issue requiring a previous model for training.
 - Fixed issue causing -p option of medaka_variant to crash.
 - Fixed issue preventing installation in a virtualenv with python <3.6.
### Added
 - R9.4.1 variant calling models for Guppy 3.6.0 and updated benchmarks.
 - Made r941_min_high_g360 the default consensus model. 


## [v1.0.1]

Minor fixes release, resolving issues introduced in v1.0.0.

### Fixed
 - Fix default model for SNP calling. 
 - Fix issue causing medaka_consensus to crash. 


## [v1.0.0]

Models, features and fixes release

### Fixed
 - Fix to methylation aggregation.
### Added
 - Consensus models for Guppy 3.6.0.
 - Add functionality for auto-download of older models.
 - VCF annotation tool.


## [v0.12.1]

Minor release.

### Changed
 - Harmonised versions of htslib/samtools dependencies.


## [v0.12.0]

Models, features and fixes release

### Fixed
 - Minor speed improvement.
 - Fix bug where force overwrite of output was always enabled.
 - Fix bug where variant calling of a region crashed if the region began with a deletion.
### Added
 - Variant calling models for R10.3 and R9.4.1 and updated benchmarks.
 - Consensus models for Guppy 3.5.1.
 - Add read group (RG) tag filtering.
 - Add option to create consensus sequence via intermediate .vcf file.
 - Update to methylation calling documentation.
 - Addition of all-context modified-base aggregation.


## [v0.11.5]

R10.3 model and small fixes

### Fixed
 - Fix index/compression issue with RLE workflow
 - Fix a rare memory error during feature generation caused by very long indels.
### Added
 - Add model for R10.3 on MinION.
 - Write and empty vcf when no variants are found in medaka_variant.


## [v0.11.4]

Bugfix release

### Fixed
 - Fix invalid specification of variant calling model.


## [v0.11.3]

Model release

### Added
 - Models for guppy 3.4.4. 


## [v0.11.2]

Minor fix release

### Fixed
 - Fix a memory error in pileup calculation.
### Added
 - Update variant calling models and benchmarks.


## [v0.11.1]

Minor fix release

### Fixed
 - Detect NaNs during training and halt early.
 - Workaround pysam interface changes (for conda package).
### Added
 - Preliminary hard-RLE model for R9.4.1
 - --regions argument can now be a .bed file.
 - Support soft-RLE network training.

This release includes an experimental consensus mode using run-length encoded
alignments. Use of this algorithm can be specified using the new "rle" model:

    medaka_consensus -m r941_min_high_g340_rle -i basecalls.fasta -d draft.fa

## [v0.11.0]

Feature release

### Added
 - Consensus models for guppy 3.3 and 3.4.
 - Aggregation of Guppy modified base probability tables.
 - Multi-thread stitching of inference chunks in `medaka_consensus`.
 - Optionally run whatshap phase at the end of `medaka_variant`.


## [v0.10.1]

Minor fix release

### Fixed
 - Fix bug where feature matrix was misaligned with coordinate system.
 - Fixed issue with `medaka_variant` failing on zero-coverage regions.
 - Rename incorrectly named diploid SNP calling model.
 - Made variant calling faster by resolving trivial bottleneck in variant classification.
### Added
 - Add missing arguments from `smolecule` command.
 - Output contig names are no longer written as samtools-style regions.


## [v0.10.0]

Feature release

### Fixed
 - Corrected parsing of region strings with multiple `:` charaters
 - Fixed bug causing larger than requested overlap in inference chunks.
 - Fixed rare consensus stitching error.
### Added
 - Added a `-f` force overwrite option to `medaka_consenses`.
 - Added *C. elegans* assembly benchmarks to documentation.
### Changed
 - Switched variant calling to an explicitely diploid calling model.
 - Refreshed *E. coli* benchmark to include effect of `racon`.
 - Refreshed variant calling benchmarks.


## [v0.9.2]

Minor fix release.

### Fixed
 - Additional fix to handling lowercase reference sequences.
 - Fix bug in creation of RLE alignments.
### Added
 - Update `update_model.py` script.
### Changed
 - Unify how LabelSchemes store training data.
### Removed
 - Remove option to select labelling scheme during training.


## [v0.9.1]

Minor fix release

### Fixed
 - Fix regression in medaka stitch and medaka snp speed.
 - Remove dill and yaml requirements.
### Added
 - Handle lowercase letters in reference sequences.


## [v0.9.0]

Bugfix and training refactor release

### Fixed
 - Fix readlink issue on MacOS
 - Fix bug where medaka_variant did not call indels by default
 - Fix bug in determining when to split contigs
 - Make network feature generation 2x faster
### Added
 - Add smolecule command
 - Log use of GPU and cuDNN, noting workaround for RTX cards
### Changed
 - Store models in git-lfs
 - Simplify medaka_variant workflow for speed
 - Refactor labelling of training data and storing of models
 - Reimplement RLE feature generation
### Removed
 - Drop support for older basecaller models (guppy<3.0.3)


## [v0.8.2]

Documentation release

### Added
 - Clarify suggested workflows in documentation.


## [v0.8.1]

Patch release

### Fixed
 - Patch import of loading of older models


## [v0.8.0]

Model release and development release

### Added
 - Add support for R10 basecaller
 - Add diploid multi-labelling
### Changed
 - Upgrade to tensorflow 1.14.0


## [v0.7.1]

Bug fix release

### Fixed
 - Fix regression in consensus stitching when chunks do not overlap.


## [v0.7.0]

Feature release

### Added
 - Indel calling for `medaka_variant`.
 - New models for MinION/GridION and PromethION paired to high accuracy an fast
  guppy basecallers.
### Changed
 - Overhaul of chunk handling and overlapping.


## [v0.6.5]

Bug fix release

### Fixed
 - Fix Makefile for parallel build.
 - Ensure medaka consensus is given absolute path to model.
### Added
 - Tidy up some parsing and sorting of regions from strings.
### Changed
 - Disable by default validation of output HDF during consensus.
 - Refactor variant handling code.


## [v0.6.4]

Bug fix release

### Fixed
 - Fix for models not specifiying data types.


## [v0.6.3]

Bug fix release

### Fixed
 - Split pileup when reads do not span rather than silently deleting region.
 - Fix error in stitching occuring with a single region.
### Changed
 - Refactor handling of short and remainder regions.
### Removed
 - Drop 3.4 support.


## [v0.6.2]

Bug fix release

### Added
 - Enhanced verification of training feature samples
### Changed
 - Pin pyyaml version

## [v0.6.1]

Bug fix release

### Fixed
 - Fixed bug in `medaka_consensus` incorrectly calling python


## [v0.6.0]

SNP calling, model, and bugfix release release

### Fixed
 - Workaround short-contig/no-coverage corner case during pileup.
### Added
 - Prototype SNP calling and phasing, [benchmarks](https://nanoporetech.github.io/medaka/snp.html)
 - Add model for improved Flip-flop model in Guppy 2.3.5
### Changed
 - Rename models to be more logical
 - Update to htslib version 1.9 for long cigars


## [v0.5.2]

Bug fix release

### Fixed
 - Fix bug leading to dropping of pileup chunks during loading


## [v0.5.1]

Development and performance release

### Changed
 - Refactor batch queuing in preparation to using keras Sequence
 - Asynchrounous feature loading during inference
 - Pin version of h5py to work around intermittent errors in saving models


## [v0.5.0]

Training and bug-fix release.

### Fixed
 - Resolve issue with contained chunks during stitching
 - Resolve hanging at the end of training
### Added
 - Improved storage and retrieval of features for better IO
 - Training speed improved >10X
### Changed
 - Switch to CuDNN for GRU layers
 - Check presence of `minimap2` and `samtools`
 - Provide more feedback on error


## [v0.4.3]

Model release.

### Added
 - Add support for R9.4.1 flip-flop basecaller


## [v0.4.1]

Development release.

### Added
 - Adds build infrastructure for source distributions and manylinux wheels.


## [v0.4.0]

Performance and bugfix release.

### Added
 - Large refactoring of feature and sample generation. Fixes many small bugs
  and edge cases
 - Resize models for small contigs
 - Faster Generation of inference features
 - Model updates
 - Ability to handle multiple read types
### Changed
 - Remove redundant samtools tview code
 - Limit CPU usage when running without a GPU


## [v0.3.0]

Model and userbility release.

### Fixed
 - Many small bug fixes
### Added
 - New non-RLE model
 - Updated documentation and benchmarks
 - Dockerfile to build a medaka Docker image
### Changed
 - `medaka_consensus` no longer needs a pomoxis installation to run
