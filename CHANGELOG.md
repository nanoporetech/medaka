v1.0.0
-------
Models, features and fixes release

* Consensus models for Guppy 3.6.0.
* Fast consensus models for Guppy 3.4.4.
* Fix to methylation aggregation.
* VCF annotation tool.


v0.12.1
-------
Minor release harmonising versions of htslib/samtools dependencies.

v0.12.0
-------
Models, features and fixes release

* Variant calling models for R10.3 and R9.4.1 and updated benchmarks.
* Consensus models for Guppy 3.5.1.
* Add read group (RG) tag filtering.
* Add option to create consensus sequence via intermediate .vcf file.
* Update to methylation calling documentation.
* Addition of all-context modified-base aggregation.
* Minor speed improvement.
* Fix bug where force overwrite of output was always enabled.
* Fix bug where variant calling of a region crashed if the region began with a deletion.

v0.11.5
-------
R10.3 model and small fixes

* Add model for R10.3 on MinION.
* Fix index/compression issue with RLE workflow
* Write and empty vcf when no variants are found in medaka_variant.
* Fix a rare memory error during feature generation caused by very long indels.

v0.11.4
-------
Bugfix

* Fix invalid specification of variant calling model.

v0.11.3
-------
Model release

* Models for guppy 3.4.4. 

v0.11.2
-------
Minor fix release

* Fix a memory error in pileup calculation.
* Update variant calling models and benchmarks.

v0.11.1
-------
Minor fix release

* Preliminary hard-RLE model for R9.4.1
* --regions argument can now be a .bed file.
* Detect NaNs during training and halt early.
* Workaround pysam interface changes (for conda package).
* Support soft-RLE network training.

This release includes an experimental consensus mode using run-length encoded
alignments. Use of this algorithm can be specified using the new "rle" model:

    medaka_consensus -m r941_min_high_g340_rle -i basecalls.fasta -d draft.fa


v0.11.0
-------
Feature release

* Consensus models for guppy 3.3 and 3.4.
* Aggregation of Guppy modified base probability tables.
* Multi-thread stitching of inference chunks in `medaka_consensus`.
* Optionally run whatshap phase at the end of `medaka_variant`.

v0.10.1
-------
Minor fix release

* Fix bug where feature matrix was misaligned with coordinate system.
* Add missing arguments from `smolecule` command.
* Output contig names are no longer written as samtools-style regions.
* Fixed issue with `medaka_variant` failing on zero-coverage regions.
* Rename incorrectly named diploid SNP calling model.
* Made variant calling faster by resolving trivial bottleneck in variant classification.

v0.10.0
-------
Feature release

* Switched variant calling to an explicitely diploid calling model.
* Added a `-f` force overwrite option to `medaka_consenses`.
* Refreshed *E. coli* benchmark to include effect of `racon`.
* Refreshed variant calling benchmarks.
* Added *C. elegans* assembly benchmarks to documentation.
* Fixed bug causing larger than requested overlap in inference chunks.
* Corrected parsing of region strings with multiple `:` charaters
* Fixed rare consensus stitching error.

v0.9.2
------
Minor fix release.

* Additional fix to handling lowercase reference sequences.
* Fix bug in creation of RLE alignments.
* Update `update_model.py` script.
* Remove option to select labelling scheme during training.
* Unify how LabelSchemes store training data.

v0.9.1
------
Minor fix release

* Fix regression in medaka stitch and medaka snp speed.
* Handle lowercase letters in reference sequences.
* Remove dill and yaml requirements.

v0.9.0
------
Bugfix and training refactor release

* Fix readlink issue on MacOS
* Fix bug where medaka_variant did not call indels by default
* Fix bug in determining when to split contigs
* Drop support for older basecaller models (guppy<3.0.3)
* Store models in git-lfs
* Simplify medaka_variant workflow for speed
* Make network feature generation 2x faster
* Add smolecule command
* Log use of GPU and cuDNN, noting workaround for RTX cards
* Refactor labelling of training data and storing of models
* Reimplement RLE feature generation


v0.8.2
------

Documentation release

* Clarify suggested workflows in documentation.


v0.8.1
------
Patch release

* Patch import of loading of older models


v0.8.0
------
Model release and development release

* Add support for R10 basecaller
* Add diploid multi-labelling
* Upgrade to tensorflow 1.14.0


v0.7.1
------
Bug fix release

* Fix regression in consensus stitching when chunks do not overlap


v0.7.0
------
Feature release

* Indel calling for `medaka_variant`.
* New models for MinION/GridION and PromethION paired to high accuracy an fast
  guppy basecallers.
* Overhaul of chunk handling and overlapping.


v0.6.5
------
Bug fix release

* Tidy up some parsing and sorting of regions from strings.
* Disable by default validation of output HDF during consensus.
* Refactor variant handling code.
* Ensure medaka consensus is given absolute path to model.
* Fix Makefile for parallel build.


v0.6.4
------
Bug fix release

* Fix for models not specifiying data types.


v0.6.3
------
Bug fix release

* Split pileup when reads do not span rather than silently deleting region.
* Refactor handling of short and remainder regions.
* Drop 3.4 support.
* Fix error in stitching occuring with a single region.


v0.6.2
------
Bug fix release

* Enhanced verification of training feature samples
* Pin pyyaml version

v0.6.1
------
Bug fix release

* Fixed bug in `medaka_consensus` incorrectly calling python


v0.6.0
------
SNP calling, model, and bugfix release release

* Prototype SNP calling and phasing, [benchmarks](https://nanoporetech.github.io/medaka/snp.html)
* Add model for improved Flip-flop model in Guppy 2.3.5
* Rename models to be more logical
* Update to htslib version 1.9 for long cigars
* Workaround short-contig/no-coverage corner case during pileup.

v0.5.2
------
Bug fix release

* Fix bug leading to dropping of pileup chunks during loading


v0.5.1
------
Development and performance release

* Refactor batch queuing in preparation to using keras Sequence
* Asynchrounous feature loading during inference
* Pin version of h5py to work around intermittent errors in saving models


v0.5.0
------
Training and bug-fix release.

* Large refactor of training code
    * Resolve hanging at the end of training #15
    * Switch to CuDNN for GRU layers
    * Improved storage and retrieval of features for better IO
    * Training speed improved >10X
* Resolve issue with contained chunks during stitching #20
* Changes to `medaka_consensus`
    * checks presence of `minimap2` and `samtools`
    * provides more feedback on error #16
    
    
v0.4.3
------
Model release.

* Add support for R9.4.1 flip-flop basecaller


v0.4.1
------
Development release.

* Adds build infrastructure for source distributions and manylinux wheels.
 
 
v0.4.0
------
Performance and bugfix release.
 
* Large refactoring of feature and sample generation #10. Fixes many small bugs
  and edge cases
* Resize models for small contigs #9
* Faster Generation of inference features
* Model updates
* Remove redundant samtools tview code
* Ability to handle multiple read types
* Limit CPU usage when running without a GPU


v0.3.0
------
Model and userbility release.

* New non-RLE model
* Updated documentation and benchmarks
* Many small bug fixes
* `medaka_consensus` no longer needs a pomoxis installation to run
* Dockerfile to build a medaka Docker image
