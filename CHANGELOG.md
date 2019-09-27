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
