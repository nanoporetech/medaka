
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
