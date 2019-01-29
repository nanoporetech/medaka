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
 
* Large refactoring of feature and sample generation #10. Fixes many small bugs and edge cases
* Resize models for small contigs #9
* Faster Generation of inference features
* Model updates
* Remove redundant samtools tview code
* Ability to handle multiple read types
* Limit CPU usage when running without a GPU


v0.5.0
------
Model and userbility release.

* New non-RLE model
* Updated documentation and benchmarks
* Many small bug fixes
* `medaka_consensus` no longer needs a pomoxis installation to run
* Dockerfile to build a medaka Docker image