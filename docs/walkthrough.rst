.. _Walkthrough:

Walkthrough
===========

In the following we will show how to train and use a simple form of consensus
model. We will start from sequencing data obtained from R9.4 flowcells and the
ligation sequencing kit.

.. note:: The below assumes that you are using a Linux environment, some
    changes may need to be made on macOS. Windows environments are not
    supported. Environment variables set in one code block are assumed to
    persist through to other code blocks. 

The below serves to demonstrate the process at a simple level. It does not
represent a best-practices or state-of-the-art workflow. To train models
which generalise well to other datasets more careful preparation of a larger
dataset is required. `Medaka` is designed for flexibility over performance.


Obtaining data
--------------

We start by downloading training data in the form of Oxford Nanopore
Technologies `.fast5` files for a set of reads.

.. code-block:: bash

    WALKTHROUGH=${PWD}/medaka_walkthrough
    mkdir ${WALKTHROUGH} && cd ${WALKTHROUGH}
    wget https://s3-eu-west-1.amazonaws.com/ont-medaka-demo/medaka_demo.tar
    tar -xvf medaka_demo.tar
    DATA=${PWD}/data


Creating a software environment
-------------------------------

We will use several open source tools from Oxford Nanopore Technologies. To
obtain these run the following:

.. code-block:: bash

    cd ${WALKTHROUGH}
    git clone https://github.com/nanoporetech/pomoxis --recursive
    git clone https://github.com/nanoporetech/medaka
    git clone https://github.com/nanoporetech/scrappie
    
    # While it is possible to install pomoxis and medaka into the same
    #   virtual environment, we will install each package into its own
    #   environment for simplicity. For more details see the readme for
    #   each of the packages.

    cd pomoxis && make install && cd ..
    cd medaka && make install && cd ..
    cd scrappie && mkdir build && cd build && cmake .. && make && cd ../../

If any of the above steps fail consult the install instructions of the
individual packages. 
 
In order to train effectively consensus models it is desirable to use a system
with a GPU. In order to use `medaka` training with a GPU, one should further
install the `tensorflow-gpu` package from pip into the medaka virtual
environment with:

.. code-block:: bash

    cd ${WALKTHROUGH}
    . medaka/venv/bin/activate
    pip install tensorflow-gpu

Depending on your environment, specifically the versions of `CUDA` and `cuDNN`
that one has installed, it may be necessary to use versions of this package other
than the latest. `medaka` has been used with all release versions of `tensorflow`
after and including version 1.0.0.


Basecalling and draft assembly
------------------------------

Our first step is to basecall the data downloaded above. We do this using `scrappie`:

.. code-block:: bash

    cd ${WALKTHROUGH}
    SCRAPPIE=${WALKTHROUGH}/scrappie/build/scrappie
    NPROC=$(nproc)
    export OMP_NUM_THREADS=${NPROC}
    export OPENBLAS_NUM_THREADS=1
    ${SCRAPPIE} raw ${DATA}/reads --model rgrgr_r94 > basecalls.fa

We can now form a draft assembly using the
`miniasm <https://github.com/lh3/miniasm>`_ based pipeline from `pomoxis`.
Alternatively one could use `canu <https://github.com/marbl/canu>`_ at this step.

.. code-block:: bash

    cd ${WALKTHROUGH}
    . pomoxis/venv/bin/activate
    cat basecalls.fa | fast_convert --mock_q 10 aq > basecalls.fq 
    mini_assemble -i basecalls.fq -o draft_assm -p assm -t ${NPROC}

This will create a draft assembly at `draft_assm/assm_final.fa`. The
`mini_assemble` script has two useful options not used here:

    * specifying `-c` will run `porechop <https://github.com/rrwick/Porechop>`_
      on the reads to first trim sequencing adapters and,
    * specifying `-e 10` will perform error correction on the longest 10% of
      reads prior to assembly (similar to the strategy of canu).

Both these steps can improve the assembly quality at the expense of speed.


Preparing training data
-----------------------

In order to correct a draft assembly, medaka currently uses a strategy of
aligning reads to a draft and learning corrections to be made taking the pileup
data as input. In order to train the network we therefore need to perform
such an alignment and also an alignment of the ground truth data to the same
draft. In this way we can learn the correct base(s), or gaps, to call for each
pileup column.

To align the basecalls and truth sequence to the draft assembly we use the
`mini_align` script in `pomoxis` which conveniently wraps 
`minimap2 <https://github.com/lh3/minimap2>`_:

.. code-block:: bash

    cd ${WALKTHROUGH}
    . pomoxis/venv/bin/activate
    DRAFT=draft_assm/assm_final.fa
    CALLS2DRAFT=calls_to_draft
    TRUTH2DRAFT=truth_to_draft
    mini_align -P -r ${DRAFT} -i basecalls.fa -a -t ${NPROC} -p ${CALLS2DRAFT}
    mini_align -P -r ${DRAFT} -i ${DATA}/truth.fa -a -t ${NPROC} -p ${TRUTH2DRAFT}

At the end of this process we have two `.bam` files which we use in the
following training step effectively to learn how to predict the contents of
one from the other.

We now use `medaka prepare` to generate training examples in the form of 
chunks of labelled pileup, and save them to HDF5.

.. code-block:: bash

    cd ${WALKTHROUGH}
    . medaka/venv/bin/activate
    LABELLEDPILEUP=labelled_pileup.hdf
    TRAINREFNAME=Consensus_utg000002l
    EVALREFNAME=Consensus_utg000001l
    medaka prepare ${DRAFT} ${CALLS2DRAFT}.bam ${LABELLEDPILEUP} --truth ${TRUTH2DRAFT}.bam --ref_name $TRAINREFNAME


Training The Consensus Network
------------------------------

We now have everything we need to train a consensus network using `medaka train`:

.. code-block:: bash

    cd ${WALKTHROUGH}
    . medaka/venv/bin/activate
    TRAINNAME=training
    medaka train ${LABELLEDPILEUP} --train_name ${TRAINNAME}

During training, models are regularly checkpointed so that training may be
easily resumed if interrupted. At the end of training, we have a number of
output models including in particular:

    * `model.best.hdf5`: model with the best accuracy over the training set  
    * `model.best.val.hdf5`: model with the best accuracy over the validation set

which can be used to calculate a consensus. Other ancilliary output are
also produced.


Polishing a Consensus 
---------------------

Having trained a model we can run `medaka consensus` to calculate consensus
using our trained model:

.. code-block:: bash

    cd ${WALKTHROUGH}
    . medaka/venv/bin/activate
    medaka consensus ${TRAINNAME}/model.best.val.hdf5 --alignments ${CALLS2DRAFT}.bam ${DRAFT} ${EVALREFNAME} 

The program outputs a single .fasta file containing consensus sequences for
overlapping chunks of the input. This stage may be parallelised trivially
by running the program on distinct sections of a draft assembly.

It is possible also to save probabilities obtained from the neural network
using the `--output_probs` option. This may be of use in assessing consensus
quality and in variant calling; future releases of `medaka` may implement the
writing of VCF files.

Next we can recombine the consensus chunks using `medaka stitch`.

.. code-block:: bash

    cd ${WALKTHROUGH}
    . medaka/venv/bin/activate
    medaka stitch consensus_chunks.fa medaka_consensus.fasta

The output file `medaka_consensus.fasta` now contains our neural network consensus.

Finally, we can run `stats_from_bam` to assess to what extent `medaka` has improved accuracy. 

.. code-block:: bash

    cd ${WALKTHROUGH}
    . pomoxis/venv/bin/activate
    EVALREGION=$(awk -F '[>:-]' '{if(NR==1){printf("%s:%i-%i\n",$2, $3, $4)}}' medaka_consensus.fasta)
    samtools faidx ${DRAFT} ${EVALREGION} > ${EVALREGION}_draft_assm.fa
    mini_align -P -r ${DATA}/truth.fa -i ${EVALREGION}_draft_assm.fa -a -t ${NPROC} -p draft_to_truth
    echo "Draft assembly"
    stats_from_bam --bam draft_to_truth.bam > draft_to_truth.stats.txt
    mini_align -P -r ${DATA}/truth.fa -i medaka_consensus.fasta -a -t ${NPROC} -p consensus_to_truth
    echo "Medaka consensus"
    stats_from_bam --bam consensus_to_truth.bam > consensus_to_truth.stats.txt

An increase in accuracy should be observed. 
