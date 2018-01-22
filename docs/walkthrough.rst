Getting Started
===============

After installing all the necessary software (see :ref:`CreatingSoftwareEnv`), 
it is possible to use the pre-trained compressed-homopolymer model released with `medaka`
to perform a consensus call using `medaka_consensus`. An assembly in .fasta format and basecalls
in .fasta/.fastq format is required, see :ref:`BasecallingAndDraftAssembly`.

.. code-block:: bash

    source ${MEDAKA}
    NPROC=$(nproc)
    medaka_consensus -i basecalls.fa -d draft_assembly.fa -o medaka_consensus -t ${NPROC} -p ${POMOXIS}

When `medaka_consensus` has finished running, the consensus will be saved to medaka_consensus/consensus.fasta.


Walkthrough
===========

The following shows how to train and use a simple form of consensus
model from sequencing data obtained from R9.4 flowcells and the
ligation sequencing kit.

.. note:: The below assumes a Linux environment, some
    changes may need to be made on macOS. Windows environments are not
    supported. Environment variables set in one code block are assumed to
    persist through to other code blocks. 

The below serves to demonstrate the process at a simple level. It does not
represent a best-practices or state-of-the-art workflow. To train models
which generalise well to other datasets more careful preparation of a larger
dataset is required. `Medaka` is designed for flexibility over performance.


Obtaining Data
--------------

Start by downloading training data in the form of Oxford Nanopore
Technologies `.fast5` files for a set of reads.

.. code-block:: bash

    WALKTHROUGH=${PWD}/medaka_walkthrough
    mkdir -p ${WALKTHROUGH} && cd ${WALKTHROUGH}
    wget https://s3-eu-west-1.amazonaws.com/ont-medaka-demo/medaka_demo.tar
    tar -xvf medaka_demo.tar
    DATA=${PWD}/data

.. _CreatingSoftwareEnv:

Creating a Software Environment
-------------------------------

This walkthrough uses several open-source tools from Oxford Nanopore Technologies. To
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

    POMOXIS=${WALKTHROUGH}/pomoxis/venv/bin/activate
    MEDAKA=${WALKTHROUGH}/medaka/venv/bin/activate

If any of the above steps fail consult the install instructions of the
individual packages. 
 
In order to train effectively consensus models it is desirable to use a system
with a GPU. In order to use `medaka` training with a GPU, one should further
install the `tensorflow-gpu` package from pip into the medaka virtual
environment with:

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${MEDAKA}
    pip install tensorflow-gpu

Depending on your environment, specifically the versions of `CUDA` and `cuDNN`
that one has installed, it may be necessary to use versions of this package other
than the latest. `medaka` has been used with all release versions of `tensorflow`
after and including version 1.0.0.

.. _BasecallingAndDraftAssembly:

Basecalling and Draft Assembly
------------------------------

First basecall the data downloaded above using `scrappie`:

.. code-block:: bash

    cd ${WALKTHROUGH}
    SCRAPPIE=${WALKTHROUGH}/scrappie/build/scrappie
    NPROC=$(nproc)
    export OMP_NUM_THREADS=${NPROC}
    export OPENBLAS_NUM_THREADS=1
    BASECALLS=basecalls.fa
    ${SCRAPPIE} raw ${DATA}/reads --model rgrgr_r94 > ${BASECALLS}

Now form a draft assembly using the
`miniasm <https://github.com/lh3/miniasm>`_ based pipeline from `pomoxis`.
Alternatively one could use `canu <https://github.com/marbl/canu>`_ at this step.

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${POMOXIS}
    mini_assemble -i ${BASECALLS} -o draft_assm -p assm -t ${NPROC}

This will create a draft assembly at `draft_assm/assm_final.fa`. The
`mini_assemble` script has two useful options not used here:

    * specifying `-c` will run `porechop <https://github.com/rrwick/Porechop>`_
      on the reads to first trim sequencing adapters and,
    * specifying `-e 10` will perform error correction on the longest 10% of
      reads prior to assembly (similar to the strategy of canu).

Both these steps can improve the assembly quality at the expense of speed.

Now check the number and length of the assembled contigs.

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${POMOXIS}
    DRAFT=draft_assm/assm_final.fa
    awk '{if(/>/){n=$1}else{print n " " length($0)}}' ${DRAFT}

The expected output is a contig 4701891 bases long (Consensus_utg000001c) and a short contig just 408 bases long (Consensus_utg000002c). 
If this is not the case, the assembly step can be skipped by using the assembly in the data directory:

.. code-block:: bash

    cd ${WALKTHROUGH}
    rm -f draft_assm/* 
    cp ${DATA}/draft_assm.fa ${DRAFT} 
    awk '{if(/>/){n=$1}else{print n " " length($0)}}' ${DRAFT}

We will work with the long contig (the short one is likely an artefact of the
assembly), so create a fasta file containing just the longer contig.  

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${POMOXIS}
    REFNAME=Consensus_utg000001c
    samtools faidx ${DRAFT} Consensus_utg000001c >draft_assm/assm_final_filt.fa
    DRAFT=draft_assm/assm_final_filt.fa
    awk '{if(/>/){n=$1}else{print n " " length($0)}}' ${DRAFT}



.. _polishing_with_compressed_hp:

Polishing a Consensus with Run-length Encoding
----------------------------------------------

An experimental feature of medaka is to compress input basecalls and draft
assembly using run-length encoding and perform alignments using these
compressed sequences.  Limited tests on Ecoli suggests this improves consensus
accuracy, providing similar results to nanopolish (with homopolymer corrections
turned on), albeit at significantly higher speed. 

After performing all steps up to `Basecalling and draft assembly`, use the
following commands to run consensus using a model released with medaka. This
model was trained on Ecoli, Yeast and Human data. In this protype model, the
maximum homopolymer length is limited to 10. 

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${MEDAKA}
    DRAFT=draft_assm/assm_final_filt
    CONSENSUS=consensus
    medaka_consensus -i ${BASECALLS} -d ${DRAFT}.fa -o ${CONSENSUS} -t ${NPROC} -p ${POMOXIS}

To polish an assembly using another model (see `Training a Consensus Network Compressed Homopolymers`), use the `-m` option to specify the model. 

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${MEDAKA}
    DRAFT=draft_assm/assm_final_filt
    CONSENSUS=consensus
    MODEL=${TRAINNAME}/model.best.val.hdf5
    medaka_consensus -m ${MODEL} -i ${BASECALLS} -d ${DRAFT}.fa -o consensus -t ${NPROC} -p ${POMOXIS}

Finally, run `stats_from_bam` to assess to what extent `medaka` with
run-length encoding has improved accuracy. 

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${POMOXIS}
    TRUTH=${DATA}/truth 
    DRAFT2TRUTH=draft_to_truth
    CONSENSUS2TRUTH=${CONSENSUS}_to_truth
    CHUNK=10000
    mini_align -P -c ${CHUNK} -r ${TRUTH}.fasta -i ${DRAFT}.fa -p $DRAFT2TRUTH -t ${NPROC} 
    echo "Draft assembly"
    stats_from_bam --bam ${DRAFT2TRUTH}.bam > ${DRAFT2TRUTH}.stats.txt
    mini_align -P -c ${CHUNK} -r ${TRUTH}.fasta -i ${CONSENSUS}/consensus.fasta -p $CONSENSUS2TRUTH -t ${NPROC} 
    echo "Medaka hompolymer compression consensus"
    stats_from_bam --bam ${CONSENSUS2TRUTH}.bam > ${CONSENSUS2TRUTH}.stats.txt
    source ${MEDAKA}
    python -c "import sys; import pandas as pd; d=pd.read_table(sys.argv[-2]); m=pd.read_table(sys.argv[-1]); d['n']='draft'; m['n']='medaka'; c=pd.concat([d,m]); print(c.groupby('n')['acc','iden'].mean().T)" ${DRAFT2TRUTH}.stats.txt ${CONSENSUS2TRUTH}.stats.txt


Training a Consensus Network with Run-length Encoding
-----------------------------------------------------

After performing all steps up to `Basecalling and draft assembly`, use the
following commands to train a model using ren-length encoded features.

First compress the draft, reference and basecalls:

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${MEDAKA}
    DRAFT=draft_assm/assm_final_filt
    TRUTH=${DATA}/truth 
    DRAFTCOMPRFQ=${DRAFT}_compr.fq 
    TRUTHCOMPRFQ=${TRUTH}_compr.fq 
    BASECALLSCOMPRFQ=basecalls_compr.fq
    hp_compress compress ${DRAFT}.fa -t ${NPROC} > ${DRAFTCOMPRFQ}
    hp_compress compress ${TRUTH}.fasta -t ${NPROC} > ${TRUTHCOMPRFQ}
    hp_compress compress ${BASECALLS} -t ${NPROC} > ${BASECALLSCOMPRFQ}

Now align the compressed basecalls and compressed truth to the compressed draft consensus. Note that we chunk the truth reference prior to aligning it. 

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${POMOXIS}
    DRAFTCOMPRFA=${DRAFT}_compr.fa 
    fast_convert qa < ${DRAFTCOMPRFQ} > ${DRAFTCOMPRFA}
    COMPRCALLS2COMPRDRAFT=compr_calls_to_compr_draft
    COMPRTRUTH2COMPRDRAFT=compr_truth_to_compr_draft
    CHUNKSIZE=100000

    mini_align -P -m -r ${DRAFTCOMPRFA} -i ${BASECALLSCOMPRFQ} -t ${NPROC} -p ${COMPRCALLS2COMPRDRAFT}
    mini_align -c ${CHUNKSIZE} -P -m -r ${DRAFTCOMPRFA} -i ${TRUTHCOMPRFQ} -t ${NPROC} -p ${COMPRTRUTH2COMPRDRAFT}
    
Now create features for training. To reduce any IO bottlenecks during training, write training data to HDF5 in batches using the --batch_size option. To train a model which is more robust to variations in depth, use the --read_fraction option to randomly subsample reads. 

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${MEDAKA}
    REFNAME=Consensus_utg000001c
    TRAINEND=3761512
    TRAINFEATURES=hp_compress_train_features.hdf
    FRACTION="0.1 1"
    BATCHSIZE=200
    hp_compress features ${COMPRCALLS2COMPRDRAFT}.bam ${DRAFTCOMPRFQ} ${TRAINFEATURES} -T ${COMPRTRUTH2COMPRDRAFT}.bam -t ${NPROC} -r ${REFNAME}:-${TRAINEND} --batch_size ${BATCHSIZE} --read_fraction ${FRACTION} --chunk_len 1000 --chunk_ovlp 0


Now everything is in place to train a consensus network with compressed homopolymers using `medaka train`:

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${MEDAKA}
    TRAINNAME=training
    medaka train ${TRAINFEATURES} --train_name ${TRAINNAME}


Once training is finished, add feature-creation information to the model:

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${MEDAKA}
    medaka fix ${TRAINNAME}/model.best.val.hdf5 ${TRAINFEATURES}.yml

You can now use your model to polish a consensus by following steps in `Polishing a Consensus with Compressed Homopolymers`.


Training Models without Run-length Encoding
===========================================

It is still possible to train models which do not rely on run-length encoded inputs with `medaka`. 

Preparing Training Data
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
    source ${POMOXIS}
    TRUTH=${DATA}/truth.fasta 
    CALLS2DRAFT=calls_to_draft
    CHUNKSIZE=100000
    TRUTH2DRAFT=truth_to_draft
    mini_align -P -r ${DRAFT} -i ${BASECALLS} -t ${NPROC} -p ${CALLS2DRAFT}
    mini_align -P -c ${CHUNKSIZE} -r ${DRAFT} -i ${TRUTH} -t ${NPROC} -p ${TRUTH2DRAFT}

At the end of this process we have two `.bam` files which we use in the
following training step effectively to learn how to predict the contents of
one from the other.

We now use `medaka prepare` to generate training examples in the form of features calculated from  
chunks of labelled pileup, and save them to HDF5. We will train on the first 80% of our assembly,
saving the remaining 20% for evaluation. 

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${MEDAKA}
    TRAINFEATURES=train_features.hdf
    TRAINEND=3761512
    medaka prepare ${DRAFT} ${CALLS2DRAFT}.bam ${TRAINFEATURES} --truth ${TRUTH2DRAFT}.bam --ref_name ${REFNAME} --end ${TRAINEND}


Training the Consensus Network
------------------------------

We now have everything we need to train a consensus network using `medaka train`:

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${MEDAKA}
    TRAINNAME=training
    medaka train ${TRAINFEATURES} --train_name ${TRAINNAME} --max_label_len 1

During training, models are regularly checkpointed so that training may be
easily resumed if interrupted. At the end of training, we have a number of
output models including in particular:

    * `model.best.hdf5`: model with the best accuracy over the training set  
    * `model.best.val.hdf5`: model with the best accuracy over the validation set

which can be used to calculate a consensus. Other ancilliary output are
also produced.


Medaka Consensus Calling
------------------------

Having trained a model we can run `medaka consensus` to calculate a consensus
using our trained model:

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${MEDAKA}
    medaka consensus ${TRAINNAME}/model.best.val.hdf5 --alignments ${CALLS2DRAFT}.bam ${DRAFT} ${REFNAME} --start ${TRAINEND} --output_probs consensus_probs.hdf

The program outputs a HDF5 file containing consensus label probabilities for
overlapping chunks of the input. This stage may be parallelised trivially
by running the program on distinct sections of a draft assembly.

The consensus label probabilities may be of use in assessing consensus quality
and in variant calling; future releases of `medaka` may implement the writing
of VCF files.

Next we can recombine the consensus chunks using `medaka stitch`.

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${MEDAKA}
    medaka stitch consensus_probs.hdf medaka_consensus.fasta --mode hdf 


The output file `medaka_consensus.fasta` now contains our neural network consensus.

Finally, we can run `stats_from_bam` to assess to what extent `medaka` has improved accuracy. 

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${POMOXIS}
    EVALREGION=$(awk -F '[>:-]' '{if(NR==1){printf("%s:%i-%i\n",$2, $3, $4)}}' medaka_consensus.fasta)
    samtools faidx ${DRAFT} ${EVALREGION} > ${EVALREGION}_draft_assm.fa
    mini_align -P -r ${TRUTH} -i ${EVALREGION}_draft_assm.fa -t ${NPROC} -p draft_to_truth
    echo "Draft assembly"
    stats_from_bam --bam draft_to_truth.bam > draft_to_truth.stats.txt
    mini_align -P -r ${TRUTH} -i medaka_consensus.fasta -t ${NPROC} -p consensus_to_truth
    echo "Medaka consensus"
    stats_from_bam --bam consensus_to_truth.bam > consensus_to_truth.stats.txt

An increase in accuracy should be observed.
