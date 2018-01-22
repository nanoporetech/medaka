
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
dataset is required. `medaka` is designed for flexibility over performance,
though see :doc:`benchmarks` for a speed comparison with other commonly used
tools.


Obtaining Data and Software
---------------------------

Start by downloading training data in the form of Oxford Nanopore
Technologies `.fast5` files for a set of reads.

.. code-block:: bash

    WALKTHROUGH=${PWD}/medaka_walkthrough
    mkdir -p ${WALKTHROUGH} && cd ${WALKTHROUGH}
    wget https://s3-eu-west-1.amazonaws.com/ont-medaka-demo/medaka_demo.tar
    tar -xvf medaka_demo.tar
    DATA=${PWD}/data

The extracted archive contains also all the intermediate output files that
are created during the process below. Any step may be skipped by simply copying
the requisite subfolder from the `${DATA}` directory into the `${WALKTHROUGH}`
directory.

The necessary software can be sourced using the same process as described in
:ref:`CreatingSoftwareEnv`, namely:

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



.. _basecalling_and_draft_assembly:

Basecalling and Draft Assembly
------------------------------

The downloaded data can be basecalled using `scrappie`:

.. code-block:: bash

    cd ${WALKTHROUGH}
    SCRAPPIE=${WALKTHROUGH}/scrappie/build/scrappie
    NPROC=$(nproc)
    export OMP_NUM_THREADS=${NPROC}
    export OPENBLAS_NUM_THREADS=1
    BASECALLS=basecalls.fa
    ${SCRAPPIE} raw ${DATA}/reads --model rgrgr_r94 > ${BASECALLS}

and a draft assembly can be formed using the 
`miniasm <https://github.com/lh3/miniasm>`_ and
`racon <https://github.com/isovic/racon>`_ based pipeline from `pomoxis`.
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

The number and length of the assembled contigs can be checked

.. code-block:: bash

    cd ${WALKTHROUGH}
    DRAFT=draft_assm/assm_final.fa
    awk '{if(/>/){n=$1}else{print n " " length($0)}}' ${DRAFT}

The expected output is a contig 4,701,891 bases long (Consensus_utg000001c) and
a short remainder sequence just 408 bases long (Consensus_utg000002c). The
following will use only the long contig, which is around the expected genome
size. To create a `.fasta` file containing just the longer contig, run:  

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${POMOXIS}
    REFNAME=Consensus_utg000001c
    samtools faidx ${DRAFT} Consensus_utg000001c > draft_assm/assm_final_filt.fa
    DRAFT=draft_assm/assm_final_filt.fa


.. _polishing_with_rle:

Polishing a Consensus with Run-length Encoding
----------------------------------------------

An feature of `medaka` is to compress input basecalls and the draft
assembly using run-length encoding and perform alignments using these
compressed sequences. Tests with E.coli data suggests this improves consensus
accuracy, providing similar results to nanopolish (with homopolymer corrections
turned on), at significantly higher speed.

After performing all steps up to :ref:`basecalling_and_draft_assembly`, the
following commands can be run to yield a consensus using `medaka`'s default
model. This model was trained using data obtained from E.coli, S.cerevisaie,
and H.sapiens samples. The maximum homopolymer length that this model will call
successfully is limited to 10.

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${MEDAKA}
    DRAFT=draft_assm/assm_final_filt
    CONSENSUS=consensus
    medaka_consensus -i ${BASECALLS} -d ${DRAFT}.fa -o ${CONSENSUS} -t ${NPROC} -p ${POMOXIS}

To polish an assembly using another model (see :ref:`training_with_rle`), use
the `-m` option to specify the filepath of the model.. 

Alignment statistics can be calculated using the `stats_from_bam` program from
pomoxis: 

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
    echo "Medaka RLE consensus"
    stats_from_bam --bam ${CONSENSUS2TRUTH}.bam > ${CONSENSUS2TRUTH}.stats.txt
    source ${MEDAKA}
    python -c "import sys; import pandas as pd; d=pd.read_table(sys.argv[-2]); m=pd.read_table(sys.argv[-1]); d['n']='draft'; m['n']='medaka'; c=pd.concat([d,m]); print(c.groupby('n')['acc','iden'].mean().T)" ${DRAFT2TRUTH}.stats.txt ${CONSENSUS2TRUTH}.stats.txt



.. _training_with_rle:

Training a Consensus Network
----------------------------

In order to train a bespoke network first perform all the steps up to and
including :ref:`basecalling_and_draft_assembly` above. Following this the first
task is to perform the run length encoding of the three inputs: the truth
sequence, the draft assembly, and the basecalls:

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


The ultimate aim of the consensus network is to predict the truth sequence from
the alignment of basecalls to the draft. This requires understanding how the
basecalls may align to the draft and how the draft much be edited to obtain the
truth. The draft acts as a common frame-of-reference between the basecalls
and the truth.

The compressed basecalls and truth sequence are aligned to the compressed
draft. For the latter, this is performed in chunks.

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

These raw alignments must now be converted into features for input into a neural
network. To reduce any IO bottlenecks during training, the training data can be
written to the `HDF5` file in batches using the `-\\-batch_size` option. The option
`-\\-read_fraction` is used to randomly subsample reads which has the effect of
making the resultant model more robust to variations in pileup depth when the
model is used to make predictions.

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${MEDAKA}
    REFNAME=Consensus_utg000001c
    TRAINEND=3761512
    TRAINFEATURES=rle_train_features.hdf
    FRACTION="0.1 1"
    BATCHSIZE=200
    hp_compress features ${COMPRCALLS2COMPRDRAFT}.bam ${DRAFTCOMPRFQ} ${TRAINFEATURES} -T ${COMPRTRUTH2COMPRDRAFT}.bam -t ${NPROC} -r ${REFNAME}:-${TRAINEND} --batch_size ${BATCHSIZE} --read_fraction ${FRACTION} --chunk_len 1000 --chunk_ovlp 0

Now everything is in place to train a consensus network with the run-length
encoded features with `medaka train`:

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${MEDAKA}
    TRAINNAME=training
    medaka train ${TRAINFEATURES} --train_name ${TRAINNAME}

Depending on the compute resources available, this step may take some time.
During training, models are regularly checkpointed so that training may be
easily resumed if interrupted. At the end of training, we have a number of
output models including in particular:

    * `model.best.hdf5`: model with the best accuracy over the training set  
    * `model.best.val.hdf5`: model with the best accuracy over the validation set

Other ancilliary output are also produced. The final model can be combined with
its meta information in order to make it ready for use:

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${MEDAKA}
    medaka fix ${TRAINNAME}/model.best.val.hdf5 ${TRAINFEATURES}.yml

To use the model run `medaka_consensus` for the default model (specifying
the model using the `-m` option):

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${MEDAKA}
    DRAFT=draft_assm/assm_final_filt
    CONSENSUS=consensus
    MODEL=${TRAINNAME}/model.best.val.hdf5
    medaka_consensus -m ${MODEL} -i ${BASECALLS} -d ${DRAFT}.fa -o consensus -t ${NPROC} -p ${POMOXIS}

