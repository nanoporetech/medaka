
.. _walkthrough:

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
dataset is required. Medaka's standard models are trained using 
`katuali <https://nanoporetech.github.io/katuali/medaka_train.html>`_.


Obtaining Data and Software
---------------------------

Start by downloading training data, containing basecalls from Oxford Nanopore
Technologies' reads:

.. code-block:: bash

    WALKTHROUGH=${PWD}/medaka_walkthrough
    mkdir -p ${WALKTHROUGH} && cd ${WALKTHROUGH}
    wget https://s3-eu-west-1.amazonaws.com/ont-research/medaka_walkthrough_no_reads.tar.gz
    tar -xvf medaka_walkthrough_no_reads.tar.gz
    DATA=${PWD}/data

The extracted archive contains also all the intermediate output files that
are created during the process below. Any step may be skipped by simply copying
the requisite subfolder from the ``${DATA}`` directory into the ``${WALKTHROUGH}``
directory.

The necessary software can be sourced using:

.. code-block:: bash

    cd ${WALKTHROUGH}
    git clone https://github.com/nanoporetech/pomoxis --recursive
    git clone https://github.com/nanoporetech/medaka
    
    # While it is possible to install pomoxis and medaka into the same
    #   virtual environment, we will install each package into its own
    #   environment for simplicity. For more details see the readme for
    #   each of the packages.

    cd pomoxis && make install && cd ..
    cd medaka && make install && cd ..

    POMOXIS=${WALKTHROUGH}/pomoxis/venv/bin/activate
    MEDAKA=${WALKTHROUGH}/medaka/venv/bin/activate


.. _basecalling_and_draft_assembly:

Creating a Draft Assembly
-------------------------

A draft assembly can be formed from the provided basecalls using the 
`miniasm <https://github.com/lh3/miniasm>`_ and
`racon <https://github.com/isovic/racon>`_ based pipeline from ``pomoxis``.
Alternatively one could use `canu <https://github.com/marbl/canu>`_ at this step.

.. code-block:: bash

    cd ${WALKTHROUGH}
    NPROC=$(nproc)
    BASECALLS=data/basecalls.fa
    source ${POMOXIS}
    mini_assemble -i ${BASECALLS} -o draft_assm -p assm -t ${NPROC}

This will create a draft assembly at ``draft_assm/assm_final.fa``. The
``mini_assemble`` script has two useful options not used here:

    * specifying ``-c`` will run `porechop <https://github.com/rrwick/Porechop>`_
      on the reads to first trim sequencing adapters and,
    * specifying ``-e 10`` will perform error correction on the longest 10% of
      reads prior to assembly (similar to the strategy of canu).

Both these steps can improve the assembly quality at the expense of speed.

The number and length of the assembled contigs can be checked

.. code-block:: bash

    DRAFT=draft_assm/assm_final.fa
    awk '{if(/>/){n=$1}else{print n " " length($0)}}' ${DRAFT}

The expected output is a contig 4,703,280 bases long (utg000001c). 

.. _polishing:

Polishing a Consensus 
----------------------

After performing all steps up to :ref:`basecalling_and_draft_assembly`, the
following commands can be run to yield a consensus using ``medaka``'s default
model. This model was trained using data obtained from E.coli, S.cerevisaie,
and H.sapiens samples. 

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${MEDAKA}
    CONSENSUS=consensus
    DRAFT=draft_assm/assm_final.fa
    medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${CONSENSUS} -t ${NPROC}

To polish an assembly using another model, use
the ``-m`` option to specify the filepath of the model.

Alignment statistics can be calculated using the ``assess_assembly`` program from
pomoxis: 

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${POMOXIS}
    TRUTH=${DATA}/truth.fasta
    DRAFT2TRUTH=draft_to_truth
    CONSENSUS2TRUTH=${CONSENSUS}_to_truth
    echo "Draft assembly"
    assess_assembly -i ${DRAFT} -r ${TRUTH} -p ${DRAFT2TRUTH} -t ${NPROC}
    echo "Medaka consensus"
    assess_assembly -i ${CONSENSUS}/consensus.fasta -r ${TRUTH} -p ${CONSENSUS2TRUTH} -t ${NPROC}

An decrease in error rate from 0.367% to 0.070% should be observed.

.. _training:

Training a Consensus Network
----------------------------

In order to train a bespoke network first perform all the steps up to and
including :ref:`basecalling_and_draft_assembly` above. 

The ultimate aim of the consensus network is to predict the truth sequence from
the alignment of basecalls to the draft. This requires understanding how the
basecalls may align to the draft and how the draft must be edited to obtain the
truth. The draft acts as a common frame-of-reference between the basecalls
and the truth.

The basecalls and truth sequence are aligned to the draft. For the latter, this
is performed in chunks.

.. code-block:: bash

    cd ${WALKTHROUGH}
    DRAFT=draft_assm/assm_final.fa
    TRUTH=${DATA}/truth.fasta
    source ${POMOXIS}
    CHUNKSIZE=100000
    CALLS2DRAFT=calls2draft
    TRUTH2DRAFT=truth2draft

    mini_align -P -m -r ${DRAFT} -i ${BASECALLS} -t ${NPROC} -p ${CALLS2DRAFT}
    mini_align -c ${CHUNKSIZE} -P -m -r ${DRAFT} -i ${TRUTH} -t ${NPROC} -p ${TRUTH2DRAFT}

These raw alignments must now be converted into features for input into a neural
network. To reduce any IO bottlenecks during training, the training data can be
written to the ``HDF5`` file in batches using the ``-\\-batch_size`` option. The option
``-\\-read_fraction`` is used to randomly subsample reads which has the effect of
making the resultant model more robust to variations in pileup depth when the
model is used to make predictions.

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${MEDAKA}
    REFNAME=utg000001c
    TRAINEND=3762624
    TRAINFEATURES=train_features.hdf
    BATCHSIZE=100
    medaka features ${CALLS2DRAFT}.bam ${TRAINFEATURES} --truth ${TRUTH2DRAFT}.bam --threads ${NPROC} --region ${REFNAME}:-${TRAINEND} --batch_size ${BATCHSIZE} --chunk_len 1000 --chunk_ovlp 0 

Now everything is in place to train a consensus network with ``medaka train``:

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${MEDAKA}
    TRAINNAME=training
    medaka train ${TRAINFEATURES} --train_name ${TRAINNAME} --epochs 10

Depending on the compute resources available, this step may take some time.
During training, models are regularly checkpointed so that training may be
easily resumed if interrupted. At the end of training, we have a number of
output models including in particular:

    * ``model.best.hdf5``: model with the best accuracy over the training set  
    * ``model.best.val.hdf5``: model with the best accuracy over the validation set

Other ancilliary output are also produced. 

To use a model run ``medaka_consensus``, specifying the full absolute path to
the model using the ``-m`` option:

.. code-block:: bash

    cd ${WALKTHROUGH}
    source ${MEDAKA}
    CONSENSUS=consensus_trained
    MODEL=${PWD}/${TRAINNAME}/model.best.val.hdf5
    medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${CONSENSUS} -t ${NPROC} -m ${MODEL}


Training with multipe data types
--------------------------------

Medaka supports creating inference networks that handle independently data from
difference sources, for example two distinct nanopores. The is enabled through
the use of alignment tags in the input ``.bam`` files.

To train a model which can handle multiple datatypes follow the same process as above
for all datatypes independently until the ``medaka feature`` step. Before this combine
the ``.bam`` files whilst adding the string tag ``DT`` to the files. For example to
combine data for the R9 and R10 pores:

.. code-block:: bash

    samtools view r9.bam | awk 'BEGIN{OFS="\t"}{print $0, "DT:Z:r9"}' > r9.sam
    samtools view r10.bam | awk 'BEGIN{OFS="\t"}{print $0, "DT:Z:r10"}' > r10.sam
    samtools view -H r9.bam | cat - r9.sam r10.sam | samtools view -bS | samtools sort > all_reads.bam
    samtools index all_reads.bam

Having combined the data types training a model simply requires an additional
argument to ``medaka features``:

.. code-block:: bash

    medaka features all_reads.bam ... --feature_encoder_args dtypes=r9,r10

In order to create consensus sequences using models trained for multiple
datatypes, the input ``.bam`` to ``medaka consensus`` should have the ``DT``
tag added appropriately to all reads. Models trained for multiple datatypes
will not work without this tag being added to the input files.


Automated training pipeline
---------------------------

With `katuali <https://nanoporetech.github.io/katuali/medaka_train.html>`_ it
is now possible to train medaka models starting from folders of fast5s in a single
command:

.. code-block:: bash

    katuali medaka_train_replicates --keep-going

Running the above will

    * basecall all multiple runs of data,
    * align all basecalls to reference sequences,
    * create subsampled sets of basecalls over reference sequences and depths,
    * assemble those sets of basecalls into draft assemblies,
    * create medaka training features for all assemblies,
    * train a medaka in multiple replicates, and
    * evaluate the models on test data.

For further information concerning settig-up adnd running katuali, refer to
its `documentation <https://nanoporetech.github.io/>`_.

