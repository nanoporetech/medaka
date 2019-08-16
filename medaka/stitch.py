from itertools import chain

import numpy as np

import medaka.common
from medaka.common import Sample, Region
import medaka.datastore


def write_fasta(filename, contigs):
    with open(filename, 'w') as fasta:
        for name, seq in contigs:
            fasta.write('>{}\n{}\n'.format(name, seq))


def stitch_from_probs(probs_hdfs, regions=None, model_yml=None):
    """Decode and join overlapping label probabilities from hdf
    to assemble complete sequence

    :param probs_hdfs: iterable of hdf filepaths.
    :ref_names: iterable of ref_names to limit stitching to.

    :returns: list of (contig_name, sequence)
    """
    logger = medaka.common.get_named_logger('Stitch')

    index = medaka.datastore.DataIndex(probs_hdfs)

    label_decoding = index.meta['medaka_label_decoding']

    logger.debug("Label decoding is:\n{}".format('\n'.join(
        '{}: {}'.format(i, x) for i, x in enumerate(label_decoding)
    )))
    if regions is None:
        ref_names = index.index.keys()
    else:
        # TODO: respect entire region specification
        ref_names = list()
        for region in (Region.from_string(r) for r in regions):
            if region.start is None or region.end is None:
                logger.warning("Ignoring start:end for '{}'.".format(region))
            ref_names.append(region.ref_name)

    def get_pos(s, i):
        return '{}.{}'.format(
            s.positions[i]['major'] + 1, s.positions[i]['minor'])
    ref_assemblies = []
    for ref_name in ref_names:
        logger.info("Processing {}.".format(ref_name))
        data_gen = index.yield_from_feature_files(ref_names=(ref_name,))
        seq_parts = list()
        s1 = next(data_gen)
        start = get_pos(s1, 0)
        start_1_ind = None
        start_2_ind = None

        for s2 in chain(data_gen, (None,)):
            s1_name = 'Unknown' if s1 is None else s1.name
            s2_name = 'Unknown' if s2 is None else s2.name

            if s2 is None:  # s1 is last chunk
                end_1_ind = None
            else:
                if s2.last_pos <= s1.last_pos:
                    logger.info('{} ends before {}, skipping.'.format(
                        s2_name, s1_name
                    ))
                    continue
                elif s2.first_pos >= s1.last_pos:
                    # trigger a break
                    end_1_ind, start_2_ind = None, None
                else:
                    try:
                        end_1_ind, start_2_ind = Sample.overlap_indices(s1, s2)
                    except medaka.common.OverlapException as e:
                        logger.info(
                            "Unhandled overlap type whilst stitching chunks.")
                        raise(e)

            best = np.argmax(s1.label_probs[start_1_ind:end_1_ind], -1)
            new_seq = ''.join([label_decoding[x] for x in best])
            new_seq.replace(medaka.common._gap_, '')
            seq_parts.append(new_seq)
            if end_1_ind is None:
                key = '{}:{}-{}'.format(s1.ref_name, start, get_pos(s1, -1))
                ref_assemblies.append((key, ''.join(seq_parts)))
                seq_parts = list()
                if s2 is not None and start_2_ind is None:
                    msg = 'There is no overlap betwen {} and {}'
                    logger.info(msg.format(s1_name, s2_name))
                    start = get_pos(s2, 0)

            s1 = s2
            start_1_ind = start_2_ind
    return ref_assemblies


def stitch(args):
    joined = stitch_from_probs(args.inputs, regions=args.regions)
    write_fasta(args.output, joined)
