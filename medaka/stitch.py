import logging
from itertools import chain

import numpy as np

from medaka.common import (_gap_, Region, get_sample_overlap,)
from medaka.datastore import DataIndex


def write_fasta(filename, contigs):
    with open(filename, 'w') as fasta:
        for name, seq in contigs:
            fasta.write('>{}\n{}\n'.format(name, seq))


def stitch_from_probs(probs_hdfs, regions=None, model_yml=None):
    """Decode and join overlapping label probabilities from hdf to assemble complete sequence

    :param probs_hdfs: iterable of hdf filepaths.
    :ref_names: iterable of ref_names to limit stitching to.
    :returns: list of (contig_name, sequence)

    """
    logger = logging.getLogger(__package__)
    logger.name = 'Stitch'

    index = DataIndex(probs_hdfs)

    label_decoding = index.meta['_label_decod_path_']

    logger.debug("Label decoding is:\n{}".format('\n'.join(
        '{}: {}'.format(i, x) for i, x in enumerate(label_decoding)
    )))
    if regions is None:
        ref_names = index.index.keys()
    else:
        #TODO: respect entire region specification
        ref_names = list()
        for region in (Region.from_string(r) for r in regions):
            if region.start is None or region.end is None:
                logger.warning("Ignoring start:end for '{}'.".format(region))
            ref_names.append(region.ref_name)

    get_pos = lambda s, i: '{}.{}'.format(s.positions[i]['major'] + 1, s.positions[i]['minor'])
    ref_assemblies = []
    for ref_name in ref_names:
        logger.info("Processing {}.".format(ref_name))
        data_gen = index.yield_from_feature_files(ref_names=(ref_name,))
        seq=''
        s1 = next(data_gen)
        start = get_pos(s1, 0)
        start_1_ind = None
        for s2 in chain(data_gen, (None,)):
            if s2 is None:  # s1 is last chunk
                end_1_ind = None
            else:
                end_1_ind, start_2_ind = get_sample_overlap(s1, s2)

            best = np.argmax(s1.label_probs[start_1_ind:end_1_ind], -1)
            seq += ''.join([label_decoding[x] for x in best]).replace(_gap_, '')
            if end_1_ind is None:
                key = '{}:{}-{}'.format(s1.ref_name, start, get_pos(s1, -1))
                ref_assemblies.append((key, seq))
                seq = ''
                if start_2_ind is None:
                    msg = 'There is no overlap betwen {} and {}'
                    logging.info(msg.format(s1.name, s2.name))
                    start = get_pos(s2, 0)

            s1 = s2
            start_1_ind = start_2_ind
    return ref_assemblies


def stitch(args):
    joined = stitch_from_probs(args.inputs, regions=args.regions)
    write_fasta(args.output, joined)
