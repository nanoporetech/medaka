import logging
import numpy as np
from itertools import chain
from medaka.common import (_gap_, encode_sample_name, get_sample_overlap, get_sample_index_from_files,
                           load_yaml_data, yield_from_feature_files, _label_decod_path_)


def write_fasta(filename, contigs):
    with open(filename, 'w') as fasta:
        for name, seq in contigs:
            fasta.write('>{}\n{}\n'.format(name, seq))


def stitch_from_probs(probs_hdfs, ref_names=None, model_yml=None):
    """Decode and join overlapping label probabilities from hdf to assemble complete sequence

    :param probs_hdfs: iterable of hdf filepaths.
    :ref_names: iterable of ref_names to limit stitching to.
    :returns: list of (contig_name, sequence)

    """
    label_decoding = load_yaml_data(probs_hdfs[0], _label_decod_path_)
    if label_decoding is None:
        if model_yml is not None:
            logging.info("Loading label encoding from {}".format(model_yml))
            label_decoding = load_yaml_data(model_yml, _label_decod_path_)
        else:
            raise ValueError('Cannot decode probabilities without label decoding')
    logging.info("Label decoding is:\n{}".format('\n'.join(
        '{}: {}'.format(i, x) for i, x in enumerate(label_decoding)
    )))
    index = get_sample_index_from_files(probs_hdfs, 'hdf')
    if ref_names is None:
        ref_names = index.keys()
    get_pos = lambda s, i: '{}.{}'.format(s.positions[i]['major'] + 1, s.positions[i]['minor'])
    ref_assemblies = []
    for ref_name in ref_names:
        data_gen = yield_from_feature_files(probs_hdfs, ref_names=(ref_name,), index=index)
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
                    logging.info(msg.format(encode_sample_name(s1),
                                            encode_sample_name(s2)))
                    start = get_pos(s2, 0)

            s1 = s2
            start_1_ind = start_2_ind
    return ref_assemblies


def stitch(args):

    joined = stitch_from_probs(args.inputs, ref_names=args.ref_names, model_yml=args.model_yml)
    write_fasta(args.output, joined)
