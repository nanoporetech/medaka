import logging
import numpy as np
from Bio import SeqIO
from Bio import pairwise2
from itertools import chain
from medaka.common import (_gap_, encode_sample_name, get_sample_overlap, get_sample_index_from_files,
                           load_yaml_data, yield_from_feature_files, _label_decod_path_)

def merge_into_sequence(existing, incoming, overlap_length):
    """Append string to another string using pairwise alignment of
    overlapping ends to select join

    :param existing: str growing sequence
    :param incoming: str new sequence to be appended
    :param overlap_length: int length of overlap between sequence ends
    :returns: str new sequence

    .. note::

        reproduced almost verbatim from
        `nanopolish_merge <https://github.com/jts/nanopolish/blob/master/scripts/nanopolish_merge.py>`_
    """

    # if first segment, no overlapping needs to be done
    if existing == '':
        return incoming

    or_con = existing[-overlap_length:]
    or_inc = incoming[0:overlap_length]

    # These parameters are designed to give us the region of highest similarity
    # between the two sequences
    alignments = pairwise2.align.globalms(or_con, or_inc, 2, -10, -10, -3)

    best = alignments[0]
    aln_con, aln_inc, score, begin, end = best

    # We merge at the midpoint between the two aligned segments
    m_con = 0
    m_inc = 0

    assert(len(aln_con) == len(aln_inc))

    for i in range(0, len(aln_con) // 2):
        a = aln_con[i]
        b = aln_inc[i]

        if a != '-':
            m_con += 1
        if b != '-':
            m_inc += 1

    m_con += len(existing) - overlap_length
    merged = existing[0:m_con] + incoming[m_inc:]

    return merged


def stitch_from_fastas(fastas, min_overlap_len=50, ref_names=None):
    """Join overlapping sequences to assemble complete sequence

    :param fastas: iterable of fasta filepaths.
    :min_overlap_len: int, minimum required overlap in reference positions at sequence ends
    :ref_names: iterable of ref_names to limit stitching to.
    :returns: list of (contig_name, sequence)

    """
    index = get_sample_index_from_files(fastas, 'fasta')
    if ref_names is None:
        ref_names = index.keys()
    fhs = {f: SeqIO.index(f, 'fasta') for f in fastas}

    ref_assemblies = []

    for ref in index:
        ref_assembly = ''
        prev_chunk_end = np.inf
        start = None
        for chunk in index[ref]:
            if start is None:
                start = chunk['start']
            # if we have a gap, add what we have joined up
            # and move onto the next chunk we can join together
            is_overlapping = float(chunk['start']) < float(prev_chunk_end)
            # check we have at least min_overlap_len major reference positions
            enough_overlap = float(prev_chunk_end) - float(chunk['start']) > min_overlap_len
            if not is_overlapping or not enough_overlap:
                msg = '{} will be fragmented: insufficient overlap between end {} and start {}'
                logging.info(msg.format(ref, prev_chunk_end, chunk['start']))
                name = '{}:{}-{}'.format(ref, start, prev_chunk_end)
                ref_assemblies.append((name, ref_assembly))
                ref_assembly = ''
                prev_chunk_end = np.inf
                start = chunk['start']
                continue

            sequence = fhs[chunk['filename']][chunk['key']].seq
            if ref_assembly == '':
                ref_assembly = sequence
            else:
                overlap_len = int(float(prev_chunk_end)) - int(float(chunk['start']))
                ref_assembly = merge_into_sequence(ref_assembly, sequence, overlap_len)

            prev_chunk_end = chunk['end']

        name = '{}:{}-{}'.format(ref, start, prev_chunk_end)
        ref_assemblies.append((name, ref_assembly))

    return ref_assemblies


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

            if end_1_ind is None and start_2_ind is None:
                msg = 'There is no overlap betwen {} and {}'
                logging.info(msg.format(encode_sample_name(s1),
                                        encode_sample_name(s2)))

            best = np.argmax(s1.label_probs[start_1_ind:end_1_ind], -1)
            seq += ''.join([label_decoding[x] for x in best]).replace(_gap_, '')
            if end_1_ind is None:
                key = '{}:{}-{}'.format(s1.ref_name, start, get_pos(s1, -1))
                ref_assemblies.append((key, seq))
                seq = ''
            s1 = s2
            start_1_ind = start_2_ind
    return ref_assemblies


def stitch(args):

    if args.mode == 'fasta':
        joined = stitch_from_fastas(args.inputs, min_overlap_len=args.min_overlap, ref_names=args.ref_names)
    else:
        joined = stitch_from_probs(args.inputs, ref_names=args.ref_names, model_yml=args.model_yml)
    write_fasta(args.output, joined)
