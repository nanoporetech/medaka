import logging
import numpy as np
from Bio import SeqIO
from Bio import pairwise2
from collections import defaultdict


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


def assemble_sequences(fastas, min_overlap_len=50):
    """Join overlapping sequences to assemble complete sequence

    :param fastas: iterable of fasta filepaths.
    :min_overlap_len: minimum required length of overlaps at sequence ends
    :returns: list of (contig_name, sequence) 

    """
    seq_records = []
    for fasta in fastas:
        seq_records += list(SeqIO.parse(fasta, 'fasta'))

    ref_assemblies = []

    segments_by_name = defaultdict(dict)
    for rec in seq_records:
        (ref, segment_range) = rec.id.split(':')
        segment_start, segment_end = [int(i) for i in segment_range.split('-')]
        segment_end -= 1
        segments_by_name[ref][(segment_start, segment_end)] = str(rec.seq)

    for ref in sorted(segments_by_name.keys()):
        ref_assembly = ''
        prev_segment_end = np.inf
        start = None
        for segment_start, segment_end in sorted(segments_by_name[ref]):
            if start is None:
                start = segment_start
            # if we have a gap, add what we have joined up 
            # and move onto the next chunk we can join together
            enough_overlap = prev_segment_end - min_overlap_len > segment_start
            if not enough_overlap:
                msg = '{} will be fragmented: overlap between end {} and start {} is {}, less than required {}'
                logging.info(msg.format(ref, prev_segment_end, segment_start, 
                                        prev_segment_end - segment_start, 
                                        min_overlap_len)
                )
                # start + len(ref_assembly) could exceed prev_segment_end if correction 
                # has filled in deletion errors, so to know if we have a gap we should use
                # segment_start and prev_segment_end
                #name = '{}:{}-{}'.format(ref, start, start + len(ref_assembly) + 1) 
                name = '{}:{}-{}'.format(ref, start, prev_segment_end)
                ref_assemblies.append((name, ref_assembly))
                ref_assembly = ''
                prev_segment_end = np.inf
                start = segment_start
                continue

            sequence = segments_by_name[ref][(segment_start, segment_end)]
            overlap_len = prev_segment_end - segment_start
            ref_assembly = merge_into_sequence(ref_assembly, sequence, overlap_len)
            prev_segment_start = segment_start
            prev_segment_end = segment_end

        name = '{}:{}-{}'.format(ref, start, prev_segment_end) 
        ref_assemblies.append((name, ref_assembly))

    return ref_assemblies


def write_fasta(filename, contigs):
    with open(filename, 'w') as fasta:
        for name, seq in contigs:
            fasta.write('>{}\n{}\n'.format(name, seq))


def stitch(args):
    joined = assemble_sequences(args.fastas, min_overlap_len=args.min_overlap)
    write_fasta(args.output, joined)
