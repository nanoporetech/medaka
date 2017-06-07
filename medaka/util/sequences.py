import itertools

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2

from medaka.util.generators import serve_sample, serve_data_batch

decode = {1: 'A', 2: 'C', 3: 'G', 4: 'T', 5: '*'}


def get_benchmark(X, y, steps, batch_size, window_size):
    """Identify proportion ref positions in which reference matches truth,
    a starting point for the correction

    :param X: feature array
    :param y: label array
    :param steps: num. of batches
    :param batch_size: int batch_size
    :param window_size: int window_size
    :returns: float proportion of correct reference positions
    """
    sample_generator = itertools.islice(
        serve_sample(X, y, window_size), None, batch_size * steps)
    benchmark = [data[window_size // 2, 0] == np.argmax(label)
                 for data, label in sample_generator]
    return np.mean(benchmark)


def model_output_to_fasta(predictions, name):
    """Convert from one hot encoding to bases, excising any deletions

    :param predictions: numpy array model output, one-hot encoding (num_classes 6)
                        0=blank, 1=A, 2=C, 3=G, 4=T, 5=Del
    :param name: str prefix for output fasta files
    :returns: SeqRecord corrected fasta sequence
    """
    corrected_sequence = ''.join(decode[np.argmax(p)]
                                for p in predictions if not np.argmax(p) == 5)
    corrected_sequence = SeqRecord(Seq(corrected_sequence), name, '', '')
    return corrected_sequence


def get_original_sequence(data, name, window_size):
    """Convert from encoded reference base column, excising any deletions
    and respecting truncations due to windowing

    :param data: feature array
    :param name: str prefix for output fasta files
    :returns: SeqRecord original fasta sequence
    """
    original_sequence = ''.join(
        decode[i] for i in data[window_size // 2: data.shape[0] - window_size + 1, 0] if not i == 5
    )
    original_sequence = SeqRecord(Seq(original_sequence), name, '', '')
    return original_sequence


def write_test_sequences(test_data, model, test_steps, batch_size, window_size):
    """Extract and write original and correctd fasta sequences from test data

    :param test_data: test data array
    :param model: keras model
    :param test_steps: int number of steps (batches)
    :param batch_size: int batch_size
    :param window_size: int window_size
    """
    test_data_generator = serve_data_batch(test_data, batch_size, window_size)
    predictions = model.predict_generator(test_data_generator, test_steps)
    corrected = model_output_to_fasta(predictions, 'test_sequence')
    SeqIO.write(corrected, 'test_corrected_sequence.fa', 'fasta')
    original = get_original_sequence(test_data, 'test_sequence', window_size)
    SeqIO.write(original, 'test_original_sequence.fa', 'fasta')


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
    if existing == "":
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


def assemble_sequences(sequences, segment_len=50000, overlap_len=1000):
    """Join overlapping sequences to assemble complete sequence

    :param sequences: list of biopython SeqRecords
    :segment_len: length of sequences
    :overlap_len: length of overlaps at sequence ends
    :returns: list of reference SeqRecords

    .. note::

        reproduced almost verbatim from
         `nanopolish <https://github.com/jts/nanopolish/blob/master/scripts/nanopolish_makerange.py>`_
    """
    ref_assemblies = []

    segments_by_name = dict()
    for rec in sequences:
        (ref, segment_range) = rec.id.split(":")
        if ref not in segments_by_name:
            segments_by_name[ref] = dict()
        segment_start, segment_end = segment_range.split("-")
        segments_by_name[ref][int(segment_start)] = str(rec.seq)

    for ref in sorted(segments_by_name.keys()):
        ref_assembly = ""
        prev_segment = None
        for segment_start in sorted(segments_by_name[ref]):
            # Ensure the segments overlap
            assert(prev_segment is None or
                   prev_segment + segment_len + overlap_len > segment_start)

            sequence = segments_by_name[ref][segment_start]

            ref_assembly = merge_into_sequence(ref_assembly, sequence,
                                                   overlap_len)
            prev_segment = segment_start

        ref_assembly = SeqRecord(Seq(ref_assembly), ref, '', '')
        ref_assemblies.append(ref_assembly)

    return ref_assemblies


def segment_limits(endpoint, segment_len=50000, overlap_len=1000):
    """Generate segments of a range [0, end_point].
   
    :param endpoint: endpoint of range.
    :param segment_len: length of resultant segments.
    :param overlap_len: length of overlap between segments.
    :yields: tuples (start, end) indicating overlapping segments of
        the input range. 
    """
    for n in range(0, endpoint, segment_len):
        if (n + segment_len) > endpoint:
            yield n, endpoint - 1
        else:
            yield n, n + segment_len + overlap_len
