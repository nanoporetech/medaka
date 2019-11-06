"""Run length encoding and realignment of reads."""
import array
import concurrent.futures
import functools
import os
from pathlib import Path
import re

import h5py
import numpy as np
import parasail
import pysam

import medaka.common
import medaka.smolecule


def rle(iterable):
    """Calculate a run length encoding (rle), of an input iterable.

    :param iterable: input iterable.

    :returns: structured array with fields `start`, `length`, and `value`.
    """
    if not isinstance(iterable, np.ndarray):
        array = np.fromiter(iterable, dtype='U1', count=len(iterable))
    else:
        array = iterable

    if len(array.shape) != 1:
        raise TypeError("Input array must be one dimensional.")
    dtype = [('length', int), ('start', int), ('value', array.dtype)]

    n = len(array)
    starts = np.r_[0, np.flatnonzero(array[1:] != array[:-1]) + 1]
    rle = np.empty(len(starts), dtype=dtype)
    rle['start'] = starts
    rle['length'] = np.diff(np.r_[starts, n])
    rle['value'] = array[starts]
    return rle


class RLEConverter(object):
    """Class to convert a basecall to RLE, with coordinate conversions."""

    def __init__(self, basecall):
        """Create an RLE conversion helper.

        :param basecall: string to be converted to RLE.
        """
        self.basecall = basecall
        self.rle_conversion = rle(basecall)
        self.compact_basecall = ''.join(self.rle_conversion['value'])
        self.homop_length = self.rle_conversion['length']
        self.inverse = np.repeat(
            np.arange(0, len(self.rle_conversion)),
            self.rle_conversion['length'])

    def trimmed_compact(self, start, end):
        """Return a trimmed compressed sequence.

        :param start: start co-ordinate in uncompressed sequence.
        :param end: end co-ordinate in uncompressed sequence (exclusive).

        :returns: the trimmed compressed sequence.
        """
        s, e = self.transform_coords(start, end)
        return self.compact_basecall[s:e]

    def transform_coords(self, start, end):
        """Transform a slice co-ordinates from the input sequence.

        :param start: start co-ordinate in uncompressed sequence.
        :param end: end co-ordinate in uncompressed sequence (exclusive).

        :returns: the trimmed compressed sequence.
        """
        # Visual explanation:
        # seq         : AATTCCGG
        # seq_i       : 01234567
        # rle_i       : 00112233
        # search right: 11223344
        # => the slice [0:8] should map to [0:4]
        # => subtract 1 from s; e is fine because we want exclusive

        # older implementation
        # s, e = np.searchsorted(
        #     self.rle_conversion['start'], [start, end - 1], 'right')
        # return s - 1, e
        return self.inverse[start], self.inverse[end - 1] + 1

    def coord_compact_to_full(self, coord):
        """Map from encoded index to full basecall index.

        :param coord: array_like, compact coordinate(s) to be converted to
            full basecall

        :returns: numpy.ndarray of coordinates
        """
        return self.rle_conversion[coord]['start']


def parasail_alignment(query, ref):
    """Run a Smith-Waterman alignment between two sequences.

    :param query: the query sequence.
    :param ref: the reference sequence.

    :returns: reference start co-ordinate, cigar string
    """
    result = parasail.sw_trace_striped_32(query, ref, 5, 3, parasail.dnafull)
    rstart, cigar = medaka.smolecule.parasail_to_sam(result, query)
    return rstart, cigar


def add_extra_clipping(cigar, start_clipped, end_clipped):
    """Add extra soft clipping to begining and end of cigar string.

    :param cigar: input cigar string
    :param start_clipped: soft clipping to be added at the start of the read
    :param end_clipped: soft clipping to be added at the end of the read

    :returns: modified cigar string
    """
    begin_cigar = re.compile(r"^(?P<len>\d+)(?P<op>\D+)")
    end_cigar = re.compile(r"(?P<len>\d+)(?P<op>\D+)$")

    # Add start and end clipping
    if start_clipped:
        m = re.search(begin_cigar, cigar)
        if m.group('op') == 'S':
            orig_clip = ''.join(m.groups())
            corrected_length = int(m.group('len')) + start_clipped
            head = '{}S'.format(corrected_length)
            body = cigar[len(orig_clip):]
        else:
            head = '{}S'.format(start_clipped)
            body = cigar
    else:
        head = ''
        body = cigar

    if end_clipped:
        m = re.search(end_cigar, body)
        if m.group('op') == 'S':
            orig_clip = ''.join(m.groups())
            corrected_length = int(m.group('len')) + end_clipped
            tail = '{}S'.format(corrected_length)
            body = body[:-len(orig_clip)]
        else:
            tail = str(end_clipped) + 'S'
    else:
        tail = ''

    return ''.join((head, body, tail))


def initialise_alignment(
        query_name, reference_id, reference_start,
        query_sequence, cigarstring, flag, mapping_quality=60,
        query_qualities=None, tags=None):
    """Create a `Pysam.AlignedSegment` object.

    :param query_name: name of the query sequence
    :param reference_id: index to the reference name
    :param reference_start: 0-based index of first leftmost reference
        coordinate
    :param query_sequence: read sequence bases, including those soft clipped
    :param cigarstring: cigar string representing the alignment of query
        and reference
    :param flag: bitwise flag representing some properties of the alignment
        (see SAM format)
    :param mapping_quality: optional quality of the mapping or query to
        reference
    :param query_qualities: optional base qualities of the query, including
        soft-clipped ones!

    :returns: `pysam.AlignedSegment` object
    """
    if tags is None:
        tags = dict()

    a = pysam.AlignedSegment()
    a.query_name = query_name
    a.reference_id = reference_id
    a.reference_start = reference_start
    a.query_sequence = query_sequence
    a.cigarstring = cigarstring
    a.flag = flag
    a.mapping_quality = mapping_quality
    if query_qualities is not None:
        a.query_qualities = query_qualities

    for tag_name, tag_value in tags.items():
        a.set_tag(tag_name, tag_value)

    return a


def get_rl_params(read_name, read_fast5):
    """Get shape and scale parameters from fast5 for read."""
    data_path = (
        'read_{}/Analyses/Basecall_1D_000/'
        'BaseCalled_template/RunlengthBasecall')

    with h5py.File(read_fast5) as h:
        data = h[data_path.format(read_name)][()]

    call = ''.join(x.decode() for x in data['base'])
    wl = array.array('f', data['shape'])
    wk = array.array('f', data['scale'])

    return call, wl, wk


def _compress_alignment(alignment, ref_rle, fast5_dir=None, file_index=None):
    logger = medaka.common.get_named_logger('Compress_bam')

    if alignment.is_unmapped or alignment.is_secondary:
        msg = 'Alignment of read {} is unmapped or secondary. Skip.'
        logger.info(msg.format(alignment.query_name))
        return None

    query_rle = RLEConverter(alignment.query_sequence)

    # Get aligned query in RLE
    qstart, qend = (
        alignment.query_alignment_start, alignment.query_alignment_end)
    q_compact_start, q_compact_end = query_rle.transform_coords(qstart, qend)
    compact_query = query_rle.compact_basecall[q_compact_start:q_compact_end]

    # Get aligned reference in RLE
    rstart, rend = alignment.reference_start, alignment.reference_end
    r_compact_start, r_compact_end = ref_rle.transform_coords(rstart, rend)
    compact_ref = ref_rle.compact_basecall[r_compact_start:r_compact_end]

    # Calculate new alignment with compressed reads
    rstart, cigar = parasail_alignment(compact_query, compact_ref)

    # Cigar is now with respect to the trimmed compact_query, we need to
    # return it to the full compact query
    extra_start_clip = q_compact_start
    extra_end_clip = len(query_rle.compact_basecall) - q_compact_end
    corrected_cigar = add_extra_clipping(
        cigar, extra_start_clip, extra_end_clip)
    rstart += r_compact_start

    # If file_index is present, retrive RLE parameters from fast5 files.
    if file_index:
        try:
            fast5_fname = file_index[alignment.query_name]
        except KeyError:
            logger.warning('Not found in summary file: {}'.format(
                alignment.query_name))
            return None

        # Get full path from a root path. There should be only ONE file found
        # with that name.
        fast5_path = list(Path(fast5_dir).glob('**/{}'.format(fast5_fname)))
        if len(fast5_path) > 1:
            logger.warning(
                'Found several fast5 with name {}'.format(fast5_fname))
            return None
        fast5_path = fast5_path[0]

        fast5_call, wl, wk = get_rl_params(alignment.query_name, fast5_path)

        # params needs flipping for reverse alignments
        if alignment.is_reverse:
            wl = wl[::-1]
            wk = wk[::-1]
            fast5_call = medaka.common.reverse_complement(fast5_call)

        # ocasional errors where bases are repeated
        # may lead to discrepancies, discard these
        if fast5_call != query_rle.compact_basecall:
            logger.warning(
                'RLE table within fast5 file is inconsistent with '
                'compressed basecall for read {}'.format(alignment.query_name))
            return None
        tags = {'WL': wl, 'WK': wk}
    else:
        tags = dict()

    # Create alignment object
    a = initialise_alignment(
        alignment.query_name, alignment.reference_id, rstart,
        query_rle.compact_basecall, corrected_cigar, alignment.flag, tags=tags,
        query_qualities=array.array(
            'B', np.minimum(query_rle.homop_length, 255)))

    return a


def _compress_bam(bam_input, bam_output, ref_fname,
                  fast5_dir=None, summary=None, regions=None, threads=1):
    """Compress a bam into run length encoding (RLE).

    :param bam_input: str, name of the bam file to be compressed
    :param bam_output: str, name of the bam to be produced
    :param ref_fname: str, reference filename, used to produce bam_input
    :param fast5_dir: str, root directory to find the fast5 files
    :param summary: str, filename of a summary name, with the columns
        to link the read id to the fast5 file containing it. Must contain
        columns 'read_id' and 'filename'
    :param regions: list, genomic regions to be extracted
    :param threads: int, number of workers to be used

    :returns: None
    """
    regions = medaka.common.get_regions(bam_input, regions)
    ref_fasta = pysam.FastaFile(ref_fname)

    # If fast_dir is passed, create an index
    if fast5_dir:
        with open(summary) as input_file:
            # Summary files can be huge, avoid loading to memory
            col_names = input_file.readline().replace('#', '').split()
            col_filename = col_names.index('filename')
            col_readid = col_names.index('read_id')

            file_index = dict(
                (line.split()[col_readid], line.split()[col_filename])
                for line in input_file)
    else:
        file_index = None

    with pysam.AlignmentFile(bam_input, 'r') as alignments_bam:
        tmp_output = '{}.tmp'.format(bam_output)
        with pysam.AlignmentFile(
                tmp_output, 'wb', header=alignments_bam.header) as output:
            for region in regions:
                bam_current = alignments_bam.fetch(
                    reference=region.ref_name,
                    start=region.start,
                    end=region.end)
                ref_sequence = ref_fasta.fetch(region.ref_name)
                ref_rle = RLEConverter(ref_sequence)
                func = functools.partial(
                    _compress_alignment,
                    ref_rle=ref_rle,
                    fast5_dir=fast5_dir,
                    file_index=file_index)
                with concurrent.futures.ThreadPoolExecutor(
                        max_workers=threads) as executor:
                    for chunk in medaka.common.grouper(bam_current, 100):
                        for new_alignment in executor.map(func, chunk):
                            if new_alignment is not None:
                                output.write(new_alignment)

        pysam.sort("-o", bam_output, tmp_output)
        os.remove(tmp_output)
        pysam.index(bam_output)


def compress_seq(read):
    """Compress homopolymers within a basecall.

    Creates and RLE sequence record, with run-lengths stored in qualities.

    :param read: `pysam FastxRecord` object.

    :returns: `pysam FastxRecord`
    """
    logger = medaka.common.get_named_logger('Compress_basecalls')

    # Phred qscores `!"#$%&`...
    scores = ''.join(chr(x) for x in range(33, 127))

    rle_compressed = RLEConverter(read.sequence)

    # we can only encode up to a homopolymer length 93
    inds = np.where(rle_compressed.homop_length >= len(scores))[0]
    if len(inds) > 0:
        logger.warning(
            "Some homopolymers in {} are longer than the longest "
            "supported length\n".format(read.name))
        rle_compressed.homop_length[inds] = len(scores) - 1

    coded_lengths = ''.join([scores[x] for x in rle_compressed.homop_length])

    compressed_record = pysam.FastxRecord(
        name=read.name,
        comment=read.comment if read.comment is not None else '',
        sequence=rle_compressed.compact_basecall,
        quality=coded_lengths)

    return compressed_record


def compress_bam(args):
    """Compress a bam alignment file into an RLE system of reference."""
    if args.use_fast5_info:
        fast5_dir, summary = args.use_fast5_info
    else:
        fast5_dir = summary = None

    _compress_bam(
        args.bam_input, args.bam_output, args.ref_fname,
        fast5_dir=fast5_dir, summary=summary,
        threads=args.threads, regions=args.regions)
