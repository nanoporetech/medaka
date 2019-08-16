import array
import concurrent.futures
import functools
import itertools
from multiprocessing import Pool
import re
import sys
from timeit import default_timer as now

import numpy as np
import parasail
import pysam

import medaka.common
import medaka.smolecule


def rle(iterable, low_mem=False):
    """Calculate a run length encoding (rle), of an input iterable.

    :param iterable: input iterable.
    :param low_mem: use a lower memory implementation

    :returns: structured array with fields `start`, `length`, and `value`.
    """

    if not isinstance(iterable, np.ndarray):
        array = np.fromiter(iterable, dtype='U1', count=len(iterable))
    else:
        array = iterable

    if len(array.shape) != 1:
        raise TypeError("Input array must be one dimensional.")
    dtype = [('length', int), ('start', int), ('value', array.dtype)]

    if not low_mem:
        pos = np.where(array[:-1] != array[1:])[0]
        pos = np.concatenate(([0], pos+1, [len(array)]))

        return np.fromiter((
            (end - start, start, array[start])
            for (start, end) in zip(pos[:-1], pos[1:])),
            dtype, count=len(pos) - 1,)
    else:
        def _gen():
            start = 0
            for key, group in itertools.groupby(array):
                length = sum(1 for x in group)
                yield length, start, key
                start += length

        return np.fromiter(_gen(), dtype=dtype)


class RLEConverter(object):
    def __init__(self, basecall):
        """Class to convert a basecall to RLE, with coordinate conversions.

        :param basecall: string to be converted to RLE.
        """
        self.basecall = basecall
        self.rle_conversion = rle(basecall, low_mem=True)
        self.compact_basecall = ''.join(self.rle_conversion['value'])
        self.homop_length = self.rle_conversion['length']

    def coord_full_to_compact(self, coord):
        """Map string index from basecall to encoded string.

        :param coord: array_like, full coordinate(s) position to be converted
            to compact

        :returns: numpy.ndarray of coordinates
        """

        return np.searchsorted(
            self.rle_conversion['start'], coord, 'right') - 1

    def coord_compact_to_full(self, coord):
        """Map from encoded index to full basecall index.

        :param coord: array_like, compact coordinate(s) to be converted to
            full basecall

        :returns: numpy.ndarray of coordinates
        """

        return self.rle_conversion[coord]['start']


def parasail_alignment(query, ref):
    result = parasail.sw_trace_striped_16(query, ref, 5, 3, parasail.dnafull)
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
        query_qualities=None):
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

    return a


def _compress_alignment(alignment, ref_rle):
    logger = medaka.common.get_named_logger('Compress_bam')
    logger.info('Processing: {}'.format(alignment.query_name))

    if alignment.is_unmapped or alignment.is_secondary:
        msg = 'Alignment of read {} is unmapped or secondary. Skip.'
        logger.info(msg.format(alignment.query_name))
        return None

    query_rle = RLEConverter(alignment.query_sequence)

    # Get aligned query in RLE
    qstart = alignment.query_alignment_start
    qend = alignment.query_alignment_end
    q_compact_start = query_rle.coord_full_to_compact(qstart)
    q_compact_end = query_rle.coord_full_to_compact(qend)
    compact_query = query_rle.compact_basecall[q_compact_start:q_compact_end]

    # Get aligned reference in RLE
    rstart, rend = alignment.reference_start, alignment.reference_end
    r_compact_start = ref_rle.coord_full_to_compact(rstart)
    r_compact_end = ref_rle.coord_full_to_compact(rend)
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

    # Create alignment object
    a = initialise_alignment(
        alignment.query_name, alignment.reference_id,
        rstart, query_rle.compact_basecall,
        corrected_cigar, alignment.flag,
        query_qualities=array.array('B', list(query_rle.homop_length)))

    return a


def _compress_bam(bam_input, bam_output, ref_fname, regions=None, threads=1):
    """Compress a bam into run length encoding (RLE)

    :param bam_input: str, name of the bam file to be compressed
    :param bam_output: str, name of the bam to be produced
    :param ref_fname: str, reference filename, used to produce bam_input
    :param regions: list, genomic regions to be extracted
    :param threads: int, number of workers to be used

    :returns: None
    """
    regions = medaka.common.get_regions(bam_input, regions)
    ref_fasta = pysam.FastaFile(ref_fname)

    with pysam.AlignmentFile(bam_input, 'r') as alignments_bam:
        with pysam.AlignmentFile(
                bam_output, 'wb', header=alignments_bam.header) as output:
            for region in regions:
                bam_current = alignments_bam.fetch(
                    reference=region.ref_name,
                    start=region.start,
                    end=region.end)
                ref_sequence = ref_fasta.fetch(region.ref_name)
                ref_rle = RLEConverter(ref_sequence)
                func = functools.partial(
                    _compress_alignment,
                    ref_rle=ref_rle)
                with concurrent.futures.ThreadPoolExecutor(
                        max_workers=threads) as executor:
                    for chunk in medaka.common.grouper(bam_current, 100):
                        for new_alignment in executor.map(func, chunk):
                            if new_alignment is not None:
                                output.write(new_alignment)


def compress_seq(read):
    """Compress homopolymers within a basecall, encoding their lengths
    in qscores.

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
        comment=read.comment,
        sequence=rle_compressed.compact_basecall,
        quality=coded_lengths)

    return compressed_record


def compress_basecalls(args):
    logger = medaka.common.get_named_logger('Compress_basecalls')

    reads = pysam.FastxFile(args.input)
    if args.threads > 1:
        pool = Pool(args.threads)
        compressed = pool.imap(compress_seq, reads)
    else:
        compressed = (compress_seq(r) for r in reads)

    t0 = now()
    if args.output is None:
        fh = sys.stdout
    else:
        fh = open(args.output, 'w')

    for read in compressed:
        fh.write('@{}\n{}\n{}\n'.format(
            read.name, read.comment, read.sequence))
        fh.write('{}\n{}\n'.format('+', read.quality))
    t1 = now()
    logger.info('Compressing {} took {:.3f}s.'.format(args.input, t1 - t0))

    if args.output is not None:
        fh.close()


def compress_bam(args):
    """Compress a bam alignment file into an RLE system of reference """
    _compress_bam(
        args.bam_input, args.bam_output, args.ref_fname,
        threads=args.threads, regions=args.regions)
