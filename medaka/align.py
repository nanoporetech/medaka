"""Alignment and CIGAR string processing."""
import re

import edlib
import pysam

from medaka import parasail

re_split_cigar = re.compile(r"(?P<len>\d+)(?P<op>\D+)")
re_split_cigar_rev = re.compile(r"(?P<op>\D+)(?P<len>\d+)")


def cigar_ops_from_start(cigar):
    """Yield cigar string operations from start of cigar (in order).

    :param cigar: cigar string.

    :yields: str, op. length, str op. type
    """
    for m in re.finditer(re_split_cigar, cigar):
        yield m.group('len'), m.group('op')


def cigar_ops_from_end(cigar):
    """Yield cigar string operations from end of cigar (in reverse order).

    :param cigar: cigar string.

    :yields: str, op. length, str op. type
    """
    for m in re.finditer(re_split_cigar_rev, cigar[::-1]):
        yield m.group('len')[::-1], m.group('op')[::-1]


def trim_cigar(cigar, start=True):
    """Trim cigar string to start or end on a match.

    :param cigar: cigar string.
    :param start: if True, trim start, else trim end.

    :returns: cigar, number bases to trim from query, rstart offset
    """
    trim_cigar, rstart_offset, q_trim = 0, 0, 0
    op_gen = cigar_ops_from_start if start else cigar_ops_from_end
    for n, op in op_gen(cigar):
        if op == '=':
            break
        trim_cigar += len(n) + len(op)
        if op in {'I', 'X'}:
            # insertion or mismatch at start/end of aln, trim cigar and query
            q_trim += int(n)
            # only offset rstart for deletion at start
            rstart_offset += int(n) if (op == 'X' and start) else 0
        elif op == 'D':  # query does not need trimming
            # deletion of ref base at start/end of alignment (global alignment)
            rstart_offset += int(n) if start else 0
        else:
            msg = 'Encountered unsupported cigar operation: {}'
            raise ValueError(msg.format(op))
    cigar = cigar[trim_cigar:] if start else cigar[:len(cigar) - trim_cigar]
    return cigar, q_trim, rstart_offset


def parasail_to_sam(result, seq):
    """Extract reference start and sam compatible cigar string.

    :param result: parasail alignment result.
    :param seq: query sequence.

    :returns: reference start coordinate, cigar string.

    """
    cigstr = result.cigar.decode.decode()

    first = next(cigar_ops_from_start(cigstr))
    prefix = ''.join(first)
    rstart = result.cigar.beg_ref
    cliplen = result.cigar.beg_query
    clip = '' if cliplen == 0 else '{}S'.format(cliplen)
    if first[1] == 'I':
        pre = '{}S'.format(int(first[0]) + cliplen)
    elif first[1] == 'D':
        pre = clip
        rstart = int(first[0])
    else:
        pre = '{}{}'.format(clip, prefix)

    mid = cigstr[len(prefix):]
    end_clip = len(seq) - result.end_query - 1
    suf = '{}S'.format(end_clip) if end_clip > 0 else ''
    new_cigstr = ''.join((pre, mid, suf))
    return rstart, new_cigstr


def parasail_alignment(query, ref):
    """Run a Smith-Waterman alignment between two sequences.

    :param query: the query sequence.
    :param ref: the reference sequence.

    :returns: reference start co-ordinate, cigar string
    """
    result = parasail.sw_trace_striped_32(query, ref, 5, 3, parasail.dnafull)
    rstart, cigar = parasail_to_sam(result, query)
    return rstart, cigar


def add_extra_clipping(cigar, start_clipped, end_clipped):
    """Add extra soft clipping to begining and end of cigar string.

    :param cigar: input cigar string
    :param start_clipped: soft clipping to be added at the start of the read
    :param end_clipped: soft clipping to be added at the end of the read

    :returns: modified cigar string
    """
    # Add start and end clipping
    if start_clipped:
        oplen, op = next(cigar_ops_from_start(cigar))
        if op == 'S':
            orig_clip = '{}{}'.format(oplen, op)
            corrected_length = int(oplen) + start_clipped
            head = '{}S'.format(corrected_length)
            body = cigar[len(orig_clip):]
        else:
            head = '{}S'.format(start_clipped)
            body = cigar
    else:
        head = ''
        body = cigar

    if end_clipped:
        oplen, op = next(cigar_ops_from_end(cigar))
        if op == 'S':
            orig_clip = '{}{}'.format(oplen, op)
            corrected_length = int(oplen) + end_clipped
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
        query_qualities=None, tags=None, header=None):
    """Create a `pysam.AlignedSegment` object.

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
    :param header: optional `pysam.AlignmentHeader` object, enabling use of the
        reference_name attr of the returned `pysam.AlignedSegment` obj.

    :returns: `pysam.AlignedSegment` object
    """
    if tags is None:
        tags = dict()
    if header is None:
        a = pysam.AlignedSegment()
    else:
        a = pysam.AlignedSegment(header)
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


def chunked_edlib_align(qseq, rseq, contig_name, chunk_size=100000,
                        pad=10000, mode='NW', header=None):
    """Align query to reference in chunks using edlib.

    .. note::

        * Chunks are aligned sequentially and consecutive alignments overlap by
          1 match (= cigar op) to facilitate decoding of variants.
        * First chunk uses infix (HW) option to allow start of the first.
          query chunk to be offset from the start of the first ref chunk.
        * Subsequent chunks use the prefix(SHW) for alignment extension.
        * In both cases edlib forces all query bases to align.
        * This method assumes query can be aligned to ref in a single
          alignment, and thus the absense of any large structural variants.
        * Chunking can affect anchoring context around long indels.
        * Increase chunk_size to avoid spurious alignment of long insertions.
        * Increase pad to avoid spurious alignment of long deletions.
        * It is assumed that query contigs have the same names as ref contigs.

    :param qseq: str, query sequence.
    :param rseq: str, reference sequence.
    :param chunk_size: int, length of consensus chunks. Runtime scales
        quadratically with chunksize.
    :param pad: int, reference chunks are `chunk_size` + `pad` bases long.
    :param mode: str, alignment mode
        HW - chunked version of edlib HW mode (global in query, local in ref).
        NW - chunked version of edlib NW mode (global in query and ref).
        HWT - same as HW, but first/last alignment starts/ends on a match.
    :param header: `pysam.AlignmentHeader` obj, if not given the reference_name
        attr of yielded `pysam.AlignedSegment` obj will not be available.

    :yields: `pysam.AlignedSegment` objs

    """
    # TODO this could break down if we have variants longer than the chunksize
    # such that a chunk could not contain two matches that we trim back to as
    # start and end. We could implement something to detect such cases and
    # increase the chunk size until we span the variant - however, is this
    # really needed for a small-variant caller?

    ends_modes = {  # edlib modes for first and last chunks
        'HW': ('HW', 'SHW'),
        'NW': ('SHW', 'NW'),
        'HWT': ('HW', 'SHW')  # same as 'HW', but we trim query.
    }
    if mode not in ends_modes:
        msg = 'Unrecognised mode option {}, use one of {}'
        raise KeyError(msg.format(mode, set(ends_modes.keys())))
    mode_first, mode_last = ends_modes[mode]

    def check_cigar_starts_with_match(cigar):
        l, op = next(cigar_ops_from_start(cigar))
        if op != '=':
            msg = 'Alignment did not start with a match: {}{}'
            raise ValueError(msg.format(l, op))

    if header is not None:
        ref_id = header.references.index(contig_name)
    else:
        ref_id = 0

    flag = 0  # not rev_comp

    qend_last = 0
    qend = 0
    aln_last = None
    trim_qend = 0

    while qend + trim_qend < len(qseq):
        qstart = max(0, qend_last - 1)  # overlap by 1 match
        qend = min(qend_last + chunk_size, len(qseq))
        is_last_chunk = qend == len(qseq)
        if qstart == 0:
            rstart = 0
            rend = min(len(rseq), qend + pad)
            if is_last_chunk and mode == 'NW':
                # unchunked global alignment
                result = edlib.align(qseq, rseq, task='path', mode='NW')
                rstart_aln, cigar = result['locations'][0][0], result['cigar']
            else:
                result = edlib.align(qseq[qstart:qend], rseq[rstart:rend],
                                     task='path', mode=mode_first)
                rstart_aln, cigar = result['locations'][0][0], result['cigar']
                if mode == 'HWT':  # trim start of alignment
                    cigar, trim_qstart, r_offset = trim_cigar(cigar,
                                                              start=True)
                    qstart += trim_qstart
                    rstart_aln += r_offset
                # don't trim end if this is the last chunk unless mode is match
                if not is_last_chunk or mode == 'HWT':
                    cigar, trim_qend, _ = trim_cigar(cigar, start=False)
                    qend -= trim_qend
        elif is_last_chunk:
            rstart = aln_last.reference_end - 1  # overlap by 1 match
            rend = len(rseq)
            # For last alignment, allow gaps only at end of qseq
            result = edlib.align(qseq[qstart:qend], rseq[rstart:rend],
                                 task='path', mode=mode_last)
            rstart_aln, cigar = result['locations'][0][0], result['cigar']
            check_cigar_starts_with_match(cigar)
            # don't trim end of alignment for last chunk unless mode is match
            if mode == 'HWT':
                cigar, trim_qend, _ = trim_cigar(cigar, start=False)
                qend -= trim_qend
        else:
            rstart = aln_last.reference_end - 1  # overlap by 1 match
            rend = min(len(rseq), rstart + chunk_size)
            # For subsequent alignments, allow gaps only at end of qseq
            result = edlib.align(qseq[qstart:qend], rseq[rstart:rend],
                                 task='path', mode='SHW')
            rstart_aln, cigar = result['locations'][0][0], result['cigar']
            check_cigar_starts_with_match(cigar)

            cigar, trim_qend, _ = trim_cigar(cigar, start=False)
            qend -= trim_qend

        qname = '{}_{}_{}'.format(contig_name, qstart, qend)
        # use edit distance to set NM tag, since alignments with no edits
        # are skipped when decoding variants. Note that the edit distnce of
        # the trimmed alignment could be less than this value.
        nm = result['editDistance']
        aln = initialise_alignment(qname, ref_id, rstart + rstart_aln,
                                   qseq[qstart: qend], cigar, flag,
                                   header=header, tags=dict(NM=nm))

        yield aln

        qend_last, aln_last = qend, aln
