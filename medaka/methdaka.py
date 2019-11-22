"""Functionality for aggregating single read methylation calls."""

from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import functools
import sys

import mappy
from ont_fast5_api.analysis_tools.basecall_1d import Basecall1DTools
from ont_fast5_api.conversion_tools.conversion_utils import get_fast5_file_list
from ont_fast5_api.fast5_interface import get_fast5_file
import pysam

import medaka.common


MODBASEPATH = 'BaseCalled_template/ModBaseProbs'
BASES = 'A 6mA C 5mC G T'.split()
MODTYPE = [(b, 'uint8') for b in BASES]
MOTIFS = {
    'dcm': {
        # seq: [fwd offset, rev offset]
        # NB: code will assume rev offset is larger!
        'motifs': {
            'CCAGG': (1, 3),
            'CCTGG': (1, 3)},
        'tag': 'MC'},
    'dam': {
        'motifs': {
            'GATC': (1, 2)},
        'tag': 'MA'},
    'cpg': {
        'motifs': {
            'CG': (0, 1)},
        'tag': 'MC'}
    }


def hdf_to_sam_worker(reference, fname):
    """Extract and align basecall and methylation data from `.fast5`.

    :param reference: `.fasta` file containing reference sequence(s).
    :param fname: `.'fast5` file containing read data.
    """
    logger = medaka.common.get_named_logger('ModExtract')
    logger.info("Processing {}.".format(fname))
    results = list()
    aligner = mappy.Aligner(reference, preset='map-ont')
    with get_fast5_file(fname, mode="r") as f5:
        reads = list(f5.get_read_ids())
        logger.info("Found {} reads for {}.".format(len(reads), fname))
        for read_id in reads:
            read = f5.get_read(read_id)
            tool = Basecall1DTools(read)
            name, sequence, qstring = tool.get_called_sequence(
                'template', fastq=False)
            try:
                align = next(aligner.map(sequence, MD=True, cs=True))
            except StopIteration:
                continue
            else:
                if align.strand == +1:
                    flag = '0'
                    seq = sequence
                else:
                    flag = '16'
                    seq = medaka.common.reverse_complement(sequence)
                rname = align.ctg
                pos = str(align.r_st + 1)
                mapq = str(align.mapq)
                clip = [
                    '' if x == 0 else '{}S'.format(x)
                    for x in (align.q_st, len(sequence) - align.q_en)]
                if align.strand == -1:
                    clip = clip[::-1]
                cigar = clip[0] + align.cigar_str + clip[1]
                NM = 'NM:i:' + str(align.NM)

            latest = read.get_latest_analysis('Basecall_1D')
            mod_base = read.get_analysis_dataset(latest, MODBASEPATH)
            mod_base = mod_base.view(dtype=MODTYPE)
            mA = 'MA:B:C,{}'.format(','.join(
                str(x) for x in mod_base['6mA'].reshape(-1)))
            mC = 'MC:B:C,{}'.format(','.join(
                str(x) for x in mod_base['5mC'].reshape(-1)))

            results.append('\t'.join((
                read_id, flag, rname, pos, mapq, cigar,
                '*', '0', '0', seq, qstring, NM, mA, mC)))
    return results


def hdf_to_sam(args):
    """Entry point for converting guppy methylcalled fast5s to sam."""
    sys.stdout.write('\t'.join(('@HD', 'VN:1.5', 'SO:unsorted')))
    sys.stdout.write('\n')
    for name, seq, _ in mappy.fastx_read(
            args.reference, read_comment=False):
        sys.stdout.write('@SQ\tSN:{}\tLN:{}\n'.format(name, len(seq)))

    fast5s = get_fast5_file_list(args.path, recursive=args.recursive)
    worker = functools.partial(hdf_to_sam_worker, args.reference)
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for res in executor.map(worker, fast5s):
            for r in res:
                sys.stdout.write('{}\n'.format(r))


def call_methylation(args):
    """Entry point for calling methylation from bam file."""
    logger = medaka.common.get_named_logger('Mextract')
    logger.info("Processing {}.".format(args.bam))
    region = medaka.common.get_regions(args.bam, [args.region])[0]
    ref = pysam.FastaFile(args.reference)
    ref = ref.fetch(region.ref_name)
    motifs = MOTIFS[args.meth]['motifs']
    tag_name = MOTIFS[args.meth]['tag']

    def _find_loci(ref, motif, rel_pos, start, end):
        index = start
        while index < end:
            index = ref.find(motif, index)
            if index == -1:
                break
            yield index, index + rel_pos[0], index + rel_pos[1]
            index += 1

    with open(args.output, 'w') as fh:
        with pysam.AlignmentFile(args.bam, "rb") as bam:
            for motif, rel_pos in motifs.items():
                logger.info("Searching for motif '{}'.".format(motif))
                loci = _find_loci(
                    ref, motif, rel_pos, region.start, region.end)
                locus = next(loci)

                # we assume the reverse strand position is after
                # the forward strand
                is_rev = False
                next_pos = locus[1]
                meth = Counter()
                not_meth = Counter()
                for pileupcolumn in bam.pileup(
                        region.ref_name, region.start, region.end):
                    if pileupcolumn.pos < next_pos:
                        continue
                    for read in pileupcolumn.pileups:
                        aln = read.alignment
                        if read.is_del \
                                or read.is_refskip \
                                or aln.is_reverse != is_rev:
                            continue
                        qpos = read.query_position
                        tag = read.alignment.get_tag(tag_name)
                        if not aln.is_reverse:
                            mscore = tag[qpos]
                        else:
                            mscore = tag[len(tag) - qpos - 1]

                        if mscore > args.filter[1]:
                            meth[aln.is_reverse] += 1
                        elif mscore < args.filter[0]:
                            not_meth[aln.is_reverse] += 1
                    if is_rev:
                        fh.write(
                            '\t'.join(str(x) for x in (
                                region.ref_name, locus[0], motif,
                                meth[False], meth[True],
                                not_meth[False], not_meth[True])))
                        fh.write('\n')
                        try:
                            locus = next(loci)
                        except StopIteration:
                            break
                        else:
                            is_rev = False
                            next_pos = locus[1]
                            meth = Counter()
                            not_meth = Counter()
                    else:
                        is_rev = True
                        next_pos = locus[2]
