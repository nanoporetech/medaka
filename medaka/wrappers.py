"""Convenience pipelines and wrappers for command line tools."""
import os
import subprocess
from timeit import default_timer as now

import medaka.common


def racon(reads_fx, scaf_fasta, paf_out, fasta_out, threads=4):
    """Polish scaffold using racon.

    Assumes minimap2, bgzip and racon are in the path.

    :param reads_fx: str, input reads filepath.
    :param scaf_fasta: str, input scaffold fasta filepath.
    :param paf_out: str, output paf filepath.
    :param fasta_out: str, output racon consensus filepath.
    :param threads: int, number of threads to use for minimap2 and racon.
    """
    threads = str(threads)
    logger = medaka.common.get_named_logger('DRAFT_NORM')

    def run(cmd, out_fp):
        with open(out_fp, 'wb') as fh:
            p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(['bgzip', '-c'], stdin=p1.stdout, stdout=fh)
            p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
            _ = p2.communicate()[0]

    aln_cmd = ['minimap2', '-x', 'map-ont', '-t', threads, scaf_fasta,
               reads_fx]
    t0 = now()
    logger.info('Aligning reads.')
    logger.info(' '.join(aln_cmd))
    run(aln_cmd, paf_out)

    racon_cmd = ['racon', '--include-unpolished', '--no-trimming', '-q', '-1',
                 '-t', threads, reads_fx, paf_out, scaf_fasta]

    t1 = now()
    logger.info('Creating normalised draft by running racon.')
    logger.info(' '.join(racon_cmd))
    run(racon_cmd, fasta_out)
    t2 = now()
    msg = 'Racon consensus took {:.2f}s ({:.2f} for PAF, {:.2f} for Racon).'
    logger.info(msg.format(t2 - t0, t1 - t0, t2 - t1))
    logger.info('Racon consensus written to {}'.format(fasta_out))


def minimap2(query, ref, out_bam, preset='map-ont', extra_args=None,
             threads=4):
    """Align query to reference with minimap2 and sort and index bam.

    :param query: str, query filepath.
    :param ref: str, reference filepath.
    :param out_bam: str, output filepath.
    :param preset: str, minimap2 -x option.
    :param extra_args: iterable of str, any other minimap2 options which are
        applied after preset.
    :param threads: int, number of threads to use for alignment and sorting.
    """
    threads = str(threads)
    extra_args = [] if extra_args is None else list(extra_args)

    logger = medaka.common.get_named_logger('ALIGN')
    logger.info('Aligning {} to {}'.format(query, ref))
    aln_cmd = ['minimap2', ref, query, '-x', preset, '-t', threads, '-a',
               '--secondary=no', '--MD', '-L'] + extra_args
    sort_cmd = ['samtools', 'sort', '--output-fmt', 'BAM', '-o', out_bam,
                '-@', threads]
    logger.info(' '.join(aln_cmd + ['|'] + sort_cmd))
    p1 = subprocess.Popen(aln_cmd, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(sort_cmd, stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
    _ = p2.communicate()[0]
    logger.info('Indexing bam')
    index_cmd = ['samtools', 'index', out_bam, '-@', threads]
    logger.info(' '.join(index_cmd))
    p3 = subprocess.Popen(index_cmd, stdout=subprocess.PIPE)
    _ = p3.communicate()[0]
    logger.info('BAM written to {}'.format(out_bam))


def haploid_variant(args):
    """Command line haploid variant calling tool."""
    logger = medaka.common.get_named_logger('HAPLOID_VARIANT')
    medaka.common.mkdir_p(args.output_dir, info='May using existing results.')

    paf_fp = os.path.join(args.output_dir, 'reads_to_ref.paf.gz')
    racon_cons_fp = os.path.join(args.output_dir, 'racon.fa.gz')
    if os.path.isfile(racon_cons_fp):
        logger.info('Racon consensus exists, skipping.')
    else:
        racon(args.reads_fastx, args.ref_fasta, paf_fp, racon_cons_fp,
              args.threads)
    reads_to_draft_bam_fp = os.path.join(args.output_dir, 'calls_to_draft.bam')
    if os.path.isfile(reads_to_draft_bam_fp):
        logger.info('Skipping alignment of reads to draft, bam exists.')
    else:
        logger.info('Aligning reads to draft.')
        # if args.regions is subset of args.ref_fasta contigs, align
        minimap2(args.reads_fastx, racon_cons_fp, reads_to_draft_bam_fp,
                 threads=args.threads)

    args.bam = reads_to_draft_bam_fp
    args.regions = None
    consensus_probs_fp = os.path.join(args.output_dir, 'consensus_probs.hdf')
    if os.path.isfile(consensus_probs_fp):
        logger.info('Skipping medaka consensus, consensus probs exist.')
    else:
        logger.info('Running medaka consensus.')
        args.output = consensus_probs_fp
        args.check_output = False
        args.save_features = False
        args.tag_name = None
        args.tag_value = None
        args.RG = None
        args.tag_keep_missing = True
        medaka.prediction.predict(args)

    consensus_fasta_fp = os.path.join(args.output_dir, 'consensus.fasta')
    if os.path.isfile(consensus_fasta_fp):
        logger.info('Skipping medaka stitch, consensus fasta exists.')
    else:
        logger.info('Running medaka stitch.')
        args.inputs = [consensus_probs_fp]
        args.output = consensus_fasta_fp
        args.draft = racon_cons_fp
        medaka.stitch.stitch(args)

    logger.info('Running medaka tools consensus2vcf.')
    args.bam = None  # name clash with reads to draft bam for consensus
    args.consensus = consensus_fasta_fp
    args.out_prefix = os.path.join(args.output_dir, 'consensus_to_ref')
    medaka.variant.vcf_from_fasta(args)
