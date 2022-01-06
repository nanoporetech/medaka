"""Convenience pipelines and wrappers for command line tools."""
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
    logger = medaka.common.get_named_logger('RunRacon')

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

    logger = medaka.common.get_named_logger('RunAlign')
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
