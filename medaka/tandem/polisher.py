"""Polisher module for polishing and stitching consensus sequences."""

import tempfile

import medaka.common
import medaka.medaka
import medaka.prediction


class Polisher:
    """Responsible for polishing and stitching consensus sequences.

    Uses subparsers to handle different sets of arguments for each step.
    """

    def __init__(self, *, input_fasta, bam, threads, output_dir, output_fasta):
        """Initializes the Polisher with the output directory.

        Args:
            input_fasta (str): Path to fasta file for seqeunces to be polished
            bam (str): Path to the BAM containing alignemnts
            to be used for poilshing
            threads (int): Number of threads
            output_dir (str): The directory where output files will be stored.
            output_fasta (str): Path to the output fasta file.

        """  # noqa: D401
        self.output_dir = output_dir
        self.input_fasta = input_fasta
        self.bam = bam
        self.threads = threads
        self.output_fasta = output_fasta
        self.logger = medaka.common.get_named_logger("Polisher")

    def polish(
        self,
        *,
        batch_size=100,
        chunk_len=2000,
        chunk_ovlp=500,
        model=None,
        save_features=False,
        bam_chunk=1000000,
        min_mapq=0,
        full_precision=True,
    ):
        """Run the polishing step on the consensus sequences.

        Args:
            batch_size (int): Batch size for polishing.
            chunk_len (int): Length of each chunk.
            chunk_ovlp (int): Overlap between chunks.
            model (str): Model name for polishing.
            save_features (bool): Whether to save features during polishing.
            bam_chunk (int): BAM chunk size.
            min_mapq (int): Minimum mapping quality.
            full_precision (bool): Use full precision during polishing.

        """
        self.logger.info("Running medaka consensus.")
        self.temporary_consensus = tempfile.NamedTemporaryFile(
            prefix="temp_consesus", dir=self.output_dir
        ).name

        parser = medaka.medaka.medaka_parser()
        args = [
            "inference",
            "--bam_workers",
            str(self.threads),
            "--threads",
            str(self.threads),
            "--batch_size",
            str(batch_size),
            "--chunk_len",
            str(chunk_len),
            "--chunk_ovlp",
            str(chunk_ovlp),
            "--bam_chunk",
            str(bam_chunk),
            "--min_mapq",
            str(min_mapq),
            "--cpu",
        ]
        if save_features:
            args.append("--save_features")
        if model:
            args.extend(["--model", model])
        if full_precision:
            args.append("--full_precision")

        args += [self.bam, self.temporary_consensus]
        args = parser.parse_args(args)
        self.logger.info("Starting polishing step.")
        medaka.prediction.predict(args)

        self._stitch()

    def _stitch(self):
        """Run the stitching step."""
        parser = medaka.medaka.medaka_parser()
        args = parser.parse_args(
            [
                "sequence",
                "--threads",
                "1",
                "--fill_char",
                None,
                "--min_depth",
                0,
                self.temporary_consensus,
                self.input_fasta,
                self.output_fasta,
            ]
        )

        self.logger.info("Starting stitching step.")
        medaka.stitch.stitch(args)
