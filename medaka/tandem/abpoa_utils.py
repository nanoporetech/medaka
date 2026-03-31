"""Shared helpers for running abPOA consensus."""

import medaka
import medaka.common
from medaka.tandem.record_name import RecordName


def run_abpoa_consensus(
    aligner,
    subreads: list,
    max_n_cons: int = 1,
    min_freq: float = 0.25,
    context_label: str = None,
) -> tuple[object, str]:
    """Run abPOA consensus with input/result validation and stderr capture."""
    if not subreads:
        raise RuntimeError("ABPOA consensus requires at least one subread.")

    records = [RecordName.from_str(s.name) for s in subreads]
    seqs = [
        s.seq if r.strand == "fwd" else medaka.common.reverse_complement(s.seq)
        for s, r in zip(subreads, records)
    ]
    empty_seq_names = [s.name for s, seq in zip(subreads, seqs) if not seq]
    if empty_seq_names:
        label = context_label or str(records[0].to_unpadded_region())
        raise RuntimeError(
            f"ABPOA input contains empty sequence(s) while processing "
            f"{label}: {len(empty_seq_names)} / {len(subreads)} reads are "
            "empty."
        )

    if medaka.wurlitzer is None:
        raise RuntimeError("wurlitzer is required for ABPOA consensus.")
    with medaka.wurlitzer.pipes(bufsize=0) as (_out, err):
        result = aligner.msa(
            seqs,
            out_cons=True,
            out_msa=False,
            max_n_cons=max_n_cons,
            min_freq=min_freq,
        )
    stderr = err.read()

    required_attrs = ("n_seq", "n_cons", "cons_seq", "clu_read_ids")
    is_valid_result = (
        result is not None
        and all(hasattr(result, attr) for attr in required_attrs)
        and result.n_seq == len(seqs)
        and 1 <= result.n_cons <= max_n_cons
        and len(result.cons_seq) == result.n_cons
    )
    if not is_valid_result:
        label = context_label or str(records[0].to_unpadded_region())
        raise RuntimeError(
            f"ABPOA encountered an error while processing {label}. "
            f"error: {stderr}"
        )

    return result, stderr
