"""
MicrobiomeDash — Rare ASV removal: compute stats and filter by prevalence/abundance.
"""
import io

import numpy as np
from biom import Table, load_table


def _biom_to_bytes(table: Table) -> bytes:
    """Serialize a biom.Table to HDF5 bytes."""
    import h5py

    buf = io.BytesIO()
    with h5py.File(buf, "w") as f:
        table.to_hdf5(f, generated_by="MicrobiomeDash-RareASV")
    return buf.getvalue()


def compute_asv_stats(biom_path: str) -> dict:
    """Load a BIOM table and return summary stats + per-ASV prevalence/abundance.

    Returns dict with keys:
        n_asvs, n_samples, total_reads,
        per_asv: list of {id, prevalence_pct, total_abundance}
    """
    table = load_table(biom_path)
    obs_ids = table.ids(axis="observation")
    sample_ids = table.ids(axis="sample")
    n_samples = len(sample_ids)
    data = table.matrix_data.toarray()  # obs x sample

    per_asv = []
    for i, oid in enumerate(obs_ids):
        row = data[i]
        n_present = int((row > 0).sum())
        prevalence_pct = round(100.0 * n_present / n_samples, 2) if n_samples > 0 else 0.0
        total_abundance = int(row.sum())
        per_asv.append({
            "id": str(oid),
            "prevalence_pct": prevalence_pct,
            "total_abundance": total_abundance,
        })

    return {
        "n_asvs": len(obs_ids),
        "n_samples": n_samples,
        "total_reads": int(data.sum()),
        "per_asv": per_asv,
    }


def filter_rare_asvs(
    biom_path: str, min_prevalence: float, min_abundance: int
) -> dict:
    """Filter ASVs below prevalence and/or abundance thresholds.

    Args:
        biom_path: Path to input BIOM file.
        min_prevalence: Minimum % of samples an ASV must appear in (0-100).
        min_abundance: Minimum total read count across all samples.

    Returns dict with keys:
        biom_bytes, n_input, n_kept, n_removed, reads_kept, reads_removed
    """
    table = load_table(biom_path)
    obs_ids = table.ids(axis="observation")
    n_samples = len(table.ids(axis="sample"))
    data = table.matrix_data.toarray()  # obs x sample

    keep_ids = []
    for i, oid in enumerate(obs_ids):
        row = data[i]
        n_present = int((row > 0).sum())
        prevalence_pct = 100.0 * n_present / n_samples if n_samples > 0 else 0.0
        total_abundance = int(row.sum())

        if prevalence_pct >= min_prevalence and total_abundance >= min_abundance:
            keep_ids.append(oid)

    if not keep_ids:
        raise ValueError(
            "All ASVs would be removed with these thresholds. "
            "Try lowering the prevalence or abundance cutoff."
        )

    filtered = table.filter(keep_ids, axis="observation", inplace=False)

    n_removed = len(obs_ids) - len(keep_ids)
    reads_kept = int(filtered.matrix_data.toarray().sum())
    reads_total = int(data.sum())

    return {
        "biom_bytes": _biom_to_bytes(filtered),
        "n_input": len(obs_ids),
        "n_kept": len(keep_ids),
        "n_removed": n_removed,
        "reads_kept": reads_kept,
        "reads_removed": reads_total - reads_kept,
    }
