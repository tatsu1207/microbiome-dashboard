"""
MicrobiomeDash — Subsampling: filter samples and optionally rarefy read depths.
"""
import io

import numpy as np
from biom import Table, load_table


def _biom_to_bytes(table: Table) -> bytes:
    """Serialize a biom.Table to HDF5 bytes."""
    import h5py

    buf = io.BytesIO()
    with h5py.File(buf, "w") as f:
        table.to_hdf5(f, generated_by="MicrobiomeDash-Subsample")
    return buf.getvalue()


def _remove_empty_obs(table: Table) -> Table:
    """Remove observations (ASVs) with all-zero counts."""
    obs_sums = table.sum(axis="observation")
    keep = [
        oid
        for oid, s in zip(table.ids(axis="observation"), obs_sums)
        if s > 0
    ]
    if len(keep) == len(table.ids(axis="observation")):
        return table
    return table.filter(keep, axis="observation", inplace=False)


def filter_samples(biom_path: str, keep_sample_ids: list[str]) -> bytes:
    """Load BIOM, keep only specified samples, remove empty ASVs.

    Returns HDF5 BIOM bytes.
    """
    table = load_table(biom_path)
    filtered = table.filter(keep_sample_ids, axis="sample", inplace=False)
    filtered = _remove_empty_obs(filtered)
    return _biom_to_bytes(filtered)


def rarefy_samples(
    biom_path: str, keep_sample_ids: list[str], depth: int
) -> bytes:
    """Load BIOM, keep only specified samples, rarefy to depth, remove empty ASVs.

    Rarefaction subsamples each sample to exactly `depth` reads without
    replacement. Samples with fewer than `depth` reads are excluded.

    Returns HDF5 BIOM bytes.
    """
    table = load_table(biom_path)
    filtered = table.filter(keep_sample_ids, axis="sample", inplace=False)

    # Drop samples below the rarefaction depth
    sample_sums = filtered.sum(axis="sample")
    sufficient = [
        sid
        for sid, s in zip(filtered.ids(axis="sample"), sample_sums)
        if s >= depth
    ]
    if not sufficient:
        raise ValueError(
            f"No samples have >= {depth} reads after filtering."
        )
    filtered = filtered.filter(sufficient, axis="sample", inplace=False)

    # Rarefy: random subsampling without replacement per sample
    rng = np.random.default_rng()
    data = filtered.matrix_data.toarray()  # obs x sample
    n_obs, n_samples = data.shape
    rarefied = np.zeros_like(data)

    for j in range(n_samples):
        col = data[:, j].astype(int)
        total = col.sum()
        if total < depth:
            continue
        # Build pool of observation indices, one per read
        pool = np.repeat(np.arange(n_obs), col)
        chosen = rng.choice(pool, size=depth, replace=False)
        obs_indices, counts = np.unique(chosen, return_counts=True)
        rarefied[obs_indices, j] = counts

    rarefied_table = Table(
        rarefied.astype(float),
        observation_ids=filtered.ids(axis="observation"),
        sample_ids=filtered.ids(axis="sample"),
        observation_metadata=[
            filtered.metadata(oid, axis="observation")
            for oid in filtered.ids(axis="observation")
        ],
        sample_metadata=[
            filtered.metadata(sid, axis="sample")
            for sid in filtered.ids(axis="sample")
        ],
        type="OTU table",
    )
    rarefied_table = _remove_empty_obs(rarefied_table)
    return _biom_to_bytes(rarefied_table)
