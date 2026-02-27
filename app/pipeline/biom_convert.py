"""
MicrobiomeDash — Convert DADA2 ASV table (TSV) to BIOM format.

Reads the TSV output from run_dada2.R and creates an HDF5 BIOM table
with ASV sequences stored as observation metadata.
"""
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from biom import Table


def tsv_to_biom(
    tsv_path: Path,
    output_dir: Path,
    logger: logging.Logger,
    taxonomy_path: Path | None = None,
) -> Path:
    """Convert a DADA2 ASV table (TSV) to BIOM format.

    The input TSV has columns: ASV_ID, sequence, Sample1, Sample2, ...
    The output BIOM stores ASV_ID as observation IDs, sequences as
    observation metadata, and sample abundances as the data matrix.

    If taxonomy_path is provided, taxonomy assignments are embedded in
    observation metadata as a list of ranks (standard BIOM format).

    Returns the path to the .biom file.
    """
    df = pd.read_csv(tsv_path, sep="\t")

    asv_ids = df["ASV_ID"].values
    sequences = df["sequence"].values
    sample_columns = [c for c in df.columns if c not in ("ASV_ID", "sequence")]
    data = df[sample_columns].values.astype(np.float64)

    # Build taxonomy lookup from taxonomy.tsv if provided
    tax_lookup: dict[str, list[str]] = {}
    if taxonomy_path and taxonomy_path.exists():
        tax_df = pd.read_csv(taxonomy_path, sep="\t")
        rank_cols = [c for c in tax_df.columns if c not in ("ASV_ID", "sequence")]
        for _, row in tax_df.iterrows():
            ranks = [str(row[c]) if pd.notna(row[c]) else "" for c in rank_cols]
            tax_lookup[row["ASV_ID"]] = ranks
        logger.info(f"Loaded taxonomy for {len(tax_lookup)} ASVs from {taxonomy_path}")

    # Observation metadata: store sequence and optionally taxonomy per row
    obs_metadata = []
    for asv_id, seq in zip(asv_ids, sequences):
        md: dict = {"sequence": seq}
        if asv_id in tax_lookup:
            md["taxonomy"] = tax_lookup[asv_id]
        obs_metadata.append(md)

    table = Table(
        data,
        observation_ids=asv_ids,
        sample_ids=sample_columns,
        observation_metadata=obs_metadata,
        type="OTU table",
    )
    import h5py

    biom_path = output_dir / "asv_table.biom"
    with h5py.File(biom_path, "w") as f:
        table.to_hdf5(f, generated_by="MicrobiomeDash-DADA2")

    tax_msg = " (with taxonomy)" if tax_lookup else ""
    logger.info(
        f"BIOM table written{tax_msg}: {biom_path} "
        f"({len(asv_ids)} ASVs x {len(sample_columns)} samples)"
    )
    return biom_path


def extract_fasta_from_biom(
    biom_path: Path,
    output_path: Path,
    logger: logging.Logger,
) -> Path:
    """Extract representative sequences from BIOM observation metadata.

    Each observation is expected to have metadata with a "sequence" key.
    Writes a FASTA file with >ASV_ID headers.

    Returns the path to the written FASTA file.
    """
    from biom import load_table

    table = load_table(str(biom_path))
    obs_ids = table.ids(axis="observation")

    lines: list[str] = []
    for obs_id in obs_ids:
        md = table.metadata(obs_id, axis="observation")
        if not md or "sequence" not in md:
            continue
        lines.append(f">{obs_id}")
        lines.append(md["sequence"])

    if not lines:
        raise ValueError(
            "No sequences found in BIOM observation metadata. "
            "The BIOM file must contain 'sequence' in observation metadata."
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n")
    n_seqs = len(lines) // 2
    logger.info(f"Extracted {n_seqs} sequences from BIOM to {output_path}")
    return output_path
