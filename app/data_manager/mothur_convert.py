"""
MicrobiomeDash — Bidirectional BIOM / MOTHUR conversion.

BIOM → MOTHUR: produces a ZIP with {name}.fasta + {name}.count_table
MOTHUR → BIOM: reads a FASTA + count_table and produces a BIOM (HDF5) file.
"""
import io
import tempfile
import zipfile
from pathlib import Path

import numpy as np
from biom import load_table
from biom.table import Table


def biom_to_mothur_zip(biom_path: str, name: str = "sequences") -> bytes:
    """Convert a BIOM file to a ZIP containing MOTHUR-compatible files.

    The BIOM file must have 'sequence' in observation metadata.
    *name* is used as the base filename for the output files
    (e.g. name="asv_table" → asv_table.fasta + asv_table.count_table).

    Returns ZIP file contents as bytes.
    """
    table = load_table(str(biom_path))
    obs_ids = table.ids(axis="observation")
    sample_ids = table.ids(axis="sample")

    # Build FASTA and count_table simultaneously
    fasta_lines: list[str] = []
    count_rows: list[str] = []

    # Count_table header
    count_rows.append(
        "Representative_Sequence\ttotal\t" + "\t".join(sample_ids)
    )

    for obs_id in obs_ids:
        md = table.metadata(obs_id, axis="observation")
        if not md or "sequence" not in md:
            continue

        seq = md["sequence"]

        # FASTA entry
        fasta_lines.append(f">{obs_id}")
        fasta_lines.append(seq)

        # Count_table row
        counts = table.data(obs_id, axis="observation", dense=True)
        int_counts = counts.astype(np.int64)
        total = int(int_counts.sum())
        row_vals = "\t".join(str(c) for c in int_counts)
        count_rows.append(f"{obs_id}\t{total}\t{row_vals}")

    if not fasta_lines:
        raise ValueError(
            "No sequences found in BIOM observation metadata. "
            "The BIOM file must contain 'sequence' in observation metadata."
        )

    fasta_content = "\n".join(fasta_lines) + "\n"
    count_content = "\n".join(count_rows) + "\n"

    # Create ZIP in memory
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(f"{name}.fasta", fasta_content)
        zf.writestr(f"{name}.count_table", count_content)

    return buf.getvalue()


def mothur_to_biom(fasta_path: str, count_table_path: str) -> bytes:
    """Convert MOTHUR FASTA + count_table into a BIOM (HDF5) file.

    Args:
        fasta_path: path to the representative-sequences FASTA file.
        count_table_path: path to the MOTHUR count_table (TSV with
            Representative_Sequence, total, sample1, sample2, ...).

    Returns:
        BIOM file contents as bytes (HDF5 format).
    """
    # ── Parse FASTA ──────────────────────────────────────────────────────
    sequences: dict[str, str] = {}
    current_id = None
    seq_parts: list[str] = []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(seq_parts)
                current_id = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if current_id is not None:
            sequences[current_id] = "".join(seq_parts)

    if not sequences:
        raise ValueError("No sequences found in FASTA file.")

    # ── Parse count_table ────────────────────────────────────────────────
    with open(count_table_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        # header: Representative_Sequence, total, sample1, sample2, ...
        sample_ids = header[2:]  # skip Representative_Sequence and total

        obs_ids: list[str] = []
        data_rows: list[list[int]] = []
        for line in f:
            parts = line.rstrip("\n").split("\t")
            obs_id = parts[0]
            counts = [int(x) for x in parts[2:]]  # skip total column
            obs_ids.append(obs_id)
            data_rows.append(counts)

    if not obs_ids:
        raise ValueError("No observations found in count_table.")
    if not sample_ids:
        raise ValueError("No samples found in count_table header.")

    # ── Build BIOM table ─────────────────────────────────────────────────
    data_matrix = np.array(data_rows, dtype=np.int64)

    observation_metadata = []
    for obs_id in obs_ids:
        md = {}
        if obs_id in sequences:
            md["sequence"] = sequences[obs_id]
        observation_metadata.append(md)

    table = Table(
        data_matrix,
        observation_ids=obs_ids,
        sample_ids=sample_ids,
        observation_metadata=observation_metadata,
    )

    # Write to HDF5 bytes
    import h5py

    tmp = tempfile.NamedTemporaryFile(suffix=".biom", delete=False)
    tmp.close()
    tmp_path = Path(tmp.name)
    try:
        with h5py.File(tmp_path, "w") as f:
            table.to_hdf5(f, generated_by="MicrobiomeDash")
        return tmp_path.read_bytes()
    finally:
        tmp_path.unlink(missing_ok=True)
