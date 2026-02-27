"""
MicrobiomeDash — Shared helpers for analysis pages.

Common BIOM/metadata input handling used across all 5 analysis pages.
"""
import base64
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from biom import load_table

from app.utils.metadata_parser import parse_metadata


def parse_uploaded_biom(contents: str, filename: str) -> tuple[str | None, str | None]:
    """Decode a base64 dcc.Upload BIOM file, save to a temp path.

    Returns (temp_path, error).
    """
    try:
        _, content_string = contents.split(",", 1)
        decoded = base64.b64decode(content_string)
    except Exception as e:
        return None, f"Could not decode file: {e}"

    suffix = ".biom" if filename.endswith(".biom") else ""
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    tmp.write(decoded)
    tmp.close()

    # Verify it's a valid BIOM table
    try:
        load_table(tmp.name)
    except Exception as e:
        Path(tmp.name).unlink(missing_ok=True)
        return None, f"Invalid BIOM file: {e}"

    return tmp.name, None


def parse_uploaded_metadata(contents: str, filename: str) -> tuple[pd.DataFrame | None, str | None, str | None]:
    """Decode a base64 dcc.Upload metadata file.

    Returns (df, sample_id_col, error).
    """
    try:
        _, content_string = contents.split(",", 1)
        decoded = base64.b64decode(content_string).decode("utf-8")
    except Exception as e:
        return None, None, f"Could not decode file: {e}"

    result = parse_metadata(decoded, filename)
    if not result["valid"]:
        return None, None, "; ".join(result["errors"])
    return result["df"], result["sample_id_column"], None


def get_pipeline_biom_options() -> list[dict]:
    """Query DB for completed pipeline datasets with BIOM files."""
    from app.db.database import SessionLocal
    from app.db.models import Dataset

    db = SessionLocal()
    try:
        datasets = (
            db.query(Dataset)
            .filter(Dataset.status == "complete", Dataset.asv_table_path.isnot(None))
            .order_by(Dataset.id.desc())
            .all()
        )
        options = []
        for ds in datasets:
            biom_path = ds.asv_table_path
            if biom_path and Path(biom_path).exists():
                label = f"#{ds.id}: {ds.name}"
                if ds.sample_count:
                    label += f" ({ds.sample_count} samples)"
                options.append({"label": label, "value": biom_path})
        return options
    finally:
        db.close()


def get_picrust2_run_options() -> list[dict]:
    """Query DB for completed PICRUSt2 runs."""
    from app.db.database import SessionLocal
    from app.db.models import Picrust2Run

    db = SessionLocal()
    try:
        runs = (
            db.query(Picrust2Run)
            .filter(Picrust2Run.status == "complete", Picrust2Run.output_dir.isnot(None))
            .order_by(Picrust2Run.id.desc())
            .all()
        )
        options = []
        for run in runs:
            if run.output_dir and Path(run.output_dir).exists():
                options.append({
                    "label": f"#{run.id}: {run.name}",
                    "value": str(run.id),
                })
        return options
    finally:
        db.close()


def validate_metadata_vs_biom(
    meta_df: pd.DataFrame, sample_id_col: str, biom_sample_ids: list[str]
) -> dict:
    """Cross-check metadata sample IDs against BIOM sample IDs."""
    meta_ids = set(meta_df[sample_id_col].astype(str).str.strip())
    biom_ids = set(biom_sample_ids)
    matched = sorted(meta_ids & biom_ids)
    return {
        "matched": matched,
        "meta_only": sorted(meta_ids - biom_ids),
        "biom_only": sorted(biom_ids - meta_ids),
    }


def get_group_columns(meta_df: pd.DataFrame, sample_id_col: str) -> list[str]:
    """Return non-SampleID columns from the metadata DataFrame."""
    return [c for c in meta_df.columns if c != sample_id_col]


def get_dataset_metadata_df(biom_path: str) -> tuple[pd.DataFrame | None, str | None]:
    """Reconstruct metadata DataFrame from DB for a pipeline dataset.

    Looks up the Dataset by its asv_table_path, then pivots SampleMetadata
    key-value rows into a wide DataFrame.  The SampleID column uses the BIOM
    sample IDs (which may include filename suffixes) so that downstream
    validation against the BIOM table matches correctly.

    Returns (df, sample_id_col) or (None, None) if no metadata found.
    """
    from app.db.database import SessionLocal
    from app.db.models import Dataset, Sample, SampleMetadata

    db = SessionLocal()
    try:
        dataset = (
            db.query(Dataset)
            .filter(Dataset.asv_table_path == biom_path)
            .first()
        )
        if not dataset:
            return None, None

        samples = (
            db.query(Sample)
            .filter(Sample.dataset_id == dataset.id)
            .all()
        )
        if not samples:
            return None, None

        # Build mapping from DB sample_name -> BIOM sample ID.
        # BIOM sample IDs often include filename suffixes (e.g. "Cow_1_F_filt.fastq.gz")
        # while the DB stores the clean name ("Cow_1").
        try:
            table = load_table(biom_path)
            biom_ids = list(table.ids(axis="sample"))
        except Exception:
            biom_ids = []

        db_to_biom: dict[str, str] = {}
        if biom_ids:
            for db_name in [s.sample_name for s in samples]:
                # Find the BIOM ID that starts with this sample name
                for bid in biom_ids:
                    if str(bid) == db_name or str(bid).startswith(db_name + "_"):
                        db_to_biom[db_name] = str(bid)
                        break

        rows = []
        has_metadata = False
        for sample in samples:
            # Use the BIOM sample ID if we have a mapping, otherwise the DB name
            sample_id = db_to_biom.get(sample.sample_name, sample.sample_name)
            row = {"SampleID": sample_id}
            for entry in sample.metadata_entries:
                row[entry.key] = entry.value
                has_metadata = True
            rows.append(row)

        if not has_metadata:
            return None, None

        df = pd.DataFrame(rows)
        return df, "SampleID"
    finally:
        db.close()


def find_metadata_for_samples(
    biom_sample_ids: list[str],
) -> tuple[pd.DataFrame | None, str | None, str | None]:
    """Search all pipeline datasets for metadata matching the given BIOM sample IDs.

    Uses the same prefix matching as get_dataset_metadata_df to map DB sample
    names to BIOM sample IDs.  Picks the dataset with the most matched samples.

    Returns (df, sample_id_col, dataset_name) or (None, None, None).
    """
    from app.db.database import SessionLocal
    from app.db.models import Dataset, Sample

    db = SessionLocal()
    try:
        datasets = (
            db.query(Dataset)
            .filter(Dataset.status == "complete")
            .order_by(Dataset.id.desc())
            .all()
        )

        biom_ids = [str(s) for s in biom_sample_ids]
        best_df = None
        best_sid_col = None
        best_name = None
        best_overlap = 0

        for ds in datasets:
            samples = (
                db.query(Sample).filter(Sample.dataset_id == ds.id).all()
            )
            if not samples:
                continue

            # Prefix-match DB sample names to BIOM sample IDs
            db_to_biom: dict[str, str] = {}
            for sample in samples:
                db_name = sample.sample_name
                for bid in biom_ids:
                    if bid == db_name or bid.startswith(db_name + "_"):
                        db_to_biom[db_name] = bid
                        break

            if not db_to_biom:
                continue

            # Build metadata DataFrame for matched samples
            rows = []
            has_metadata = False
            for sample in samples:
                mapped_id = db_to_biom.get(sample.sample_name)
                if not mapped_id:
                    continue
                row = {"SampleID": mapped_id}
                for entry in sample.metadata_entries:
                    row[entry.key] = entry.value
                    has_metadata = True
                rows.append(row)

            if not has_metadata or not rows:
                continue

            if len(rows) > best_overlap:
                best_overlap = len(rows)
                best_df = pd.DataFrame(rows)
                best_sid_col = "SampleID"
                best_name = ds.name

        if best_df is not None:
            return best_df, best_sid_col, best_name
        return None, None, None
    finally:
        db.close()


def biom_to_count_df(biom_path: str) -> pd.DataFrame:
    """Load a BIOM table and return a features x samples DataFrame (integer counts)."""
    table = load_table(biom_path)
    df = pd.DataFrame(
        table.matrix_data.toarray().astype(int),
        index=table.ids(axis="observation"),
        columns=table.ids(axis="sample"),
    )
    return df
