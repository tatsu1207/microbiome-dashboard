"""
MicrobiomeDash — Metadata CSV/TSV parser and validator.
"""
from io import StringIO
from pathlib import Path

import pandas as pd

# Accepted sample-id column names (lowercased, stripped of separators)
_NORMALIZED_SAMPLE_ID = "sampleid"


def _normalize_col(name: str) -> str:
    return name.lower().strip().replace("-", "").replace("_", "").replace("#", "")


def parse_metadata(
    content: str | Path, filename: str = ""
) -> dict:
    """Parse and validate a metadata CSV/TSV file.

    Args:
        content: Raw file text or a Path to read from.
        filename: Original filename (helps with delimiter detection).

    Returns:
        dict with keys:
            valid (bool), errors (list[str]), warnings (list[str]),
            df (DataFrame | None), sample_ids (list[str]),
            sample_id_column (str | None)
    """
    result = {
        "valid": True,
        "errors": [],
        "warnings": [],
        "df": None,
        "sample_ids": [],
        "sample_id_column": None,
    }

    # Read content
    if isinstance(content, Path):
        try:
            text = content.read_text(encoding="utf-8")
        except Exception as e:
            result["valid"] = False
            result["errors"].append(f"Could not read file: {e}")
            return result
    else:
        text = content

    if not text.strip():
        result["valid"] = False
        result["errors"].append("File is empty.")
        return result

    # Auto-detect delimiter
    delimiter = "\t" if ("\t" in text.split("\n")[0] or filename.endswith(".tsv")) else ","

    try:
        df = pd.read_csv(StringIO(text), sep=delimiter, dtype=str, keep_default_na=False)
    except Exception as e:
        result["valid"] = False
        result["errors"].append(f"Could not parse file: {e}")
        return result

    if df.empty:
        result["valid"] = False
        result["errors"].append("File contains no data rows.")
        return result

    # Find sample-id column
    sample_col = None
    for col in df.columns:
        if _normalize_col(col) == _NORMALIZED_SAMPLE_ID:
            sample_col = col
            break

    if sample_col is None:
        result["valid"] = False
        result["errors"].append(
            "Missing required column: 'sample-id' "
            "(also accepts 'sample_id', 'SampleID', '#SampleID')"
        )
        result["df"] = df
        return result

    # Check for empty sample IDs
    empty_mask = df[sample_col].str.strip() == ""
    if empty_mask.any():
        n = empty_mask.sum()
        result["valid"] = False
        result["errors"].append(f"{n} sample ID(s) are empty.")

    # Check for duplicate sample IDs
    duplicates = df[sample_col][df[sample_col].duplicated(keep=False)].unique().tolist()
    if duplicates:
        result["valid"] = False
        result["errors"].append(f"Duplicate sample IDs: {duplicates}")

    # Warnings for potential issues
    if len(df.columns) < 2:
        result["warnings"].append("Only one column found — metadata typically has additional columns.")

    result["df"] = df
    result["sample_ids"] = df[sample_col].str.strip().tolist()
    result["sample_id_column"] = sample_col
    return result


def validate_sample_ids_match(
    metadata_ids: list[str], fastq_sample_names: list[str]
) -> dict:
    """Cross-validate metadata sample IDs against FASTQ-derived sample names.

    Returns:
        dict with keys: matched (list), in_metadata_only (list), in_fastq_only (list)
    """
    meta_set = set(metadata_ids)
    fastq_set = set(fastq_sample_names)
    return {
        "matched": sorted(meta_set & fastq_set),
        "in_metadata_only": sorted(meta_set - fastq_set),
        "in_fastq_only": sorted(fastq_set - meta_set),
    }
