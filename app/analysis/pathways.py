"""
MicrobiomeDash — PICRUSt2 pathway analysis helpers.
"""
import base64
import json
import re
import tempfile
import threading
import zipfile
from pathlib import Path

import pandas as pd

from app.analysis.r_runner import run_r_script
from app.config import PICRUST2_RUNS_DIR


# Map prediction type to expected filenames
PREDICTION_FILES = {
    "metacyc": [
        "pathways_out/path_abun_unstrat_described.tsv.gz",
        "pathways_out/path_abun_unstrat.tsv.gz",
    ],
    "ko": [
        "KO_metagenome_out/pred_metagenome_unstrat_described.tsv.gz",
        "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",
    ],
    "ec": [
        "EC_metagenome_out/pred_metagenome_unstrat_described.tsv.gz",
        "EC_metagenome_out/pred_metagenome_unstrat.tsv.gz",
    ],
}


def load_picrust2_table(
    run_id_or_dir: str, prediction_type: str
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """Load a PICRUSt2 prediction table.

    Args:
        run_id_or_dir: Either a run ID (to look up in PICRUST2_RUNS_DIR) or a directory path.
        prediction_type: 'metacyc', 'ko', or 'ec'.

    Returns:
        (counts_df, desc_df) where counts_df has numeric columns and desc_df has descriptions.
        desc_df is None if no description column found.
    """
    # Determine base directory
    try:
        run_id = int(run_id_or_dir)
        # Look up actual output_dir from DB
        from app.db.database import SessionLocal
        from app.db.models import Picrust2Run
        db = SessionLocal()
        try:
            run = db.query(Picrust2Run).filter(Picrust2Run.id == run_id).first()
            if run and run.output_dir:
                base_dir = Path(run.output_dir)
            else:
                base_dir = PICRUST2_RUNS_DIR / str(run_id)
        finally:
            db.close()
    except (ValueError, TypeError):
        base_dir = Path(run_id_or_dir)

    # Find the prediction file — first check expected relative paths,
    # then search recursively (handles ZIPs with extra nesting).
    candidates = PREDICTION_FILES.get(prediction_type, [])
    found_path = None
    for candidate in candidates:
        p = base_dir / candidate
        if p.exists():
            found_path = p
            break

    if found_path is None:
        # Recursive search by filename
        for candidate in candidates:
            fname = Path(candidate).name
            matches = list(base_dir.rglob(fname))
            if matches:
                found_path = matches[0]
                break

    if found_path is None:
        raise FileNotFoundError(
            f"No {prediction_type} prediction file found in {base_dir}. "
            f"Looked for: {', '.join(candidates)}"
        )

    return load_prediction_file(str(found_path))


def load_prediction_file(path: str) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """Load a single PICRUSt2 prediction TSV/TSV.gz file.

    Returns (counts_df, desc_df) where counts_df has features as rows and
    samples as columns.  desc_df is None if no description column found.
    """
    df = pd.read_csv(path, sep="\t")
    id_col = df.columns[0]

    desc_df = None
    desc_col = None
    for col in df.columns:
        if col.lower() == "description":
            desc_col = col
            break

    if desc_col:
        desc_df = df[[id_col, desc_col]].copy()
        desc_df.columns = ["feature", "description"]
        desc_df = desc_df.set_index("feature")

    skip_cols = {id_col}
    if desc_col:
        skip_cols.add(desc_col)
    numeric_cols = [c for c in df.columns if c not in skip_cols]

    counts_df = df[[id_col] + numeric_cols].copy()
    counts_df = counts_df.set_index(id_col)
    counts_df.index.name = "feature_id"

    return counts_df, desc_df


def parse_uploaded_prediction_file(
    contents: str, filename: str
) -> tuple[str | None, str | None]:
    """Decode a base64 uploaded TSV/TSV.gz prediction file and save to temp.

    Returns (temp_path, error).
    """
    try:
        _, content_string = contents.split(",", 1)
        decoded = base64.b64decode(content_string)
    except Exception as e:
        return None, f"Could not decode file: {e}"

    suffix = ".tsv.gz" if filename.lower().endswith(".gz") else ".tsv"
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    tmp.write(decoded)
    tmp.close()

    try:
        load_prediction_file(tmp.name)
    except Exception as e:
        Path(tmp.name).unlink(missing_ok=True)
        return None, f"Invalid prediction file: {e}"

    return tmp.name, None


def detect_prediction_label(counts_df: pd.DataFrame) -> str:
    """Auto-detect prediction type from feature IDs and return a display label."""
    ids = [str(i) for i in counts_df.index[:50]]
    ko_count = sum(1 for i in ids if re.match(r"^K\d{5}$", i))
    ec_count = sum(1 for i in ids if re.match(r"^(EC:)?\d+\.\d+\.\d+", i))
    if ko_count > len(ids) * 0.5:
        return "KO"
    if ec_count > len(ids) * 0.5:
        return "EC"
    return "MetaCyc"


def parse_picrust2_zip(contents: str, filename: str) -> tuple[str | None, str | None]:
    """Decode a base64 uploaded ZIP and extract to tempdir.

    Returns (output_dir, error).
    """
    try:
        _, content_string = contents.split(",", 1)
        decoded = base64.b64decode(content_string)
    except Exception as e:
        return None, f"Could not decode file: {e}"

    tmp_zip = tempfile.NamedTemporaryFile(delete=False, suffix=".zip")
    tmp_zip.write(decoded)
    tmp_zip.close()

    output_dir = tempfile.mkdtemp()
    try:
        with zipfile.ZipFile(tmp_zip.name, "r") as z:
            z.extractall(output_dir)
    except zipfile.BadZipFile:
        return None, "Invalid ZIP file"
    finally:
        Path(tmp_zip.name).unlink(missing_ok=True)

    # If the ZIP contained a single top-level directory, unwrap it so that
    # prediction files (e.g. pathways_out/) are found directly.
    entries = list(Path(output_dir).iterdir())
    if len(entries) == 1 and entries[0].is_dir():
        output_dir = str(entries[0])

    return output_dir, None


def run_pathway_da(
    counts_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    sample_id_col: str,
    matched_samples: list[str],
    group_col: str,
    ref_group: str,
    test_group: str,
    tool: str = "aldex2",
) -> pd.DataFrame:
    """Run differential abundance on pathway predictions.

    Args:
        tool: DA method — 'aldex2', 'deseq2', 'ancombc', 'linda', or 'maaslin2'.
    """
    from app.analysis.diff_abundance import TOOL_SCRIPTS, TOOL_LABELS
    from app.analysis.r_runner import prepare_da_inputs

    script = TOOL_SCRIPTS.get(tool)
    if not script:
        raise ValueError(f"Unknown DA tool: {tool}")

    # Round float counts to integers only for count-based methods
    if tool in ("aldex2", "deseq2"):
        work_counts = counts_df.round(0).astype(int)
    else:
        work_counts = counts_df

    output_dir = tempfile.mkdtemp()
    output_path = str(Path(output_dir) / "pathway_da_results.tsv")

    counts_path, meta_path = prepare_da_inputs(
        work_counts, meta_df, sample_id_col, matched_samples,
        group_col, ref_group, test_group, output_dir,
    )

    script_args = {
        "counts": counts_path,
        "metadata": meta_path,
        "group_col": group_col,
        "ref_group": ref_group,
        "test_group": test_group,
        "output": output_path,
    }

    result = run_r_script(script, script_args)

    if not result["success"]:
        tool_label = TOOL_LABELS.get(tool, tool)
        raise RuntimeError(f"Pathway DA ({tool_label}) failed: {result.get('error', 'unknown')}")

    return pd.read_csv(output_path, sep="\t")


def merge_descriptions(results_df: pd.DataFrame, desc_df: pd.DataFrame | None) -> pd.DataFrame:
    """Join description column onto results."""
    if desc_df is None or results_df.empty:
        return results_df

    merged = results_df.merge(
        desc_df, left_on="feature", right_index=True, how="left"
    )
    # Reorder so description is after feature
    cols = list(merged.columns)
    if "description" in cols:
        cols.remove("description")
        cols.insert(1, "description")
        merged = merged[cols]
    return merged


# ── Async pathway DA ─────────────────────────────────────────────────────────


def _pathway_da_progress_path(job_id: str) -> Path:
    return Path(f"/tmp/pw_da_{job_id}.json")


def _write_pathway_progress(job_id: str, data: dict) -> None:
    _pathway_da_progress_path(job_id).write_text(json.dumps(data))


def read_pathway_da_progress(job_id: str) -> dict | None:
    path = _pathway_da_progress_path(job_id)
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        return None


def run_pathway_da_background(
    counts_df: pd.DataFrame,
    desc_df: pd.DataFrame | None,
    meta_df: pd.DataFrame,
    sid_col: str,
    matched_samples: list[str],
    group_col: str,
    ref_group: str,
    test_group: str,
    pred_label: str,
    job_id: str,
    tool: str = "aldex2",
) -> None:
    """Spawn a background thread to run pathway DA."""
    from app.analysis.diff_abundance import TOOL_LABELS

    tool_label = TOOL_LABELS.get(tool, tool)

    def _worker():
        try:
            _write_pathway_progress(job_id, {
                "status": "running",
                "log": [f"Running {tool_label} on {pred_label} features..."],
            })

            results_df = run_pathway_da(
                counts_df, meta_df, sid_col, matched_samples,
                group_col, ref_group, test_group, tool=tool,
            )
            results_df = merge_descriptions(results_df, desc_df)

            n_total = len(results_df)
            n_sig = int((results_df["qvalue"] < 0.05).sum()) if "qvalue" in results_df.columns else 0

            _write_pathway_progress(job_id, {
                "status": "complete",
                "log": [f"{tool_label} complete: {n_sig} significant out of {n_total} (q < 0.05)"],
                "results_csv": results_df.to_csv(sep="\t", index=False),
                "pred_label": pred_label,
                "tool": tool,
                "n_total": n_total,
                "n_sig": n_sig,
            })
        except Exception as e:
            _write_pathway_progress(job_id, {
                "status": "error",
                "log": [f"Error: {e}"],
            })

    t = threading.Thread(target=_worker, daemon=True)
    t.start()
