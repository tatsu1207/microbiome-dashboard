"""
MicrobiomeDash — R subprocess wrapper for analysis pages.
"""
import json
import logging
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from app.config import R_SCRIPTS_DIR, conda_cmd

logger = logging.getLogger(__name__)


def run_r_script(script_name: str, args_dict: dict, on_line=None) -> dict:
    """Call an R script via conda run, capture stdout for JSON status.

    Args:
        on_line: optional callback invoked with each stderr line for real-time logging.

    Returns {success, output_path, error, stdout, stderr}.
    """
    script_path = R_SCRIPTS_DIR / script_name
    if not script_path.exists():
        return {"success": False, "error": f"R script not found: {script_path}"}

    cmd_args = ["Rscript", str(script_path)]
    for key, value in args_dict.items():
        cmd_args.extend([f"--{key}", str(value)])

    full_cmd = conda_cmd(cmd_args)
    logger.info(f"Running R: {' '.join(full_cmd)}")

    if on_line is not None:
        return _run_r_streaming(full_cmd, args_dict, on_line)

    try:
        result = subprocess.run(
            full_cmd,
            capture_output=True,
            text=True,
            timeout=1800,  # 30 min timeout
        )
    except subprocess.TimeoutExpired:
        return {"success": False, "error": "R script timed out after 30 minutes"}
    except Exception as e:
        return {"success": False, "error": str(e)}

    return _parse_r_result(result.returncode, result.stdout, result.stderr,
                           args_dict.get("output"))


def _run_r_streaming(full_cmd, args_dict, on_line):
    """Run R with line-by-line stderr streaming via on_line callback."""
    try:
        proc = subprocess.Popen(
            full_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        stderr_lines = []
        for line in proc.stderr:
            line = line.rstrip("\n")
            if line:
                stderr_lines.append(line)
                on_line(line)

        stdout = proc.stdout.read()
        proc.wait(timeout=1800)

    except subprocess.TimeoutExpired:
        proc.kill()
        return {"success": False, "error": "R script timed out after 30 minutes"}
    except Exception as e:
        return {"success": False, "error": str(e)}

    stderr_text = "\n".join(stderr_lines)
    return _parse_r_result(proc.returncode, stdout, stderr_text,
                           args_dict.get("output"))


def _parse_r_result(returncode, stdout, stderr, output_path):
    """Parse R script output into standard result dict."""
    status_json = {}
    if stdout:
        for line in reversed(stdout.strip().splitlines()):
            line = line.strip()
            if line.startswith("{"):
                try:
                    status_json = json.loads(line)
                    break
                except json.JSONDecodeError:
                    pass

    if returncode != 0:
        error_msg = stderr.strip() if stderr else "R script failed"
        return {
            "success": False,
            "error": error_msg,
            "stdout": stdout,
            "stderr": stderr,
        }

    return {
        "success": True,
        "status": status_json,
        "stdout": stdout,
        "stderr": stderr,
        "output_path": output_path,
    }


def prepare_da_inputs(
    count_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    sample_id_col: str,
    matched_samples: list[str],
    group_col: str,
    ref_group: str,
    test_group: str,
    output_dir: str | None = None,
) -> tuple[str, str]:
    """Write count_matrix.tsv and metadata.tsv subsetted to two groups.

    Returns (counts_path, meta_path).
    """
    if output_dir is None:
        output_dir = tempfile.mkdtemp()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Subset metadata to matched samples in the two groups
    meta_sub = meta_df[meta_df[sample_id_col].isin(matched_samples)].copy()
    meta_sub = meta_sub[meta_sub[group_col].isin([ref_group, test_group])]
    keep_samples = meta_sub[sample_id_col].tolist()

    # Subset counts to those samples
    count_sub = count_df[[s for s in count_df.columns if s in keep_samples]]

    # Remove features with all zeros
    count_sub = count_sub[count_sub.sum(axis=1) > 0]

    # Write counts TSV (feature_id + sample columns)
    counts_path = str(output_dir / "count_matrix.tsv")
    count_sub.index.name = "feature_id"
    count_sub.to_csv(counts_path, sep="\t")

    # Write metadata TSV
    meta_path = str(output_dir / "metadata.tsv")
    meta_out = meta_sub[[sample_id_col, group_col]].copy()
    meta_out.columns = ["SampleID", group_col]
    meta_out.to_csv(meta_path, sep="\t", index=False)

    return counts_path, meta_path
