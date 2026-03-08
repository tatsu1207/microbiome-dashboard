"""
MicrobiomeDash — SILVA taxonomy assignment R script wrapper.
"""
import json
import logging
import subprocess
from collections.abc import Callable
from pathlib import Path

from app.config import DADA2_ENV_NAME, R_SCRIPTS_DIR, SILVA_SPECIES, SILVA_TRAIN_SET, conda_cmd


def run_taxonomy(
    rep_seqs_path: Path,
    output_dir: Path,
    threads: int,
    logger: logging.Logger,
    proc_callback: Callable[[subprocess.Popen], None] | None = None,
    skip_species: bool = True,
) -> dict:
    """Run SILVA 138.1 taxonomy assignment via the R script.

    Returns:
        dict with: success (bool), taxonomy_path (str)
    """
    taxonomy_path = output_dir / "taxonomy.tsv"

    cmd_args = [
        "Rscript", str(R_SCRIPTS_DIR / "run_taxonomy.R"),
        "--rep_seqs", str(rep_seqs_path),
        "--output", str(taxonomy_path),
        "--silva_train", str(SILVA_TRAIN_SET),
        "--threads", str(threads),
    ]
    if skip_species:
        cmd_args.append("--skip_species")
    else:
        cmd_args.extend(["--silva_species", str(SILVA_SPECIES)])

    cmd = conda_cmd(cmd_args, env_name=DADA2_ENV_NAME)

    logger.info("Running taxonomy assignment (SILVA 138.1)...")

    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
        start_new_session=True,
    )
    if proc_callback:
        proc_callback(proc)

    last_line = ""
    for line in proc.stdout:
        stripped = line.rstrip()
        if stripped:
            logger.info(f"[Taxonomy] {stripped}")
            last_line = stripped

    proc.wait()

    if proc.returncode != 0:
        raise RuntimeError(f"Taxonomy R script failed (exit code {proc.returncode})")

    # Parse JSON status
    try:
        r_status = json.loads(last_line)
        if r_status.get("status") == "error":
            raise RuntimeError(f"Taxonomy error: {r_status.get('message', 'unknown')}")
    except json.JSONDecodeError:
        pass

    if not taxonomy_path.exists():
        raise FileNotFoundError(f"Taxonomy output not found: {taxonomy_path}")

    logger.info(f"Taxonomy assignment complete: {taxonomy_path}")

    return {"success": True, "taxonomy_path": str(taxonomy_path)}
