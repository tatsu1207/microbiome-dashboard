"""
MicrobiomeDash — Phylogenetic tree construction (MAFFT + FastTree).
"""
import logging
import os
import subprocess
from pathlib import Path

from app.config import conda_cmd


def run_phylogeny(
    rep_seqs_path: Path,
    output_dir: Path,
    threads: int,
    logger: logging.Logger,
) -> dict:
    """Build a phylogenetic tree from representative sequences.

    Steps:
        1. MAFFT multiple sequence alignment
        2. FastTree approximate maximum-likelihood tree

    Returns:
        dict with: success (bool), aligned_path (str), tree_path (str)
    """
    aligned_path = output_dir / "rep_seqs_aligned.fasta"
    tree_path = output_dir / "tree.nwk"

    # Step 1: MAFFT alignment
    logger.info("Running MAFFT alignment...")
    with open(aligned_path, "w") as out_f:
        result = subprocess.run(
            conda_cmd(["mafft", "--auto", "--thread", str(threads), str(rep_seqs_path)]),
            stdout=out_f,
            stderr=subprocess.PIPE,
            text=True,
            timeout=3600,
        )

    if result.returncode != 0:
        raise RuntimeError(f"MAFFT failed: {result.stderr[:500]}")

    logger.info(f"MAFFT alignment complete: {aligned_path}")

    # Step 2: FastTree
    logger.info("Running FastTree...")
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)

    with open(tree_path, "w") as out_f:
        result = subprocess.run(
            conda_cmd(["FastTree", "-gtr", "-nt", str(aligned_path)]),
            stdout=out_f,
            stderr=subprocess.PIPE,
            text=True,
            timeout=3600,
            env=env,
        )

    if result.returncode != 0:
        raise RuntimeError(f"FastTree failed: {result.stderr[:500]}")

    logger.info(f"Phylogenetic tree complete: {tree_path}")

    return {
        "success": True,
        "aligned_path": str(aligned_path),
        "tree_path": str(tree_path),
    }
