"""
MicrobiomeDash — PICRUSt2 functional prediction wrapper.

This step is optional — if PICRUSt2 is not installed or fails,
the pipeline continues and the dataset is still marked as complete.
"""
import logging
import shutil
import subprocess
import threading
import time
from pathlib import Path

from app.config import PICRUST2_ENV_NAME, conda_cmd

# PICRUSt2 sub-steps and their progress percentages.
# Two variants: full (EC+KO) and EC-only (skip KO).
_SUBSTEPS_FULL = [
    ("place_seqs",       "Placing sequences",          5),
    ("EC_predicted",     "EC hidden-state prediction", 20),
    ("KO_predicted",     "KO hidden-state prediction", 50),
    ("EC_metagenome_out", "EC metagenome prediction",  65),
    ("KO_metagenome_out", "KO metagenome prediction",  80),
    ("pathways_out",      "Pathway prediction",        95),
]
_SUBSTEPS_EC = [
    ("place_seqs",       "Placing sequences",          10),
    ("EC_predicted",     "EC hidden-state prediction",  35),
    ("EC_metagenome_out", "EC metagenome prediction",   65),
    ("pathways_out",      "Pathway prediction",         95),
]


def run_picrust2(
    asv_table_path: Path,
    rep_seqs_path: Path,
    output_dir: Path,
    threads: int,
    logger: logging.Logger,
    skip_ko: bool = False,
    progress_callback=None,
    biom_ready: bool = False,
    proc_callback=None,
) -> dict:
    """Run PICRUSt2 functional prediction.

    Args:
        progress_callback: optional callable(substep: str, pct: int) called
            as PICRUSt2 progresses through sub-steps.
        biom_ready: if True, asv_table_path is already a BIOM file (skip conversion).
        proc_callback: optional callable(proc: Popen) called after starting
            the subprocess, allowing callers to track or cancel it.

    Returns:
        dict with: success (bool), picrust_dir (str | None), error (str | None)
    """
    picrust_dir = output_dir / "picrust2"

    # Check if PICRUSt2 is available
    picrust_cmd = _find_picrust2()
    if picrust_cmd is None:
        msg = "PICRUSt2 not found — skipping functional prediction"
        logger.warning(msg)
        return {"success": False, "picrust_dir": None, "error": msg}

    try:
        # PICRUSt2 needs a BIOM table
        if biom_ready:
            biom_path = asv_table_path
        else:
            biom_path = output_dir / "asv_table.biom"
            _convert_to_biom(asv_table_path, biom_path, logger)

        # Run PICRUSt2 with EPA-ng (lower memory than SEPP) and mp HSP method
        cmd = picrust_cmd + [
            "-s", str(rep_seqs_path),
            "-i", str(biom_path),
            "-o", str(picrust_dir),
            "-p", str(_picrust2_processes(threads)),
            "--hsp_method", "mp",
        ]
        if skip_ko:
            cmd += ["--in_traits", "EC"]
            logger.info("KO prediction skipped by user request")

        if progress_callback:
            progress_callback("Starting PICRUSt2", 0)

        logger.info("Running PICRUSt2...")
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
            start_new_session=True,
        )
        if proc_callback:
            proc_callback(proc)

        # Monitor output directory for sub-step milestones
        stop_event = threading.Event()
        monitor = threading.Thread(
            target=_monitor_progress,
            args=(picrust_dir, skip_ko, progress_callback, logger, stop_event),
            daemon=True,
        )
        if progress_callback:
            monitor.start()

        for line in proc.stdout:
            stripped = line.rstrip()
            if stripped:
                logger.info(f"[PICRUSt2] {stripped}")
        proc.wait()
        stop_event.set()

        if proc.returncode != 0:
            msg = f"PICRUSt2 exited with code {proc.returncode}"
            logger.warning(msg)
            return {"success": False, "picrust_dir": None, "error": msg}

        # Add descriptions to output tables
        _add_descriptions(picrust_dir, skip_ko, logger)

        if progress_callback:
            progress_callback("Complete", 100)
        logger.info(f"PICRUSt2 complete: {picrust_dir}")
        return {"success": True, "picrust_dir": str(picrust_dir), "error": None}

    except Exception as e:
        msg = f"PICRUSt2 failed: {e}"
        logger.warning(msg)
        return {"success": False, "picrust_dir": None, "error": msg}


def _monitor_progress(picrust_dir, skip_ko, callback, logger, stop_event):
    """Poll the PICRUSt2 output directory to detect sub-step completion."""
    substeps = _SUBSTEPS_EC if skip_ko else _SUBSTEPS_FULL
    reached = -1

    while not stop_event.is_set():
        for i, (marker, label, pct) in enumerate(substeps):
            if i <= reached:
                continue
            # Check for directory or file (with .tsv.gz extension)
            if (picrust_dir / marker).exists() or (
                picrust_dir / f"{marker}.tsv.gz"
            ).exists() or (
                picrust_dir / "intermediate" / marker
            ).exists():
                reached = i
                logger.info(f"[PICRUSt2] Step: {label} ({pct}%)")
                if callback:
                    callback(label, pct)
        stop_event.wait(3)


def _picrust2_processes(threads: int) -> int:
    """Determine safe number of PICRUSt2 processes based on available RAM.

    Each HSP worker uses ~2 GB. We reserve 4 GB for the OS and other tasks.
    """
    import os
    try:
        total_ram_gb = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / (1024 ** 3)
    except (ValueError, OSError):
        total_ram_gb = 8  # conservative fallback

    max_by_ram = max(1, int((total_ram_gb - 4) / 2))
    return max(1, min(threads, max_by_ram))


def _find_picrust2() -> list[str] | None:
    """Locate the PICRUSt2 pipeline command via the conda environment."""
    try:
        result = subprocess.run(
            ["conda", "run", "-n", PICRUST2_ENV_NAME, "which", "picrust2_pipeline.py"],
            capture_output=True, text=True, timeout=10,
        )
        if result.returncode == 0:
            return ["conda", "run", "-n", PICRUST2_ENV_NAME, "picrust2_pipeline.py"]
    except Exception:
        pass

    return None


def _add_descriptions(picrust_dir: Path, skip_ko: bool, logger: logging.Logger):
    """Run add_descriptions.py on PICRUSt2 output tables.

    Adds a 'description' column to each output file, mapping IDs to
    human-readable names (e.g. MetaCyc pathway names, EC names, KO names).
    """
    add_desc_cmd = [
        "conda", "run", "-n", PICRUST2_ENV_NAME, "add_descriptions.py",
    ]

    # (input_file, map_type, output_file) — relative to picrust_dir
    targets = [
        ("pathways_out/path_abun_unstrat.tsv.gz", "METACYC",
         "pathways_out/path_abun_unstrat_described.tsv.gz"),
        ("EC_metagenome_out/pred_metagenome_unstrat.tsv.gz", "EC",
         "EC_metagenome_out/pred_metagenome_unstrat_described.tsv.gz"),
    ]
    if not skip_ko:
        targets.append(
            ("KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", "KO",
             "KO_metagenome_out/pred_metagenome_unstrat_described.tsv.gz"),
        )

    for infile, map_type, outfile in targets:
        in_path = picrust_dir / infile
        out_path = picrust_dir / outfile
        if not in_path.exists():
            logger.warning(f"add_descriptions: {in_path} not found, skipping")
            continue
        cmd = add_desc_cmd + ["-i", str(in_path), "-m", map_type, "-o", str(out_path)]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if result.returncode == 0:
                logger.info(f"Added {map_type} descriptions → {out_path.name}")
            else:
                logger.warning(f"add_descriptions ({map_type}) failed: {result.stderr}")
        except Exception as e:
            logger.warning(f"add_descriptions ({map_type}) error: {e}")


def _convert_to_biom(tsv_path: Path, biom_path: Path, logger: logging.Logger):
    """Convert an ASV table TSV to BIOM format for PICRUSt2.

    DADA2 outputs samples-as-rows, ASVs-as-columns.
    PICRUSt2/BIOM expects ASVs-as-rows (observations), samples-as-columns.
    We transpose the table before converting.
    """
    import csv

    # Read and transpose the TSV (samples×ASVs → ASVs×samples)
    transposed_path = tsv_path.parent / "asv_table_transposed.tsv"
    with open(tsv_path) as f:
        reader = csv.reader(f, delimiter="\t")
        rows = list(reader)

    # Transpose
    t_rows = list(zip(*rows))
    # Replace "sample" header with "#OTU ID" for BIOM compatibility
    header = list(t_rows[0])
    header[0] = "#OTU ID"
    t_rows[0] = header

    with open(transposed_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(t_rows)

    result = subprocess.run(
        conda_cmd([
            "biom", "convert",
            "-i", str(transposed_path),
            "-o", str(biom_path),
            "--table-type", "OTU table",
            "--to-hdf5",
        ]),
        capture_output=True, text=True, timeout=300,
    )
    if result.returncode == 0:
        transposed_path.unlink(missing_ok=True)
        return

    msg = f"biom convert failed: {result.stderr}"
    logger.warning(msg)
    raise RuntimeError(msg)
