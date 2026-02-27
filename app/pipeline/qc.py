"""
MicrobiomeDash — FastQC wrapper.
"""
import logging
import subprocess
import zipfile
from pathlib import Path

from app.config import DADA2_DEFAULTS, conda_cmd


def run_fastqc(
    fastq_dir: Path,
    output_dir: Path,
    logger: logging.Logger,
) -> dict:
    """Run FastQC on all FASTQ files in a directory.

    Returns:
        dict with: success (bool), report_dir (str), sample_metrics (dict)
    """
    qc_dir = output_dir / "fastqc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    # Find FASTQ files
    fastq_files = sorted(
        list(fastq_dir.glob("*.fastq.gz")) + list(fastq_dir.glob("*.fq.gz"))
    )
    if not fastq_files:
        raise FileNotFoundError(f"No FASTQ files in {fastq_dir}")

    threads = DADA2_DEFAULTS.get("threads", 1)

    logger.info(f"Running FastQC on {len(fastq_files)} files...")
    result = subprocess.run(
        conda_cmd(["fastqc", "--outdir", str(qc_dir), "--threads", str(threads)]
        + [str(f) for f in fastq_files]),
        capture_output=True,
        text=True,
        timeout=3600,
    )

    if result.returncode != 0:
        logger.error(f"FastQC stderr: {result.stderr}")
        raise RuntimeError(f"FastQC failed: {result.stderr[:500]}")

    logger.info("FastQC completed")

    # Parse metrics from FastQC zip files
    sample_metrics = {}
    for zpath in sorted(qc_dir.glob("*_fastqc.zip")):
        sample_name = zpath.stem.replace("_fastqc", "")
        metrics = _parse_fastqc_zip(zpath)
        if metrics:
            sample_metrics[sample_name] = metrics

    return {
        "success": True,
        "report_dir": str(qc_dir),
        "sample_metrics": sample_metrics,
    }


def _parse_fastqc_zip(zip_path: Path) -> dict:
    """Extract key metrics from a FastQC zip archive."""
    metrics = {}
    try:
        with zipfile.ZipFile(zip_path) as zf:
            # Find the fastqc_data.txt inside the zip
            data_files = [n for n in zf.namelist() if n.endswith("fastqc_data.txt")]
            if not data_files:
                return metrics

            with zf.open(data_files[0]) as f:
                for line in f:
                    line = line.decode("utf-8", errors="replace").strip()
                    if line.startswith("Total Sequences"):
                        metrics["total_sequences"] = int(line.split("\t")[1])
                    elif line.startswith("%GC"):
                        metrics["gc_percent"] = float(line.split("\t")[1])
                    elif line.startswith("Sequence length"):
                        metrics["sequence_length"] = line.split("\t")[1]
                    elif line.startswith("Sequences flagged as poor quality"):
                        metrics["poor_quality"] = int(line.split("\t")[1])

            # Parse summary.txt for pass/warn/fail
            summary_files = [n for n in zf.namelist() if n.endswith("summary.txt")]
            if summary_files:
                with zf.open(summary_files[0]) as f:
                    for line in f:
                        parts = line.decode("utf-8", errors="replace").strip().split("\t")
                        if len(parts) >= 2:
                            metrics[f"module_{parts[1]}"] = parts[0]

    except Exception:
        pass
    return metrics
