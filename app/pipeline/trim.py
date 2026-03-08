"""
MicrobiomeDash — Cutadapt primer/adapter trimming wrapper.
"""
import logging
import os
import subprocess
from pathlib import Path

from app.config import conda_cmd
from app.pipeline.detect import (
    PRIMERS,
    _primer_matches,
    _read_fastq_sequences,
    detect_sequencing_type,
    extract_sample_name,
)


def run_cutadapt(
    fastq_dir: Path,
    output_dir: Path,
    sequencing_type: str,
    variable_region: str | None,
    logger: logging.Logger,
    threads: int = 1,
    custom_fwd_primer: str | None = None,
    custom_rev_primer: str | None = None,
) -> dict:
    """Run Cutadapt to remove primers from FASTQ files.

    Cutadapt threads are capped at 32 to avoid hitting OS open-file limits,
    regardless of how many threads the caller requests.

    If custom primers are provided, they override region-based lookup and
    Cutadapt always runs (no auto-skip).
    If variable_region is None and no custom primers, files are symlinked.

    Returns:
        dict with: success (bool), trimmed_dir (str), trim_stats (dict)
    """
    trimmed_dir = output_dir / "trimmed"
    trimmed_dir.mkdir(parents=True, exist_ok=True)

    # Determine primer sequences: custom primers take precedence
    using_custom = bool(custom_fwd_primer)
    if using_custom:
        fwd_primer = custom_fwd_primer
        rev_primer = custom_rev_primer or ""
        logger.info(f"Using custom primers: F={fwd_primer}, R={rev_primer or '(none)'}")
    elif variable_region and variable_region in PRIMERS:
        fwd_primer = PRIMERS[variable_region]["forward"]
        rev_primer = PRIMERS[variable_region]["reverse"]
        logger.info(f"Trimming primers for {variable_region}: F={fwd_primer}, R={rev_primer}")
    else:
        # No primers to trim — symlink files through
        logger.warning(
            f"Variable region '{variable_region}' not recognized or not set. "
            "Skipping primer trimming."
        )
        _symlink_files(fastq_dir, trimmed_dir)
        return {
            "success": True,
            "trimmed_dir": str(trimmed_dir),
            "trim_stats": {"skipped": True},
        }

    is_paired = sequencing_type == "paired-end"
    fastq_files = sorted(
        list(fastq_dir.glob("*.fastq.gz")) + list(fastq_dir.glob("*.fq.gz"))
    )

    # Skip auto-detection when custom primers are provided — always trim
    if not using_custom and not _primers_present(fastq_files, fwd_primer, logger):
        logger.info(
            "Primers not detected in reads — files appear already trimmed. "
            "Skipping Cutadapt."
        )
        _symlink_files(fastq_dir, trimmed_dir)
        return {
            "success": True,
            "trimmed_dir": str(trimmed_dir),
            "trim_stats": {"skipped": True, "reason": "primers_absent"},
        }

    trim_stats = {"total_reads": 0, "trimmed_reads": 0, "samples_processed": 0}

    # Cap Cutadapt threads to avoid "Too many open files" errors
    threads = min(threads, 32)

    if is_paired:
        _trim_paired(
            fastq_files, trimmed_dir, fwd_primer, rev_primer, logger, trim_stats, threads
        )
    else:
        _trim_single(
            fastq_files, trimmed_dir, fwd_primer, logger, trim_stats, threads
        )

    logger.info(
        f"Cutadapt done: {trim_stats['samples_processed']} samples, "
        f"{trim_stats['trimmed_reads']}/{trim_stats['total_reads']} reads trimmed"
    )

    return {
        "success": True,
        "trimmed_dir": str(trimmed_dir),
        "trim_stats": trim_stats,
    }


def _primers_present(
    fastq_files: list[Path],
    fwd_primer: str,
    logger: logging.Logger,
    n_reads: int = 100,
    threshold: float = 0.3,
) -> bool:
    """Check whether the forward primer is present in the first FASTQ file.

    Reads *n_reads* sequences from the first file and returns True if at least
    *threshold* fraction of them start with the primer (allowing degenerate
    bases and up to 3 mismatches).
    """
    if not fastq_files:
        return False

    test_file = fastq_files[0]
    sequences = _read_fastq_sequences(test_file, n_reads)
    if not sequences:
        logger.warning(f"Could not read sequences from {test_file.name} for primer check")
        return True  # assume present if we can't check

    matches = sum(1 for seq in sequences if _primer_matches(seq, fwd_primer))
    ratio = matches / len(sequences)
    logger.info(
        f"Primer check: {matches}/{len(sequences)} reads ({ratio:.0%}) "
        f"start with forward primer in {test_file.name}"
    )
    return ratio >= threshold


def _trim_paired(
    fastq_files: list[Path],
    trimmed_dir: Path,
    fwd_primer: str,
    rev_primer: str,
    logger: logging.Logger,
    stats: dict,
    threads: int = 1,
):
    """Trim paired-end files sample by sample."""
    filenames = [f.name for f in fastq_files]
    detection = detect_sequencing_type(filenames)

    for sample_name, info in detection["samples"].items():
        r1_name = info.get("R1")
        r2_name = info.get("R2")
        if not r1_name or not r2_name:
            logger.warning(f"Skipping sample {sample_name}: missing R1 or R2")
            continue

        r1_in = fastq_files[0].parent / r1_name
        r2_in = fastq_files[0].parent / r2_name
        r1_out = trimmed_dir / r1_name
        r2_out = trimmed_dir / r2_name

        primer_args = ["-g", f"^{fwd_primer}"]
        if rev_primer:
            primer_args += ["-G", f"^{rev_primer}"]
        cmd = conda_cmd([
            "cutadapt",
            "-j", str(threads),
            *primer_args,
            "--minimum-length", "50",
            "-o", str(r1_out),
            "-p", str(r2_out),
            str(r1_in),
            str(r2_in),
        ])

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if result.returncode != 0:
            raise RuntimeError(f"Cutadapt failed for {sample_name}: {result.stderr[:300]}")

        stats["samples_processed"] += 1
        _parse_cutadapt_stats(result.stderr, stats)
        logger.info(f"  Trimmed: {sample_name}")


def _trim_single(
    fastq_files: list[Path],
    trimmed_dir: Path,
    fwd_primer: str,
    logger: logging.Logger,
    stats: dict,
    threads: int = 1,
):
    """Trim single-end files."""
    for fq in fastq_files:
        out = trimmed_dir / fq.name
        cmd = conda_cmd([
            "cutadapt",
            "-j", str(threads),
            "-g", f"^{fwd_primer}",
            "--minimum-length", "50",
            "-o", str(out),
            str(fq),
        ])

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if result.returncode != 0:
            raise RuntimeError(f"Cutadapt failed for {fq.name}: {result.stderr[:300]}")

        sample = extract_sample_name(fq.name)
        stats["samples_processed"] += 1
        _parse_cutadapt_stats(result.stderr, stats)
        logger.info(f"  Trimmed: {sample}")


def _parse_cutadapt_stats(stderr: str, stats: dict):
    """Parse cutadapt stderr for read counts."""
    for line in stderr.splitlines():
        if "Total reads processed:" in line:
            try:
                n = int(line.split(":")[-1].strip().replace(",", ""))
                stats["total_reads"] += n
            except ValueError:
                pass
        elif "Reads with adapters:" in line:
            try:
                n = int(line.split(":")[-1].strip().split("(")[0].strip().replace(",", ""))
                stats["trimmed_reads"] += n
            except ValueError:
                pass


def _process_one_fastq(args: tuple) -> tuple[str, int]:
    """Process a single FASTQ file: trim leading bases and fix N bases.

    Designed to be called via multiprocessing.Pool.
    Returns (filename, n_fixed).
    """
    import gzip

    fq_path, out_path, trim_n = args
    n_fixed = 0
    with gzip.open(fq_path, "rt") as fin, gzip.open(out_path, "wt") as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline().rstrip()[trim_n:]
            plus = fin.readline()
            qual = fin.readline().rstrip()[trim_n:]

            # Replace N bases with A at quality Q0
            if "N" in seq:
                seq_list = list(seq)
                qual_list = list(qual)
                for j, base in enumerate(seq_list):
                    if base == "N" and j < len(qual_list):
                        seq_list[j] = "A"
                        qual_list[j] = "#"  # Phred 2 (matches original N quality)
                        n_fixed += 1
                seq = "".join(seq_list)
                qual = "".join(qual_list)

            fout.write(header)
            fout.write(seq + "\n")
            fout.write(plus)
            fout.write(qual + "\n")

    return (os.path.basename(fq_path), n_fixed)


def trim_leading_bases(
    trimmed_dir: Path,
    trim_left_f: int,
    trim_left_r: int,
    sequencing_type: str,
    logger: logging.Logger,
    threads: int = 1,
) -> Path:
    """Physically trim leading bases and replace N bases in FASTQ files.

    Two fixes for DADA2's maxN=0 requirement:
    1. Trim leading bases (works around maxN applying before trimLeft).
    2. Replace remaining N bases with 'A' at quality Q0 (Phred 0, '!'),
       so DADA2 doesn't reject reads but maxEE properly penalizes them.

    Files are processed in parallel using multiprocessing.

    Writes processed files to a 'lead_trimmed' subdirectory alongside trimmed_dir.
    Returns the path to the directory with processed files.
    """
    from multiprocessing import Pool

    out_dir = trimmed_dir.parent / "lead_trimmed"
    out_dir.mkdir(parents=True, exist_ok=True)

    fastq_files = sorted(
        list(trimmed_dir.glob("*.fastq.gz")) + list(trimmed_dir.glob("*.fq.gz"))
    )

    if not fastq_files:
        logger.warning("No FASTQ files found for leading-base trimming")
        return trimmed_dir

    is_paired = sequencing_type == "paired-end"
    filenames = [f.name for f in fastq_files]
    detection = detect_sequencing_type(filenames)

    # Build set of R2 filenames for paired-end
    r2_names = set()
    if is_paired:
        for info in detection["samples"].values():
            if info.get("R2"):
                r2_names.add(info["R2"])

    # Build work items for parallel processing
    work_items = []
    for fq in fastq_files:
        is_r2 = fq.name in r2_names
        trim_n = trim_left_r if is_r2 else trim_left_f
        out_path = out_dir / fq.name
        work_items.append((str(fq), str(out_path), trim_n))

    n_workers = min(threads, len(work_items))
    logger.info(f"Processing {len(work_items)} files with {n_workers} workers")

    with Pool(processes=n_workers) as pool:
        results = pool.map(_process_one_fastq, work_items)

    for fname, n_fixed in results:
        if n_fixed:
            logger.info(f"  {fname}: fixed {n_fixed} N bases → A at Q0")

    logger.info(f"Leading-base trimming complete → {out_dir}")
    return out_dir


def _symlink_files(src_dir: Path, dst_dir: Path):
    """Symlink all FASTQ files from src to dst."""
    for ext in ("*.fastq.gz", "*.fq.gz"):
        for f in src_dir.glob(ext):
            link = dst_dir / f.name
            if not link.exists():
                os.symlink(f, link)
