"""
MicrobiomeDash — Quality profiling and auto-detection of DADA2 truncation parameters.

Reads quality scores from FASTQ files and determines optimal
trunc_len_f / trunc_len_r so that >=70% of reads pass DADA2's maxEE<=2 filter,
while maintaining sufficient overlap for paired-end merging.

For paired-end data, DADA2's filterAndTrim discards the entire pair if either
direction fails, so we use a higher per-direction threshold (sqrt of target)
to achieve the desired combined retention.
"""
import gzip
import logging
import math
from pathlib import Path

from app.pipeline.detect import PRIMERS, detect_sequencing_type

# Minimum truncation length we'll consider (shorter than this loses too much data)
_MIN_TRUNC = 100


def detect_truncation_params(
    trimmed_dir: Path,
    sequencing_type: str,
    variable_region: str | None,
    min_overlap: int = 12,
    target_retention: float = 0.70,
    n_reads: int = 5000,
    primer_offset_f: int = 0,
    primer_offset_r: int = 0,
    logger: logging.Logger | None = None,
) -> dict:
    """Auto-detect optimal DADA2 truncation lengths from quality profiles.

    Args:
        trimmed_dir: Directory with FASTQ(.gz) files.
        sequencing_type: "paired-end" or "single-end".
        variable_region: Detected 16S region (e.g., "V3-V4") for amplicon length.
        min_overlap: Minimum bases of overlap required for PE merging.
        target_retention: Fraction of reads that must pass EE<=2 (default 0.70).
        n_reads: Number of reads to sample for profiling.
        primer_offset_f: Bases to skip at start of R1 (primer length, for raw reads).
        primer_offset_r: Bases to skip at start of R2 (primer length, for raw reads).
        logger: Optional logger.

    Returns:
        dict with: trunc_len_f, trunc_len_r, details (str)
    """
    log = logger or logging.getLogger(__name__)

    # Find FASTQ files and use detect_sequencing_type to properly identify R1/R2
    all_fastq = sorted(trimmed_dir.glob("*.fastq.gz")) + sorted(
        trimmed_dir.glob("*.fq.gz")
    )
    if not all_fastq:
        log.warning("No FASTQ files found for quality profiling")
        return {"trunc_len_f": 0, "trunc_len_r": 0, "details": "No files found"}

    filenames = [f.name for f in all_fastq]
    file_map = {f.name: f for f in all_fastq}
    detection = detect_sequencing_type(filenames)

    r1_files = []
    r2_files = []
    for sample_info in detection["samples"].values():
        if sample_info.get("R1") and sample_info["R1"] in file_map:
            r1_files.append(file_map[sample_info["R1"]])
        if sample_info.get("R2") and sample_info["R2"] in file_map:
            r2_files.append(file_map[sample_info["R2"]])

    if not r1_files:
        r1_files = all_fastq

    is_paired = sequencing_type == "paired-end" and r2_files

    # For paired-end, DADA2 discards the pair if EITHER direction fails EE≤2.
    # To achieve the desired combined retention, raise the per-direction target:
    #   per_dir_target ≈ sqrt(target) so that per_dir² ≈ target
    per_dir_target = math.sqrt(target_retention) if is_paired else target_retention

    # Sample quality scores from R1 (skip primer bases)
    log.info("Profiling R1 quality scores...")
    r1_scores = _read_quality_scores(r1_files[0], n_reads, skip_bases=primer_offset_f)
    trunc_f = _find_optimal_trunc(r1_scores, target_retention=per_dir_target)
    r1_pass = _pass_rate_at(r1_scores, trunc_f) if trunc_f else 0.0

    trunc_r = 0
    r2_pass = 0.0

    if is_paired:
        log.info("Profiling R2 quality scores...")
        r2_scores = _read_quality_scores(r2_files[0], n_reads, skip_bases=primer_offset_r)
        trunc_r = _find_optimal_trunc(r2_scores, target_retention=per_dir_target)
        r2_pass = _pass_rate_at(r2_scores, trunc_r) if trunc_r else 0.0

        # Enforce overlap constraint.
        # expected_amplicon_length includes primers, but truncation operates
        # on post-cutadapt reads where primers are removed.  Subtract primer
        # lengths to get the insert length that needs overlapping.
        amplicon_len = _expected_amplicon_length(variable_region)
        if amplicon_len:
            primer_info = PRIMERS.get(variable_region, {})
            fwd_len = len(primer_info.get("forward", ""))
            rev_len = len(primer_info.get("reverse", ""))
            insert_len = amplicon_len - fwd_len - rev_len
            log.info(
                f"Amplicon ~{amplicon_len}bp, primers {fwd_len}+{rev_len}bp, "
                f"insert ~{insert_len}bp"
            )
            amplicon_len = insert_len
        if amplicon_len and trunc_f and trunc_r:
            overlap = trunc_f + trunc_r - amplicon_len
            if overlap < min_overlap:
                shortfall = min_overlap - overlap
                log.warning(
                    f"Overlap too small ({overlap}bp < {min_overlap}bp), "
                    f"need {shortfall}bp more"
                )
                # Try increasing trunc_r first (R2 usually the bottleneck)
                # then trunc_f if needed — within limits of read length
                max_r2 = max(len(s) for s in r2_scores) if r2_scores else trunc_r
                max_r1 = max(len(s) for s in r1_scores) if r1_scores else trunc_f
                if trunc_r + shortfall <= max_r2:
                    trunc_r += shortfall
                    r2_pass = _pass_rate_at(r2_scores, trunc_r)
                elif trunc_f + shortfall <= max_r1:
                    trunc_f += shortfall
                    r1_pass = _pass_rate_at(r1_scores, trunc_f)
                else:
                    # Split the increase
                    half = shortfall // 2
                    trunc_f = min(trunc_f + half + shortfall % 2, max_r1)
                    trunc_r = min(trunc_r + half, max_r2)
                    r1_pass = _pass_rate_at(r1_scores, trunc_f)
                    r2_pass = _pass_rate_at(r2_scores, trunc_r)
                overlap = trunc_f + trunc_r - amplicon_len
                log.info(f"Adjusted for overlap: {overlap}bp")

    details_parts = [
        f"trunc_len_f={trunc_f} ({r1_pass:.0%} pass EE≤2)",
    ]
    if is_paired:
        combined_pass = r1_pass * r2_pass
        details_parts.append(f"trunc_len_r={trunc_r} ({r2_pass:.0%} pass EE≤2)")
        details_parts.append(f"~{combined_pass:.0%} combined retention")
        if amplicon_len:
            overlap = trunc_f + trunc_r - amplicon_len
            details_parts.append(f"overlap={overlap}bp")

    details = ", ".join(details_parts)
    log.info(f"Auto-detected: {details}")

    return {
        "trunc_len_f": trunc_f,
        "trunc_len_r": trunc_r,
        "details": details,
    }


def _read_quality_scores(
    fastq_path: Path, n_reads: int = 5000, skip_bases: int = 0
) -> list[list[int]]:
    """Read Phred+33 quality scores from first n_reads of a FASTQ(.gz) file.

    Args:
        skip_bases: Number of bases to skip from the start of each read
                    (e.g., primer length when profiling raw reads).

    Returns:
        List of lists, each inner list is per-base Phred scores for one read.
    """
    scores = []
    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    try:
        with opener(fastq_path, "rt") as f:
            line_num = 0
            for line in f:
                line_num += 1
                if line_num % 4 == 0:  # Quality line (4th line of each record)
                    qual = line.rstrip()
                    phred = [ord(c) - 33 for c in qual]
                    if skip_bases:
                        phred = phred[skip_bases:]
                    if phred:
                        scores.append(phred)
                    if len(scores) >= n_reads:
                        break
    except Exception:
        pass
    return scores


def _find_optimal_trunc(
    quality_scores: list[list[int]],
    max_ee: float = 2.0,
    target_retention: float = 0.70,
) -> int:
    """Find the longest truncation length where >=target_retention of reads pass EE<=max_ee.

    Returns 0 if no suitable length found.
    """
    if not quality_scores:
        return 0

    max_len = max(len(s) for s in quality_scores)

    # Precompute cumulative EE for each read
    cum_ee = []
    for scores in quality_scores:
        ee = []
        running = 0.0
        for q in scores:
            running += 10 ** (-q / 10)
            ee.append(running)
        cum_ee.append(ee)

    # Search from max length down to _MIN_TRUNC
    for trunc_len in range(max_len, _MIN_TRUNC - 1, -1):
        passing = 0
        total = 0
        for ee_list in cum_ee:
            if len(ee_list) >= trunc_len:
                total += 1
                if ee_list[trunc_len - 1] <= max_ee:
                    passing += 1
        if total == 0:
            continue
        rate = passing / total
        if rate >= target_retention:
            return trunc_len

    return _MIN_TRUNC


def _pass_rate_at(quality_scores: list[list[int]], trunc_len: int) -> float:
    """Calculate fraction of reads passing EE<=2 at a given truncation length."""
    if not quality_scores or trunc_len <= 0:
        return 0.0
    passing = 0
    total = 0
    for scores in quality_scores:
        if len(scores) >= trunc_len:
            total += 1
            ee = sum(10 ** (-q / 10) for q in scores[:trunc_len])
            if ee <= 2.0:
                passing += 1
    return passing / total if total > 0 else 0.0


def _expected_amplicon_length(variable_region: str | None) -> int | None:
    """Get the expected amplicon length (midpoint) for a variable region."""
    if not variable_region or variable_region not in PRIMERS:
        return None
    lo, hi = PRIMERS[variable_region]["expected_len"]
    return (lo + hi) // 2
