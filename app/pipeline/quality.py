"""
16S Analyzer — Quality profiling and auto-detection of DADA2 truncation parameters.

Strategy:
  1. Determine insert length from the detected variable region and primer lengths.
  2. Profile per-base quality and truncate from the tail where a 10bp sliding
     window mean drops below Q20.
  3. Ensure truncation lengths satisfy the overlap constraint for PE merging
     (trunc_f + trunc_r >= insert_len + min_overlap). Extend R1 first (better
     quality), then R2 if needed.
  4. For single-end reads, truncate at the quality drop-off only.
"""
import gzip
import logging
from pathlib import Path

from app.pipeline.detect import PRIMERS, detect_sequencing_type

# Absolute floor — never truncate shorter than this
_MIN_TRUNC = 100

# Quality threshold: truncate where per-base mean Phred drops below this
_QUAL_FLOOR = 20


def detect_truncation_params(
    trimmed_dir: Path,
    sequencing_type: str,
    variable_region: str | None,
    min_overlap: int = 50,
    n_reads: int = 5000,
    primer_offset_f: int = 0,
    primer_offset_r: int = 0,
    logger: logging.Logger | None = None,
    platform: str | None = None,
    **_kwargs,
) -> dict:
    """Auto-detect optimal DADA2 truncation lengths.

    For paired-end:
      - Computes insert length from variable region and primers.
      - Finds quality cutoff where 10bp sliding window mean drops below Q20.
      - Ensures trunc_f + trunc_r >= insert_len + min_overlap.
      - Extends R1 first (better quality), then R2 if overlap is insufficient.

    For single-end:
      - Truncates at quality drop-off only.

    Returns dict with: trunc_len_f, trunc_len_r, details (str)
    """
    log = logger or logging.getLogger(__name__)

    # Long-read platforms don't use truncation (they use length filtering instead)
    if platform in ("pacbio", "nanopore"):
        log.info(f"Skipping truncation detection for {platform} long-read data")
        return {"trunc_len_f": 0, "trunc_len_r": 0, "details": f"Not applicable for {platform}"}

    # Find FASTQ files
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

    # Detect leading N bases that DADA2's maxN=0 would reject
    trim_left_f = _detect_leading_n(r1_files[0], n_reads)
    if trim_left_f:
        log.info(f"R1: detected {trim_left_f} leading N bases, will trim")

    # Profile R1 quality
    log.info("Profiling R1 quality scores...")
    r1_scores = _read_quality_scores(r1_files[0], n_reads, skip_bases=primer_offset_f + trim_left_f)
    r1_median_len = _median_read_length(r1_scores)

    if not is_paired:
        trunc_f = _quality_cutoff(r1_scores)
        r1_pass = _pass_rate_at(r1_scores, trunc_f)
        details = f"trunc_len_f={trunc_f} ({r1_pass:.0%} pass EE≤5)"
        if trim_left_f:
            details += f", trim_left_f={trim_left_f} (leading N)"
        log.info(f"Auto-detected: {details}")
        return {"trunc_len_f": trunc_f, "trunc_len_r": 0,
                "trim_left_f": trim_left_f, "trim_left_r": 0,
                "details": details}

    # Paired-end: profile R2 and compute insert length
    trim_left_r = _detect_leading_n(r2_files[0], n_reads)
    if trim_left_r:
        log.info(f"R2: detected {trim_left_r} leading N bases, will trim")

    log.info("Profiling R2 quality scores...")
    r2_scores = _read_quality_scores(r2_files[0], n_reads, skip_bases=primer_offset_r + trim_left_r)
    r2_median_len = _median_read_length(r2_scores)

    insert_len = None
    amplicon_range = _expected_amplicon_range(variable_region)
    if amplicon_range:
        primer_info = PRIMERS.get(variable_region, {})
        fwd_len = len(primer_info.get("forward", ""))
        rev_len = len(primer_info.get("reverse", ""))
        # Use the upper bound so even the longest amplicons get enough overlap
        insert_len = amplicon_range[1] - fwd_len - rev_len
        insert_mid = (amplicon_range[0] + amplicon_range[1]) // 2 - fwd_len - rev_len
        log.info(
            f"Amplicon {amplicon_range[0]}-{amplicon_range[1]}bp, "
            f"primers {fwd_len}+{rev_len}bp, "
            f"insert {insert_mid}-{insert_len}bp (using upper bound)"
        )

    # Quality-based truncation (where 10bp window mean drops below Q20,
    # accounting for truncQ=2 per-read truncation)
    qual_trunc_f = _quality_cutoff(r1_scores)
    qual_trunc_r = _quality_cutoff(r2_scores)

    # Start from quality cutoffs, clamped to read length
    trunc_f = max(_MIN_TRUNC, min(qual_trunc_f, r1_median_len))
    trunc_r = max(_MIN_TRUNC, min(qual_trunc_r, r2_median_len))

    log.info(f"Quality cutoff: R1={qual_trunc_f}, R2={qual_trunc_r}")

    if insert_len:
        # Cap truncation at insert length — reading past the insert means
        # sequencing into adapter/primer, wasting the maxEE error budget.
        # Allow a small buffer (10bp) for insert length variation.
        max_useful = insert_len + 10
        if trunc_f > max_useful:
            log.info(f"Capping R1 trunc from {trunc_f} to {max_useful} (insert={insert_len})")
            trunc_f = max_useful
        if trunc_r > max_useful:
            log.info(f"Capping R2 trunc from {trunc_r} to {max_useful} (insert={insert_len})")
            trunc_r = max_useful

        # Optimize: search for the truncation pair with best combined pass rate
        # while maintaining minimum overlap.
        min_total = insert_len + min_overlap
        best_combined = _pass_rate_at(r1_scores, trunc_f) * _pass_rate_at(r2_scores, trunc_r)
        best_tf, best_tr = trunc_f, trunc_r

        # Try shorter R2 lengths (R2 usually has worse quality)
        for try_tr in range(trunc_r, _MIN_TRUNC - 1, -10):
            # Ensure overlap: adjust R1 up if needed
            try_tf = max(trunc_f, min_total - try_tr)
            if try_tf > r1_median_len:
                continue
            r1p = _pass_rate_at(r1_scores, try_tf)
            r2p = _pass_rate_at(r2_scores, try_tr)
            combined = r1p * r2p
            if combined > best_combined:
                best_combined = combined
                best_tf, best_tr = try_tf, try_tr

        trunc_f, trunc_r = best_tf, best_tr

        # Final overlap check — extend if needed
        if trunc_f + trunc_r < min_total:
            shortfall = min_total - (trunc_f + trunc_r)
            extend_f = max(0, min(shortfall, r1_median_len - trunc_f))
            trunc_f += extend_f
            shortfall -= extend_f
            if shortfall > 0:
                extend_r = max(0, min(shortfall, r2_median_len - trunc_r))
                trunc_r += extend_r

        overlap = trunc_f + trunc_r - insert_len
        log.info(f"Overlap-adjusted: R1={trunc_f}, R2={trunc_r}, overlap={overlap}bp")
    else:
        log.info(f"No region info, using quality cutoff: R1={trunc_f}, R2={trunc_r}")

    r1_pass = _pass_rate_at(r1_scores, trunc_f)
    r2_pass = _pass_rate_at(r2_scores, trunc_r)
    combined = r1_pass * r2_pass

    details_parts = [
        f"trunc_len_f={trunc_f} ({r1_pass:.0%} pass EE≤5)",
        f"trunc_len_r={trunc_r} ({r2_pass:.0%} pass EE≤5)",
        f"~{combined:.0%} combined retention",
    ]
    if insert_len:
        overlap = trunc_f + trunc_r - insert_len
        details_parts.append(f"overlap={overlap}bp")
    if trim_left_f:
        details_parts.append(f"trim_left_f={trim_left_f} (leading N)")
    if trim_left_r:
        details_parts.append(f"trim_left_r={trim_left_r} (leading N)")

    details = ", ".join(details_parts)
    log.info(f"Auto-detected: {details}")

    return {
        "trunc_len_f": trunc_f,
        "trunc_len_r": trunc_r,
        "trim_left_f": trim_left_f,
        "trim_left_r": trim_left_r,
        "details": details,
    }


def _detect_leading_n(fastq_path: Path, n_reads: int = 500, check_first: int = 10) -> int:
    """Detect N bases in the first few positions that would cause DADA2 maxN=0 to reject reads.

    Checks the first `check_first` positions. If >50% of reads have an N
    at any position within that window, returns trim length to remove it.

    Returns the number of leading bases to trim (0 if no N issue).
    """
    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    # Count N occurrences at each of the first `check_first` positions
    n_at_pos = [0] * check_first
    total = 0
    try:
        with opener(fastq_path, "rt") as f:
            for line_num, line in enumerate(f, 1):
                if line_num % 4 == 2:  # sequence line
                    seq = line.rstrip()
                    for i in range(min(check_first, len(seq))):
                        if seq[i] == "N":
                            n_at_pos[i] += 1
                    total += 1
                    if total >= n_reads:
                        break
    except Exception:
        return 0

    if total == 0:
        return 0

    # Find the rightmost position where >50% of reads have N
    trim_to = 0
    for i in range(check_first):
        if n_at_pos[i] / total >= 0.5:
            trim_to = i + 1  # trim up to and including this position

    return trim_to


def _read_quality_scores(
    fastq_path: Path, n_reads: int = 5000, skip_bases: int = 0
) -> list[list[int]]:
    """Read Phred+33 quality scores sampled evenly across a FASTQ(.gz) file.

    Reads the entire file to count total reads, then samples n_reads evenly
    spaced across the file. This avoids bias from Illumina's tendency to have
    better quality at the start of a run.
    """
    opener = gzip.open if str(fastq_path).endswith(".gz") else open

    # First pass: collect ALL quality lines
    all_quals = []
    try:
        with opener(fastq_path, "rt") as f:
            line_num = 0
            for line in f:
                line_num += 1
                if line_num % 4 == 0:
                    all_quals.append(line.rstrip())
    except Exception:
        pass

    if not all_quals:
        return []

    # Subsample evenly across the file
    total = len(all_quals)
    if total <= n_reads:
        indices = range(total)
    else:
        step = total / n_reads
        indices = [int(i * step) for i in range(n_reads)]

    scores = []
    for idx in indices:
        qual = all_quals[idx]
        phred = [ord(c) - 33 for c in qual]
        if skip_bases:
            phred = phred[skip_bases:]
        if phred:
            scores.append(phred)

    return scores


def _median_read_length(quality_scores: list[list[int]]) -> int:
    """Return the median read length from quality score lists."""
    if not quality_scores:
        return 0
    lengths = sorted(len(s) for s in quality_scores)
    return lengths[len(lengths) // 2]


def _quality_cutoff(quality_scores: list[list[int]], q_floor: int = _QUAL_FLOOR) -> int:
    """Find the position where per-base mean quality drops below q_floor.

    Scans from the end of the read backward. Returns the last position
    where a sliding window of 10bp has mean quality >= q_floor.
    """
    if not quality_scores:
        return _MIN_TRUNC

    # Compute per-position mean quality
    max_len = max(len(s) for s in quality_scores)
    mean_qual = []
    for pos in range(max_len):
        quals = [s[pos] for s in quality_scores if pos < len(s)]
        if quals:
            mean_qual.append(sum(quals) / len(quals))
        else:
            mean_qual.append(0)

    # Use a 10bp sliding window from the end to find the quality drop-off
    window = 10
    cutoff = max_len
    for pos in range(max_len - 1, window - 2, -1):
        start = max(0, pos - window + 1)
        win_mean = sum(mean_qual[start:pos + 1]) / (pos + 1 - start)
        if win_mean >= q_floor:
            cutoff = pos + 1  # truncate at this length (1-indexed)
            break
    else:
        cutoff = _MIN_TRUNC

    cutoff = max(_MIN_TRUNC, cutoff)

    return cutoff


def _pass_rate_at(quality_scores: list[list[int]], trunc_len: int, max_ee: float = 5.0) -> float:
    """Calculate fraction of ALL reads passing EE<=max_ee at a given truncation length.

    Reads shorter than trunc_len count as failures (DADA2 discards them).
    Default max_ee=5.0 matches DADA2's maxEE=c(5,5).
    Uses truncQ=0 (disabled) to match DADA2 settings.
    """
    if not quality_scores or trunc_len <= 0:
        return 0.0
    passing = 0
    for scores in quality_scores:
        if len(scores) >= trunc_len:
            ee = sum(10 ** (-q / 10) for q in scores[:trunc_len])
            if ee <= max_ee:
                passing += 1
    return passing / len(quality_scores)


def _expected_amplicon_range(variable_region: str | None) -> tuple[int, int] | None:
    """Get the expected amplicon length range (lo, hi) for a variable region."""
    if not variable_region or variable_region not in PRIMERS:
        return None
    return PRIMERS[variable_region]["expected_len"]
