"""
MicrobiomeDash — Auto-detection of sequencing type (SE/PE) and variable region.
"""
import gzip
import logging
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

# ── Filename patterns for paired-end detection ────────────────────────────────

# Each tuple: (R1 regex, R2 regex) applied to filenames
_PE_PATTERNS = [
    # Illumina standard: sample_R1_001.fastq.gz / sample_R2_001.fastq.gz
    (r"^(.+?)_R1(_\d+)?\.fastq\.gz$", r"^(.+?)_R2(_\d+)?\.fastq\.gz$"),
    (r"^(.+?)_R1(_\d+)?\.fq\.gz$", r"^(.+?)_R2(_\d+)?\.fq\.gz$"),
    # Simpler: sample_R1.fastq.gz / sample_R2.fastq.gz
    (r"^(.+?)_R1\.fastq\.gz$", r"^(.+?)_R2\.fastq\.gz$"),
    (r"^(.+?)_R1\.fq\.gz$", r"^(.+?)_R2\.fq\.gz$"),
    # Numeric: sample_1.fastq.gz / sample_2.fastq.gz
    (r"^(.+?)_1\.fastq\.gz$", r"^(.+?)_2\.fastq\.gz$"),
    (r"^(.+?)_1\.fq\.gz$", r"^(.+?)_2\.fq\.gz$"),
]


def extract_sample_name(filename: str) -> str:
    """Extract sample name from a FASTQ filename by stripping read direction and extension."""
    name = filename
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    # Strip read direction suffixes — try R1/R2 first, fall back to _1/_2
    stripped = re.sub(r"_R[12](_\d+)?$", "", name)
    if stripped != name:
        return stripped
    return re.sub(r"_[12]$", "", name)


def detect_sequencing_type(filenames: list[str]) -> dict:
    """Detect single-end vs paired-end from a list of FASTQ filenames.

    Returns:
        dict with keys:
            type: 'single-end' | 'paired-end' | 'mixed'
            samples: dict mapping sample_name -> {'R1': filename, 'R2': filename | None}
            errors: list of issues found
    """
    if not filenames:
        return {"type": "unknown", "samples": {}, "errors": ["No files provided."]}

    # Try to match R1/R2 pairs
    r1_files = {}
    r2_files = {}
    unmatched = []

    for fn in filenames:
        matched = False
        for r1_pat, r2_pat in _PE_PATTERNS:
            m1 = re.match(r1_pat, fn)
            if m1:
                sample = m1.group(1)
                r1_files[sample] = fn
                matched = True
                break
            m2 = re.match(r2_pat, fn)
            if m2:
                sample = m2.group(1)
                r2_files[sample] = fn
                matched = True
                break
        if not matched:
            unmatched.append(fn)

    # Build sample map
    all_samples = sorted(set(list(r1_files.keys()) + list(r2_files.keys())))
    samples = {}
    errors = []

    if all_samples:
        # We found R1/R2 patterns
        for sample in all_samples:
            r1 = r1_files.get(sample)
            r2 = r2_files.get(sample)
            if r1 and r2:
                samples[sample] = {"R1": r1, "R2": r2}
            elif r1 and not r2:
                errors.append(f"Sample '{sample}' has R1 but no R2: {r1}")
                samples[sample] = {"R1": r1, "R2": None}
            else:
                errors.append(f"Sample '{sample}' has R2 but no R1: {r2}")
                samples[sample] = {"R1": None, "R2": r2}

        if unmatched:
            errors.append(
                f"{len(unmatched)} file(s) don't match paired-end patterns: "
                f"{unmatched[:3]}{'...' if len(unmatched) > 3 else ''}"
            )

        # Determine overall type
        has_pairs = any(s["R1"] and s["R2"] for s in samples.values())
        has_singles = bool(unmatched) or any(
            (s["R1"] is None) != (s["R2"] is None) for s in samples.values()
        )
        if has_pairs and not has_singles:
            seq_type = "paired-end"
        elif has_singles and not has_pairs:
            seq_type = "single-end"
        else:
            seq_type = "mixed"
    else:
        # No R1/R2 patterns found → single-end
        seq_type = "single-end"
        for fn in filenames:
            sample = extract_sample_name(fn)
            samples[sample] = {"R1": fn, "R2": None}

    return {"type": seq_type, "samples": samples, "errors": errors}


# ── Variable region detection ─────────────────────────────────────────────────

# Known 16S primer sequences (IUPAC degenerate bases)
PRIMERS = {
    "V1-V2": {
        "forward": "AGAGTTTGATCMTGGCTCAG",    # 27F
        "reverse": "TGCTGCCTCCCGTAGGAGT",      # 338R
        "expected_len": (300, 350),
    },
    "V1-V3": {
        "forward": "AGAGTTTGATCMTGGCTCAG",    # 27F
        "reverse": "ATTACCGCGGCTGCTGG",         # 534R
        "expected_len": (490, 540),
    },
    "V3-V4": {
        "forward": "CCTACGGGNGGCWGCAG",        # 341F
        "reverse": "GACTACHVGGGTATCTAATCC",     # 806R
        "expected_len": (420, 480),
    },
    "V4": {
        "forward": "GTGYCAGCMGCCGCGGTAA",      # 515F
        "reverse": "GGACTACNVGGGTWTCTAAT",      # 806R
        "expected_len": (240, 260),
    },
    "V4-V5": {
        "forward": "GTGYCAGCMGCCGCGGTAA",      # 515F
        "reverse": "CCGYCAATTYMTTTRAGTTT",      # 926R
        "expected_len": (370, 430),
    },
    "V5-V6": {
        "forward": "AYTGGGYDTAAAGNG",          # 784F
        "reverse": "CCGTCAATTCMTTTRAGT",        # 1061R
        "expected_len": (280, 350),
    },
    "V1-V9": {
        "forward": "AGRGTTYGATYMTGGCTCAG",    # 27F (degenerate)
        "reverse": "RGYTACCTTGTTACGACTT",      # 1492R (degenerate)
        "expected_len": (1400, 1500),
    },
}

# IUPAC degenerate base lookup
_IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
    "H": "ACT", "V": "ACG", "N": "ACGT",
}


def _primer_matches(sequence: str, primer: str, max_mismatches: int = 3) -> bool:
    """Check if a sequence starts with a primer, allowing degenerate bases and mismatches."""
    if len(sequence) < len(primer):
        return False
    mismatches = 0
    for seq_base, primer_base in zip(sequence.upper(), primer.upper()):
        allowed = _IUPAC.get(primer_base, primer_base)
        if seq_base not in allowed:
            mismatches += 1
            if mismatches > max_mismatches:
                return False
    return True


def _read_fastq_sequences(fastq_path: Path, n_reads: int = 100) -> list[str]:
    """Read the first n sequences from a (gzipped) FASTQ file."""
    sequences = []
    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    try:
        with opener(fastq_path, "rt") as f:
            line_num = 0
            for line in f:
                line_num += 1
                if line_num % 4 == 2:  # Sequence line
                    sequences.append(line.strip())
                    if len(sequences) >= n_reads:
                        break
    except Exception:
        pass
    return sequences


def detect_variable_region(
    fastq_path: Path, n_reads: int = 100
) -> dict:
    """Detect the 16S variable region from primer sequences in FASTQ reads.

    Args:
        fastq_path: Path to a FASTQ(.gz) file.
        n_reads: Number of reads to check.

    Returns:
        dict with keys:
            region (str | None), confidence (float 0-1),
            method (str), details (str)
    """
    sequences = _read_fastq_sequences(fastq_path, n_reads)
    if not sequences:
        return {
            "region": None,
            "confidence": 0.0,
            "method": "none",
            "details": "Could not read sequences from file.",
        }

    # Count primer matches for each region
    region_scores = defaultdict(int)
    for seq in sequences:
        for region_name, primers in PRIMERS.items():
            if _primer_matches(seq, primers["forward"]):
                region_scores[region_name] += 1

    if not region_scores:
        # Fallback: align reads to E. coli 16S with bbmap
        bbmap_result = _detect_region_bbmap(fastq_path)
        if bbmap_result["region"]:
            return bbmap_result

        return {
            "region": None,
            "confidence": 0.0,
            "method": "none",
            "details": "No primer matches found and bbmap alignment inconclusive.",
        }

    # Pick the region with the most matches
    best_region = max(region_scores, key=region_scores.get)
    best_count = region_scores[best_region]
    confidence = best_count / len(sequences)

    return {
        "region": best_region,
        "confidence": round(confidence, 2),
        "method": "primer_match",
        "details": f"{best_count}/{len(sequences)} reads matched {best_region} forward primer.",
    }


# ── BBMap-based variable region detection ─────────────────────────────────────

# E. coli 16S rRNA coordinate ranges for each variable region (J01859-based, 1-indexed)
_REGION_COORDS = {
    "V1-V2": (69, 337),
    "V3-V4": (341, 805),
    "V4":    (515, 806),
    "V4-V5": (515, 926),
    "V5-V6": (784, 1061),
    "V1-V9": (69, 1492),
}


def _detect_region_bbmap(fastq_path: Path, n_reads: int = 500) -> dict:
    """Detect variable region by aligning reads to E. coli 16S with bbmap.

    Subsamples reads, aligns with bbmap, parses alignment positions,
    and determines which variable region the reads cover.
    """
    from app.config import REFERENCE_DIR, conda_cmd

    ref_path = REFERENCE_DIR / "ecoli_16S.fasta"
    if not ref_path.exists():
        return {
            "region": None, "confidence": 0.0,
            "method": "bbmap", "details": "E. coli 16S reference not found.",
        }

    logger = logging.getLogger("detect")

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Subsample reads to a temp file
            subset_path = tmpdir / "subset.fastq"
            _subsample_fastq(fastq_path, subset_path, n_reads)

            # Run bbmap
            sam_path = tmpdir / "mapped.sam"
            cmd = conda_cmd([
                "bbmap.sh",
                f"ref={ref_path}",
                f"in={subset_path}",
                f"out={sam_path}",
                "nodisk=t",
                "ambiguous=best",
                "minid=0.5",
                "threads=2",
                "overwrite=t",
                "-Xmx1g",
            ])

            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=120,
            )

            if result.returncode != 0:
                logger.debug(f"bbmap failed: {result.stderr[:200]}")
                return {
                    "region": None, "confidence": 0.0,
                    "method": "bbmap", "details": f"bbmap alignment failed.",
                }

            # Parse SAM to get alignment start/end positions
            positions = _parse_sam_positions(sam_path)

            if len(positions) < 10:
                return {
                    "region": None, "confidence": 0.0,
                    "method": "bbmap",
                    "details": f"Too few reads aligned ({len(positions)}).",
                }

            # Determine variable region from alignment coordinates
            return _coords_to_region(positions)

    except Exception as e:
        logger.debug(f"bbmap detection failed: {e}")
        return {
            "region": None, "confidence": 0.0,
            "method": "bbmap", "details": f"bbmap detection error: {e}",
        }


def _subsample_fastq(fastq_path: Path, output_path: Path, n_reads: int):
    """Extract the first n_reads from a FASTQ file."""
    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    count = 0
    with opener(fastq_path, "rt") as fin, open(output_path, "w") as fout:
        for line in fin:
            fout.write(line)
            count += 1
            if count >= n_reads * 4:  # 4 lines per FASTQ record
                break


def _parse_sam_positions(sam_path: Path) -> list[tuple[int, int]]:
    """Parse a SAM file and return (start, end) positions of mapped reads."""
    positions = []
    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.split("\t")
            if len(fields) < 10:
                continue
            flag = int(fields[1])
            if flag & 4:  # unmapped
                continue
            pos = int(fields[3])  # 1-based leftmost position
            seq = fields[9]
            # Calculate end position from CIGAR
            cigar = fields[5]
            ref_len = _cigar_ref_length(cigar)
            end = pos + ref_len - 1
            positions.append((pos, end))
    return positions


def _cigar_ref_length(cigar: str) -> int:
    """Calculate reference-consuming length from a CIGAR string."""
    length = 0
    num = ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            if ch in "MDN=X":  # reference-consuming operations
                length += int(num)
            num = ""
    return length if length > 0 else 1


def _coords_to_region(positions: list[tuple[int, int]]) -> dict:
    """Determine variable region from alignment coordinates."""
    # Get median start and end
    starts = sorted(p[0] for p in positions)
    ends = sorted(p[1] for p in positions)
    median_start = starts[len(starts) // 2]
    median_end = ends[len(ends) // 2]

    # Score each region by how well the reads overlap it.
    # Prefer the smallest (most specific) region with high coverage.
    tolerance = 50
    candidates = []
    for region_name, (reg_start, reg_end) in _REGION_COORDS.items():
        hits = sum(
            1 for s, e in positions
            if s >= reg_start - tolerance and e <= reg_end + tolerance
        )
        score = hits / len(positions)
        region_span = reg_end - reg_start
        if score >= 0.3:
            candidates.append((region_name, score, region_span))

    # Among regions with >= 90% of the best score, pick the smallest span
    best_region = None
    best_score = 0.0
    if candidates:
        top_score = max(c[1] for c in candidates)
        good = [c for c in candidates if c[1] >= top_score * 0.9]
        best = min(good, key=lambda c: c[2])
        best_region = best[0]
        best_score = best[1]

    if best_region:
        return {
            "region": best_region,
            "confidence": round(best_score, 2),
            "method": "bbmap_alignment",
            "details": (
                f"Reads align to E. coli 16S positions {median_start}-{median_end}, "
                f"matching {best_region} ({best_score:.0%} of reads)."
            ),
        }

    return {
        "region": None,
        "confidence": 0.0,
        "method": "bbmap_alignment",
        "details": (
            f"Reads align to E. coli 16S positions {median_start}-{median_end}, "
            f"but no clear variable region match."
        ),
    }


# ── Platform detection (Illumina / PacBio / Nanopore) ─────────────────────────


def _read_fastq_qualities(fastq_path: Path, n_reads: int = 500) -> list[list[int]]:
    """Read Phred+33 quality scores from the first n_reads of a FASTQ(.gz) file."""
    qualities = []
    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    try:
        with opener(fastq_path, "rt") as f:
            line_num = 0
            for line in f:
                line_num += 1
                if line_num % 4 == 0:  # Quality line
                    qualities.append([ord(c) - 33 for c in line.strip()])
                    if len(qualities) >= n_reads:
                        break
    except Exception:
        pass
    return qualities


def detect_platform(fastq_path: Path, n_reads: int = 500) -> dict:
    """Detect sequencing platform from read lengths and quality scores.

    Returns:
        dict with keys:
            platform: "illumina" | "pacbio" | "nanopore"
            median_len: int
            median_qual: float
            details: str
    """
    sequences = _read_fastq_sequences(fastq_path, n_reads)
    qualities = _read_fastq_qualities(fastq_path, n_reads)

    if not sequences:
        return {
            "platform": "illumina",
            "median_len": 0,
            "median_qual": 0.0,
            "details": "Could not read sequences; defaulting to Illumina.",
        }

    lengths = sorted(len(s) for s in sequences)
    median_len = lengths[len(lengths) // 2]

    # Compute median per-read mean quality
    if qualities:
        mean_quals = sorted(sum(q) / len(q) for q in qualities if q)
        median_qual = mean_quals[len(mean_quals) // 2] if mean_quals else 0.0
    else:
        median_qual = 0.0

    if median_len > 1000:
        # Long read — distinguish PacBio HiFi (high quality) from Nanopore
        if median_qual >= 25:
            platform = "pacbio"
            details = (
                f"Long reads detected (median {median_len}bp, Q{median_qual:.0f}). "
                f"High quality suggests PacBio HiFi."
            )
        else:
            platform = "nanopore"
            details = (
                f"Long reads detected (median {median_len}bp, Q{median_qual:.0f}). "
                f"Lower quality suggests Nanopore."
            )
    else:
        platform = "illumina"
        details = f"Short reads detected (median {median_len}bp). Illumina platform."

    return {
        "platform": platform,
        "median_len": median_len,
        "median_qual": round(median_qual, 1),
        "details": details,
    }
