"""
MicrobiomeDash — BIOM operations: region extraction, detection, and combining.

Provides tools for:
- Auto-detecting variable region from BIOM observation metadata
- Extracting sub-regions (e.g., V4 from V3-V4) via internal primer search
- Combining multiple BIOM files by sequence (same region) or taxonomy (cross-region)
"""
import io
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from biom import Table, load_table

from app.pipeline.detect import PRIMERS, _IUPAC

# ── Constants ────────────────────────────────────────────────────────────────

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

TAX_LEVELS = {
    "Kingdom": 0,
    "Phylum": 1,
    "Class": 2,
    "Order": 3,
    "Family": 4,
    "Genus": 5,
    "Species": 6,
}

# E. coli 16S rRNA positions (0-based on J01859) for each amplicon region.
# These are the spans that amplicon sequences typically cover, including the
# primer binding sites at the boundaries.
REGION_ECOLI_BOUNDS: dict[str, tuple[int, int]] = {
    # Standard amplicon regions
    "V1-V2": (8, 338),
    "V1-V3": (8, 534),
    "V3-V4": (341, 806),
    "V4": (515, 806),
    "V4-V5": (515, 926),
    "V5-V6": (784, 1061),
    "V1-V9": (8, 1492),
    # Sub-regions (from extraction or single-region amplicons)
    "V3": (341, 534),
    "V5": (806, 926),
    "V6": (926, 1061),
}

_ECOLI_16S: str | None = None  # lazy-loaded


def _load_ecoli_ref() -> str:
    """Load the E. coli 16S reference, caching on first call."""
    global _ECOLI_16S
    if _ECOLI_16S is not None:
        return _ECOLI_16S
    fasta_path = Path(__file__).resolve().parent.parent.parent / "data" / "references" / "ecoli_16S.fasta"
    lines = fasta_path.read_text().strip().splitlines()
    _ECOLI_16S = "".join(l.strip() for l in lines if not l.startswith(">")).upper()
    return _ECOLI_16S


# ── Helpers ──────────────────────────────────────────────────────────────────


def _reverse_complement_str(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.upper().translate(_COMPLEMENT)[::-1]


def _find_primer_internal(
    sequence: str, primer: str, max_mismatches: int = 2
) -> tuple[int, int] | None:
    """Scan *sequence* for an internal match to *primer* (IUPAC-aware).

    Searches all positions (not just position 0) for the primer, allowing
    degenerate bases and up to *max_mismatches*.

    Returns (start, end) of the match or None if not found.
    """
    seq = sequence.upper()
    plen = len(primer)
    if len(seq) < plen:
        return None

    primer = primer.upper()
    for i in range(len(seq) - plen + 1):
        mismatches = 0
        ok = True
        for j in range(plen):
            allowed = _IUPAC.get(primer[j], primer[j])
            if seq[i + j] not in allowed:
                mismatches += 1
                if mismatches > max_mismatches:
                    ok = False
                    break
        if ok:
            return (i, i + plen)
    return None


def _align_to_ecoli(sequence: str, ref: str, kmer_size: int = 15) -> tuple[int, int] | None:
    """Find the approximate position of *sequence* on the E. coli 16S reference.

    Uses kmer seed matching: builds an index of all kmers in the reference,
    finds matching kmers from the query, and picks the consensus offset.
    Returns (start, end) positions on the reference, or None if no alignment.
    """
    seq = sequence.upper()
    if len(seq) < kmer_size:
        return None

    # Build kmer index of reference (done per-call; ref is ~1.5 kb so this is fast)
    ref_kmers: dict[str, list[int]] = {}
    for i in range(len(ref) - kmer_size + 1):
        kmer = ref[i : i + kmer_size]
        if "N" in kmer:
            continue
        ref_kmers.setdefault(kmer, []).append(i)

    # Collect offset votes
    offsets: list[int] = []
    for i in range(len(seq) - kmer_size + 1):
        kmer = seq[i : i + kmer_size]
        if kmer in ref_kmers:
            for ref_pos in ref_kmers[kmer]:
                offsets.append(ref_pos - i)

    if not offsets:
        return None

    best_offset = Counter(offsets).most_common(1)[0][0]
    start = best_offset
    end = best_offset + len(seq)
    return (start, end)


def _check_primers_in_seqs(sequences: list[str], max_check: int = 50) -> bool:
    """Check whether ASV sequences still contain primer sequences at their ends.

    Returns True if primers are detected (i.e. NOT trimmed).
    """
    check = sequences[:max_check]
    # Check forward primers at the start of sequences
    for region_name, info in PRIMERS.items():
        fwd = info["forward"]
        matches = 0
        for seq in check:
            if len(seq) < len(fwd):
                continue
            mismatches = 0
            ok = True
            for s_base, p_base in zip(seq.upper(), fwd.upper()):
                allowed = _IUPAC.get(p_base, p_base)
                if s_base not in allowed:
                    mismatches += 1
                    if mismatches > 3:
                        ok = False
                        break
            if ok:
                matches += 1
        if matches / len(check) >= 0.3:
            return True
    return False


def _biom_to_bytes(table: Table) -> bytes:
    """Serialize a biom.Table to HDF5 bytes (same pattern as biom_convert.py)."""
    import h5py

    buf = io.BytesIO()
    with h5py.File(buf, "w") as f:
        table.to_hdf5(f, generated_by="MicrobiomeDash-DataManager")
    return buf.getvalue()


def _load_biom(biom_path: str) -> Table:
    """Load a BIOM table from a file path."""
    return load_table(str(biom_path))


# ── Extraction rules ──────────────────────────────────────────────────────────
#
# Boundary primers in forward orientation at the junctions between regions.
# In a forward-oriented amplicon:
#   - Forward primers appear as-is at region starts
#   - Reverse primers appear as reverse-complements at region ends

_CUT_V2_V3 = PRIMERS["V3-V4"]["forward"]                       # 341F
_CUT_V3_V4 = PRIMERS["V4"]["forward"]                          # 515F
_CUT_V4_V5 = _reverse_complement_str(PRIMERS["V4"]["reverse"]) # RC(806R)
_CUT_V5_V6 = _reverse_complement_str(PRIMERS["V4-V5"]["reverse"])  # RC(926R)

# (source_region, target_region) -> (left_primer | None, right_primer | None)
# left_primer:  find internally, keep sequence AFTER it
# right_primer: find internally, keep sequence BEFORE it
# Both None is invalid; both set = two-cut (keep between)
EXTRACTION_RULES: dict[tuple[str, str], tuple[str | None, str | None]] = {
    # From V1-V3
    ("V1-V3", "V1-V2"): (None, _CUT_V2_V3),
    ("V1-V3", "V3"):    (_CUT_V2_V3, None),
    # From V3-V4
    ("V3-V4", "V3"): (None, _CUT_V3_V4),
    ("V3-V4", "V4"): (_CUT_V3_V4, None),
    # From V4-V5
    ("V4-V5", "V4"): (None, _CUT_V4_V5),
    ("V4-V5", "V5"): (_CUT_V4_V5, None),
    # From V5-V6
    ("V5-V6", "V5"): (None, _CUT_V5_V6),
    ("V5-V6", "V6"): (_CUT_V5_V6, None),
    # From V1-V9 (full-length) — single-cut targets
    ("V1-V9", "V1-V2"): (None, _CUT_V2_V3),
    ("V1-V9", "V1-V3"): (None, _CUT_V3_V4),
    # From V1-V9 (full-length) — two-cut targets
    ("V1-V9", "V3"):    (_CUT_V2_V3, _CUT_V3_V4),
    ("V1-V9", "V3-V4"): (_CUT_V2_V3, _CUT_V4_V5),
    ("V1-V9", "V4"):    (_CUT_V3_V4, _CUT_V4_V5),
    ("V1-V9", "V4-V5"): (_CUT_V3_V4, _CUT_V5_V6),
    ("V1-V9", "V5"):    (_CUT_V4_V5, _CUT_V5_V6),
    ("V1-V9", "V5-V6"): (_CUT_V4_V5, None),
}


# ── Public API ───────────────────────────────────────────────────────────────


def get_valid_extractions(source_region: str) -> list[str]:
    """Return available target regions for a given source region."""
    return [
        target
        for (src, target) in EXTRACTION_RULES
        if src == source_region
    ]


def detect_region_from_biom(biom_path: str) -> dict:
    """Auto-detect the variable region from a BIOM file.

    Detection strategy:
    1. If the BIOM table_id contains a region tag (set by extract_region), use it.
    2. Otherwise, align a sample of ASV sequences to the E. coli 16S reference
       and determine which variable region they cover based on alignment coordinates.
    3. Also checks whether primer sequences are still present in the ASVs.

    Returns dict with: region, confidence, n_asvs, n_samples, median_len, primers_present
    """
    table = _load_biom(biom_path)
    obs_ids = table.ids(axis="observation")
    sample_ids = table.ids(axis="sample")

    # Extract sequences from observation metadata
    sequences = []
    for obs_id in obs_ids:
        md = table.metadata(obs_id, axis="observation")
        if md and "sequence" in md:
            sequences.append(md["sequence"])

    median_len = None
    if sequences:
        lengths = [len(s) for s in sequences]
        median_len = sorted(lengths)[len(lengths) // 2]

    # Check for primer contamination
    primers_present = _check_primers_in_seqs(sequences) if sequences else False

    base_result = {
        "n_asvs": len(obs_ids),
        "n_samples": len(sample_ids),
        "median_len": median_len,
        "primers_present": primers_present,
    }

    if not sequences:
        return {**base_result, "region": None, "confidence": 0.0}

    # Align a sample of ASV sequences to E. coli 16S reference
    ref = _load_ecoli_ref()
    check_seqs = sequences[:50]
    aligned_starts: list[int] = []
    aligned_ends: list[int] = []

    for seq in check_seqs:
        hit = _align_to_ecoli(seq, ref)
        if hit:
            aligned_starts.append(hit[0])
            aligned_ends.append(hit[1])

    if not aligned_starts:
        return {**base_result, "region": None, "confidence": 0.0}

    # Consensus alignment span
    med_start = sorted(aligned_starts)[len(aligned_starts) // 2]
    med_end = sorted(aligned_ends)[len(aligned_ends) // 2]

    # Score each region by overlap (Jaccard-like: intersection / union)
    best_region = None
    best_score = 0.0
    for region_name, (rstart, rend) in REGION_ECOLI_BOUNDS.items():
        overlap = max(0, min(med_end, rend) - max(med_start, rstart))
        union = max(med_end, rend) - min(med_start, rstart)
        if union > 0:
            score = overlap / union
            if score > best_score:
                best_score = score
                best_region = region_name

    return {
        **base_result,
        "region": best_region,
        "confidence": round(best_score, 2),
    }


def extract_region(
    biom_path: str, source_region: str, target_region: str
) -> dict:
    """Extract a sub-region from a BIOM file.

    Finds boundary primer(s) internally in each ASV sequence, trims to
    the target region, collapses identical trimmed sequences, and sums
    their abundances.

    Supports single-cut (one boundary primer) and two-cut (left + right
    boundary primers) extraction.

    Returns dict with: biom_bytes, n_input, n_output, n_failed, n_collapsed, target_region
    """
    key = (source_region, target_region)
    if key not in EXTRACTION_RULES:
        raise ValueError(
            f"No extraction rule for {source_region} -> {target_region}. "
            f"Valid: {list(EXTRACTION_RULES.keys())}"
        )

    left_primer, right_primer = EXTRACTION_RULES[key]
    table = _load_biom(biom_path)
    obs_ids = table.ids(axis="observation")
    sample_ids = table.ids(axis="sample")
    data = table.matrix_data.toarray()  # obs x sample

    # Trim each ASV
    trimmed: dict[str, list[tuple[int, dict | None]]] = defaultdict(list)
    n_failed = 0

    for idx, obs_id in enumerate(obs_ids):
        md = table.metadata(obs_id, axis="observation")
        if not md or "sequence" not in md:
            n_failed += 1
            continue

        seq = md["sequence"]

        # Left cut: find primer, keep everything AFTER it
        if left_primer:
            hit = _find_primer_internal(seq, left_primer, max_mismatches=2)
            if hit is None:
                n_failed += 1
                continue
            seq = seq[hit[1]:]

        # Right cut: find primer, keep everything BEFORE it
        if right_primer:
            hit = _find_primer_internal(seq, right_primer, max_mismatches=2)
            if hit is None:
                n_failed += 1
                continue
            seq = seq[:hit[0]]

        if len(seq) < 50:  # Too short after trimming
            n_failed += 1
            continue

        trimmed[seq].append((idx, md))

    if not trimmed:
        raise ValueError(
            f"No ASVs could be trimmed for {source_region} -> {target_region}. "
            f"Boundary primer(s) not found internally in any of {len(obs_ids)} ASV sequences."
        )

    # Build new table: collapse identical trimmed sequences
    new_obs_ids = []
    new_obs_metadata = []
    new_data = []
    n_collapsed = 0

    for i, (seq, group) in enumerate(trimmed.items()):
        asv_id = f"ASV_{i+1}"
        new_obs_ids.append(asv_id)

        # Sum abundances across all ASVs mapping to this trimmed sequence
        row = np.zeros(len(sample_ids), dtype=np.float64)
        for orig_idx, _ in group:
            row += data[orig_idx]
        new_data.append(row)

        # Carry forward taxonomy from first ASV
        first_md = group[0][1]
        md_new: dict = {"sequence": seq}
        if first_md and "taxonomy" in first_md:
            md_new["taxonomy"] = first_md["taxonomy"]
        new_obs_metadata.append(md_new)

        if len(group) > 1:
            n_collapsed += len(group) - 1

    new_table = Table(
        np.array(new_data),
        observation_ids=new_obs_ids,
        sample_ids=sample_ids,
        observation_metadata=new_obs_metadata,
        type="OTU table",
    )

    return {
        "biom_bytes": _biom_to_bytes(new_table),
        "n_input": len(obs_ids),
        "n_output": len(new_obs_ids),
        "n_failed": n_failed,
        "n_collapsed": n_collapsed,
        "target_region": target_region,
    }


def combine_biom_same_region(
    biom_paths: list[str], source_names: list[str]
) -> dict:
    """Combine multiple BIOM files from the same region by merging identical sequences.

    Sample names are prefixed with source name if collisions exist.

    Returns dict with: biom_bytes, total_asvs, unique_asvs, total_samples, collisions
    """
    tables = [_load_biom(p) for p in biom_paths]

    # Detect sample name collisions
    all_sample_ids = []
    for t in tables:
        all_sample_ids.extend(t.ids(axis="sample"))
    has_collisions = len(all_sample_ids) != len(set(all_sample_ids))

    # Collect all sequences -> {seq: [(table_idx, obs_id, metadata)]}
    seq_groups: dict[str, list[tuple[int, str, dict | None]]] = defaultdict(list)
    total_asvs = 0

    for t_idx, table in enumerate(tables):
        for obs_id in table.ids(axis="observation"):
            md = table.metadata(obs_id, axis="observation")
            seq = md.get("sequence", obs_id) if md else obs_id
            seq_groups[seq].append((t_idx, obs_id, md))
            total_asvs += 1

    # Build unified sample list (prefix if collisions)
    unified_samples = []
    sample_map: list[dict[str, int]] = []

    for t_idx, table in enumerate(tables):
        mapping = {}
        for sid in table.ids(axis="sample"):
            new_name = f"{source_names[t_idx]}_{sid}" if has_collisions else sid
            mapping[sid] = len(unified_samples)
            unified_samples.append(new_name)
        sample_map.append(mapping)

    n_unified = len(unified_samples)

    # Build merged matrix
    new_obs_ids = []
    new_obs_metadata = []
    new_data = []

    for i, (seq, group) in enumerate(seq_groups.items()):
        asv_id = f"ASV_{i+1}"
        new_obs_ids.append(asv_id)

        row = np.zeros(n_unified, dtype=np.float64)
        for t_idx, obs_id, _ in group:
            table = tables[t_idx]
            obs_data = table.data(obs_id, axis="observation", dense=True)
            for j, sid in enumerate(table.ids(axis="sample")):
                row[sample_map[t_idx][sid]] += obs_data[j]
        new_data.append(row)

        first_md = group[0][2]
        md_new: dict = {"sequence": seq}
        if first_md and "taxonomy" in first_md:
            md_new["taxonomy"] = first_md["taxonomy"]
        new_obs_metadata.append(md_new)

    new_table = Table(
        np.array(new_data),
        observation_ids=new_obs_ids,
        sample_ids=unified_samples,
        observation_metadata=new_obs_metadata,
        type="OTU table",
    )

    return {
        "biom_bytes": _biom_to_bytes(new_table),
        "total_asvs": total_asvs,
        "unique_asvs": len(new_obs_ids),
        "total_samples": n_unified,
        "collisions": has_collisions,
    }


def combine_biom_by_taxonomy(
    biom_paths: list[str],
    source_names: list[str],
    tax_level: str = "Genus",
) -> dict:
    """Combine BIOM files from different regions by taxonomy.

    Aggregates ASV counts to the chosen taxonomic level, then merges
    across tables. Output is a TSV with taxa as rows and samples as columns.

    Returns dict with: tsv_bytes, total_taxa, total_samples, tax_level
    """
    if tax_level not in TAX_LEVELS:
        raise ValueError(f"Invalid tax_level '{tax_level}'. Valid: {list(TAX_LEVELS.keys())}")

    level_idx = TAX_LEVELS[tax_level]
    tables = [_load_biom(p) for p in biom_paths]

    # Detect sample collisions
    all_sample_ids = []
    for t in tables:
        all_sample_ids.extend(t.ids(axis="sample"))
    has_collisions = len(all_sample_ids) != len(set(all_sample_ids))

    # Build unified sample list
    unified_samples = []
    sample_map: list[dict[str, int]] = []

    for t_idx, table in enumerate(tables):
        mapping = {}
        for sid in table.ids(axis="sample"):
            new_name = f"{source_names[t_idx]}_{sid}" if has_collisions else sid
            mapping[sid] = len(unified_samples)
            unified_samples.append(new_name)
        sample_map.append(mapping)

    n_unified = len(unified_samples)

    # Aggregate counts per taxon
    taxon_counts: dict[str, np.ndarray] = defaultdict(
        lambda: np.zeros(n_unified, dtype=np.float64)
    )

    for t_idx, table in enumerate(tables):
        for obs_id in table.ids(axis="observation"):
            md = table.metadata(obs_id, axis="observation")
            if not md or "taxonomy" not in md:
                taxon_label = "Unassigned"
            else:
                tax = md["taxonomy"]
                if isinstance(tax, str):
                    tax = [t.strip() for t in tax.split(";")]
                if len(tax) > level_idx:
                    taxon_label = "; ".join(tax[: level_idx + 1])
                else:
                    taxon_label = "; ".join(tax)

            obs_data = table.data(obs_id, axis="observation", dense=True)
            for j, sid in enumerate(table.ids(axis="sample")):
                taxon_counts[taxon_label][sample_map[t_idx][sid]] += obs_data[j]

    # Build TSV
    lines = ["#Taxon\t" + "\t".join(unified_samples)]
    for taxon in sorted(taxon_counts.keys()):
        counts = taxon_counts[taxon]
        vals = "\t".join(str(int(c)) if c == int(c) else f"{c:.1f}" for c in counts)
        lines.append(f"{taxon}\t{vals}")

    tsv_bytes = ("\n".join(lines) + "\n").encode("utf-8")

    return {
        "tsv_bytes": tsv_bytes,
        "total_taxa": len(taxon_counts),
        "total_samples": n_unified,
        "tax_level": tax_level,
    }
