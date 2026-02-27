"""
MicrobiomeDash — KEGG pathway map helpers.

Fetch pathway members from KEGG REST API, cross-reference with PICRUSt2
predictions, and compute per-sample pathway activity.
"""
import re
import threading
import time
import urllib.request
import urllib.parse

import pandas as pd


# ── Thread-safe cache for KEGG API results ───────────────────────────────────

_cache_lock = threading.Lock()
_cache: dict[str, tuple[float, object]] = {}
_CACHE_TTL = 3600  # 1 hour


def _cache_get(key: str):
    with _cache_lock:
        entry = _cache.get(key)
        if entry and (time.time() - entry[0]) < _CACHE_TTL:
            return entry[1]
    return None


def _cache_set(key: str, value):
    with _cache_lock:
        _cache[key] = (time.time(), value)


# ── KEGG REST helpers ────────────────────────────────────────────────────────

def _kegg_get(url: str) -> str:
    """Fetch a KEGG REST URL and return the response text."""
    req = urllib.request.Request(url, headers={"User-Agent": "MicrobiomeDash/1.0"})
    with urllib.request.urlopen(req, timeout=15) as resp:
        return resp.read().decode("utf-8")


def normalize_pathway_id(raw: str) -> str:
    """Normalize user input to 'map00680' format."""
    raw = raw.strip().lower()
    raw = re.sub(r"^map", "", raw)
    digits = re.sub(r"\D", "", raw)
    if not digits:
        raise ValueError(f"Invalid pathway ID: {raw!r}")
    return f"map{digits.zfill(5)}"


def fetch_pathway_kos(pathway_id: str) -> list[str]:
    """Fetch KO IDs belonging to a KEGG pathway. Returns e.g. ['K00018', ...]."""
    key = f"kos:{pathway_id}"
    cached = _cache_get(key)
    if cached is not None:
        return cached

    url = f"https://rest.kegg.jp/link/ko/{pathway_id}"
    text = _kegg_get(url)
    kos = []
    for line in text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            ko = parts[1].replace("ko:", "")
            kos.append(ko)
    kos = sorted(set(kos))
    _cache_set(key, kos)
    return kos


def fetch_pathway_ecs(pathway_id: str) -> list[str]:
    """Fetch EC numbers belonging to a KEGG pathway. Returns e.g. ['1.1.1.1', ...]."""
    key = f"ecs:{pathway_id}"
    cached = _cache_get(key)
    if cached is not None:
        return cached

    url = f"https://rest.kegg.jp/link/ec/{pathway_id}"
    text = _kegg_get(url)
    ecs = []
    for line in text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            ec = parts[1].replace("ec:", "")
            ecs.append(ec)
    ecs = sorted(set(ecs))
    _cache_set(key, ecs)
    return ecs


def fetch_pathway_name(pathway_id: str) -> str:
    """Fetch the human-readable name for a KEGG pathway."""
    key = f"name:{pathway_id}"
    cached = _cache_get(key)
    if cached is not None:
        return cached

    url = f"https://rest.kegg.jp/list/{pathway_id}"
    text = _kegg_get(url)
    name = pathway_id
    for line in text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            name = parts[1].strip()
            break
    _cache_set(key, name)
    return name


def build_kegg_color_url(
    pathway_id: str,
    color_map: dict[str, str],
) -> str:
    """Build a KEGG Mapper Color Pathway URL.

    Uses the KEGG show_pathway endpoint with slash-delimited dataset.
    See: https://www.kegg.jp/kegg/webapp/color_url.html

    KEGG rejects URLs longer than ~2000 characters with 403 Forbidden,
    so we cap at 60 entries (prioritizing colored/significant ones).

    Args:
        pathway_id: e.g. "map00680"
        color_map: Mapping of normalized IDs to hex colors,
            e.g. {"K00399": "#ff0000", "1.1.1.1": "#00ff00"}

    Returns:
        URL string for KEGG Mapper Color tool.
    """
    if not color_map:
        return f"https://www.kegg.jp/pathway/{pathway_id}"

    # KEGG rejects URLs > ~2000 chars. Cap at 60 entries, prioritizing
    # non-yellow (significant) entries over yellow (detected-only).
    MAX_ENTRIES = 60
    if len(color_map) > MAX_ENTRIES:
        significant = {k: v for k, v in color_map.items() if v != "#ffff00"}
        neutral = {k: v for k, v in color_map.items() if v == "#ffff00"}
        items = list(significant.items())[:MAX_ENTRIES]
        remaining = MAX_ENTRIES - len(items)
        if remaining > 0:
            items.extend(list(neutral.items())[:remaining])
    else:
        items = list(color_map.items())

    # Build slash-delimited dataset for KEGG show_pathway endpoint.
    # Format: show_pathway?{mapid}/{kid1}%09{bgcolor},{fgcolor}/{kid2}%09...
    # %09 = tab, %23 = #
    parts = []
    for oid, hex_color in items:
        # Strip leading '#' and encode as %23
        color_code = hex_color.lstrip("#")
        parts.append(f"{oid}%09%23{color_code},black")

    return f"https://www.kegg.jp/kegg-bin/show_pathway?{pathway_id}/{'/'.join(parts)}"


def cross_reference_ids(
    pathway_ids: list[str], counts_df: pd.DataFrame, id_type: str
) -> dict:
    """Match pathway KOs/ECs against PICRUSt2 counts.

    Args:
        pathway_ids: List of KO or EC IDs from the KEGG pathway.
        counts_df: PICRUSt2 counts DataFrame (features x samples).
        id_type: "ko" or "ec".

    Returns:
        Dict with 'detected', 'missing', 'total', 'coverage_pct'.
    """
    # Normalize PICRUSt2 index for matching
    if id_type == "ec":
        # PICRUSt2 EC format: "EC:1.1.1.1" -> "1.1.1.1"
        data_ids = {idx.replace("EC:", "") for idx in counts_df.index}
        data_id_map = {}
        for idx in counts_df.index:
            normalized = idx.replace("EC:", "")
            data_id_map[normalized] = idx
    else:
        # KO format: "K00018" — direct match
        data_ids = set(counts_df.index)
        data_id_map = {idx: idx for idx in counts_df.index}

    pathway_set = set(pathway_ids)
    detected = sorted(pathway_set & data_ids)
    missing = sorted(pathway_set - data_ids)

    return {
        "detected": detected,
        "missing": missing,
        "total": len(pathway_ids),
        "n_detected": len(detected),
        "coverage_pct": (len(detected) / len(pathway_ids) * 100) if pathway_ids else 0,
        "id_map": data_id_map,  # normalized -> original index
    }


def compute_pathway_activity(
    detected_ids: list[str],
    counts_df: pd.DataFrame,
    id_map: dict[str, str],
) -> pd.Series:
    """Sum detected pathway feature abundances per sample.

    Args:
        detected_ids: Normalized IDs (e.g. 'K00018' or '1.1.1.1') that are in the pathway.
        counts_df: PICRUSt2 counts (features x samples).
        id_map: Mapping from normalized ID -> original counts_df index label.

    Returns:
        pd.Series indexed by sample ID with summed pathway activity.
    """
    original_ids = [id_map[d] for d in detected_ids if d in id_map]
    if not original_ids:
        return pd.Series(dtype=float)
    subset = counts_df.loc[original_ids]
    return subset.sum(axis=0)


def pathway_activity_stats(
    activity: pd.Series,
    meta_df: pd.DataFrame,
    sid_col: str,
    group_col: str,
    ref_group: str,
    test_group: str,
) -> dict:
    """Mann-Whitney U test on pathway activity between two groups.

    Returns dict with 'U', 'pvalue', 'ref_median', 'test_median', 'ref_n', 'test_n'.
    """
    meta_df = meta_df.copy()
    meta_df[sid_col] = meta_df[sid_col].astype(str)

    ref_ids = meta_df.loc[meta_df[group_col].astype(str) == str(ref_group), sid_col].tolist()
    test_ids = meta_df.loc[meta_df[group_col].astype(str) == str(test_group), sid_col].tolist()

    ref_vals = activity.reindex(ref_ids).dropna()
    test_vals = activity.reindex(test_ids).dropna()

    result = {
        "ref_group": ref_group,
        "test_group": test_group,
        "ref_n": len(ref_vals),
        "test_n": len(test_vals),
        "ref_median": float(ref_vals.median()) if len(ref_vals) > 0 else None,
        "test_median": float(test_vals.median()) if len(test_vals) > 0 else None,
        "U": None,
        "pvalue": None,
    }

    if len(ref_vals) >= 1 and len(test_vals) >= 1:
        try:
            from scipy.stats import mannwhitneyu
            u_stat, p_val = mannwhitneyu(ref_vals, test_vals, alternative="two-sided")
            result["U"] = float(u_stat)
            result["pvalue"] = float(p_val)
        except Exception:
            pass

    return result


def ko_coverage_by_group(
    pathway_ids: list[str],
    counts_df: pd.DataFrame,
    id_map: dict[str, str],
    meta_df: pd.DataFrame,
    sid_col: str,
    group_col: str,
    ref_group: str,
    test_group: str,
    desc_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Per-KO/EC table with detected status and mean abundance per group.

    Returns DataFrame with columns: feature_id, detected, description,
    mean_{ref_group}, mean_{test_group}.
    """
    meta_df = meta_df.copy()
    meta_df[sid_col] = meta_df[sid_col].astype(str)

    ref_ids = meta_df.loc[meta_df[group_col].astype(str) == str(ref_group), sid_col].tolist()
    test_ids = meta_df.loc[meta_df[group_col].astype(str) == str(test_group), sid_col].tolist()

    # Intersect with available sample columns
    ref_ids = [s for s in ref_ids if s in counts_df.columns]
    test_ids = [s for s in test_ids if s in counts_df.columns]

    rows = []
    for pid in pathway_ids:
        orig = id_map.get(pid)
        detected = orig is not None
        ref_mean = None
        test_mean = None

        if detected and orig in counts_df.index:
            if ref_ids:
                ref_mean = float(counts_df.loc[orig, ref_ids].mean())
            if test_ids:
                test_mean = float(counts_df.loc[orig, test_ids].mean())

        row = {
            "feature_id": pid,
            "detected": detected,
            "mean_ref": ref_mean,
            "mean_test": test_mean,
        }

        # Try to get description
        if desc_df is not None and orig is not None and orig in desc_df.index:
            row["description"] = desc_df.loc[orig, "description"]
        else:
            row["description"] = ""

        rows.append(row)

    df = pd.DataFrame(rows)
    # Rename columns for display
    df = df.rename(columns={
        "mean_ref": f"mean_{ref_group}",
        "mean_test": f"mean_{test_group}",
    })
    # Sort: detected first, then by abundance descending
    ref_col = f"mean_{ref_group}"
    if ref_col in df.columns:
        df = df.sort_values(
            by=["detected", ref_col],
            ascending=[False, False],
        )
    return df
