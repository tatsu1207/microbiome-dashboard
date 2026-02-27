"""
MicrobiomeDash — KEGG KO-to-pathway aggregation and annotation.

Provides disk-cached KEGG API data (24h TTL) and functions to aggregate
KO-level PICRUSt2 predictions to KEGG pathway level.
"""
import json
import time
import urllib.request
from pathlib import Path

import pandas as pd

from app.config import DATA_DIR

KEGG_CACHE_DIR = DATA_DIR / "kegg_cache"
CACHE_TTL = 86400  # 24 hours


def _ensure_cache_dir():
    KEGG_CACHE_DIR.mkdir(parents=True, exist_ok=True)


def _cache_read(filename: str) -> object | None:
    """Read a JSON cache file if it exists and is fresh (< TTL)."""
    path = KEGG_CACHE_DIR / filename
    if not path.exists():
        return None
    try:
        age = time.time() - path.stat().st_mtime
        if age > CACHE_TTL:
            return None
        return json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        return None


def _cache_write(filename: str, data: object):
    """Write data to a JSON cache file."""
    _ensure_cache_dir()
    path = KEGG_CACHE_DIR / filename
    path.write_text(json.dumps(data))


def _kegg_get(url: str) -> str:
    """Fetch a KEGG REST URL and return the response text."""
    req = urllib.request.Request(url, headers={"User-Agent": "MicrobiomeDash/1.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        return resp.read().decode("utf-8")


# ── KO → Pathway mapping ────────────────────────────────────────────────────


def fetch_ko_pathway_mapping() -> dict[str, list[str]]:
    """Fetch full KO → pathway mapping from KEGG.

    Returns dict mapping KO IDs (e.g. 'K00001') to lists of pathway IDs
    (e.g. ['map00010', 'map00071']).
    """
    cached = _cache_read("ko_pathway_map.json")
    if cached is not None:
        return cached

    text = _kegg_get("https://rest.kegg.jp/link/pathway/ko")
    mapping: dict[str, list[str]] = {}
    for line in text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        ko_id = parts[0].replace("ko:", "")
        pathway_id = parts[1].replace("path:", "")
        # Only keep reference pathways (map*), skip organism-specific
        if not pathway_id.startswith("map"):
            continue
        mapping.setdefault(ko_id, []).append(pathway_id)

    _cache_write("ko_pathway_map.json", mapping)
    return mapping


# ── Pathway names ────────────────────────────────────────────────────────────


def fetch_pathway_names() -> dict[str, str]:
    """Fetch pathway ID → name mapping from KEGG.

    Returns dict like {'map00010': 'Glycolysis / Gluconeogenesis', ...}.
    """
    cached = _cache_read("pathway_names.json")
    if cached is not None:
        return cached

    text = _kegg_get("https://rest.kegg.jp/list/pathway/map")
    names: dict[str, str] = {}
    for line in text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            pathway_id = parts[0].replace("path:", "")
            names[pathway_id] = parts[1].strip()

    _cache_write("pathway_names.json", names)
    return names


# ── BRITE hierarchy ──────────────────────────────────────────────────────────


def fetch_brite_hierarchy() -> dict[str, dict]:
    """Fetch KEGG BRITE pathway hierarchy (br08901).

    Returns dict mapping pathway IDs to
    {'class_a': str, 'class_b': str, 'name': str}.
    """
    cached = _cache_read("brite_hierarchy.json")
    if cached is not None:
        return cached

    text = _kegg_get("https://rest.kegg.jp/get/br:br08901")
    hierarchy: dict[str, dict] = {}
    current_a = ""
    current_b = ""

    for line in text.splitlines():
        if not line:
            continue
        first_char = line[0] if line else ""
        content = line[1:].strip()

        if first_char == "A":
            # Top-level class (e.g. "A09100 Metabolism")
            current_a = content
        elif first_char == "B":
            # Sub-class (e.g. "B  09101 Carbohydrate metabolism")
            current_b = content
        elif first_char == "C":
            # Pathway entry (e.g. "C    00010 Glycolysis / Gluconeogenesis [PATH:map00010]")
            # Extract pathway number from the line
            parts = content.strip().split(None, 1)
            if parts:
                pathway_num = parts[0]
                pathway_name = parts[1] if len(parts) > 1 else ""
                # Remove [PATH:mapXXXXX] suffix
                if "[PATH:" in pathway_name:
                    pathway_name = pathway_name[:pathway_name.index("[PATH:")].strip()
                pathway_id = f"map{pathway_num}"
                hierarchy[pathway_id] = {
                    "class_a": current_a,
                    "class_b": current_b,
                    "name": pathway_name,
                }

    _cache_write("brite_hierarchy.json", hierarchy)
    return hierarchy


# ── KO → Pathway aggregation ────────────────────────────────────────────────


def aggregate_ko_to_pathways(
    ko_counts_df: pd.DataFrame,
    ko_pathway_map: dict[str, list[str]] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate KO-level counts to KEGG pathway level.

    Sums constituent KO abundances per pathway per sample.

    Args:
        ko_counts_df: Features (KO IDs) x samples count DataFrame.
        ko_pathway_map: Optional pre-fetched mapping. Fetched if None.

    Returns:
        (pathway_counts_df, pathway_desc_df)
        - pathway_counts_df: Pathways x samples counts
        - pathway_desc_df: DataFrame with columns [feature, description, pathway_class]
    """
    if ko_pathway_map is None:
        ko_pathway_map = fetch_ko_pathway_mapping()

    # Build reverse mapping: pathway -> [KO IDs present in counts_df]
    pathway_kos: dict[str, list[str]] = {}
    available_kos = set(ko_counts_df.index)
    n_mapped = 0

    for ko_id, pathways in ko_pathway_map.items():
        if ko_id in available_kos:
            n_mapped += 1
            for pw in pathways:
                pathway_kos.setdefault(pw, []).append(ko_id)

    if not pathway_kos:
        raise ValueError(
            f"No KO IDs could be mapped to KEGG pathways. "
            f"Found {len(available_kos)} KO IDs in data, "
            f"{len(ko_pathway_map)} in KEGG mapping."
        )

    # Sum KO abundances per pathway
    pathway_data: dict[str, pd.Series] = {}
    for pw_id, kos in pathway_kos.items():
        present_kos = [k for k in kos if k in ko_counts_df.index]
        if present_kos:
            pathway_data[pw_id] = ko_counts_df.loc[present_kos].sum(axis=0)

    pathway_counts_df = pd.DataFrame(pathway_data).T
    pathway_counts_df.index.name = "feature_id"

    # Build description DataFrame
    try:
        pathway_names = fetch_pathway_names()
    except Exception:
        pathway_names = {}

    try:
        brite = fetch_brite_hierarchy()
    except Exception:
        brite = {}

    desc_rows = []
    for pw_id in pathway_counts_df.index:
        name = pathway_names.get(pw_id, "")
        brite_info = brite.get(pw_id, {})
        pathway_class = brite_info.get("class_b", brite_info.get("class_a", ""))
        desc_rows.append({
            "feature": pw_id,
            "description": name,
            "pathway_class": pathway_class,
        })

    pathway_desc_df = pd.DataFrame(desc_rows)
    if not pathway_desc_df.empty:
        pathway_desc_df = pathway_desc_df.set_index("feature")

    return pathway_counts_df, pathway_desc_df


# ── Annotation helper ────────────────────────────────────────────────────────


def annotate_pathway_results(
    results_df: pd.DataFrame,
    pred_type: str,
    desc_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Add annotation columns (description, pathway_class) to DA results.

    For KEGG pathways (aggregated): uses BRITE hierarchy.
    For MetaCyc/KO/EC: merges desc_df descriptions.

    Args:
        results_df: DA results with 'feature' column.
        pred_type: 'kegg_pathway', 'metacyc', 'ko', or 'ec'.
        desc_df: Optional DataFrame with description/pathway_class columns.

    Returns:
        Annotated results DataFrame.
    """
    if results_df.empty or "feature" not in results_df.columns:
        return results_df

    df = results_df.copy()

    if pred_type == "kegg_pathway":
        # Annotate from BRITE hierarchy
        try:
            brite = fetch_brite_hierarchy()
            pathway_names = fetch_pathway_names()
        except Exception:
            brite = {}
            pathway_names = {}

        df["description"] = df["feature"].map(
            lambda f: pathway_names.get(f, "")
        )
        df["pathway_class"] = df["feature"].map(
            lambda f: brite.get(f, {}).get("class_b", "")
        )
    elif desc_df is not None and not desc_df.empty:
        # Merge from provided descriptions
        merge_cols = [c for c in desc_df.columns if c in ("description", "pathway_class")]
        if merge_cols:
            df = df.merge(
                desc_df[merge_cols],
                left_on="feature",
                right_index=True,
                how="left",
            )

    # Reorder: feature, description, pathway_class, then rest
    priority_cols = ["feature"]
    for col in ["description", "pathway_class"]:
        if col in df.columns:
            priority_cols.append(col)
    other_cols = [c for c in df.columns if c not in priority_cols]
    df = df[priority_cols + other_cols]

    return df
