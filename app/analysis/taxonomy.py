"""
MicrobiomeDash — Taxonomy aggregation and visualization helpers.
"""
import numpy as np
import pandas as pd
from biom import load_table


LEVEL_MAP = {
    "ASV": -1,
    "Kingdom": 0,
    "Phylum": 1,
    "Class": 2,
    "Order": 3,
    "Family": 4,
    "Genus": 5,
    "Species": 6,
}


def aggregate_taxonomy(
    biom_path: str, level: str, top_n: int = 20
) -> pd.DataFrame:
    """Aggregate counts by taxonomy level, keep top N, convert to relative abundance.

    Returns DataFrame with rows=taxa, columns=samples (relative abundance 0-1).
    Returns empty DataFrame with error message as column if no taxonomy found.
    """
    table = load_table(biom_path)
    obs_ids = table.ids(axis="observation")
    sample_ids = table.ids(axis="sample")

    level_idx = LEVEL_MAP.get(level, 5)

    # Extract taxonomy assignments
    tax_map = {}
    has_taxonomy = False
    for obs_id in obs_ids:
        md = table.metadata(obs_id, axis="observation")
        if md and "taxonomy" in md:
            has_taxonomy = True
            ranks = md["taxonomy"]
            if isinstance(ranks, str):
                ranks = [r.strip() for r in ranks.split(";")]
            if level_idx < len(ranks) and ranks[level_idx]:
                tax_map[obs_id] = ranks[level_idx]
            else:
                tax_map[obs_id] = "Unassigned"
        else:
            tax_map[obs_id] = "Unassigned"

    if not has_taxonomy:
        return pd.DataFrame()

    # Build counts matrix
    data = pd.DataFrame(
        table.matrix_data.toarray(),
        index=obs_ids,
        columns=sample_ids,
    )

    # Map features to taxa and aggregate
    data["taxon"] = data.index.map(tax_map)
    agg = data.groupby("taxon").sum()

    # Keep top N by total abundance
    totals = agg.sum(axis=1).sort_values(ascending=False)
    top_taxa = totals.head(top_n).index.tolist()

    other = agg.loc[~agg.index.isin(top_taxa)].sum(axis=0)
    result = agg.loc[agg.index.isin(top_taxa)].copy()
    if other.sum() > 0:
        result.loc["Other"] = other

    # Convert to relative abundance
    col_sums = result.sum(axis=0)
    col_sums = col_sums.replace(0, 1)
    result = result.div(col_sums, axis=1)

    # Sort rows by total abundance (descending)
    result["_total"] = result.sum(axis=1)
    result = result.sort_values("_total", ascending=False).drop(columns="_total")

    return result


def build_heatmap_data(
    tax_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    sample_id_col: str,
    group_col: str,
) -> tuple[pd.DataFrame, list[str]]:
    """Reorder samples by group membership for heatmap display.

    Returns (reordered DataFrame, group labels list matching column order).
    """
    meta = meta_df.set_index(meta_df[sample_id_col].astype(str))
    common = [s for s in tax_df.columns if s in meta.index]

    if not common:
        return tax_df, []

    # Sort samples by group
    group_series = meta.loc[common, group_col]
    sorted_samples = group_series.sort_values().index.tolist()

    reordered = tax_df[sorted_samples]
    group_labels = group_series.loc[sorted_samples].tolist()

    return reordered, group_labels


def aggregate_counts_by_level(biom_path: str, level: str) -> pd.DataFrame:
    """Aggregate raw integer counts by taxonomic level.

    Returns DataFrame with rows=taxa (or ASV IDs), columns=samples.
    Unlike aggregate_taxonomy(), this returns raw counts (not relative abundance)
    and does NOT apply top-N filtering — suitable for DA tools.

    If level is "ASV", returns the original count matrix unchanged.
    """
    table = load_table(biom_path)
    obs_ids = table.ids(axis="observation")
    sample_ids = table.ids(axis="sample")

    data = pd.DataFrame(
        table.matrix_data.toarray().astype(int),
        index=obs_ids,
        columns=sample_ids,
    )

    level_idx = LEVEL_MAP.get(level, -1)
    if level_idx < 0:
        return data

    # Map each ASV to its taxon at the requested level
    tax_map = {}
    for obs_id in obs_ids:
        md = table.metadata(obs_id, axis="observation")
        if md and "taxonomy" in md:
            ranks = md["taxonomy"]
            if isinstance(ranks, str):
                ranks = [r.strip() for r in ranks.split(";")]
            if level_idx < len(ranks) and ranks[level_idx]:
                tax_map[obs_id] = ranks[level_idx]
            else:
                tax_map[obs_id] = "Unassigned"
        else:
            tax_map[obs_id] = "Unassigned"

    data.index = data.index.map(tax_map)
    return data.groupby(data.index).sum()
