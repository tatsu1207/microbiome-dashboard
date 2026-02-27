"""
MicrobiomeDash — Alpha diversity computation and statistics.
"""
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
from skbio.diversity import alpha_diversity


def compute_alpha(count_df: pd.DataFrame, metrics: list[str]) -> pd.DataFrame:
    """Compute alpha diversity metrics for each sample.

    Args:
        count_df: Features (rows) x samples (columns) integer count DataFrame.
        metrics: List of metric names (shannon, simpson, observed_otus, chao1, pielou_e).

    Returns:
        DataFrame with index=sample_id, columns=metric names.
    """
    # Transpose: samples as rows, features as columns
    counts = count_df.T
    sample_ids = counts.index.tolist()
    counts_array = counts.values.astype(int)

    results = {}

    skbio_metric_map = {
        "shannon": "shannon",
        "simpson": "simpson",
        "observed_otus": "observed_otus",
        "chao1": "chao1",
    }

    for metric in metrics:
        if metric == "pielou_e":
            # Pielou's evenness = H / ln(S), where H=Shannon, S=observed OTUs
            if "shannon" not in results:
                h = alpha_diversity("shannon", counts_array, ids=sample_ids)
                results["shannon"] = h
            if "observed_otus" not in results:
                s = alpha_diversity("observed_otus", counts_array, ids=sample_ids)
                results["observed_otus"] = s
            h_vals = results["shannon"].values
            s_vals = results["observed_otus"].values
            with np.errstate(divide="ignore", invalid="ignore"):
                pielou = np.where(s_vals > 1, h_vals / np.log(s_vals), 0.0)
            results["pielou_e"] = pd.Series(pielou, index=sample_ids)
        elif metric in skbio_metric_map:
            if metric not in results:
                results[metric] = alpha_diversity(
                    skbio_metric_map[metric], counts_array, ids=sample_ids
                )

    # Build output DataFrame with only requested metrics
    out = pd.DataFrame({m: results[m] for m in metrics if m in results})
    out.index.name = "sample_id"
    return out


def run_alpha_stats(
    diversity_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    sample_id_col: str,
    group_col: str,
    metric: str,
) -> dict:
    """Run Kruskal-Wallis and pairwise Mann-Whitney U tests.

    Returns {kruskal_H, kruskal_p, pairwise: [{group1, group2, U, pvalue, qvalue}]}.
    """
    # Merge diversity with metadata
    merged = diversity_df[[metric]].copy()
    merged.index = merged.index.astype(str)
    meta_map = meta_df.set_index(meta_df[sample_id_col].astype(str))[group_col]
    merged["group"] = merged.index.map(meta_map)
    merged = merged.dropna(subset=["group"])

    groups = merged["group"].unique().tolist()
    group_values = [merged[merged["group"] == g][metric].values for g in groups]

    # Kruskal-Wallis
    if len(groups) < 2 or any(len(gv) < 1 for gv in group_values):
        return {"kruskal_H": None, "kruskal_p": None, "pairwise": []}

    try:
        h_stat, kw_p = kruskal(*group_values)
    except ValueError:
        return {"kruskal_H": None, "kruskal_p": None, "pairwise": []}

    result = {"kruskal_H": float(h_stat), "kruskal_p": float(kw_p), "pairwise": []}

    # Pairwise Mann-Whitney U with BH correction if significant
    if kw_p < 0.05 and len(groups) >= 2:
        from statsmodels.stats.multitest import multipletests

        pairs = []
        p_values = []
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                g1_vals = merged[merged["group"] == groups[i]][metric].values
                g2_vals = merged[merged["group"] == groups[j]][metric].values
                if len(g1_vals) > 0 and len(g2_vals) > 0:
                    try:
                        u_stat, p_val = mannwhitneyu(g1_vals, g2_vals, alternative="two-sided")
                        pairs.append({
                            "group1": groups[i],
                            "group2": groups[j],
                            "U": float(u_stat),
                            "pvalue": float(p_val),
                        })
                        p_values.append(p_val)
                    except ValueError:
                        pass

        if p_values:
            _, q_values, _, _ = multipletests(p_values, method="fdr_bh")
            for pair, qval in zip(pairs, q_values):
                pair["qvalue"] = float(qval)

        result["pairwise"] = pairs

    return result
