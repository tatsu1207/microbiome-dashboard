"""
Core Microbiome Analysis — Identify taxa present in a given proportion of samples.
"""
import pandas as pd

from app.analysis.taxonomy import aggregate_counts_by_level


def compute_core_microbiome(
    biom_path: str,
    meta_df: pd.DataFrame,
    sample_id_col: str,
    group_col: str | None,
    level: str,
    threshold: float = 0.8,
) -> dict:
    """
    Compute core microbiome at a given taxonomic level.

    Args:
        biom_path: Path to BIOM file
        meta_df: Metadata DataFrame
        sample_id_col: Column name for sample IDs
        group_col: Column to group by (None for all samples)
        level: Taxonomic level (Phylum, Class, Order, Family, Genus, Species, ASV)
        threshold: Minimum prevalence (0.0–1.0) to be considered "core"

    Returns:
        dict with keys:
            "overall": DataFrame of core taxa across all samples (taxon, prevalence)
            "per_group": {group_name: DataFrame} if group_col provided
            "all_prevalences": Full prevalence DataFrame for plotting
    """
    # Get counts aggregated by taxonomic level
    counts = aggregate_counts_by_level(biom_path, level)

    if counts.empty:
        return {"overall": pd.DataFrame(), "per_group": {}, "all_prevalences": pd.DataFrame()}

    # Compute prevalence across all samples
    presence = (counts > 0).astype(int)
    overall_prev = presence.sum(axis=1) / presence.shape[1]
    overall_prev.name = "prevalence"

    overall_core = overall_prev[overall_prev >= threshold].sort_values(ascending=False)
    overall_df = pd.DataFrame({"taxon": overall_core.index, "prevalence": overall_core.values})

    # Per-group analysis
    per_group = {}
    all_prev_data = {"taxon": counts.index.tolist()}

    if group_col and group_col in meta_df.columns:
        groups = meta_df[group_col].dropna().unique()

        for group in sorted(groups):
            group_samples = meta_df[meta_df[group_col] == group][sample_id_col].tolist()
            # Match to available columns in counts
            matched = [s for s in group_samples if s in counts.columns]
            if not matched:
                continue

            group_counts = counts[matched]
            group_presence = (group_counts > 0).astype(int)
            group_prev = group_presence.sum(axis=1) / group_presence.shape[1]

            group_core = group_prev[group_prev >= threshold].sort_values(ascending=False)
            per_group[str(group)] = pd.DataFrame({
                "taxon": group_core.index,
                "prevalence": group_core.values,
            })

            all_prev_data[str(group)] = group_prev.values

    all_prev_data["overall"] = overall_prev.values
    all_prevalences = pd.DataFrame(all_prev_data)

    return {
        "overall": overall_df,
        "per_group": per_group,
        "all_prevalences": all_prevalences,
    }
