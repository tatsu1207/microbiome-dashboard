"""
Venn Diagram Analysis — Compute shared/unique taxa sets between groups.
"""
import base64
import io

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from app.analysis.taxonomy import aggregate_counts_by_level


def compute_venn_sets(
    biom_path: str,
    meta_df: pd.DataFrame,
    sample_id_col: str,
    group_col: str,
    groups: list[str],
    level: str,
) -> dict[str, set]:
    """
    Compute sets of taxa present in each group.

    Returns {group_name: set_of_taxa_names}.
    """
    counts = aggregate_counts_by_level(biom_path, level)
    if counts.empty:
        return {}

    sets = {}
    for group in groups:
        group_samples = meta_df[meta_df[group_col] == group][sample_id_col].tolist()
        matched = [s for s in group_samples if s in counts.columns]
        if not matched:
            sets[str(group)] = set()
            continue

        group_counts = counts[matched]
        # A taxon is "present" in a group if it's detected in at least one sample
        present = group_counts.sum(axis=1) > 0
        sets[str(group)] = set(present[present].index.tolist())

    return sets


def render_venn_figure(sets: dict[str, set]) -> str:
    """
    Render a Venn diagram as a base64-encoded PNG string.
    Supports 2 or 3 groups. For 4+ groups, renders an intersection bar chart.
    Returns a data URI string for use in html.Img(src=...).
    """
    n_groups = len(sets)
    group_names = list(sets.keys())

    fig, ax = plt.subplots(figsize=(8, 6))

    if n_groups == 2:
        try:
            from matplotlib_venn import venn2
            venn2(
                [sets[group_names[0]], sets[group_names[1]]],
                set_labels=group_names,
                ax=ax,
            )
        except ImportError:
            return _fallback_bar_chart(sets)
    elif n_groups == 3:
        try:
            from matplotlib_venn import venn3
            venn3(
                [sets[group_names[0]], sets[group_names[1]], sets[group_names[2]]],
                set_labels=group_names,
                ax=ax,
            )
        except ImportError:
            return _fallback_bar_chart(sets)
    else:
        plt.close(fig)
        return _fallback_bar_chart(sets)

    ax.set_title(f"Shared Taxa ({n_groups} groups)")

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")
    return f"data:image/png;base64,{b64}"


def _fallback_bar_chart(sets: dict[str, set]) -> str:
    """Render an UpSet-style intersection bar chart for 4+ groups."""
    from itertools import combinations

    group_names = list(sets.keys())

    # Compute intersections
    bars = []
    for r in range(1, len(group_names) + 1):
        for combo in combinations(group_names, r):
            intersection = sets[combo[0]]
            for g in combo[1:]:
                intersection = intersection & sets[g]
            # Exclusive to this combination
            exclusive = intersection.copy()
            for g in group_names:
                if g not in combo:
                    exclusive -= sets[g]
            if exclusive:
                bars.append((" & ".join(combo), len(exclusive)))

    if not bars:
        return ""

    bars.sort(key=lambda x: x[1], reverse=True)
    labels, values = zip(*bars[:20])  # Top 20

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(range(len(labels)), values, color="#4C78A8")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_xlabel("Number of unique taxa")
    ax.set_title("Shared and Unique Taxa")
    ax.invert_yaxis()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")
    return f"data:image/png;base64,{b64}"


def compute_intersection_table(sets: dict[str, set]) -> pd.DataFrame:
    """
    Compute an intersection summary table.
    Returns DataFrame with columns: region, groups, count, taxa.
    """
    from itertools import combinations

    group_names = list(sets.keys())
    rows = []

    # Individual groups
    for g in group_names:
        unique = sets[g].copy()
        for other in group_names:
            if other != g:
                unique -= sets[other]
        rows.append({
            "region": f"Only {g}",
            "groups": g,
            "count": len(unique),
            "taxa": ", ".join(sorted(unique)[:20]) + ("..." if len(unique) > 20 else ""),
        })

    # Pairwise and higher intersections
    for r in range(2, len(group_names) + 1):
        for combo in combinations(group_names, r):
            shared = sets[combo[0]]
            for g in combo[1:]:
                shared = shared & sets[g]
            rows.append({
                "region": f"Shared: {' ∩ '.join(combo)}",
                "groups": " & ".join(combo),
                "count": len(shared),
                "taxa": ", ".join(sorted(shared)[:20]) + ("..." if len(shared) > 20 else ""),
            })

    return pd.DataFrame(rows)
