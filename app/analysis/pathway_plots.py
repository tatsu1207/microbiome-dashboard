"""
MicrobiomeDash — ggpicrust2-inspired pathway visualizations.

Provides errorbar, heatmap, and PCA plots for pathway DA results.
"""
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.stats import chi2


# ── Errorbar plot ────────────────────────────────────────────────────────────


def build_pathway_errorbar(
    results_df: pd.DataFrame,
    counts_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    sid_col: str,
    group_col: str,
    ref_group: str,
    test_group: str,
    top_n: int = 20,
    group_by_class: bool = False,
) -> go.Figure:
    """Build a grouped errorbar + log2FC direction plot.

    Left panel: horizontal grouped bars (mean relative abundance +/- SD per group).
    Right panel: log2FC direction bars (red = up in test, blue = up in ref).

    Args:
        results_df: DA results with 'feature', 'qvalue', 'log2fc' columns.
        counts_df: Features x samples count DataFrame.
        meta_df: Metadata DataFrame.
        sid_col: Sample ID column name in meta_df.
        group_col: Group column name in meta_df.
        ref_group: Reference group value.
        test_group: Test group value.
        top_n: Number of top significant features to show.
        group_by_class: If True, sort by pathway_class column.
    """
    if results_df.empty or "qvalue" not in results_df.columns:
        return go.Figure().update_layout(
            title="No significant results for errorbar plot",
            template="plotly_dark",
        )

    # Select top N significant features
    sig_df = results_df.dropna(subset=["qvalue"]).sort_values("qvalue")
    if group_by_class and "pathway_class" in sig_df.columns:
        sig_df = sig_df.sort_values(["pathway_class", "qvalue"])
    sig_df = sig_df.head(top_n)

    if sig_df.empty:
        return go.Figure().update_layout(
            title="No results to display",
            template="plotly_dark",
        )

    features = sig_df["feature"].tolist()

    # Get sample IDs per group
    meta = meta_df.copy()
    meta[sid_col] = meta[sid_col].astype(str)
    ref_ids = meta.loc[meta[group_col].astype(str) == str(ref_group), sid_col].tolist()
    test_ids = meta.loc[meta[group_col].astype(str) == str(test_group), sid_col].tolist()
    ref_ids = [s for s in ref_ids if s in counts_df.columns]
    test_ids = [s for s in test_ids if s in counts_df.columns]

    # Compute relative abundance (TSS normalization)
    col_sums = counts_df[ref_ids + test_ids].sum(axis=0)
    col_sums = col_sums.replace(0, 1)
    rel_abund = counts_df[ref_ids + test_ids].div(col_sums, axis=1) * 100

    # Compute mean and SD per group for selected features
    present_features = [f for f in features if f in rel_abund.index]
    if not present_features:
        return go.Figure().update_layout(
            title="No matching features in counts",
            template="plotly_dark",
        )

    # Build display labels
    labels = []
    for f in present_features:
        row = sig_df[sig_df["feature"] == f].iloc[0]
        desc = row.get("description", "")
        if desc and pd.notna(desc) and str(desc).strip():
            label = f"{f}: {str(desc)[:50]}"
        else:
            label = str(f)
        labels.append(label)

    ref_means = rel_abund.loc[present_features, ref_ids].mean(axis=1).values
    ref_sds = rel_abund.loc[present_features, ref_ids].std(axis=1).fillna(0).values
    test_means = rel_abund.loc[present_features, test_ids].mean(axis=1).values
    test_sds = rel_abund.loc[present_features, test_ids].std(axis=1).fillna(0).values

    # Get log2FC values
    log2fc_vals = []
    for f in present_features:
        match = sig_df[sig_df["feature"] == f]
        if not match.empty and "log2fc" in match.columns:
            log2fc_vals.append(float(match.iloc[0]["log2fc"]))
        else:
            log2fc_vals.append(0.0)

    fig = make_subplots(
        rows=1, cols=2,
        shared_yaxes=True,
        column_widths=[0.65, 0.35],
        horizontal_spacing=0.02,
        subplot_titles=["Mean Relative Abundance (%)", "log2 Fold Change"],
    )

    # Left panel: grouped bars with error bars
    fig.add_trace(
        go.Bar(
            y=labels,
            x=ref_means,
            error_x=dict(type="data", array=ref_sds, visible=True),
            name=str(ref_group),
            orientation="h",
            marker_color="#3498db",
            opacity=0.85,
        ),
        row=1, col=1,
    )
    fig.add_trace(
        go.Bar(
            y=labels,
            x=test_means,
            error_x=dict(type="data", array=test_sds, visible=True),
            name=str(test_group),
            orientation="h",
            marker_color="#e74c3c",
            opacity=0.85,
        ),
        row=1, col=1,
    )

    # Right panel: log2FC direction bars
    bar_colors = ["#e74c3c" if v > 0 else "#3498db" for v in log2fc_vals]
    fig.add_trace(
        go.Bar(
            y=labels,
            x=log2fc_vals,
            orientation="h",
            marker_color=bar_colors,
            showlegend=False,
            hovertemplate="log2FC: %{x:.2f}<extra></extra>",
        ),
        row=1, col=2,
    )

    # Add zero line to right panel
    fig.add_vline(x=0, line_dash="dash", line_color="gray", opacity=0.5, row=1, col=2)

    n_features = len(present_features)
    height = max(500, 50 + n_features * 28)

    fig.update_layout(
        template="plotly_dark",
        height=height,
        barmode="group",
        legend=dict(orientation="h", yanchor="bottom", y=1.02),
        margin=dict(l=10, r=10, t=60, b=40),
    )
    fig.update_yaxes(autorange="reversed", row=1, col=1)
    fig.update_yaxes(autorange="reversed", row=1, col=2)

    return fig


# ── Heatmap ──────────────────────────────────────────────────────────────────


def build_pathway_heatmap(
    results_df: pd.DataFrame,
    counts_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    sid_col: str,
    group_col: str,
    top_n: int = 30,
    cluster_rows: bool = True,
) -> go.Figure:
    """Build a z-score normalized heatmap of top significant pathways.

    Samples are ordered by group. Rows optionally clustered via
    scipy.cluster.hierarchy.linkage.

    Args:
        results_df: DA results with 'feature' and 'qvalue' columns.
        counts_df: Features x samples count DataFrame.
        meta_df: Metadata DataFrame.
        sid_col: Sample ID column in meta_df.
        group_col: Group column in meta_df.
        top_n: Number of top significant features.
        cluster_rows: Whether to cluster rows hierarchically.
    """
    if results_df.empty or "qvalue" not in results_df.columns:
        return go.Figure().update_layout(
            title="No significant results for heatmap",
            template="plotly_dark",
        )

    # Select top N significant features
    sig_df = results_df.dropna(subset=["qvalue"]).sort_values("qvalue").head(top_n)
    features = [f for f in sig_df["feature"].tolist() if f in counts_df.index]

    if not features:
        return go.Figure().update_layout(
            title="No matching features in counts",
            template="plotly_dark",
        )

    # Order samples by group
    meta = meta_df.copy()
    meta[sid_col] = meta[sid_col].astype(str)
    available_samples = [s for s in counts_df.columns if s in meta[sid_col].values]
    meta_sub = meta[meta[sid_col].isin(available_samples)].sort_values(group_col)
    ordered_samples = meta_sub[sid_col].tolist()

    # Subset and z-score normalize
    heat_data = counts_df.loc[features, ordered_samples].copy()

    # Z-score per row (feature)
    row_means = heat_data.mean(axis=1)
    row_stds = heat_data.std(axis=1).replace(0, 1)
    z_data = heat_data.sub(row_means, axis=0).div(row_stds, axis=0)

    # Cluster rows if requested
    if cluster_rows and len(features) > 2:
        try:
            Z = linkage(z_data.values, method="ward")
            order = leaves_list(Z)
            z_data = z_data.iloc[order]
        except Exception:
            pass

    # Build labels
    y_labels = []
    for f in z_data.index:
        match = sig_df[sig_df["feature"] == f]
        if not match.empty:
            desc = match.iloc[0].get("description", "")
            if desc and pd.notna(desc) and str(desc).strip():
                y_labels.append(f"{f}: {str(desc)[:40]}")
            else:
                y_labels.append(str(f))
        else:
            y_labels.append(str(f))

    # Sample group annotations for x-axis
    sample_groups = meta_sub.set_index(sid_col)[group_col]
    x_labels = [f"{s} ({sample_groups.get(s, '')})" for s in ordered_samples]

    fig = go.Figure(data=go.Heatmap(
        z=z_data.values,
        x=x_labels,
        y=y_labels,
        colorscale="RdBu_r",
        zmid=0,
        colorbar=dict(title="Z-score"),
        hovertemplate="Sample: %{x}<br>Feature: %{y}<br>Z-score: %{z:.2f}<extra></extra>",
    ))

    n_features = len(features)
    height = max(500, 80 + n_features * 22)

    fig.update_layout(
        title="Pathway Abundance Heatmap (Z-score normalized)",
        template="plotly_dark",
        height=height,
        xaxis=dict(tickangle=45, tickfont=dict(size=9)),
        yaxis=dict(tickfont=dict(size=10)),
        margin=dict(l=10, r=10, t=60, b=120),
    )

    return fig


# ── PCA ──────────────────────────────────────────────────────────────────────


def _compute_pca_svd(data_matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """PCA via numpy SVD on centered/scaled data.

    Returns (scores, proportion_explained) for first 2 components.
    """
    # Center
    centered = data_matrix - data_matrix.mean(axis=0)
    # Scale (unit variance)
    stds = centered.std(axis=0)
    stds[stds == 0] = 1
    scaled = centered / stds

    U, S, Vt = np.linalg.svd(scaled, full_matrices=False)
    # Scores for first 2 PCs
    scores = U[:, :2] * S[:2]
    # Proportion explained
    var_explained = (S ** 2) / (S ** 2).sum()
    return scores, var_explained[:2]


def _confidence_ellipse(
    x: np.ndarray, y: np.ndarray, confidence: float = 0.95
) -> tuple[np.ndarray, np.ndarray] | None:
    """Compute 95% confidence ellipse for a set of 2D points.

    Reuses the same math as beta.py:compute_confidence_ellipse.
    Returns (ellipse_x, ellipse_y) arrays or None if < 3 points.
    """
    if len(x) < 3:
        return None
    pts = np.column_stack([x, y])
    center = pts.mean(axis=0)
    cov = np.cov(pts, rowvar=False)
    scale = chi2.ppf(confidence, df=2)
    eigvals, eigvecs = np.linalg.eigh(cov)
    radii = np.sqrt(eigvals * scale)
    theta = np.linspace(0, 2 * np.pi, 100)
    circle = np.column_stack([np.cos(theta), np.sin(theta)])
    ellipse_pts = circle * radii @ eigvecs.T + center
    return ellipse_pts[:, 0], ellipse_pts[:, 1]


def build_pathway_pca(
    counts_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    sid_col: str,
    group_col: str,
    ref_group: str,
    test_group: str,
) -> go.Figure:
    """Build a PCA scatter plot (PC1 vs PC2) colored by group.

    PCA computed via numpy SVD (no sklearn dependency).
    95% confidence ellipses drawn for each group.

    Args:
        counts_df: Features x samples count DataFrame.
        meta_df: Metadata DataFrame.
        sid_col: Sample ID column.
        group_col: Group column.
        ref_group: Reference group value.
        test_group: Test group value.
    """
    meta = meta_df.copy()
    meta[sid_col] = meta[sid_col].astype(str)

    # Filter to the two groups
    meta_sub = meta[meta[group_col].astype(str).isin([str(ref_group), str(test_group)])]
    sample_ids = [s for s in meta_sub[sid_col].tolist() if s in counts_df.columns]

    if len(sample_ids) < 3:
        return go.Figure().update_layout(
            title="Need at least 3 samples for PCA",
            template="plotly_dark",
        )

    # Samples as rows, features as columns
    data_matrix = counts_df[sample_ids].T.values.astype(float)

    # Remove zero-variance features
    col_var = data_matrix.var(axis=0)
    keep = col_var > 0
    if keep.sum() < 2:
        return go.Figure().update_layout(
            title="Not enough variable features for PCA",
            template="plotly_dark",
        )
    data_matrix = data_matrix[:, keep]

    scores, prop_explained = _compute_pca_svd(data_matrix)

    # Build mapping from sample to group
    sample_group = meta_sub.set_index(sid_col)[group_col].to_dict()

    colors = {"ref": "#3498db", "test": "#e74c3c"}
    group_colors = {str(ref_group): colors["ref"], str(test_group): colors["test"]}

    fig = go.Figure()

    for grp, color in group_colors.items():
        mask = [sample_group.get(s, "") == grp for s in sample_ids]
        grp_scores = scores[mask]
        grp_names = [s for s, m in zip(sample_ids, mask) if m]

        fig.add_trace(go.Scatter(
            x=grp_scores[:, 0],
            y=grp_scores[:, 1],
            mode="markers",
            name=str(grp),
            marker=dict(color=color, size=10, opacity=0.8),
            text=grp_names,
            hovertemplate="<b>%{text}</b><br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<extra></extra>",
        ))

        # 95% confidence ellipse
        ellipse = _confidence_ellipse(grp_scores[:, 0], grp_scores[:, 1])
        if ellipse is not None:
            ex, ey = ellipse
            fig.add_trace(go.Scatter(
                x=np.append(ex, ex[0]),
                y=np.append(ey, ey[0]),
                mode="lines",
                line=dict(color=color, dash="dash", width=1.5),
                opacity=0.5,
                showlegend=False,
                hoverinfo="skip",
            ))

    pct1 = prop_explained[0] * 100
    pct2 = prop_explained[1] * 100

    fig.update_layout(
        title="PCA of Pathway Abundances",
        xaxis_title=f"PC1 ({pct1:.1f}%)",
        yaxis_title=f"PC2 ({pct2:.1f}%)",
        template="plotly_dark",
        height=500,
        legend=dict(orientation="h", yanchor="bottom", y=1.02),
    )

    return fig
