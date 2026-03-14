"""
16S Pipeline — Generate a comprehensive analysis report PDF.

Sections: dataset summary, methods text, alpha diversity, beta diversity (PCoA),
taxonomy bar plot, and read tracking.  Uses matplotlib (Agg) for charts and
fpdf2 for PDF assembly.
"""
import logging
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from fpdf import FPDF

from app.analysis.shared import biom_to_count_df, get_dataset_metadata_df
from app.report.methods_text import generate_methods_text

logger = logging.getLogger(__name__)


# ── Public API ───────────────────────────────────────────────────────────────

def generate_report_pdf(
    dataset_id: int,
    group_col: str | None = None,
    sections: list[str] | None = None,
) -> str:
    """Generate a multi-section analysis report PDF.

    Args:
        dataset_id: Database dataset ID.
        group_col: Metadata column for grouping (required for alpha/beta).
        sections: List of section names to include. Default: all.
            Options: "summary", "methods", "alpha", "beta", "taxonomy", "reads"

    Returns:
        Path to the generated PDF file.
    """
    from app.config import EXPORT_DIR
    from app.db.database import get_session
    from app.db.models import Dataset, Sample

    if sections is None:
        sections = ["summary", "methods", "alpha", "beta", "taxonomy", "reads"]

    # Load dataset info
    with get_session() as db:
        dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
        if not dataset:
            raise ValueError(f"Dataset {dataset_id} not found")

        ds_info = {
            "id": dataset.id,
            "name": dataset.name,
            "sequencing_type": dataset.sequencing_type or "unknown",
            "variable_region": dataset.variable_region or "unknown",
            "platform": dataset.platform or "illumina",
            "sample_count": dataset.sample_count or 0,
            "asv_count": dataset.asv_count or 0,
            "trunc_len_f": dataset.trunc_len_f,
            "trunc_len_r": dataset.trunc_len_r,
            "biom_path": dataset.asv_table_path,
        }

        samples = (
            db.query(Sample)
            .filter(Sample.dataset_id == dataset_id)
            .order_by(Sample.sample_name)
            .all()
        )
        sample_data = []
        for s in samples:
            sample_data.append({
                "name": s.sample_name,
                "raw": s.read_count_raw or 0,
                "filtered": s.read_count_filtered or 0,
                "denoised": s.read_count_denoised or 0,
                "nonchim": s.read_count_nonchimeric or 0,
            })

    biom_path = ds_info["biom_path"]
    if not biom_path or not Path(biom_path).exists():
        raise ValueError("BIOM file not found for this dataset")

    # Load BIOM data and metadata
    count_df = biom_to_count_df(biom_path)
    meta_df, sid_col = get_dataset_metadata_df(biom_path)
    has_meta = meta_df is not None and group_col and group_col in (meta_df.columns if meta_df is not None else [])

    # ── Build PDF ────────────────────────────────────────────────────────
    pdf = FPDF(orientation="L", format="A4")
    pdf.set_auto_page_break(auto=True, margin=15)

    # 1. Summary page
    if "summary" in sections:
        _add_summary_page(pdf, ds_info, sample_data)

    # 2. Methods text
    if "methods" in sections:
        _add_methods_page(pdf, dataset_id)

    # 3. Alpha diversity
    if "alpha" in sections and has_meta:
        _add_alpha_section(pdf, count_df, meta_df, sid_col, group_col)

    # 4. Beta diversity (PCoA)
    if "beta" in sections and has_meta:
        _add_beta_section(pdf, count_df, meta_df, sid_col, group_col)

    # 5. Taxonomy bar plot
    if "taxonomy" in sections:
        _add_taxonomy_section(pdf, biom_path, meta_df, sid_col, group_col if has_meta else None)

    # 6. Read tracking
    if "reads" in sections and sample_data:
        _add_reads_section(pdf, sample_data)

    # Write PDF
    EXPORT_DIR.mkdir(parents=True, exist_ok=True)
    pdf_path = EXPORT_DIR / f"report_dataset_{dataset_id}.pdf"
    pdf.output(str(pdf_path))
    logger.info(f"Analysis report PDF written: {pdf_path}")
    return str(pdf_path)


# ── Section builders ─────────────────────────────────────────────────────────

def _add_summary_page(pdf: FPDF, ds_info: dict, sample_data: list[dict]):
    """Title page with dataset summary statistics."""
    pdf.add_page()
    pdf.set_font("Helvetica", "B", 22)
    pdf.cell(text="16S Analysis Report", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(8)

    pdf.set_font("Helvetica", "", 12)
    pdf.cell(text=f"Dataset: {ds_info['name']}  (ID: {ds_info['id']})",
             new_x="LMARGIN", new_y="NEXT")
    pdf.ln(2)

    info_lines = [
        f"Sequencing: {ds_info['sequencing_type'].upper()}  |  "
        f"Region: {ds_info['variable_region']}  |  "
        f"Platform: {ds_info['platform']}",
        f"Samples: {ds_info['sample_count']}  |  ASVs: {ds_info['asv_count']}",
    ]
    if ds_info["trunc_len_f"]:
        trunc = f"Truncation: F={ds_info['trunc_len_f']}"
        if ds_info["trunc_len_r"]:
            trunc += f", R={ds_info['trunc_len_r']}"
        info_lines.append(trunc)

    for line in info_lines:
        pdf.cell(text=line, new_x="LMARGIN", new_y="NEXT")
        pdf.ln(1)

    # Read count summary
    if sample_data:
        nonchim = [s["nonchim"] for s in sample_data if s["nonchim"] > 0]
        if nonchim:
            pdf.ln(4)
            pdf.set_font("Helvetica", "B", 12)
            pdf.cell(text="Read Count Summary", new_x="LMARGIN", new_y="NEXT")
            pdf.ln(2)
            pdf.set_font("Helvetica", "", 11)
            pdf.cell(
                text=f"Mean non-chimeric reads: {int(np.mean(nonchim)):,}  "
                     f"(min: {min(nonchim):,}, max: {max(nonchim):,})",
                new_x="LMARGIN", new_y="NEXT",
            )
            total_raw = sum(s["raw"] for s in sample_data)
            total_nonchim = sum(nonchim)
            if total_raw > 0:
                pdf.cell(
                    text=f"Overall retention: {total_nonchim / total_raw * 100:.1f}%  "
                         f"({total_nonchim:,} / {total_raw:,} reads)",
                    new_x="LMARGIN", new_y="NEXT",
                )


def _add_methods_page(pdf: FPDF, dataset_id: int):
    """Materials & Methods paragraph."""
    pdf.add_page()
    pdf.set_font("Helvetica", "B", 16)
    pdf.cell(text="Materials & Methods", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(6)

    text = generate_methods_text(dataset_id)
    pdf.set_font("Helvetica", "", 10)
    pdf.multi_cell(w=0, h=5, text=text)


def _add_alpha_section(pdf: FPDF, count_df, meta_df, sid_col, group_col):
    """Alpha diversity box plots and Kruskal-Wallis statistics."""
    from app.analysis.alpha import compute_alpha, run_alpha_stats

    metrics = ["shannon", "observed_otus", "chao1", "simpson"]
    try:
        alpha_df = compute_alpha(count_df, metrics)
    except Exception as e:
        logger.warning(f"Alpha diversity computation failed: {e}")
        return

    # Merge with metadata
    alpha_df.index = alpha_df.index.astype(str)
    meta_map = meta_df.set_index(meta_df[sid_col].astype(str))[group_col]
    alpha_df["group"] = alpha_df.index.map(meta_map)
    alpha_df = alpha_df.dropna(subset=["group"])

    if alpha_df.empty:
        return

    groups = sorted(alpha_df["group"].unique())
    colors = _get_group_colors(len(groups))

    # Box plots — 2x2 grid
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    metric_labels = {
        "shannon": "Shannon Index",
        "observed_otus": "Observed ASVs",
        "chao1": "Chao1",
        "simpson": "Simpson Index",
    }

    for ax, metric in zip(axes.flat, metrics):
        data_by_group = [alpha_df[alpha_df["group"] == g][metric].values for g in groups]
        bp = ax.boxplot(data_by_group, labels=groups, patch_artist=True, widths=0.6)
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        ax.set_title(metric_labels.get(metric, metric), fontsize=11)
        ax.tick_params(axis="x", rotation=30)

    fig.suptitle("Alpha Diversity", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    chart_path = _save_temp_chart(fig)
    pdf.add_page()
    pdf.set_font("Helvetica", "B", 16)
    pdf.cell(text="Alpha Diversity", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(4)
    pdf.image(chart_path, w=260)
    Path(chart_path).unlink(missing_ok=True)

    # Statistics table
    pdf.ln(4)
    pdf.set_font("Helvetica", "B", 11)
    pdf.cell(text="Kruskal-Wallis Test Results", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(2)

    headers = ["Metric", "H statistic", "p-value"]
    col_widths = [50, 40, 40]
    pdf.set_font("Helvetica", "B", 9)
    for h, w in zip(headers, col_widths):
        pdf.cell(w, 7, h, border=1)
    pdf.ln()

    pdf.set_font("Helvetica", "", 9)
    for metric in metrics:
        try:
            stats = run_alpha_stats(
                compute_alpha(count_df, [metric]),
                meta_df, sid_col, group_col, metric,
            )
            h_val = f"{stats['kruskal_H']:.3f}" if stats["kruskal_H"] is not None else "-"
            p_val = _format_pvalue(stats["kruskal_p"]) if stats["kruskal_p"] is not None else "-"
        except Exception:
            h_val, p_val = "-", "-"

        row = [metric_labels.get(metric, metric), h_val, p_val]
        for val, w in zip(row, col_widths):
            pdf.cell(w, 6, val, border=1)
        pdf.ln()


def _add_beta_section(pdf: FPDF, count_df, meta_df, sid_col, group_col):
    """PCoA scatter plot and PERMANOVA results."""
    from app.analysis.beta import compute_distance, run_pcoa, run_permanova_global

    try:
        dm = compute_distance(count_df, "braycurtis")
        coords, prop_exp = run_pcoa(dm)
    except Exception as e:
        logger.warning(f"Beta diversity computation failed: {e}")
        return

    # Merge with metadata
    coords.index = coords.index.astype(str)
    meta_map = meta_df.set_index(meta_df[sid_col].astype(str))[group_col]
    coords["group"] = coords.index.map(meta_map)
    coords = coords.dropna(subset=["group"])

    if coords.empty:
        return

    groups = sorted(coords["group"].unique())
    colors = _get_group_colors(len(groups))
    color_map = dict(zip(groups, colors))

    fig, ax = plt.subplots(figsize=(9, 7))
    for grp in groups:
        mask = coords["group"] == grp
        ax.scatter(
            coords.loc[mask, "PC1"], coords.loc[mask, "PC2"],
            label=grp, color=color_map[grp], s=50, alpha=0.8, edgecolors="white",
        )

    pc1_pct = prop_exp.get("PC1", 0) * 100
    pc2_pct = prop_exp.get("PC2", 0) * 100
    ax.set_xlabel(f"PC1 ({pc1_pct:.1f}%)", fontsize=11)
    ax.set_ylabel(f"PC2 ({pc2_pct:.1f}%)", fontsize=11)
    ax.set_title("PCoA (Bray-Curtis)", fontsize=14, fontweight="bold")
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()

    chart_path = _save_temp_chart(fig)
    pdf.add_page()
    pdf.set_font("Helvetica", "B", 16)
    pdf.cell(text="Beta Diversity", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(4)
    pdf.image(chart_path, w=200)
    Path(chart_path).unlink(missing_ok=True)

    # PERMANOVA
    try:
        perm = run_permanova_global(dm, meta_df, sid_col, group_col)
        if perm:
            pdf.ln(4)
            pdf.set_font("Helvetica", "B", 11)
            pdf.cell(text="PERMANOVA (999 permutations)", new_x="LMARGIN", new_y="NEXT")
            pdf.ln(2)
            pdf.set_font("Helvetica", "", 10)
            pdf.cell(
                text=f"pseudo-F = {perm['pseudo_F']:.3f},  "
                     f"p = {_format_pvalue(perm['pvalue'])},  "
                     f"n = {perm['n_samples']} samples, "
                     f"{perm['n_groups']} groups",
                new_x="LMARGIN", new_y="NEXT",
            )
    except Exception as e:
        logger.warning(f"PERMANOVA failed: {e}")


def _add_taxonomy_section(pdf: FPDF, biom_path, meta_df, sid_col, group_col):
    """Stacked bar plot at Phylum level."""
    from app.analysis.taxonomy import aggregate_taxonomy

    try:
        tax_df = aggregate_taxonomy(biom_path, "Phylum", top_n=15)
    except Exception as e:
        logger.warning(f"Taxonomy aggregation failed: {e}")
        return

    if tax_df.empty:
        return

    # Reorder by group if metadata available
    sample_order = list(tax_df.columns)
    group_labels = None
    if meta_df is not None and sid_col and group_col and group_col in meta_df.columns:
        meta_map = meta_df.set_index(meta_df[sid_col].astype(str))[group_col]
        mapped = {s: meta_map.get(s, "") for s in sample_order}
        sample_order = sorted(sample_order, key=lambda s: (mapped.get(s, ""), s))
        group_labels = [mapped.get(s, "") for s in sample_order]

    tax_df = tax_df[sample_order]

    # Stacked bar chart
    n_samples = len(sample_order)
    fig_width = max(10, n_samples * 0.45)
    fig, ax = plt.subplots(figsize=(min(fig_width, 22), 7))

    cmap = plt.cm.get_cmap("tab20", len(tax_df))
    bottom = np.zeros(n_samples)
    x = np.arange(n_samples)

    for i, (taxon, row) in enumerate(tax_df.iterrows()):
        vals = row.values
        ax.bar(x, vals, bottom=bottom, label=taxon, color=cmap(i), width=0.85)
        bottom += vals

    ax.set_xticks(x)
    x_labels = sample_order
    if group_labels:
        x_labels = [f"{s}\n({g})" for s, g in zip(sample_order, group_labels)]
    ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Relative Abundance", fontsize=11)
    ax.set_title("Taxonomy Composition (Phylum Level)", fontsize=14, fontweight="bold")
    ax.set_ylim(0, 1)
    ax.legend(
        bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7,
        ncol=1, frameon=False,
    )
    fig.tight_layout()

    chart_path = _save_temp_chart(fig)
    pdf.add_page()
    pdf.set_font("Helvetica", "B", 16)
    pdf.cell(text="Taxonomy Composition", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(4)
    pdf.image(chart_path, w=270)
    Path(chart_path).unlink(missing_ok=True)


def _add_reads_section(pdf: FPDF, sample_data: list[dict]):
    """Read tracking bar chart and retention table."""
    # Chart — split if many samples
    chart_paths = _make_read_tracking_charts(sample_data)
    first = True
    for chart_path in chart_paths:
        pdf.add_page()
        if first:
            pdf.set_font("Helvetica", "B", 16)
            pdf.cell(text="Read Tracking", new_x="LMARGIN", new_y="NEXT")
            pdf.ln(4)
            first = False
        pdf.image(chart_path, w=260)
        Path(chart_path).unlink(missing_ok=True)

    # Retention table
    pdf.add_page()
    pdf.set_font("Helvetica", "B", 14)
    pdf.cell(text="Read Retention Summary", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(4)

    headers = ["Sample", "Raw", "Filtered", "% Filtered",
               "Denoised", "Non-chimeric", "% Retained"]
    col_widths = [55, 30, 30, 28, 30, 35, 28]

    pdf.set_font("Helvetica", "B", 9)
    for h, w in zip(headers, col_widths):
        pdf.cell(w, 7, h, border=1)
    pdf.ln()

    pdf.set_font("Helvetica", "", 9)
    for s in sample_data:
        raw = s["raw"]
        nonchim = s["nonchim"]
        pct_filt = f"{s['filtered'] / raw * 100:.1f}%" if raw > 0 else "-"
        pct_ret = f"{nonchim / raw * 100:.1f}%" if raw > 0 else "-"
        row = [
            s["name"], f"{raw:,}", f"{s['filtered']:,}",
            pct_filt, f"{s['denoised']:,}", f"{nonchim:,}", pct_ret,
        ]
        for val, w in zip(row, col_widths):
            pdf.cell(w, 6, val, border=1)
        pdf.ln()


# ── Helpers ──────────────────────────────────────────────────────────────────

def _save_temp_chart(fig) -> str:
    """Save a matplotlib figure to a temp PNG and close it."""
    tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
    fig.savefig(tmp.name, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return tmp.name


def _make_read_tracking_charts(
    sample_data: list[dict], samples_per_page: int = 25,
) -> list[str]:
    """Create grouped bar chart PNGs, splitting across pages if needed."""
    total = len(sample_data)
    chunks = (
        [sample_data]
        if total <= samples_per_page
        else [sample_data[i:i + samples_per_page]
              for i in range(0, total, samples_per_page)]
    )

    paths = []
    for page_idx, chunk in enumerate(chunks):
        names = [s["name"] for s in chunk]
        n = len(names)
        x = range(n)
        width = 0.2

        fig, ax = plt.subplots(figsize=(max(8, n * 0.55), 5))
        ax.bar([i - 1.5 * width for i in x], [s["raw"] for s in chunk],
               width, label="Raw", color="#636EFA")
        ax.bar([i - 0.5 * width for i in x], [s["filtered"] for s in chunk],
               width, label="Filtered", color="#EF553B")
        ax.bar([i + 0.5 * width for i in x], [s["denoised"] for s in chunk],
               width, label="Denoised", color="#00CC96")
        ax.bar([i + 1.5 * width for i in x], [s["nonchim"] for s in chunk],
               width, label="Non-chimeric", color="#AB63FA")

        ax.set_xticks(list(x))
        ax.set_xticklabels(names, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Read Count")
        title = "Read Tracking Through Pipeline"
        if len(chunks) > 1:
            title += f"  (page {page_idx + 1}/{len(chunks)})"
        ax.set_title(title)
        ax.legend(loc="upper right", fontsize=8)
        fig.tight_layout()
        paths.append(_save_temp_chart(fig))

    return paths


def _get_group_colors(n: int) -> list[str]:
    """Return n distinct colors for groups."""
    palette = [
        "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
        "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52",
    ]
    if n <= len(palette):
        return palette[:n]
    cmap = plt.cm.get_cmap("tab20", n)
    return [matplotlib.colors.rgb2hex(cmap(i)) for i in range(n)]


def _format_pvalue(p: float) -> str:
    """Format p-value for display."""
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.4f}"
