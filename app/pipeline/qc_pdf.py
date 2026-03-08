"""
MicrobiomeDash — Generate a QC report PDF for a completed pipeline run.

Contains: read tracking bar chart, retention summary table, and dataset info.
Uses matplotlib (Agg backend) for charts and fpdf2 for PDF assembly.
"""
import logging
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fpdf import FPDF


def generate_qc_pdf(
    dataset_id: int,
    output_dir: Path,
    logger: logging.Logger,
) -> Path:
    """Generate a QC report PDF for the given dataset.

    Returns the path to the generated PDF file.
    """
    from app.db.database import get_session
    from app.db.models import Dataset, Sample

    with get_session() as db:
        dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
        if not dataset:
            raise ValueError(f"Dataset {dataset_id} not found")

        samples = (
            db.query(Sample)
            .filter(Sample.dataset_id == dataset_id)
            .order_by(Sample.sample_name)
            .all()
        )

        # Extract data while session is open
        ds_info = {
            "name": dataset.name,
            "status": dataset.status,
            "sequencing_type": dataset.sequencing_type or "unknown",
            "variable_region": dataset.variable_region or "unknown",
            "sample_count": dataset.sample_count or 0,
            "asv_count": dataset.asv_count or 0,
            "trunc_len_f": dataset.trunc_len_f,
            "trunc_len_r": dataset.trunc_len_r,
        }

        sample_data = []
        for s in samples:
            sample_data.append({
                "name": s.sample_name,
                "raw": s.read_count_raw or 0,
                "filtered": s.read_count_filtered or 0,
                "denoised": s.read_count_denoised or 0,
                "nonchim": s.read_count_nonchimeric or 0,
            })

    # Build PDF
    pdf = FPDF(orientation="L", format="A4")
    pdf.set_auto_page_break(auto=True, margin=15)

    # ── Page 1: Summary + Read Tracking Chart ─────────────────────────
    pdf.add_page()
    pdf.set_font("Helvetica", "B", 18)
    pdf.cell(text="QC Report", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(4)

    pdf.set_font("Helvetica", "", 11)
    pdf.cell(text=f"Dataset: {ds_info['name']}  (ID: {dataset_id})",
             new_x="LMARGIN", new_y="NEXT")
    pdf.cell(text=(
        f"Type: {ds_info['sequencing_type'].upper()}  |  "
        f"Region: {ds_info['variable_region']}  |  "
        f"Samples: {ds_info['sample_count']}  |  "
        f"ASVs: {ds_info['asv_count']}"
    ), new_x="LMARGIN", new_y="NEXT")

    if ds_info["trunc_len_f"]:
        trunc_info = f"Truncation: F={ds_info['trunc_len_f']}"
        if ds_info["trunc_len_r"]:
            trunc_info += f", R={ds_info['trunc_len_r']}"
        pdf.cell(text=trunc_info, new_x="LMARGIN", new_y="NEXT")

    pdf.ln(6)

    # Generate chart(s) as PNG and embed — split across pages if many samples
    if sample_data:
        chart_paths = _make_read_tracking_charts(sample_data)
        for i, chart_path in enumerate(chart_paths):
            if i > 0:
                pdf.add_page()
            pdf.image(chart_path, w=260)
            Path(chart_path).unlink(missing_ok=True)

    # ── Page 2: Retention Table ───────────────────────────────────────
    if sample_data:
        pdf.add_page()
        pdf.set_font("Helvetica", "B", 14)
        pdf.cell(text="Read Retention Summary", new_x="LMARGIN", new_y="NEXT")
        pdf.ln(4)
        _add_retention_table(pdf, sample_data)

    # Write PDF
    pdf_path = output_dir / "qc_report.pdf"
    pdf.output(str(pdf_path))
    logger.info(f"QC report PDF written: {pdf_path}")
    return pdf_path


def _make_read_tracking_charts(sample_data: list[dict],
                                samples_per_page: int = 25) -> list[str]:
    """Create grouped bar chart PNGs, splitting across pages if needed.

    Returns a list of temp file paths (one per page).
    """
    total = len(sample_data)
    if total <= samples_per_page:
        chunks = [sample_data]
    else:
        chunks = [
            sample_data[i:i + samples_per_page]
            for i in range(0, total, samples_per_page)
        ]

    paths = []
    for page_idx, chunk in enumerate(chunks):
        names = [s["name"] for s in chunk]
        raw = [s["raw"] for s in chunk]
        filtered = [s["filtered"] for s in chunk]
        denoised = [s["denoised"] for s in chunk]
        nonchim = [s["nonchim"] for s in chunk]

        n = len(names)
        x = range(n)
        width = 0.2

        fig, ax = plt.subplots(figsize=(max(8, n * 0.55), 5))
        ax.bar([i - 1.5 * width for i in x], raw, width,
               label="Raw", color="#636EFA")
        ax.bar([i - 0.5 * width for i in x], filtered, width,
               label="Filtered", color="#EF553B")
        ax.bar([i + 0.5 * width for i in x], denoised, width,
               label="Denoised", color="#00CC96")
        ax.bar([i + 1.5 * width for i in x], nonchim, width,
               label="Non-chimeric", color="#AB63FA")

        ax.set_xticks(list(x))
        ax.set_xticklabels(names, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Read Count")
        if len(chunks) > 1:
            ax.set_title(
                f"Read Tracking Through Pipeline  "
                f"(page {page_idx + 1}/{len(chunks)})")
        else:
            ax.set_title("Read Tracking Through Pipeline")
        ax.legend(loc="upper right", fontsize=8)
        fig.tight_layout()

        tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
        fig.savefig(tmp.name, dpi=150)
        plt.close(fig)
        paths.append(tmp.name)

    return paths


def _add_retention_table(pdf: FPDF, sample_data: list[dict]):
    """Add a read retention table to the PDF."""
    headers = ["Sample", "Raw", "Filtered", "% Filtered",
               "Denoised", "Non-chimeric", "% Retained"]
    col_widths = [55, 30, 30, 28, 30, 35, 28]

    # Header row
    pdf.set_font("Helvetica", "B", 9)
    for header, w in zip(headers, col_widths):
        pdf.cell(w, 7, header, border=1)
    pdf.ln()

    # Data rows
    pdf.set_font("Helvetica", "", 9)
    for s in sample_data:
        raw = s["raw"]
        filtered = s["filtered"]
        nonchim = s["nonchim"]
        pct_filt = f"{filtered / raw * 100:.1f}%" if raw > 0 else "-"
        pct_ret = f"{nonchim / raw * 100:.1f}%" if raw > 0 else "-"

        row = [
            s["name"],
            f"{raw:,}",
            f"{filtered:,}",
            pct_filt,
            f"{s['denoised']:,}",
            f"{nonchim:,}",
            pct_ret,
        ]
        for val, w in zip(row, col_widths):
            pdf.cell(w, 6, val, border=1)
        pdf.ln()
