"""
MicrobiomeDash — Rare ASV Removal page.

Remove low-prevalence and low-abundance ASVs from a BIOM file before
downstream analysis. Users set prevalence (% of samples) and minimum
total abundance thresholds, preview the filtering effect, and download
a cleaned BIOM file.
"""
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, no_update

from app.analysis.shared import get_pipeline_biom_options, parse_uploaded_biom
from app.dashboard.app import app as dash_app


def get_layout():
    pipeline_opts = get_pipeline_biom_options()

    return dbc.Container(
        [
            html.H3("Rare ASV Removal", className="mb-2"),
            html.P(
                "Remove rare or low-abundance ASVs from a BIOM file to reduce "
                "noise from sequencing artifacts and low-confidence taxa.",
                className="text-muted mb-4",
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Label("Pipeline Dataset", className="fw-bold"),
                            dbc.Select(
                                id="rar-pipeline-select",
                                options=[{"label": "— none —", "value": ""}]
                                + pipeline_opts,
                                value="",
                            ),
                        ],
                        md=6,
                    ),
                    dbc.Col(
                        [
                            dbc.Label("Or Upload BIOM", className="fw-bold"),
                            dcc.Upload(
                                id="rar-biom-upload",
                                children=html.Div(
                                    [
                                        "Drag & drop or ",
                                        html.A("select .biom", className="text-info"),
                                    ],
                                    className="text-center py-2",
                                ),
                                style={
                                    "borderWidth": "2px",
                                    "borderStyle": "dashed",
                                    "borderRadius": "5px",
                                    "borderColor": "#555",
                                },
                                accept=".biom,application/octet-stream,application/x-hdf5",
                            ),
                        ],
                        md=6,
                    ),
                ],
                className="mb-3",
            ),
            html.Div(id="rar-source-status", className="mb-2"),
            dbc.Button(
                "Load BIOM",
                id="rar-btn-load",
                color="primary",
                disabled=True,
                className="mb-4",
            ),
            html.Div(id="rar-main-area"),
            dcc.Download(id="rar-download-biom"),
            dcc.Store(id="rar-biom-path"),
            dcc.Store(id="rar-biom-filename"),
            dcc.Store(id="rar-asv-stats"),
        ],
        fluid=True,
    )


# ── Callback 1: resolve BIOM source ─────────────────────────────────────────


@dash_app.callback(
    Output("rar-biom-path", "data"),
    Output("rar-biom-filename", "data"),
    Output("rar-source-status", "children"),
    Output("rar-btn-load", "disabled"),
    Output("rar-pipeline-select", "value"),
    Input("rar-pipeline-select", "value"),
    Input("rar-biom-upload", "contents"),
    State("rar-biom-upload", "filename"),
    prevent_initial_call=True,
)
def rar_on_biom_input(pipeline_val, upload_contents, upload_filename):
    from pathlib import Path

    from dash import ctx

    triggered = ctx.triggered_id

    if triggered == "rar-biom-upload" and upload_contents:
        tmp_path, err = parse_uploaded_biom(upload_contents, upload_filename)
        if err:
            return no_update, no_update, dbc.Alert(err, color="danger"), True, ""
        # Strip .biom extension to get the base name for output
        base = Path(upload_filename).stem
        return (
            tmp_path, base,
            html.Span(f"Uploaded: {upload_filename}", className="text-success"),
            False, "",
        )

    if triggered == "rar-pipeline-select" and pipeline_val:
        # Look up the dataset name from DB
        from app.db.database import SessionLocal
        from app.db.models import Dataset

        db = SessionLocal()
        try:
            ds = db.query(Dataset).filter(Dataset.asv_table_path == pipeline_val).first()
            base = ds.name if ds else Path(pipeline_val).stem
        finally:
            db.close()
        return (
            pipeline_val, base,
            html.Span("Pipeline dataset selected", className="text-success"),
            False, no_update,
        )

    return None, no_update, "", True, no_update


# ── Callback 2: load BIOM → compute stats → build controls ──────────────────


@dash_app.callback(
    Output("rar-asv-stats", "data"),
    Output("rar-main-area", "children"),
    Input("rar-btn-load", "n_clicks"),
    State("rar-biom-path", "data"),
    prevent_initial_call=True,
)
def rar_on_load(n_clicks, biom_path):
    if not n_clicks or not biom_path:
        return no_update, no_update

    try:
        from app.data_manager.rare_asv import compute_asv_stats

        stats = compute_asv_stats(biom_path)

        controls = html.Div([
            dbc.Alert(
                f"Loaded: {stats['n_asvs']:,} ASVs, {stats['n_samples']:,} samples, "
                f"{stats['total_reads']:,} total reads",
                color="info",
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Label("Min Prevalence (%)"),
                            html.Small(
                                "ASV must appear in at least this % of samples",
                                className="text-muted d-block mb-1",
                            ),
                            dbc.Input(
                                id="rar-min-prevalence",
                                type="number",
                                min=0, max=100, step=0.1,
                                value=5,
                                placeholder="e.g. 5",
                            ),
                        ],
                        md=4,
                    ),
                    dbc.Col(
                        [
                            dbc.Label("Min Total Abundance"),
                            html.Small(
                                "Minimum total read count across all samples",
                                className="text-muted d-block mb-1",
                            ),
                            dbc.Input(
                                id="rar-min-abundance",
                                type="number",
                                min=0, step=1,
                                value=2,
                                placeholder="e.g. 2",
                            ),
                        ],
                        md=4,
                    ),
                ],
                className="mb-3",
            ),
            html.Div(id="rar-preview", className="mb-3"),
            dbc.Button(
                "Filter & Download",
                id="rar-btn-filter",
                color="success",
                size="lg",
                className="w-100 mb-3",
            ),
            html.Div(id="rar-result"),
        ])

        return stats, controls

    except Exception as e:
        return (
            no_update,
            dbc.Alert(f"Error loading BIOM: {e}", color="danger"),
        )


# ── Callback 3: preview filtering effect ────────────────────────────────────


@dash_app.callback(
    Output("rar-preview", "children"),
    Input("rar-min-prevalence", "value"),
    Input("rar-min-abundance", "value"),
    State("rar-asv-stats", "data"),
    prevent_initial_call=True,
)
def rar_on_preview(min_prev, min_abund, stats):
    if not stats:
        return no_update

    min_prev = float(min_prev) if min_prev is not None else 0
    min_abund = int(min_abund) if min_abund is not None else 0

    per_asv = stats["per_asv"]
    n_total = stats["n_asvs"]
    total_reads = stats["total_reads"]

    n_removed = 0
    reads_removed = 0
    for asv in per_asv:
        if asv["prevalence_pct"] < min_prev or asv["total_abundance"] < min_abund:
            n_removed += 1
            reads_removed += asv["total_abundance"]

    n_kept = n_total - n_removed
    reads_kept = total_reads - reads_removed
    pct_removed = 100.0 * n_removed / n_total if n_total > 0 else 0
    pct_reads = 100.0 * reads_removed / total_reads if total_reads > 0 else 0

    if n_removed == n_total:
        return dbc.Alert(
            "All ASVs would be removed! Lower the thresholds.",
            color="danger",
        )

    return dbc.Card(
        dbc.CardBody([
            html.Span(
                f"Removing {n_removed:,} of {n_total:,} ASVs ({pct_removed:.1f}%)",
                className="d-block",
            ),
            html.Span(
                f"Remaining reads: {reads_kept:,} of {total_reads:,} "
                f"({pct_reads:.1f}% reads lost)",
                className="d-block text-muted",
            ),
            html.Span(
                f"Kept: {n_kept:,} ASVs",
                className="d-block text-success",
            ),
        ]),
        color="dark",
        className="border-secondary",
    )


# ── Callback 4: filter & download ───────────────────────────────────────────


@dash_app.callback(
    Output("rar-result", "children"),
    Output("rar-download-biom", "data"),
    Input("rar-btn-filter", "n_clicks"),
    State("rar-biom-path", "data"),
    State("rar-biom-filename", "data"),
    State("rar-min-prevalence", "value"),
    State("rar-min-abundance", "value"),
    prevent_initial_call=True,
)
def rar_on_filter(n_clicks, biom_path, biom_filename, min_prev, min_abund):
    if not n_clicks or not biom_path:
        return no_update, no_update

    min_prev = float(min_prev) if min_prev is not None else 0
    min_abund = int(min_abund) if min_abund is not None else 0

    try:
        from app.data_manager.rare_asv import filter_rare_asvs

        result = filter_rare_asvs(biom_path, min_prev, min_abund)

        base = biom_filename or "filtered"
        p_tag = str(min_prev).rstrip("0").rstrip(".")
        a_tag = str(min_abund)
        filename = f"{base}_p{p_tag}_a{a_tag}.biom"
        msg = (
            f"Removed {result['n_removed']:,} ASVs "
            f"({result['n_kept']:,} kept, "
            f"{result['reads_kept']:,} reads remaining)."
        )
        return (
            dbc.Alert(msg, color="success"),
            dcc.send_bytes(result["biom_bytes"], filename=filename),
        )
    except Exception as e:
        return dbc.Alert(f"Filter failed: {e}", color="danger"), no_update
