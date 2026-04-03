"""
MicrobiomeDash — Differential Abundance analysis page.
"""
import io
import itertools
import traceback
import uuid

import dash_bootstrap_components as dbc
import pandas as pd
from biom import load_table
from dash import Input, Output, State, callback_context, dcc, html, no_update

from app.analysis.diff_abundance import (
    TOOL_LABELS,
    build_volcano,
    cancel_da_pairwise,
    read_da_pairwise_progress,
    run_pairwise_da_background,
)
from app.analysis.taxonomy import LEVEL_MAP
from app.analysis.shared import (
    find_metadata_for_samples,
    get_dataset_metadata_df,
    get_group_columns,
    get_pipeline_biom_options,
    parse_uploaded_biom,
    parse_uploaded_metadata,
    validate_metadata_vs_biom,
)
from app.dashboard.app import app as dash_app

import os as _os
_MAX_CPUS = min(32, _os.cpu_count() or 1)

TOOL_OPTIONS = [
    {"label": "ANCOM-BC2", "value": "ancombc"},
    {"label": "ALDEx2", "value": "aldex2"},
    {"label": "LinDA", "value": "linda"},
    {"label": "DESeq2", "value": "deseq2"},
    {"label": "MaAsLin2", "value": "maaslin2"},
]

DA_LEVEL_OPTIONS = [{"label": k, "value": k} for k in LEVEL_MAP.keys()]


def get_layout():
    pipeline_opts = get_pipeline_biom_options()
    return dbc.Container([
        html.H3("Differential Abundance", className="mb-2"),
        html.P(
            "Identify features that differ significantly between two groups.",
            className="text-muted mb-4",
        ),
        dcc.Store(id="da-biom-path"),
        dcc.Store(id="da-biom-name"),
        dcc.Store(id="da-meta-store"),
        dcc.Store(id="da-sample-id-col"),
        dcc.Store(id="da-results-csv"),
        dcc.Store(id="da-pairwise-job-id", storage_type="session"),
        dcc.Store(id="da-pairwise-csv"),
        dcc.Store(id="da-tool-store"),
        dcc.Interval(id="da-poll", interval=2000, disabled=True),

        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("Input Data"),
                    dbc.CardBody([
                        dbc.Label("BIOM Table", className="fw-bold"),
                        dbc.Select(
                            id="da-select-pipeline",
                            options=[{"label": "— none —", "value": ""}]
                            + pipeline_opts,
                            value="",
                            className="mb-2",
                        ),
                        html.Div("— or upload —", className="text-center text-muted small mb-2"),
                        dcc.Upload(
                            id="da-upload-biom",
                            children=html.Div(["Drag & drop or ", html.A("select .biom")]),
                            style={"borderWidth": "2px", "borderStyle": "dashed",
                                   "borderRadius": "5px", "borderColor": "#555",
                                   "textAlign": "center", "padding": "8px"},
                            multiple=False,
                        ),
                        html.Div(id="da-biom-status", className="mt-1 small mb-3"),

                        dbc.Label("Metadata (CSV/TSV)", className="fw-bold"),
                        html.Div(
                            "Auto-loaded from pipeline dataset if available. "
                            "Upload to override.",
                            className="text-muted small mb-1",
                        ),
                        dcc.Upload(
                            id="da-upload-meta",
                            children=html.Div(["Drag & drop or ", html.A("select metadata file")]),
                            style={"borderWidth": "2px", "borderStyle": "dashed",
                                   "borderRadius": "5px", "borderColor": "#555",
                                   "textAlign": "center", "padding": "8px"},
                            multiple=False,
                        ),
                        html.Div(id="da-meta-status", className="mt-1 small mb-3"),

                        dbc.Label("DA Method", className="fw-bold"),
                        dbc.RadioItems(
                            id="da-tool",
                            options=TOOL_OPTIONS,
                            value="aldex2",
                            className="mb-3",
                        ),

                        html.Div(id="da-cores-section", style={"display": "none"}, children=[
                            dbc.Label("CPU Cores", className="fw-bold"),
                            dbc.Input(
                                id="da-cores",
                                type="number",
                                value=_MAX_CPUS,
                                min=1,
                                max=_MAX_CPUS,
                                className="mb-3",
                            ),
                        ]),

                        dbc.Label("Taxonomic Level", className="fw-bold"),
                        dbc.Select(
                            id="da-tax-level",
                            options=DA_LEVEL_OPTIONS,
                            value="ASV",
                            className="mb-3",
                        ),

                        dbc.Label("Group Column", className="fw-bold"),
                        dbc.Select(id="da-group-col", placeholder="Select or upload metadata...",
                                   className="mb-2"),
                        html.Div(id="da-group-info", className="mb-2"),

                        # Single-pair selectors (hidden in pairwise mode)
                        html.Div(id="da-single-group-section", children=[
                            dbc.Label("Reference Group", className="fw-bold"),
                            dbc.Select(id="da-ref-group",
                                       placeholder="Select group column first...",
                                       className="mb-2"),
                            dbc.Label("Test Group", className="fw-bold"),
                            dbc.Select(id="da-test-group",
                                       placeholder="Select group column first...",
                                       className="mb-3"),
                        ]),

                        html.Hr(),
                        dbc.Switch(
                            id="da-pairwise-toggle",
                            label="Run all pairwise comparisons",
                            value=False,
                            className="mb-3",
                        ),

                        dbc.Button("Run Analysis", id="da-btn-run", color="primary",
                                   className="w-100", disabled=True),
                        dbc.Button("Cancel", id="da-btn-cancel", color="danger",
                                   className="w-100 mt-2",
                                   style={"display": "none"}),
                    ]),
                ], className="mb-3"),

                dbc.Card([
                    dbc.CardHeader("Filter Criteria"),
                    dbc.CardBody([
                        dbc.Label("p-value", className="fw-bold"),
                        dbc.Select(
                            id="da-filter-pval-raw",
                            options=[
                                {"label": "— None —", "value": ""},
                                {"label": "< 0.001", "value": "0.001"},
                                {"label": "< 0.01", "value": "0.01"},
                                {"label": "< 0.05", "value": "0.05"},
                                {"label": "< 0.1", "value": "0.1"},
                            ],
                            value="",
                            className="mb-2",
                        ),
                        dbc.Label("q-value (adj. p-value)", className="fw-bold"),
                        dbc.Select(
                            id="da-filter-pval",
                            options=[
                                {"label": "— None —", "value": ""},
                                {"label": "< 0.001", "value": "0.001"},
                                {"label": "< 0.01", "value": "0.01"},
                                {"label": "< 0.05", "value": "0.05"},
                                {"label": "< 0.1", "value": "0.1"},
                            ],
                            value="",
                            className="mb-2",
                        ),
                        dbc.Label("Abundance Diff (|log2FC|)", className="fw-bold"),
                        dbc.Select(
                            id="da-filter-lfc",
                            options=[
                                {"label": "— None —", "value": ""},
                                {"label": "> 0.5", "value": "0.5"},
                                {"label": "> 1.0", "value": "1.0"},
                                {"label": "> 2.0", "value": "2.0"},
                                {"label": "> 3.0", "value": "3.0"},
                            ],
                            value="",
                            className="mb-2",
                        ),
                        html.Div(id="da-filter-effect-section", style={"display": "none"}, children=[
                            dbc.Label("Effect Size (|effect|, ALDEx2 only)", className="fw-bold"),
                            dbc.Select(
                                id="da-filter-effect",
                                options=[
                                    {"label": "— None —", "value": ""},
                                    {"label": "> 0.5", "value": "0.5"},
                                    {"label": "> 1.0", "value": "1.0"},
                                    {"label": "> 1.5", "value": "1.5"},
                                    {"label": "> 2.0", "value": "2.0"},
                                ],
                                value="",
                                className="mb-2",
                            ),
                        ]),
                    ]),
                ], className="mb-3"),
            ], md=4),

            dbc.Col([
                # Single-pair results (existing)
                dbc.Spinner([
                    html.Div(id="da-error", className="mb-2"),
                    html.Div(id="da-summary", className="mb-2"),
                    dcc.Graph(id="da-volcano", style={"display": "none"},
                             config={"toImageButtonOptions": {"format": "svg", "scale": 2}}),
                    html.Div(id="da-results-table"),
                    html.Div(id="da-download-area"),
                ], color="primary"),

                # Pairwise progress section
                html.Div(id="da-progress-section", style={"display": "none"}, children=[
                    dbc.Progress(id="da-progress-bar", style={"display": "none"}),
                    html.Pre(id="da-progress-log",
                             style={"maxHeight": "350px", "overflowY": "auto",
                                    "backgroundColor": "#1e1e1e", "color": "#ccc",
                                    "padding": "10px", "borderRadius": "5px",
                                    "fontSize": "0.85rem"}),
                ]),

                # Pairwise combined results
                html.Div(id="da-pairwise-results"),
            ], md=8),
        ]),

        # Download components for pairwise
        dcc.Download(id="da-pairwise-download"),
    ], fluid=True)


# ── Helpers ──────────────────────────────────────────────────────────────────


def _apply_da_filters(df, pval_raw_thresh, pval_thresh, lfc_thresh, effect_thresh):
    """Apply filter criteria to a DA results DataFrame. Returns filtered copy."""
    out = df.copy()
    if pval_raw_thresh and "pvalue" in out.columns:
        out = out[out["pvalue"] < float(pval_raw_thresh)]
    if pval_thresh and "qvalue" in out.columns:
        out = out[out["qvalue"] < float(pval_thresh)]
    if lfc_thresh and "log2fc" in out.columns:
        out = out[out["log2fc"].abs() >= float(lfc_thresh)]
    if effect_thresh and "effect" in out.columns:
        out = out[out["effect"].abs() >= float(effect_thresh)]
    return out


def _build_pairwise_results_filtered(records_df, tool_name, header_alert,
                                     pval_raw_thresh, pval_thresh, lfc_thresh,
                                     effect_thresh):
    """Build the combined results UI for pairwise comparisons with filters applied."""
    if records_df.empty:
        return header_alert

    comparisons = records_df["comparison"].unique().tolist()

    MAX_VOLCANOS = 10
    summary_rows = []
    accordion_items = []
    for i, comp in enumerate(comparisons):
        sub = records_df[records_df["comparison"] == comp]
        filtered = _apply_da_filters(sub, pval_raw_thresh, pval_thresh,
                                     lfc_thresh, effect_thresh)
        n_feat = len(filtered)
        n_sig = int((filtered["qvalue"] < 0.05).sum()) if "qvalue" in filtered.columns else 0
        summary_rows.append({"Comparison": comp, "Shown": n_feat,
                             "Total": len(sub),
                             "Significant (q<0.05)": n_sig})
        if i < MAX_VOLCANOS:
            fig = build_volcano(filtered, f"{tool_name} — {comp}")
            accordion_items.append(dbc.AccordionItem(dcc.Graph(figure=fig), title=comp))

    summary_df = pd.DataFrame(summary_rows)
    summary_table = dbc.Table.from_dataframe(
        summary_df, striped=True, bordered=True, hover=True,
        color="dark", size="sm",
    )

    volcano_note = ""
    if len(comparisons) > MAX_VOLCANOS:
        volcano_note = html.Small(
            f"Showing volcano plots for first {MAX_VOLCANOS} of "
            f"{len(comparisons)} comparisons.",
            className="text-muted",
        )

    accordion = dbc.Accordion(
        accordion_items,
        start_collapsed=len(accordion_items) > 5,
        active_item="item-0" if accordion_items and len(accordion_items) <= 5 else None,
        className="mt-3",
    )

    # Build filtered results table (cap display at 500 rows for performance)
    all_filtered = pd.concat(
        [_apply_da_filters(records_df[records_df["comparison"] == c],
                           pval_raw_thresh, pval_thresh, lfc_thresh, effect_thresh)
         for c in comparisons],
        ignore_index=True,
    )
    n_total = len(all_filtered)
    MAX_DISPLAY = 500
    # Sort by qvalue so most significant are shown first
    if "qvalue" in all_filtered.columns:
        all_filtered = all_filtered.sort_values("qvalue")
    display_truncated = n_total > MAX_DISPLAY
    display_slice = all_filtered.head(MAX_DISPLAY)
    # Select display columns (skip internal/bulky ones)
    display_cols = [c for c in display_slice.columns
                    if c not in ("ref_group", "test_group", "neg_log10_q", "sig")]
    # Round numeric columns for readability
    display_df = display_slice[display_cols].copy()
    for col in display_df.select_dtypes(include="number").columns:
        display_df[col] = display_df[col].apply(
            lambda x: f"{x:.2e}" if abs(x) < 0.001 or abs(x) > 1000 else f"{x:.4f}"
            if pd.notna(x) else ""
        )
    results_table = dbc.Table.from_dataframe(
        display_df, striped=True, bordered=True, hover=True,
        color="dark", size="sm",
        style={"fontSize": "0.8rem"},
    )
    truncation_note = ""
    if display_truncated:
        truncation_note = f" (showing top {MAX_DISPLAY} by q-value; download TSV for all)"
    results_section = html.Div([
        html.H5(f"Results ({n_total} features{truncation_note})", className="mt-3"),
        html.Div(results_table, style={"maxHeight": "500px", "overflowY": "auto"}),
    ])

    download_btn = dbc.Button(
        "Download All Results (TSV)", id="da-btn-pairwise-download",
        color="secondary", size="sm", className="mt-3",
    )

    children = [
        header_alert,
        html.H5("Summary", className="mt-3"),
        summary_table,
        accordion,
    ]
    if volcano_note:
        children.append(volcano_note)
    children.extend([results_section, download_btn])
    return html.Div(children)


# ── Callbacks ────────────────────────────────────────────────────────────────


@dash_app.callback(
    Output("da-biom-path", "data"),
    Output("da-biom-name", "data"),
    Output("da-biom-status", "children"),
    Output("da-meta-store", "data"),
    Output("da-sample-id-col", "data"),
    Output("da-meta-status", "children"),
    Output("da-group-col", "options"),
    Output("da-group-info", "children", allow_duplicate=True),
    Output("da-btn-run", "disabled"),
    Input("da-upload-biom", "contents"),
    Input("da-select-pipeline", "value"),
    Input("da-upload-meta", "contents"),
    State("da-upload-biom", "filename"),
    State("da-upload-meta", "filename"),
    State("da-biom-path", "data"),
    State("da-meta-store", "data"),
    State("da-sample-id-col", "data"),
    prevent_initial_call=True,
)
def on_input_change(biom_contents, pipeline_value, meta_contents,
                    biom_filename, meta_filename,
                    prev_biom_path, prev_meta_json, prev_sid_col):
    trigger = callback_context.triggered_id

    biom_path = prev_biom_path
    biom_name = no_update
    biom_status = no_update
    meta_json = prev_meta_json
    sid_col = prev_sid_col
    meta_status = no_update
    group_opts = no_update
    group_info = no_update
    btn_disabled = no_update

    if trigger == "da-select-pipeline" and not pipeline_value:
        biom_path = None
        biom_name = None
        biom_status = ""
        meta_json = None
        sid_col = None
        meta_status = ""
        group_opts = []
        group_info = ""
        btn_disabled = True
        return (biom_path, biom_name, biom_status, meta_json, sid_col,
                meta_status, group_opts, group_info, btn_disabled)

    from pathlib import Path as _P

    if trigger == "da-select-pipeline" and pipeline_value:
        try:
            table = load_table(pipeline_value)
            n = len(table.ids(axis="sample"))
            biom_path = pipeline_value
            biom_name = _P(pipeline_value).stem
            biom_status = dbc.Alert(f"Pipeline dataset: {n} samples", color="success")

            db_meta, db_sid = get_dataset_metadata_df(pipeline_value)
            if db_meta is not None:
                meta_json = db_meta.to_json(date_format="iso", orient="split")
                sid_col = db_sid
                gcols = get_group_columns(db_meta, db_sid)
                group_opts = [{"label": c, "value": c} for c in gcols]
                meta_status = dbc.Alert(
                    f"Metadata auto-loaded: {len(db_meta)} samples, "
                    f"{len(gcols)} group columns",
                    color="success",
                )
                btn_disabled = len(gcols) == 0
            else:
                meta_json = None
                sid_col = None
                group_opts = []
                meta_status = dbc.Alert(
                    "No metadata found for this dataset. Upload a metadata file.",
                    color="info",
                )
                btn_disabled = True
        except Exception as e:
            biom_path = None
            biom_status = dbc.Alert(f"Error: {e}", color="danger")
            btn_disabled = True

    elif trigger == "da-upload-biom" and biom_contents:
        path, error = parse_uploaded_biom(biom_contents, biom_filename or "")
        if error:
            biom_path = None
            biom_status = dbc.Alert(error, color="danger")
            btn_disabled = True
        else:
            table = load_table(path)
            sample_ids = list(table.ids(axis="sample"))
            n = len(sample_ids)
            biom_path = path
            biom_name = _P(biom_filename).stem if biom_filename else "uploaded"
            biom_status = dbc.Alert(f"Uploaded BIOM: {n} samples", color="success")

            # Try to find matching metadata from pipeline datasets
            match_df, match_sid, match_name = find_metadata_for_samples(sample_ids)
            if match_df is not None:
                meta_json = match_df.to_json(date_format="iso", orient="split")
                sid_col = match_sid
                gcols = get_group_columns(match_df, match_sid)
                group_opts = [{"label": c, "value": c} for c in gcols]
                meta_status = dbc.Alert(
                    f"Metadata auto-matched from \"{match_name}\": "
                    f"{len(match_df)} samples, {len(gcols)} group columns",
                    color="success",
                )
                btn_disabled = len(gcols) == 0
            else:
                btn_disabled = meta_json is None

    elif trigger == "da-upload-meta" and meta_contents:
        df, s_col, error = parse_uploaded_metadata(meta_contents, meta_filename or "")
        if error:
            meta_json = None
            sid_col = None
            meta_status = dbc.Alert(error, color="danger")
            group_opts = []
            btn_disabled = True
        else:
            meta_json = df.to_json(date_format="iso", orient="split")
            sid_col = s_col
            gcols = get_group_columns(df, s_col)
            group_opts = [{"label": c, "value": c} for c in gcols]
            meta_status = dbc.Alert(f"Metadata loaded: {len(df)} samples", color="success")
            btn_disabled = biom_path is None or len(gcols) == 0

    return (biom_path, biom_name, biom_status, meta_json, sid_col,
            meta_status, group_opts, group_info, btn_disabled)


@dash_app.callback(
    Output("da-ref-group", "options"),
    Output("da-test-group", "options"),
    Output("da-group-info", "children"),
    Input("da-group-col", "value"),
    State("da-meta-store", "data"),
    State("da-sample-id-col", "data"),
    State("da-biom-path", "data"),
    prevent_initial_call=True,
)
def on_group_col_change(group_col, meta_json, sid_col, biom_path):
    if not group_col or not meta_json:
        return [], [], ""

    meta_df = pd.read_json(io.StringIO(meta_json), orient="split")

    # Intersect with BIOM samples if available
    if biom_path:
        try:
            table = load_table(biom_path)
            biom_ids = set(str(s) for s in table.ids(axis="sample"))
            meta_df = meta_df[meta_df[sid_col].astype(str).isin(biom_ids)]
        except Exception:
            pass

    counts = meta_df[group_col].dropna().value_counts()
    unique_vals = sorted(counts.index.tolist())
    options = [{"label": str(v), "value": str(v)} for v in unique_vals]

    # Build per-group sample count summary
    lines = [html.Span(f"{v}: {counts[v]} samples", className="d-block small")
             for v in unique_vals]

    low_groups = [str(v) for v in unique_vals if counts[v] < 3]
    if low_groups:
        info = dbc.Alert(
            [
                html.Strong("Warning: "),
                f"Groups with < 3 samples ({', '.join(low_groups)}) "
                "may cause DA tools to fail.",
                html.Div(lines, className="mt-1"),
            ],
            color="warning",
            className="py-2 small",
        )
    else:
        info = dbc.Alert(
            lines,
            color="info",
            className="py-2 small",
        )

    return options, options, info


# ── Pairwise toggle ─────────────────────────────────────────────────────────


@dash_app.callback(
    Output("da-single-group-section", "style"),
    Output("da-btn-run", "children"),
    Input("da-pairwise-toggle", "value"),
)
def on_pairwise_toggle(toggle_on):
    if toggle_on:
        return {"display": "none"}, "Run All Pairwise"
    return {"display": "block"}, "Run Analysis"


# ── Show/hide effect size filter based on tool ────────────────────────────────


@dash_app.callback(
    Output("da-filter-effect-section", "style"),
    Output("da-cores-section", "style"),
    Input("da-tool", "value"),
)
def on_tool_change_filter(tool):
    effect_style = {"display": "block"} if tool == "aldex2" else {"display": "none"}
    cores_style = {"display": "block"} if tool in ("ancombc", "aldex2", "deseq2", "maaslin2") else {"display": "none"}
    return effect_style, cores_style


# ── Resume pairwise job on page load ─────────────────────────────────────────


@dash_app.callback(
    Output("da-poll", "disabled", allow_duplicate=True),
    Output("da-progress-section", "style", allow_duplicate=True),
    Output("da-btn-cancel", "style", allow_duplicate=True),
    Output("da-btn-cancel", "disabled", allow_duplicate=True),
    Output("da-btn-run", "disabled", allow_duplicate=True),
    Output("da-progress-bar", "value", allow_duplicate=True),
    Output("da-progress-bar", "label", allow_duplicate=True),
    Output("da-progress-log", "children", allow_duplicate=True),
    Input("da-pairwise-job-id", "data"),
    prevent_initial_call="initial_duplicate",
)
def on_job_restore(job_id):
    """Restore polling UI when returning to a page with an active pairwise job."""
    no_all = (True, {"display": "none"}, {"display": "none"},
              False, no_update, no_update, no_update, no_update)
    if not job_id:
        return no_all

    prog = read_da_pairwise_progress(job_id)
    if prog is None:
        return no_all

    status = prog.get("status", "running")
    completed = prog.get("completed", 0)
    total = prog.get("total", 1)
    pct = int(completed / total * 100) if total else 0
    bar_label = f"{completed} of {total} comparisons"
    log_text = "\n".join(prog.get("log", []))

    if status == "running":
        # Job still running — enable polling, show progress, disable run
        return (False, {"display": "block"}, {"display": "block"},
                False, True, pct, bar_label, log_text)

    # Terminal state — enable one poll cycle to render final results
    # Keep progress section hidden; poll will render results directly
    return (False, {"display": "none"}, {"display": "none"},
            False, no_update, pct, bar_label, log_text)


# ── Run button (single-pair OR pairwise launch) ─────────────────────────────


@dash_app.callback(
    Output("da-volcano", "figure"),
    Output("da-volcano", "style"),
    Output("da-results-table", "children"),
    Output("da-summary", "children"),
    Output("da-error", "children"),
    Output("da-download-area", "children"),
    Output("da-results-csv", "data"),
    Output("da-pairwise-job-id", "data"),
    Output("da-poll", "disabled"),
    Output("da-progress-section", "style"),
    Output("da-btn-cancel", "style"),
    Output("da-progress-bar", "value"),
    Output("da-progress-bar", "label"),
    Output("da-progress-log", "children"),
    Output("da-pairwise-results", "children"),
    Output("da-btn-run", "disabled", allow_duplicate=True),
    Input("da-btn-run", "n_clicks"),
    State("da-biom-path", "data"),
    State("da-meta-store", "data"),
    State("da-sample-id-col", "data"),
    State("da-tool", "value"),
    State("da-group-col", "value"),
    State("da-ref-group", "value"),
    State("da-test-group", "value"),
    State("da-pairwise-toggle", "value"),
    State("da-tax-level", "value"),
    State("da-cores", "value"),
    State("da-pairwise-job-id", "data"),
    prevent_initial_call=True,
)
def on_run(n_clicks, biom_path, meta_json, sid_col, tool, group_col,
           ref_group, test_group, pairwise_on, tax_level, cores, prev_job_id):
    # 16 outputs: fig, fig_style, table, summary, error, dl_area, csv,
    #   job_id, poll_disabled, progress_style, cancel_style,
    #   bar_value, bar_label, log_text, pairwise_results, btn_disabled
    no_change = (no_update,) * 16

    # Cancel any previously running job before starting a new one
    if prev_job_id:
        cancel_da_pairwise(prev_job_id)

    if not all([biom_path, meta_json, tool, group_col]):
        return (*no_change[:4],
                dbc.Alert("Please fill all inputs.", color="warning"),
                *no_change[5:])

    # ── Pairwise mode ────────────────────────────────────────────────────
    if pairwise_on:
        try:
            meta_df = pd.read_json(io.StringIO(meta_json), orient="split")
            table = load_table(biom_path)
            biom_ids = list(table.ids(axis="sample"))
            match_info = validate_metadata_vs_biom(meta_df, sid_col, biom_ids)

            if len(match_info["matched"]) < 3:
                return (*no_change[:4],
                        dbc.Alert("Need at least 3 matched samples.", color="danger"),
                        *no_change[5:])

            groups = sorted(
                meta_df.loc[
                    meta_df[sid_col].isin(match_info["matched"]), group_col
                ].dropna().unique().tolist()
            )
            groups = [str(g) for g in groups]
            if len(groups) < 2:
                return (*no_change[:4],
                        dbc.Alert("Need at least 2 groups for pairwise comparisons.",
                                  color="danger"),
                        *no_change[5:])

            job_id = str(uuid.uuid4())
            level = tax_level or "ASV"
            n_threads = int(cores) if cores else None
            run_pairwise_da_background(
                tool, biom_path, meta_df, sid_col,
                match_info["matched"], group_col, groups, job_id,
                level=level, threads=n_threads,
            )

            n_pairs = len(list(itertools.combinations(groups, 2)))

            return (
                no_update,                          # fig
                {"display": "none"},                 # fig style (hide single-pair)
                "",                                  # table
                "",                                  # summary
                "",                                  # error
                "",                                  # dl area
                no_update,                           # csv
                job_id,                              # job id
                False,                               # poll enabled
                {"display": "block"},                # progress visible
                {"display": "block"},                # cancel visible
                0,                                   # bar value
                f"0 of {n_pairs} comparisons",       # bar label
                "",                                  # log
                "",                                  # pairwise results (clear)
                True,                                # run btn disabled
            )
        except Exception as e:
            return (*no_change[:4],
                    dbc.Alert(f"Error: {e}\n{traceback.format_exc()}", color="danger"),
                    *no_change[5:])

    # ── Single-pair mode — run asynchronously via background thread ──────
    if not all([ref_group, test_group]):
        return (*no_change[:4],
                dbc.Alert("Please select reference and test groups.", color="warning"),
                *no_change[5:])

    if ref_group == test_group:
        return (*no_change[:4],
                dbc.Alert("Reference and test groups must differ.", color="warning"),
                *no_change[5:])

    try:
        meta_df = pd.read_json(io.StringIO(meta_json), orient="split")
        table = load_table(biom_path)
        biom_ids = list(table.ids(axis="sample"))
        match_info = validate_metadata_vs_biom(meta_df, sid_col, biom_ids)

        if len(match_info["matched"]) < 3:
            return (*no_change[:4],
                    dbc.Alert("Need at least 3 matched samples.", color="danger"),
                    *no_change[5:])

        job_id = str(uuid.uuid4())
        level = tax_level or "ASV"
        n_threads = int(cores) if cores else None
        run_pairwise_da_background(
            tool, biom_path, meta_df, sid_col,
            match_info["matched"], group_col,
            [ref_group, test_group], job_id,
            level=level, threads=n_threads,
        )

        return (
            no_update,                          # fig
            {"display": "none"},                 # fig style (hide single-pair)
            "",                                  # table
            "",                                  # summary
            "",                                  # error
            "",                                  # dl area
            no_update,                           # csv
            job_id,                              # job id
            False,                               # poll enabled
            {"display": "block"},                # progress visible
            {"display": "block"},                # cancel visible
            0,                                   # bar value
            "0 of 1 comparisons",                # bar label
            "",                                  # log
            "",                                  # pairwise results (clear)
            True,                                # run btn disabled
        )

    except Exception as e:
        return (*no_change[:4],
                dbc.Alert(f"Error: {e}\n{traceback.format_exc()}", color="danger"),
                *no_change[5:])


# ── Pairwise poll ────────────────────────────────────────────────────────────


@dash_app.callback(
    Output("da-progress-bar", "value", allow_duplicate=True),
    Output("da-progress-bar", "label", allow_duplicate=True),
    Output("da-progress-log", "children", allow_duplicate=True),
    Output("da-poll", "disabled", allow_duplicate=True),
    Output("da-pairwise-results", "children", allow_duplicate=True),
    Output("da-pairwise-csv", "data"),
    Output("da-tool-store", "data", allow_duplicate=True),
    Output("da-btn-run", "disabled", allow_duplicate=True),
    Output("da-btn-cancel", "style", allow_duplicate=True),
    Output("da-progress-section", "style", allow_duplicate=True),
    Output("da-pairwise-job-id", "data", allow_duplicate=True),
    Input("da-poll", "n_intervals"),
    State("da-pairwise-job-id", "data"),
    State("da-tool", "value"),
    State("da-filter-pval-raw", "value"),
    State("da-filter-pval", "value"),
    State("da-filter-lfc", "value"),
    State("da-filter-effect", "value"),
    prevent_initial_call=True,
)
def on_poll(n_intervals, job_id, tool, filt_pval_raw, filt_pval, filt_lfc,
            filt_effect):
    # 11 outputs
    no_all = (no_update,) * 11
    if not job_id:
        return (no_update, no_update, no_update, True,
                *no_all[4:])

    prog = read_da_pairwise_progress(job_id)
    if prog is None:
        return no_all

    completed = prog.get("completed", 0)
    total = prog.get("total", 1)
    status = prog.get("status", "running")
    log_text = "\n".join(prog.get("log", []))
    pct = int(completed / total * 100) if total else 0
    bar_label = f"{completed} of {total} comparisons"

    if status == "running":
        return (pct, bar_label, log_text, False,
                no_update, no_update, no_update, no_update, no_update, no_update, no_update)

    # Terminal states: complete, cancelled, error — clear job_id from session
    # Use tool name from job progress (survives page reload) with radio button as fallback
    tool_name = prog.get("tool_label") or TOOL_LABELS.get(tool, tool)
    results_records = prog.get("results", [])
    records_df = pd.DataFrame(results_records) if results_records else pd.DataFrame()
    csv_data = _pairwise_csv(results_records)

    if status == "cancelled":
        msg = dbc.Alert(
            f"Cancelled after {completed} of {total} comparisons. "
            "Partial results shown below." if results_records else "No results.",
            color="warning",
        )
        results_ui = _build_pairwise_results_filtered(
            records_df, tool_name, msg, filt_pval_raw, filt_pval, filt_lfc,
            filt_effect)
        return (pct, bar_label, log_text, True,
                results_ui, csv_data, tool_name,
                False, {"display": "none"}, {"display": "none"}, None)

    if status == "error":
        msg = dbc.Alert("Error during pairwise run. Check log for details.", color="danger")
        results_ui = _build_pairwise_results_filtered(
            records_df, tool_name, msg, filt_pval_raw, filt_pval, filt_lfc,
            filt_effect)
        return (100, bar_label, log_text, True,
                results_ui, csv_data, tool_name,
                False, {"display": "none"}, {"display": "none"}, None)

    # Complete
    if records_df.empty:
        msg = dbc.Alert("All comparisons failed. Check log for details.", color="danger")
        return (100, bar_label, log_text, True,
                msg, None, tool_name,
                False, {"display": "none"}, {"display": "none"}, None)

    msg = dbc.Alert(
        f"{tool_name}: {total} pairwise comparisons complete.", color="success",
    )
    results_ui = _build_pairwise_results_filtered(
        records_df, tool_name, msg, filt_pval_raw, filt_pval, filt_lfc,
        filt_effect)
    return (100, bar_label, log_text, True,
            results_ui, csv_data, tool_name,
            False, {"display": "none"}, {"display": "none"}, None)


def _pairwise_csv(records):
    if not records:
        return None
    return pd.DataFrame(records).to_csv(sep="\t", index=False)


# ── Live filter on stored results ─────────────────────────────────────────────


@dash_app.callback(
    Output("da-pairwise-results", "children", allow_duplicate=True),
    Input("da-filter-pval-raw", "value"),
    Input("da-filter-pval", "value"),
    Input("da-filter-lfc", "value"),
    Input("da-filter-effect", "value"),
    State("da-pairwise-csv", "data"),
    State("da-tool-store", "data"),
    prevent_initial_call=True,
)
def on_da_filter_change(filt_pval_raw, filt_pval, filt_lfc, filt_effect,
                        csv_data, tool_name):
    if not csv_data or not tool_name:
        return no_update

    records_df = pd.read_csv(io.StringIO(csv_data), sep="\t")
    if records_df.empty:
        return no_update

    msg = dbc.Alert(f"{tool_name}: results filtered.", color="success")
    return _build_pairwise_results_filtered(
        records_df, tool_name, msg, filt_pval_raw, filt_pval, filt_lfc,
        filt_effect)


# ── Cancel ───────────────────────────────────────────────────────────────────


@dash_app.callback(
    Output("da-btn-cancel", "disabled"),
    Input("da-btn-cancel", "n_clicks"),
    State("da-pairwise-job-id", "data"),
    prevent_initial_call=True,
)
def on_cancel(n_clicks, job_id):
    if job_id:
        cancel_da_pairwise(job_id)
    return True


# ── Downloads ────────────────────────────────────────────────────────────────


def _da_download_filename(biom_name, tool_name):
    """Build download filename from BIOM name and tool name."""
    stem = biom_name or "results"
    tool_tag = (tool_name or "DA").replace(" ", "_").replace("-", "")
    return f"{stem}_{tool_tag}.tsv"


@dash_app.callback(
    Output("da-download", "data"),
    Input("da-btn-download", "n_clicks"),
    State("da-results-csv", "data"),
    State("da-biom-name", "data"),
    State("da-tool-store", "data"),
    prevent_initial_call=True,
)
def on_download(n_clicks, csv_data, biom_name, tool_name):
    if not n_clicks or not csv_data:
        return no_update
    return dcc.send_string(csv_data, _da_download_filename(biom_name, tool_name))


@dash_app.callback(
    Output("da-pairwise-download", "data"),
    Input("da-btn-pairwise-download", "n_clicks"),
    State("da-pairwise-csv", "data"),
    State("da-biom-name", "data"),
    State("da-tool-store", "data"),
    prevent_initial_call=True,
)
def on_pairwise_download(n_clicks, csv_data, biom_name, tool_name):
    if not n_clicks or not csv_data:
        return no_update
    return dcc.send_string(csv_data, _da_download_filename(biom_name, tool_name))
