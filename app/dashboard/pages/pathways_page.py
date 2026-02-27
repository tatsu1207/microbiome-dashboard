"""
MicrobiomeDash — Pathways (PICRUSt2) analysis page.

Inspired by ggpicrust2: multi-tool DA, KO-to-KEGG aggregation, errorbar,
heatmap, and PCA visualizations.
"""
import io
import traceback
import uuid

import dash_bootstrap_components as dbc
import pandas as pd
from dash import Input, Output, State, callback_context, dcc, html, no_update

from app.analysis.diff_abundance import TOOL_LABELS, build_volcano
from app.analysis.kegg_aggregation import aggregate_ko_to_pathways, annotate_pathway_results
from app.analysis.pathway_plots import (
    build_pathway_errorbar,
    build_pathway_heatmap,
    build_pathway_pca,
)
from app.analysis.pathways import (
    detect_prediction_label,
    load_picrust2_table,
    load_prediction_file,
    merge_descriptions,
    parse_uploaded_prediction_file,
    parse_picrust2_zip,
    read_pathway_da_progress,
    run_pathway_da_background,
)
from app.analysis.shared import (
    find_metadata_for_samples,
    get_group_columns,
    get_picrust2_run_options,
    parse_uploaded_metadata,
    validate_metadata_vs_biom,
)
from app.dashboard.app import app as dash_app

HIDDEN = {"display": "none"}

DA_TOOL_OPTIONS = [
    {"label": "ALDEx2", "value": "aldex2"},
    {"label": "DESeq2", "value": "deseq2"},
    {"label": "ANCOM-BC2", "value": "ancombc"},
    {"label": "LinDA", "value": "linda"},
    {"label": "MaAsLin2", "value": "maaslin2"},
]


def _try_auto_load_meta_for_run(run_id: str):
    """Try to auto-load metadata from the pipeline dataset that sourced this PICRUSt2 run."""
    from app.analysis.shared import get_dataset_metadata_df
    from app.db.database import SessionLocal
    from app.db.models import Picrust2Run

    try:
        run_id_int = int(run_id)
    except (ValueError, TypeError):
        return None, None

    db = SessionLocal()
    try:
        run = db.query(Picrust2Run).filter(Picrust2Run.id == run_id_int).first()
        if not run or not run.biom_path:
            return None, None
        return get_dataset_metadata_df(run.biom_path)
    finally:
        db.close()


def get_layout():
    run_opts = [{"label": "— None —", "value": ""}] + get_picrust2_run_options()
    return dbc.Container([
        html.H3("Pathway Analysis", className="mb-2"),
        html.P(
            "Analyze PICRUSt2 functional predictions. Compare pathway abundances between groups.",
            className="text-muted mb-4",
        ),
        # Stores
        dcc.Store(id="pw-picrust-dir"),
        dcc.Store(id="pw-uploaded-file"),
        dcc.Store(id="pw-meta-store"),
        dcc.Store(id="pw-sample-id-col"),
        dcc.Store(id="pw-results-csv"),
        dcc.Store(id="pw-pred-label"),
        dcc.Store(id="pw-job-id", storage_type="session"),
        dcc.Store(id="pw-counts-json"),
        dcc.Store(id="pw-group-col-store"),
        dcc.Store(id="pw-ref-store"),
        dcc.Store(id="pw-test-store"),
        dcc.Store(id="pw-da-tool-store"),
        dcc.Interval(id="pw-poll", interval=2000, disabled=True),

        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("Input Data"),
                    dbc.CardBody([
                        dbc.Label("PICRUSt2 Results", className="fw-bold"),
                        dbc.Select(
                            id="pw-select-run",
                            options=run_opts,
                            placeholder="Select a completed PICRUSt2 run...",
                            className="mb-2",
                        ),
                        html.Div("— or upload prediction file —", className="text-center text-muted small mb-2"),
                        dcc.Upload(
                            id="pw-upload-file",
                            children=html.Div(["Drag & drop or ", html.A("select prediction file (.tsv / .tsv.gz)")]),
                            style={"borderWidth": "2px", "borderStyle": "dashed",
                                   "borderRadius": "5px", "borderColor": "#555",
                                   "textAlign": "center", "padding": "8px"},
                            multiple=False,
                        ),
                        html.Div(id="pw-input-status", className="mt-1 small mb-3"),

                        dbc.Label("Metadata (CSV/TSV)", className="fw-bold"),
                        html.Div(
                            "Auto-loaded from matching pipeline datasets. "
                            "Upload to override.",
                            className="text-muted small mb-1",
                        ),
                        dcc.Upload(
                            id="pw-upload-meta",
                            children=html.Div(["Drag & drop or ", html.A("select metadata file")]),
                            style={"borderWidth": "2px", "borderStyle": "dashed",
                                   "borderRadius": "5px", "borderColor": "#555",
                                   "textAlign": "center", "padding": "8px"},
                            multiple=False,
                        ),
                        html.Div(id="pw-meta-status", className="mt-1 small mb-3"),

                        html.Div(id="pw-pred-type-section", children=[
                            dbc.Label("Prediction Type", className="fw-bold"),
                            dbc.RadioItems(
                                id="pw-pred-type",
                                options=[
                                    {"label": "MetaCyc Pathways", "value": "metacyc"},
                                    {"label": "KEGG Orthologs (KO)", "value": "ko"},
                                    {"label": "Enzyme Commission (EC)", "value": "ec"},
                                ],
                                value="metacyc",
                                className="mb-3",
                            ),
                        ]),

                        # KO aggregation toggle (visible only when pred_type=ko)
                        html.Div(id="pw-aggregate-section", style=HIDDEN, children=[
                            dbc.Switch(
                                id="pw-aggregate-ko",
                                label="Aggregate KOs to KEGG Pathways",
                                value=False,
                                className="mb-1",
                            ),
                            html.Div(id="pw-aggregate-status", className="small mb-2"),
                        ]),

                        dbc.Label("DA Method", className="fw-bold"),
                        dbc.RadioItems(
                            id="pw-da-tool",
                            options=DA_TOOL_OPTIONS,
                            value="aldex2",
                            className="mb-3",
                            inline=True,
                        ),

                        dbc.Label("Group Column", className="fw-bold"),
                        dbc.Select(id="pw-group-col", placeholder="Select or upload metadata...",
                                   className="mb-2"),

                        dbc.Label("Reference Group", className="fw-bold"),
                        dbc.Select(id="pw-ref-group", placeholder="Select group column first...",
                                   className="mb-2"),

                        dbc.Label("Test Group", className="fw-bold"),
                        dbc.Select(id="pw-test-group", placeholder="Select group column first...",
                                   className="mb-3"),

                        dbc.Button("Run Analysis", id="pw-btn-run", color="primary",
                                   className="w-100", disabled=True),
                    ]),
                ], className="mb-3"),

                dbc.Card([
                    dbc.CardHeader("Filter Criteria"),
                    dbc.CardBody([
                        dbc.Label("q-value (adj. p-value)", className="fw-bold"),
                        dbc.Select(
                            id="pw-filter-pval",
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
                        html.Div(id="pw-effect-filter-section", children=[
                            dbc.Label("Effect Size (|effect|)", className="fw-bold"),
                            dbc.Select(
                                id="pw-filter-effect",
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
                        dbc.Label("Abundance Diff (|log2FC|)", className="fw-bold"),
                        dbc.Select(
                            id="pw-filter-abund",
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
                    ]),
                ], className="mb-3"),
            ], md=4),

            dbc.Col([
                # Progress section (shown during background run)
                html.Div(id="pw-progress-section", style=HIDDEN, children=[
                    html.H5("Running Pathway DA...", className="mt-2"),
                    dbc.Progress(id="pw-progress-bar", value=100, striped=True,
                                 animated=True, className="mb-2",
                                 style={"height": "24px"}),
                    html.Pre(id="pw-progress-log",
                             style={"maxHeight": "200px", "overflowY": "auto",
                                    "fontSize": "0.8rem", "whiteSpace": "pre-wrap",
                                    "backgroundColor": "#1a1a2e", "color": "#ccc",
                                    "padding": "8px", "borderRadius": "4px"}),
                ]),

                # Results section
                html.Div(id="pw-error", className="mb-2"),
                html.Div(id="pw-summary", className="mb-2"),

                # Tabbed results
                html.Div(id="pw-tabs-section", style=HIDDEN, children=[
                    dbc.Tabs(id="pw-result-tabs", active_tab="tab-volcano", children=[
                        dbc.Tab(label="Volcano", tab_id="tab-volcano"),
                        dbc.Tab(label="Errorbar", tab_id="tab-errorbar"),
                        dbc.Tab(label="Heatmap", tab_id="tab-heatmap"),
                        dbc.Tab(label="PCA", tab_id="tab-pca"),
                    ], className="mb-2"),
                    html.Div(id="pw-tab-content"),
                ]),

                html.Div(id="pw-results-table"),
                html.Div(id="pw-download-area"),
            ], md=8),
        ]),
    ], fluid=True)


# ── Callbacks ────────────────────────────────────────────────────────────────


def _apply_filters(results_df, pval_thresh, effect_thresh, abund_thresh):
    """Apply filter criteria to results DataFrame. Returns filtered copy."""
    df = results_df.copy()
    if pval_thresh and "qvalue" in df.columns:
        df = df[df["qvalue"] < float(pval_thresh)]
    if effect_thresh and "effect" in df.columns:
        df = df[df["effect"].abs() >= float(effect_thresh)]
    if abund_thresh and "log2fc" in df.columns:
        df = df[df["log2fc"].abs() >= float(abund_thresh)]
    return df


def _build_results_table(filtered_df, n_total, pred_label):
    """Build summary alert, results table, and download button."""
    n_filtered = len(filtered_df)
    n_sig = int((filtered_df["qvalue"] < 0.05).sum()) if "qvalue" in filtered_df.columns else 0
    parts = [f"{pred_label} features: {n_filtered} shown"]
    if n_filtered < n_total:
        parts.append(f"(filtered from {n_total})")
    parts.append(f"— {n_sig} with q < 0.05")
    summary = dbc.Alert(" ".join(parts), color="info")

    display_df = filtered_df.head(50).copy()
    for col in ["log2fc", "effect", "pvalue", "qvalue"]:
        if col in display_df.columns:
            display_df[col] = display_df[col].map(
                lambda x: f"{x:.4f}" if pd.notna(x) else "NA"
            )
    table_component = html.Div([
        html.H5(f"Top Results ({min(50, n_filtered)} of {n_filtered})", className="mt-3"),
        dbc.Table.from_dataframe(
            display_df, striped=True, bordered=True, hover=True,
            color="dark", size="sm",
        ),
    ])

    dl = html.Div([
        dcc.Download(id="pw-download"),
        dbc.Button("Download Full Results (TSV)", id="pw-btn-download",
                   color="secondary", size="sm", className="mt-2"),
    ])

    return summary, table_component, dl


def _try_find_meta_by_samples(picrust_dir: str, pred_type: str):
    """Load a PICRUSt2 table and search DB for metadata matching its sample names."""
    try:
        counts_df, _ = load_picrust2_table(picrust_dir, pred_type)
        sample_ids = counts_df.columns.tolist()
        if not sample_ids:
            return None, None, None
        return find_metadata_for_samples(sample_ids)
    except Exception:
        return None, None, None


# ── Prediction type toggle → show/hide KO aggregation ────────────────────────


@dash_app.callback(
    Output("pw-aggregate-section", "style"),
    Input("pw-pred-type", "value"),
    prevent_initial_call=True,
)
def on_pred_type_toggle(pred_type):
    if pred_type == "ko":
        return {}
    return HIDDEN


# ── Effect size filter visibility (ALDEx2 only) ─────────────────────────────


@dash_app.callback(
    Output("pw-effect-filter-section", "style"),
    Input("pw-da-tool", "value"),
    prevent_initial_call=True,
)
def on_tool_toggle_effect_filter(tool):
    if tool == "aldex2":
        return {}
    return HIDDEN


# ── Input data change ────────────────────────────────────────────────────────


@dash_app.callback(
    Output("pw-picrust-dir", "data"),
    Output("pw-input-status", "children"),
    Output("pw-meta-store", "data"),
    Output("pw-sample-id-col", "data"),
    Output("pw-meta-status", "children"),
    Output("pw-group-col", "options"),
    Output("pw-btn-run", "disabled"),
    Output("pw-uploaded-file", "data"),
    Output("pw-pred-type-section", "style"),
    Input("pw-select-run", "value"),
    Input("pw-upload-file", "contents"),
    Input("pw-upload-meta", "contents"),
    State("pw-upload-file", "filename"),
    State("pw-upload-meta", "filename"),
    State("pw-picrust-dir", "data"),
    State("pw-meta-store", "data"),
    State("pw-sample-id-col", "data"),
    State("pw-uploaded-file", "data"),
    prevent_initial_call=True,
)
def on_input_change(run_id, file_contents, meta_contents,
                    file_filename, meta_filename,
                    prev_picrust_dir, prev_meta_json, prev_sid_col,
                    prev_uploaded_file):
    trigger = callback_context.triggered_id

    picrust_dir = prev_picrust_dir
    input_status = no_update
    meta_json = prev_meta_json
    sid_col = prev_sid_col
    meta_status = no_update
    group_opts = no_update
    btn_disabled = no_update
    uploaded_file = prev_uploaded_file
    pred_type_style = no_update

    if trigger == "pw-select-run":
        uploaded_file = None
        pred_type_style = {}  # show prediction type
        if not run_id:
            return None, "", None, None, "", [], True, None, {}
        picrust_dir = run_id
        input_status = dbc.Alert(f"PICRUSt2 run #{run_id} selected", color="success")

        db_meta, db_sid = _try_auto_load_meta_for_run(run_id)
        ds_name = None
        if db_meta is None:
            for pt in ["metacyc", "ko", "ec"]:
                db_meta, db_sid, ds_name = _try_find_meta_by_samples(run_id, pt)
                if db_meta is not None:
                    break

        if db_meta is not None:
            meta_json = db_meta.to_json(date_format="iso", orient="split")
            sid_col = db_sid
            gcols = get_group_columns(db_meta, db_sid)
            group_opts = [{"label": c, "value": c} for c in gcols]
            source = f" (from dataset '{ds_name}')" if ds_name else ""
            meta_status = dbc.Alert(
                f"Metadata auto-loaded{source}: {len(db_meta)} samples, "
                f"{len(gcols)} group columns",
                color="success",
            )
            btn_disabled = len(gcols) == 0
        else:
            meta_json = None
            sid_col = None
            group_opts = []
            meta_status = dbc.Alert(
                "No linked metadata found. Upload a metadata file.",
                color="info",
            )
            btn_disabled = True

    elif trigger == "pw-upload-file" and file_contents:
        picrust_dir = None  # clear run selection
        pred_type_style = HIDDEN  # hide prediction type

        temp_path, error = parse_uploaded_prediction_file(
            file_contents, file_filename or ""
        )
        if error:
            uploaded_file = None
            input_status = dbc.Alert(error, color="danger")
            btn_disabled = True
        else:
            uploaded_file = temp_path
            input_status = dbc.Alert(
                f"Prediction file loaded: {file_filename}", color="success"
            )

            # Auto-load metadata by matching sample names
            try:
                counts_df, _ = load_prediction_file(temp_path)
                sample_ids = counts_df.columns.tolist()
                db_meta, db_sid, ds_name = (
                    find_metadata_for_samples(sample_ids) if sample_ids
                    else (None, None, None)
                )
            except Exception:
                db_meta, db_sid, ds_name = None, None, None

            if db_meta is not None:
                meta_json = db_meta.to_json(date_format="iso", orient="split")
                sid_col = db_sid
                gcols = get_group_columns(db_meta, db_sid)
                group_opts = [{"label": c, "value": c} for c in gcols]
                source = f" (from dataset '{ds_name}')" if ds_name else ""
                meta_status = dbc.Alert(
                    f"Metadata auto-loaded{source}: {len(db_meta)} samples, "
                    f"{len(gcols)} group columns",
                    color="success",
                )
                btn_disabled = len(gcols) == 0
            else:
                meta_json = None
                sid_col = None
                group_opts = []
                meta_status = dbc.Alert(
                    "No matching metadata found. Upload a metadata file.",
                    color="info",
                )
                btn_disabled = True

    elif trigger == "pw-upload-meta" and meta_contents:
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
            has_data = picrust_dir is not None or uploaded_file is not None
            btn_disabled = (not has_data) or len(gcols) == 0

    return (picrust_dir, input_status, meta_json, sid_col, meta_status, group_opts,
            btn_disabled, uploaded_file, pred_type_style)


@dash_app.callback(
    Output("pw-ref-group", "options"),
    Output("pw-test-group", "options"),
    Input("pw-group-col", "value"),
    State("pw-meta-store", "data"),
    prevent_initial_call=True,
)
def on_group_col_change(group_col, meta_json):
    if not group_col or not meta_json:
        return [], []

    meta_df = pd.read_json(io.StringIO(meta_json), orient="split")
    unique_vals = sorted(meta_df[group_col].dropna().unique().tolist())
    options = [{"label": str(v), "value": str(v)} for v in unique_vals]
    return options, options


# ── Run button (launches background DA) ──────────────────────────────────────


@dash_app.callback(
    Output("pw-error", "children"),
    Output("pw-summary", "children"),
    Output("pw-tabs-section", "style"),
    Output("pw-results-table", "children"),
    Output("pw-download-area", "children"),
    Output("pw-results-csv", "data"),
    Output("pw-job-id", "data"),
    Output("pw-poll", "disabled"),
    Output("pw-progress-section", "style"),
    Output("pw-btn-run", "disabled", allow_duplicate=True),
    Output("pw-progress-log", "children"),
    Output("pw-counts-json", "data"),
    Output("pw-group-col-store", "data"),
    Output("pw-ref-store", "data"),
    Output("pw-test-store", "data"),
    Output("pw-da-tool-store", "data"),
    Output("pw-aggregate-status", "children"),
    Input("pw-btn-run", "n_clicks"),
    State("pw-picrust-dir", "data"),
    State("pw-uploaded-file", "data"),
    State("pw-meta-store", "data"),
    State("pw-sample-id-col", "data"),
    State("pw-pred-type", "value"),
    State("pw-group-col", "value"),
    State("pw-ref-group", "value"),
    State("pw-test-group", "value"),
    State("pw-da-tool", "value"),
    State("pw-aggregate-ko", "value"),
    prevent_initial_call=True,
)
def on_run(n_clicks, picrust_dir, uploaded_file, meta_json, sid_col, pred_type,
           group_col, ref_group, test_group, da_tool, aggregate_ko):
    # 17 outputs
    n_out = 17
    no_all = (no_update,) * n_out

    if not all([meta_json, group_col, ref_group, test_group]):
        return (dbc.Alert("Please fill all inputs.", color="warning"),
                *no_all[1:])

    if not picrust_dir and not uploaded_file:
        return (dbc.Alert("Please select a PICRUSt2 run or upload a prediction file.",
                          color="warning"), *no_all[1:])

    if ref_group == test_group:
        return (dbc.Alert("Reference and test groups must differ.", color="warning"),
                *no_all[1:])

    try:
        meta_df = pd.read_json(io.StringIO(meta_json), orient="split")

        if uploaded_file:
            counts_df, desc_df = load_prediction_file(uploaded_file)
            pred_label = detect_prediction_label(counts_df)
        else:
            counts_df, desc_df = load_picrust2_table(picrust_dir, pred_type)
            pred_label = {"metacyc": "MetaCyc", "ko": "KO", "ec": "EC"}.get(pred_type, pred_type)

        # KO → KEGG pathway aggregation
        agg_status = no_update
        if pred_type == "ko" and aggregate_ko and not uploaded_file:
            try:
                pw_counts, pw_desc = aggregate_ko_to_pathways(counts_df)
                n_kos = len(counts_df)
                n_pathways = len(pw_counts)
                agg_status = dbc.Alert(
                    f"{n_kos} KOs mapped to {n_pathways} KEGG pathways",
                    color="info",
                )
                counts_df = pw_counts
                desc_df = pw_desc
                pred_label = "KEGG Pathway"
            except Exception as e:
                agg_status = dbc.Alert(f"KO aggregation failed: {e}", color="warning")

        biom_sample_ids = counts_df.columns.tolist()
        match_info = validate_metadata_vs_biom(meta_df, sid_col, biom_sample_ids)
        if len(match_info["matched"]) < 3:
            return (dbc.Alert("Need at least 3 matched samples.", color="danger"),
                    *no_all[1:])

        job_id = str(uuid.uuid4())
        tool = da_tool or "aldex2"
        tool_label = TOOL_LABELS.get(tool, tool)

        # Store counts for plot tabs
        # Limit size: only keep matched samples
        matched = match_info["matched"]
        counts_for_store = counts_df[[s for s in matched if s in counts_df.columns]]
        counts_json = counts_for_store.to_json(orient="split")

        run_pathway_da_background(
            counts_df, desc_df, meta_df, sid_col, matched,
            group_col, ref_group, test_group, pred_label, job_id,
            tool=tool,
        )

        return (
            "",                                  # error (clear)
            "",                                  # summary (clear)
            HIDDEN,                              # tabs hidden
            "",                                  # table (clear)
            "",                                  # dl area (clear)
            no_update,                           # csv
            job_id,                              # job id
            False,                               # poll enabled
            {"display": "block"},                # progress visible
            True,                                # run btn disabled
            f"Running {tool_label} on {pred_label} features...",  # progress log
            counts_json,                         # counts store
            group_col,                           # group col store
            ref_group,                           # ref store
            test_group,                          # test store
            tool,                                # da tool store
            agg_status,                          # aggregation status
        )

    except Exception as e:
        return (dbc.Alert(f"Error: {e}\n{traceback.format_exc()}", color="danger"),
                *no_all[1:])


# ── Poll for background results ──────────────────────────────────────────────


@dash_app.callback(
    Output("pw-tab-content", "children", allow_duplicate=True),
    Output("pw-tabs-section", "style", allow_duplicate=True),
    Output("pw-results-table", "children", allow_duplicate=True),
    Output("pw-summary", "children", allow_duplicate=True),
    Output("pw-error", "children", allow_duplicate=True),
    Output("pw-download-area", "children", allow_duplicate=True),
    Output("pw-results-csv", "data", allow_duplicate=True),
    Output("pw-pred-label", "data", allow_duplicate=True),
    Output("pw-poll", "disabled", allow_duplicate=True),
    Output("pw-progress-section", "style", allow_duplicate=True),
    Output("pw-btn-run", "disabled", allow_duplicate=True),
    Output("pw-job-id", "data", allow_duplicate=True),
    Output("pw-progress-log", "children", allow_duplicate=True),
    Input("pw-poll", "n_intervals"),
    State("pw-job-id", "data"),
    State("pw-filter-pval", "value"),
    State("pw-filter-effect", "value"),
    State("pw-filter-abund", "value"),
    prevent_initial_call=True,
)
def on_poll(n_intervals, job_id, filt_pval, filt_effect, filt_abund):
    # 13 outputs
    no_all = (no_update,) * 13
    if not job_id:
        return no_all

    prog = read_pathway_da_progress(job_id)
    if prog is None:
        return no_all

    status = prog.get("status", "running")
    log_text = "\n".join(prog.get("log", []))

    if status == "running":
        return (*no_all[:12], log_text)

    if status == "error":
        return (
            no_update,                           # tab content
            no_update,                           # tabs style
            no_update,                           # table
            no_update,                           # summary
            dbc.Alert(f"Pathway DA failed: {log_text}", color="danger"),
            no_update,                           # dl area
            no_update,                           # csv
            no_update,                           # pred_label
            True,                                # poll disabled
            HIDDEN,                              # progress hidden
            False,                               # btn enabled
            None,                                # clear job id
            log_text,                            # log
        )

    # Complete — build results
    csv_data = prog.get("results_csv", "")
    pred_label = prog.get("pred_label", "Pathway")

    results_df = pd.read_csv(io.StringIO(csv_data), sep="\t") if csv_data else pd.DataFrame()
    filtered_df = _apply_filters(results_df, filt_pval, filt_effect, filt_abund)

    # Build volcano for initial tab
    fig = build_volcano(filtered_df, f"Pathway DA ({pred_label})")
    tab_content = dcc.Graph(figure=fig)

    summary, table_component, dl = _build_results_table(
        filtered_df, len(results_df), pred_label
    )

    return (
        tab_content,                         # tab content (volcano)
        {"display": "block"},                # tabs visible
        table_component,                     # table
        summary,                             # summary
        "",                                  # error (clear)
        dl,                                  # dl area
        csv_data,                            # csv
        pred_label,                          # pred_label
        True,                                # poll disabled
        HIDDEN,                              # progress hidden
        False,                               # btn enabled
        None,                                # clear job id
        log_text,                            # log
    )


# ── Tab switching (lazy rendering) ───────────────────────────────────────────


@dash_app.callback(
    Output("pw-tab-content", "children", allow_duplicate=True),
    Input("pw-result-tabs", "active_tab"),
    State("pw-results-csv", "data"),
    State("pw-pred-label", "data"),
    State("pw-counts-json", "data"),
    State("pw-meta-store", "data"),
    State("pw-sample-id-col", "data"),
    State("pw-group-col-store", "data"),
    State("pw-ref-store", "data"),
    State("pw-test-store", "data"),
    State("pw-filter-pval", "value"),
    State("pw-filter-effect", "value"),
    State("pw-filter-abund", "value"),
    prevent_initial_call=True,
)
def on_tab_switch(active_tab, csv_data, pred_label, counts_json, meta_json,
                  sid_col, group_col, ref_group, test_group,
                  filt_pval, filt_effect, filt_abund):
    if not csv_data or not pred_label:
        return no_update

    results_df = pd.read_csv(io.StringIO(csv_data), sep="\t")
    filtered_df = _apply_filters(results_df, filt_pval, filt_effect, filt_abund)

    if active_tab == "tab-volcano":
        fig = build_volcano(filtered_df, f"Pathway DA ({pred_label})")
        return dcc.Graph(figure=fig)

    # For other tabs we need counts + metadata
    if not counts_json or not meta_json or not sid_col or not group_col:
        return dbc.Alert("Counts or metadata not available for this visualization.", color="warning")

    try:
        counts_df = pd.read_json(io.StringIO(counts_json), orient="split")
        meta_df = pd.read_json(io.StringIO(meta_json), orient="split")
    except Exception as e:
        return dbc.Alert(f"Error loading data: {e}", color="danger")

    try:
        if active_tab == "tab-errorbar":
            fig = build_pathway_errorbar(
                filtered_df, counts_df, meta_df, sid_col, group_col,
                ref_group, test_group, top_n=20,
            )
            return dcc.Graph(figure=fig)

        elif active_tab == "tab-heatmap":
            fig = build_pathway_heatmap(
                filtered_df, counts_df, meta_df, sid_col, group_col,
                top_n=30, cluster_rows=True,
            )
            return dcc.Graph(figure=fig)

        elif active_tab == "tab-pca":
            fig = build_pathway_pca(
                counts_df, meta_df, sid_col, group_col,
                ref_group, test_group,
            )
            return dcc.Graph(figure=fig)
    except Exception as e:
        return dbc.Alert(f"Error building plot: {e}\n{traceback.format_exc()}", color="danger")

    return no_update


# ── Live filter on stored results ─────────────────────────────────────────────


@dash_app.callback(
    Output("pw-tab-content", "children", allow_duplicate=True),
    Output("pw-results-table", "children", allow_duplicate=True),
    Output("pw-summary", "children", allow_duplicate=True),
    Output("pw-download-area", "children", allow_duplicate=True),
    Input("pw-filter-pval", "value"),
    Input("pw-filter-effect", "value"),
    Input("pw-filter-abund", "value"),
    State("pw-results-csv", "data"),
    State("pw-pred-label", "data"),
    State("pw-result-tabs", "active_tab"),
    State("pw-counts-json", "data"),
    State("pw-meta-store", "data"),
    State("pw-sample-id-col", "data"),
    State("pw-group-col-store", "data"),
    State("pw-ref-store", "data"),
    State("pw-test-store", "data"),
    prevent_initial_call=True,
)
def on_filter_change(filt_pval, filt_effect, filt_abund, csv_data, pred_label,
                     active_tab, counts_json, meta_json, sid_col,
                     group_col, ref_group, test_group):
    if not csv_data or not pred_label:
        return no_update, no_update, no_update, no_update

    results_df = pd.read_csv(io.StringIO(csv_data), sep="\t")
    filtered_df = _apply_filters(results_df, filt_pval, filt_effect, filt_abund)

    summary, table_component, dl = _build_results_table(
        filtered_df, len(results_df), pred_label
    )

    # Re-render active tab with new filters
    tab_content = no_update
    try:
        if active_tab == "tab-volcano":
            fig = build_volcano(filtered_df, f"Pathway DA ({pred_label})")
            tab_content = dcc.Graph(figure=fig)
        elif counts_json and meta_json and sid_col and group_col:
            counts_df = pd.read_json(io.StringIO(counts_json), orient="split")
            meta_df = pd.read_json(io.StringIO(meta_json), orient="split")

            if active_tab == "tab-errorbar":
                fig = build_pathway_errorbar(
                    filtered_df, counts_df, meta_df, sid_col, group_col,
                    ref_group, test_group, top_n=20,
                )
                tab_content = dcc.Graph(figure=fig)
            elif active_tab == "tab-heatmap":
                fig = build_pathway_heatmap(
                    filtered_df, counts_df, meta_df, sid_col, group_col,
                    top_n=30, cluster_rows=True,
                )
                tab_content = dcc.Graph(figure=fig)
            elif active_tab == "tab-pca":
                fig = build_pathway_pca(
                    counts_df, meta_df, sid_col, group_col,
                    ref_group, test_group,
                )
                tab_content = dcc.Graph(figure=fig)
    except Exception:
        pass  # Keep existing tab content on error

    return tab_content, table_component, summary, dl


# ── Restore job on page load ─────────────────────────────────────────────────


@dash_app.callback(
    Output("pw-poll", "disabled", allow_duplicate=True),
    Output("pw-progress-section", "style", allow_duplicate=True),
    Output("pw-progress-log", "children", allow_duplicate=True),
    Output("pw-btn-run", "disabled", allow_duplicate=True),
    Input("pw-job-id", "data"),
    prevent_initial_call="initial_duplicate",
)
def on_job_restore(job_id):
    """Restore polling UI when returning to a page with an active pathway DA job."""
    defaults = (True, HIDDEN, no_update, no_update)
    if not job_id:
        return defaults

    prog = read_pathway_da_progress(job_id)
    if prog is None:
        return defaults

    status = prog.get("status", "running")
    log_text = "\n".join(prog.get("log", []))

    if status == "running":
        return (False, {"display": "block"}, log_text, True)

    # Terminal state — enable one poll cycle to render final results
    return (False, {"display": "block"}, log_text, no_update)


# ── Download ──────────────────────────────────────────────────────────────────


@dash_app.callback(
    Output("pw-download", "data"),
    Input("pw-btn-download", "n_clicks"),
    State("pw-results-csv", "data"),
    prevent_initial_call=True,
)
def on_download(n_clicks, csv_data):
    if not n_clicks or not csv_data:
        return no_update
    return dcc.send_string(csv_data, "pathway_results.tsv")
