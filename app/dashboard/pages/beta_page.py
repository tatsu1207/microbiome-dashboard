"""
MicrobiomeDash — Beta Diversity analysis page.
"""
import io
import traceback
import uuid

import dash_bootstrap_components as dbc
import pandas as pd
import plotly.graph_objects as go
from biom import load_table
from dash import Input, Output, State, callback_context, dcc, html, no_update

from app.analysis.beta import (
    cleanup_permanova_progress,
    compute_confidence_ellipse,
    compute_distance,
    read_permanova_progress,
    run_nmds,
    run_pcoa,
    run_permanova_global,
    run_permanova_pairwise_background,
)
from app.analysis.shared import (
    biom_to_count_df,
    find_metadata_for_samples,
    get_dataset_metadata_df,
    get_group_columns,
    get_pipeline_biom_options,
    parse_uploaded_biom,
    parse_uploaded_metadata,
    validate_metadata_vs_biom,
)
from app.dashboard.app import app as dash_app


COLORS = ["#3498db", "#e74c3c", "#2ecc71", "#f39c12", "#9b59b6",
          "#1abc9c", "#e67e22", "#34495e"]
SYMBOLS = ["circle", "square", "diamond", "cross", "triangle-up",
           "triangle-down", "star", "hexagon"]

HIDDEN = {"display": "none"}


def get_layout():
    pipeline_opts = get_pipeline_biom_options()
    return dbc.Container([
        html.H3("Beta Diversity", className="mb-2"),
        html.P(
            "Compute between-sample distances, perform ordination, and run PERMANOVA.",
            className="text-muted mb-4",
        ),
        dcc.Store(id="bd-biom-path"),
        dcc.Store(id="bd-meta-store"),
        dcc.Store(id="bd-sample-id-col"),
        dcc.Store(id="bd-permanova-job-id", storage_type="session"),
        dcc.Interval(id="bd-poll", interval=2000, disabled=True),

        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("Input Data"),
                    dbc.CardBody([
                        dbc.Label("BIOM Table", className="fw-bold"),
                        dbc.Select(
                            id="bd-select-pipeline",
                            options=[{"label": "— none —", "value": ""}]
                            + pipeline_opts,
                            value="",
                            className="mb-2",
                        ),
                        html.Div("— or upload —", className="text-center text-muted small mb-2"),
                        dcc.Upload(
                            id="bd-upload-biom",
                            children=html.Div(["Drag & drop or ", html.A("select .biom")]),
                            style={"borderWidth": "2px", "borderStyle": "dashed",
                                   "borderRadius": "5px", "borderColor": "#555",
                                   "textAlign": "center", "padding": "8px"},
                            multiple=False,
                        ),
                        html.Div(id="bd-biom-status", className="mt-1 small mb-3"),

                        dbc.Label("Metadata (CSV/TSV)", className="fw-bold"),
                        html.Div(
                            "Auto-loaded from pipeline dataset if available. "
                            "Upload to override.",
                            className="text-muted small mb-1",
                        ),
                        dcc.Upload(
                            id="bd-upload-meta",
                            children=html.Div(["Drag & drop or ", html.A("select metadata file")]),
                            style={"borderWidth": "2px", "borderStyle": "dashed",
                                   "borderRadius": "5px", "borderColor": "#555",
                                   "textAlign": "center", "padding": "8px"},
                            multiple=False,
                        ),
                        html.Div(id="bd-meta-status", className="mt-1 small mb-3"),

                        dbc.Label("Distance Metric", className="fw-bold"),
                        dbc.RadioItems(
                            id="bd-distance-metric",
                            options=[
                                {"label": "Bray-Curtis", "value": "braycurtis"},
                                {"label": "Jaccard", "value": "jaccard"},
                            ],
                            value="braycurtis",
                            inline=True,
                            className="mb-3",
                        ),

                        dbc.Label("Ordination Method", className="fw-bold"),
                        dbc.RadioItems(
                            id="bd-ordination",
                            options=[
                                {"label": "PCoA", "value": "pcoa"},
                                {"label": "NMDS", "value": "nmds"},
                            ],
                            value="pcoa",
                            inline=True,
                            className="mb-3",
                        ),

                        dbc.Label("Color by", className="fw-bold"),
                        dbc.Select(id="bd-color-col", placeholder="Select column...",
                                   className="mb-2"),

                        dbc.Label("Shape by (optional)", className="fw-bold"),
                        dbc.Select(id="bd-shape-col", placeholder="None",
                                   className="mb-2"),

                        dbc.Label("PERMANOVA Group Column", className="fw-bold"),
                        dbc.Select(id="bd-group-col", placeholder="Select column...",
                                   className="mb-3"),

                        dbc.Button("Run Analysis", id="bd-btn-run", color="primary",
                                   className="w-100", disabled=True),
                    ]),
                ], className="mb-3"),
            ], md=4),

            dbc.Col([
                dcc.Loading(
                    id="bd-loading",
                    type="default",
                    children=html.Div(id="bd-loading-target"),
                    className="mb-3",
                ),
                html.Div(id="bd-error", className="mb-2"),
                html.Div(id="bd-progress-section", style=HIDDEN, children=[
                    html.H5("Pairwise PERMANOVA", className="mt-3 mb-2"),
                    dbc.Progress(id="bd-progress-bar", value=0, striped=True,
                                 animated=True, className="mb-2",
                                 style={"height": "24px"}),
                    html.Pre(id="bd-progress-log",
                             style={"maxHeight": "200px", "overflowY": "auto",
                                    "fontSize": "0.8rem", "whiteSpace": "pre-wrap",
                                    "backgroundColor": "#1a1a2e", "color": "#ccc",
                                    "padding": "8px", "borderRadius": "4px"}),
                ]),
                html.Div(id="bd-pairwise-table"),
                html.Div(id="bd-permanova-table"),
                dcc.Graph(id="bd-scatter", style=HIDDEN),
            ], md=8),
        ]),
    ], fluid=True)


# ── Callbacks ────────────────────────────────────────────────────────────────


def _meta_outputs(df, sid_col):
    """Build common metadata outputs for beta page."""
    gcols = get_group_columns(df, sid_col)
    options = [{"label": c, "value": c} for c in gcols]
    shape_opts = [{"label": "None", "value": ""}] + options
    return options, shape_opts, options


@dash_app.callback(
    Output("bd-biom-path", "data"),
    Output("bd-biom-status", "children"),
    Output("bd-meta-store", "data"),
    Output("bd-sample-id-col", "data"),
    Output("bd-meta-status", "children"),
    Output("bd-color-col", "options"),
    Output("bd-shape-col", "options"),
    Output("bd-group-col", "options"),
    Output("bd-btn-run", "disabled"),
    Input("bd-upload-biom", "contents"),
    Input("bd-select-pipeline", "value"),
    Input("bd-upload-meta", "contents"),
    State("bd-upload-biom", "filename"),
    State("bd-upload-meta", "filename"),
    State("bd-biom-path", "data"),
    State("bd-meta-store", "data"),
    State("bd-sample-id-col", "data"),
    prevent_initial_call=True,
)
def on_input_change(biom_contents, pipeline_value, meta_contents,
                    biom_filename, meta_filename,
                    prev_biom_path, prev_meta_json, prev_sid_col):
    trigger = callback_context.triggered_id

    biom_path = prev_biom_path
    biom_status = no_update
    meta_json = prev_meta_json
    sid_col = prev_sid_col
    meta_status = no_update
    color_opts = no_update
    shape_opts = no_update
    group_opts = no_update
    btn_disabled = no_update

    if trigger == "bd-select-pipeline" and not pipeline_value:
        biom_path = None
        biom_status = ""
        meta_json = None
        sid_col = None
        meta_status = ""
        color_opts = []
        shape_opts = [{"label": "None", "value": ""}]
        group_opts = []
        btn_disabled = True
        return (biom_path, biom_status, meta_json, sid_col, meta_status,
                color_opts, shape_opts, group_opts, btn_disabled)

    if trigger == "bd-select-pipeline" and pipeline_value:
        try:
            table = load_table(pipeline_value)
            n = len(table.ids(axis="sample"))
            biom_path = pipeline_value
            biom_status = dbc.Alert(f"Pipeline dataset: {n} samples", color="success")

            db_meta, db_sid = get_dataset_metadata_df(pipeline_value)
            if db_meta is not None:
                meta_json = db_meta.to_json(date_format="iso", orient="split")
                sid_col = db_sid
                color_opts, shape_opts, group_opts = _meta_outputs(db_meta, db_sid)
                gcols = get_group_columns(db_meta, db_sid)
                meta_status = dbc.Alert(
                    f"Metadata auto-loaded: {len(db_meta)} samples, "
                    f"{len(gcols)} group columns",
                    color="success",
                )
                btn_disabled = len(gcols) == 0
            else:
                meta_json = None
                sid_col = None
                color_opts, shape_opts, group_opts = [], [{"label": "None", "value": ""}], []
                meta_status = dbc.Alert(
                    "No metadata found for this dataset. Upload a metadata file.",
                    color="info",
                )
                btn_disabled = True
        except Exception as e:
            biom_path = None
            biom_status = dbc.Alert(f"Error: {e}", color="danger")
            btn_disabled = True

    elif trigger == "bd-upload-biom" and biom_contents:
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
            biom_status = dbc.Alert(f"Uploaded BIOM: {n} samples", color="success")

            # Try to find matching metadata from pipeline datasets
            match_df, match_sid, match_name = find_metadata_for_samples(sample_ids)
            if match_df is not None:
                meta_json = match_df.to_json(date_format="iso", orient="split")
                sid_col = match_sid
                color_opts, shape_opts, group_opts = _meta_outputs(match_df, match_sid)
                gcols = get_group_columns(match_df, match_sid)
                meta_status = dbc.Alert(
                    f"Metadata auto-matched from \"{match_name}\": "
                    f"{len(match_df)} samples, {len(gcols)} group columns",
                    color="success",
                )
                btn_disabled = len(gcols) == 0
            else:
                btn_disabled = meta_json is None

    elif trigger == "bd-upload-meta" and meta_contents:
        df, s_col, error = parse_uploaded_metadata(meta_contents, meta_filename or "")
        if error:
            meta_json = None
            sid_col = None
            meta_status = dbc.Alert(error, color="danger")
            color_opts, shape_opts, group_opts = [], [{"label": "None", "value": ""}], []
            btn_disabled = True
        else:
            meta_json = df.to_json(date_format="iso", orient="split")
            sid_col = s_col
            color_opts, shape_opts, group_opts = _meta_outputs(df, s_col)
            gcols = get_group_columns(df, s_col)
            meta_status = dbc.Alert(f"Metadata loaded: {len(df)} samples, {len(gcols)} columns", color="success")
            btn_disabled = biom_path is None or len(gcols) == 0

    return (biom_path, biom_status, meta_json, sid_col, meta_status,
            color_opts, shape_opts, group_opts, btn_disabled)


# ── Callback A: on_run — immediate results + start background pairwise ──────

@dash_app.callback(
    Output("bd-loading-target", "children"),
    Output("bd-scatter", "figure"),
    Output("bd-scatter", "style"),
    Output("bd-permanova-table", "children"),
    Output("bd-error", "children"),
    Output("bd-permanova-job-id", "data"),
    Output("bd-poll", "disabled"),
    Output("bd-progress-section", "style", allow_duplicate=True),
    Output("bd-btn-run", "disabled", allow_duplicate=True),
    Output("bd-pairwise-table", "children", allow_duplicate=True),
    Output("bd-progress-bar", "value", allow_duplicate=True),
    Output("bd-progress-bar", "label", allow_duplicate=True),
    Output("bd-progress-log", "children", allow_duplicate=True),
    Input("bd-btn-run", "n_clicks"),
    State("bd-biom-path", "data"),
    State("bd-meta-store", "data"),
    State("bd-sample-id-col", "data"),
    State("bd-distance-metric", "value"),
    State("bd-ordination", "value"),
    State("bd-color-col", "value"),
    State("bd-shape-col", "value"),
    State("bd-group-col", "value"),
    prevent_initial_call=True,
)
def on_run(n_clicks, biom_path, meta_json, sid_col, dist_metric, ord_method,
           color_col, shape_col, group_col):
    # Default outputs for early returns (13 outputs)
    loading_done = ""
    no_fig = no_update
    no_style = no_update
    no_perm = no_update
    no_err = no_update
    no_job = None
    poll_off = True
    progress_hide = HIDDEN
    btn_enabled = False
    pairwise_clear = ""
    bar_zero = 0
    bar_label = ""
    log_clear = ""

    if not all([biom_path, meta_json, dist_metric]):
        return (loading_done, no_fig, no_style, no_perm,
                dbc.Alert("Please fill required inputs.", color="warning"),
                no_job, poll_off, progress_hide, btn_enabled,
                pairwise_clear, bar_zero, bar_label, log_clear)

    try:
        meta_df = pd.read_json(io.StringIO(meta_json), orient="split")
        table = load_table(biom_path)
        biom_ids = list(table.ids(axis="sample"))
        match_info = validate_metadata_vs_biom(meta_df, sid_col, biom_ids)

        if len(match_info["matched"]) < 3:
            return (loading_done, no_fig, no_style, no_perm,
                    dbc.Alert("Need at least 3 matched samples.", color="danger"),
                    no_job, poll_off, progress_hide, btn_enabled,
                    pairwise_clear, bar_zero, bar_label, log_clear)

        count_df = biom_to_count_df(biom_path)
        matched = match_info["matched"]
        count_sub = count_df[[s for s in count_df.columns if s in matched]]

        dm = compute_distance(count_sub, dist_metric)

        if ord_method == "pcoa":
            coords, prop_exp = run_pcoa(dm)
            ax1_label = f"PC1 ({prop_exp.get('PC1', 0):.1%})"
            ax2_label = f"PC2 ({prop_exp.get('PC2', 0):.1%})"
            coords.columns = ["Axis1", "Axis2"]
            title_suffix = "PCoA"
        else:
            coords, stress = run_nmds(dm)
            ax1_label = "NMDS1"
            ax2_label = "NMDS2"
            coords.columns = ["Axis1", "Axis2"]
            title_suffix = f"NMDS (stress={stress:.3f})"

        coords.index = coords.index.astype(str)
        meta_indexed = meta_df.set_index(meta_df[sid_col].astype(str))

        if color_col:
            coords["color"] = coords.index.map(meta_indexed[color_col])
        else:
            coords["color"] = "All"

        if shape_col and shape_col != "":
            coords["shape"] = coords.index.map(meta_indexed[shape_col])
        else:
            coords["shape"] = "All"

        coords["sample_id"] = coords.index

        # ── Build scatter plot ──
        fig = go.Figure()
        color_groups = sorted(coords["color"].unique())
        shape_groups = sorted(coords["shape"].unique())

        for c_idx, cg in enumerate(color_groups):
            for s_idx, sg in enumerate(shape_groups):
                mask = (coords["color"] == cg) & (coords["shape"] == sg)
                sub = coords[mask]
                if sub.empty:
                    continue
                name = str(cg)
                if sg != "All":
                    name += f" / {sg}"
                fig.add_trace(go.Scatter(
                    x=sub["Axis1"], y=sub["Axis2"],
                    mode="markers",
                    name=name,
                    marker=dict(
                        color=COLORS[c_idx % len(COLORS)],
                        symbol=SYMBOLS[s_idx % len(SYMBOLS)],
                        size=10, opacity=0.8,
                    ),
                    text=sub["sample_id"],
                    hovertemplate="<b>%{text}</b><br>%{x:.3f}, %{y:.3f}<extra></extra>",
                ))

        # ── Confidence ellipses ──
        if color_col:
            group_values = coords.index.map(meta_indexed[color_col])
            ellipses = compute_confidence_ellipse(coords, group_values)
            for ell in ellipses:
                grp = ell["group"]
                c_idx = color_groups.index(grp) if grp in color_groups else 0
                fig.add_trace(go.Scatter(
                    x=ell["x"], y=ell["y"],
                    mode="lines",
                    fill="toself",
                    fillcolor=f"rgba({_hex_to_rgb(COLORS[c_idx % len(COLORS)])}, 0.10)",
                    line=dict(color=COLORS[c_idx % len(COLORS)], width=1.5, dash="dot"),
                    name=f"{grp} (95% CI)",
                    showlegend=False,
                    hoverinfo="skip",
                ))

        metric_label = "Bray-Curtis" if dist_metric == "braycurtis" else "Jaccard"
        fig.update_layout(
            title=f"{metric_label} — {title_suffix}",
            xaxis_title=ax1_label,
            yaxis_title=ax2_label,
            template="plotly_dark",
            height=550,
        )

        # ── Global PERMANOVA (synchronous) ──
        perm_table = ""
        job_id = None
        poll_disabled = True
        progress_style = HIDDEN
        btn_disabled = False

        if group_col:
            perm_children = []
            global_res = run_permanova_global(dm, meta_df, sid_col, group_col,
                                              n_permutations=99)
            if global_res:
                perm_children.append(html.H5("Global PERMANOVA", className="mt-3"))
                global_df = pd.DataFrame([{
                    "Groups": global_res["n_groups"],
                    "Samples": global_res["n_samples"],
                    "pseudo_F": f"{global_res['pseudo_F']:.4f}",
                    "p-value": f"{global_res['pvalue']:.4f}",
                    "Significant": "Yes" if global_res["pvalue"] < 0.05 else "No",
                }])
                perm_children.append(dbc.Table.from_dataframe(
                    global_df, striped=True, bordered=True, hover=True,
                    color="dark", size="sm",
                ))

            if perm_children:
                perm_table = html.Div(perm_children)

            # ── Launch background pairwise PERMANOVA ──
            meta_indexed_dm = meta_df.set_index(meta_df[sid_col].astype(str))
            n_groups = len(meta_indexed_dm.loc[
                meta_indexed_dm.index.intersection(dm.ids), group_col
            ].unique())

            if n_groups >= 2:
                job_id = str(uuid.uuid4())
                run_permanova_pairwise_background(
                    dm, meta_df, sid_col, group_col,
                    n_permutations=99, job_id=job_id,
                )
                poll_disabled = False
                progress_style = {"display": "block"}
                btn_disabled = True

        return ("", fig, {"display": "block"}, perm_table, "",
                job_id, poll_disabled, progress_style, btn_disabled,
                "", 0, "", "Starting pairwise PERMANOVA..." if job_id else "")

    except Exception as e:
        return (loading_done, no_fig, no_style, no_perm,
                dbc.Alert(f"Error: {e}\n{traceback.format_exc()}", color="danger"),
                no_job, poll_off, progress_hide, btn_enabled,
                pairwise_clear, bar_zero, bar_label, log_clear)


# ── Callback B: on_poll — progress updates for background pairwise ──────────

@dash_app.callback(
    Output("bd-progress-bar", "value"),
    Output("bd-progress-bar", "label"),
    Output("bd-progress-log", "children"),
    Output("bd-poll", "disabled", allow_duplicate=True),
    Output("bd-progress-section", "style"),
    Output("bd-pairwise-table", "children"),
    Output("bd-btn-run", "disabled", allow_duplicate=True),
    Output("bd-permanova-job-id", "data", allow_duplicate=True),
    Input("bd-poll", "n_intervals"),
    State("bd-permanova-job-id", "data"),
    prevent_initial_call=True,
)
def on_poll(n_intervals, job_id):
    if not job_id:
        return 0, "", "", True, HIDDEN, no_update, False, no_update

    progress = read_permanova_progress(job_id)
    if progress is None:
        return 0, "Starting...", "", False, {"display": "block"}, no_update, True, no_update

    total = progress.get("total", 1) or 1
    completed = progress.get("completed", 0)
    pct = int(100 * completed / total) if total else 0
    label = f"{completed}/{total}"
    log_text = "\n".join(progress.get("log", []))
    status = progress.get("status", "running")

    if status == "complete":
        # Build the pairwise results table
        results = progress.get("results", [])
        cleanup_permanova_progress(job_id)

        if results:
            pw_df = pd.DataFrame(results)
            display_df = pw_df.copy()
            for col in ["pseudo_F", "pvalue", "qvalue"]:
                if col in display_df.columns:
                    display_df[col] = display_df[col].map(lambda x: f"{x:.4f}")
            pairwise_table = html.Div([
                html.H5("Pairwise PERMANOVA (BH-corrected)", className="mt-3"),
                dbc.Table.from_dataframe(
                    display_df, striped=True, bordered=True, hover=True,
                    color="dark", size="sm",
                ),
            ])
        else:
            pairwise_table = html.Small(
                "No valid pairwise comparisons.", className="text-muted mt-2",
            )

        return 100, f"{total}/{total}", log_text, True, HIDDEN, pairwise_table, False, None

    if status == "error":
        error_log = "\n".join(progress.get("log", ["Unknown error"]))
        cleanup_permanova_progress(job_id)
        return (0, "Error", error_log, True, {"display": "block"},
                dbc.Alert("Pairwise PERMANOVA failed. See log above.", color="danger"),
                False, None)

    # Still running
    return pct, label, log_text, False, {"display": "block"}, no_update, True, no_update


# ── Callback: restore pairwise job on page load ─────────────────────────────


@dash_app.callback(
    Output("bd-poll", "disabled", allow_duplicate=True),
    Output("bd-progress-section", "style", allow_duplicate=True),
    Output("bd-progress-bar", "value", allow_duplicate=True),
    Output("bd-progress-bar", "label", allow_duplicate=True),
    Output("bd-progress-log", "children", allow_duplicate=True),
    Output("bd-btn-run", "disabled", allow_duplicate=True),
    Input("bd-permanova-job-id", "data"),
    prevent_initial_call="initial_duplicate",
)
def on_job_restore(job_id):
    """Restore polling UI when returning to a page with an active pairwise job."""
    defaults = (True, HIDDEN, no_update, no_update, no_update, no_update)
    if not job_id:
        return defaults

    prog = read_permanova_progress(job_id)
    if prog is None:
        return defaults

    status = prog.get("status", "running")
    completed = prog.get("completed", 0)
    total = prog.get("total", 1) or 1
    pct = int(completed / total * 100) if total else 0
    label = f"{completed}/{total}"
    log_text = "\n".join(prog.get("log", []))

    if status == "running":
        # Job still running — enable polling, show progress, disable run
        return (False, {"display": "block"}, pct, label, log_text, True)

    # Terminal state — enable one poll cycle to render final results
    return (False, {"display": "block"}, pct, label, log_text, no_update)


# ── Callback C: re-enable Run button when analysis settings change ────────

@dash_app.callback(
    Output("bd-btn-run", "disabled", allow_duplicate=True),
    Input("bd-color-col", "value"),
    Input("bd-shape-col", "value"),
    Input("bd-group-col", "value"),
    Input("bd-distance-metric", "value"),
    Input("bd-ordination", "value"),
    State("bd-biom-path", "data"),
    State("bd-meta-store", "data"),
    prevent_initial_call=True,
)
def on_settings_change(_color, _shape, _group, _dist, _ord, biom_path, meta_json):
    """Re-enable Run button when the user changes any analysis option."""
    if biom_path and meta_json:
        return False
    return True


# ── Helper ──────────────────────────────────────────────────────────────────

def _hex_to_rgb(hex_color: str) -> str:
    """Convert '#3498db' to '52, 152, 219'."""
    h = hex_color.lstrip("#")
    return ", ".join(str(int(h[i:i+2], 16)) for i in (0, 2, 4))
