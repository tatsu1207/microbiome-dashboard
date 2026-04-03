"""
MicrobiomeDash — Alpha Diversity analysis page.
"""
import io
import traceback

import dash_bootstrap_components as dbc
import pandas as pd
import plotly.graph_objects as go
from biom import load_table
from dash import Input, Output, State, callback_context, dcc, html, no_update
from plotly.subplots import make_subplots

from app.analysis.alpha import compute_alpha, run_alpha_stats
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

METRIC_OPTIONS = [
    {"label": "Shannon", "value": "shannon"},
    {"label": "Simpson", "value": "simpson"},
    {"label": "Observed OTUs", "value": "observed_otus"},
    {"label": "Chao1", "value": "chao1"},
    {"label": "Pielou's Evenness", "value": "pielou_e"},
]


def get_layout():
    pipeline_opts = get_pipeline_biom_options()
    return dbc.Container([
        html.H3("Alpha Diversity", className="mb-2"),
        html.P(
            "Compute within-sample diversity metrics and compare across groups.",
            className="text-muted mb-4",
        ),
        dcc.Store(id="ad-biom-path"),
        dcc.Store(id="ad-meta-store"),
        dcc.Store(id="ad-sample-id-col"),

        dbc.Row([
            # Left column: inputs
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("Input Data"),
                    dbc.CardBody([
                        dbc.Label("BIOM Table", className="fw-bold"),
                        dbc.Select(
                            id="ad-select-pipeline",
                            options=[{"label": "— none —", "value": ""}]
                            + pipeline_opts,
                            value="",
                            className="mb-2",
                        ),
                        html.Div("— or upload —", className="text-center text-muted small mb-2"),
                        dcc.Upload(
                            id="ad-upload-biom",
                            children=html.Div(["Drag & drop or ", html.A("select .biom")]),
                            style={"borderWidth": "2px", "borderStyle": "dashed",
                                   "borderRadius": "5px", "borderColor": "#555",
                                   "textAlign": "center", "padding": "8px"},
                            multiple=False,
                        ),
                        html.Div(id="ad-biom-status", className="mt-1 small mb-3"),

                        dbc.Label("Metadata (CSV/TSV)", className="fw-bold"),
                        html.Div(
                            "Auto-loaded from pipeline dataset if available. "
                            "Upload to override.",
                            className="text-muted small mb-1",
                        ),
                        dcc.Upload(
                            id="ad-upload-meta",
                            children=html.Div(["Drag & drop or ", html.A("select metadata file")]),
                            style={"borderWidth": "2px", "borderStyle": "dashed",
                                   "borderRadius": "5px", "borderColor": "#555",
                                   "textAlign": "center", "padding": "8px"},
                            multiple=False,
                        ),
                        html.Div(id="ad-meta-status", className="mt-1 small mb-3"),

                        dbc.Label("Metrics", className="fw-bold"),
                        dbc.Checklist(
                            id="ad-metrics",
                            options=METRIC_OPTIONS,
                            value=[o["value"] for o in METRIC_OPTIONS],
                            className="mb-3",
                        ),

                        dbc.Label("Group Column", className="fw-bold"),
                        dbc.Select(id="ad-group-col", placeholder="Select or upload metadata...",
                                   className="mb-3"),

                        dbc.Button("Run Analysis", id="ad-btn-run", color="primary",
                                   className="w-100", disabled=True),
                    ]),
                ], className="mb-3"),
            ], md=4),

            # Right column: results
            dbc.Col([
                dbc.Spinner([
                    html.Div(id="ad-error", className="mb-2"),
                    dcc.Graph(id="ad-boxplot", style={"display": "none"},
                             config={"toImageButtonOptions": {"format": "svg", "scale": 2}}),
                    html.Div(id="ad-stats-table"),
                    html.Div(id="ad-download-area"),
                ], color="primary"),
            ], md=8),
        ]),
    ], fluid=True)


# ── Callbacks ────────────────────────────────────────────────────────────────


@dash_app.callback(
    Output("ad-biom-path", "data"),
    Output("ad-biom-status", "children"),
    Output("ad-meta-store", "data"),
    Output("ad-sample-id-col", "data"),
    Output("ad-meta-status", "children"),
    Output("ad-group-col", "options"),
    Output("ad-btn-run", "disabled"),
    Input("ad-upload-biom", "contents"),
    Input("ad-select-pipeline", "value"),
    Input("ad-upload-meta", "contents"),
    State("ad-upload-biom", "filename"),
    State("ad-upload-meta", "filename"),
    State("ad-biom-path", "data"),
    State("ad-meta-store", "data"),
    State("ad-sample-id-col", "data"),
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
    group_opts = no_update
    btn_disabled = no_update

    # ── Handle BIOM input ──
    if trigger == "ad-select-pipeline" and not pipeline_value:
        biom_path = None
        biom_status = ""
        meta_json = None
        sid_col = None
        meta_status = ""
        group_opts = []
        btn_disabled = True
        return biom_path, biom_status, meta_json, sid_col, meta_status, group_opts, btn_disabled

    if trigger == "ad-select-pipeline" and pipeline_value:
        try:
            table = load_table(pipeline_value)
            n = len(table.ids(axis="sample"))
            biom_path = pipeline_value
            biom_status = dbc.Alert(f"Pipeline dataset loaded: {n} samples", color="success")

            # Auto-load metadata from DB
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

    elif trigger == "ad-upload-biom" and biom_contents:
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

    # ── Handle metadata upload (always overrides auto-loaded) ──
    elif trigger == "ad-upload-meta" and meta_contents:
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

            msg = f"Metadata loaded: {len(df)} samples, {len(gcols)} group columns"
            if biom_path:
                try:
                    table = load_table(biom_path)
                    biom_ids = list(table.ids(axis="sample"))
                    match = validate_metadata_vs_biom(df, s_col, biom_ids)
                    msg += f" | {len(match['matched'])} matched with BIOM"
                except Exception:
                    pass
            meta_status = dbc.Alert(msg, color="success")
            btn_disabled = biom_path is None or len(gcols) == 0

    return biom_path, biom_status, meta_json, sid_col, meta_status, group_opts, btn_disabled


@dash_app.callback(
    Output("ad-boxplot", "figure"),
    Output("ad-boxplot", "style"),
    Output("ad-stats-table", "children"),
    Output("ad-error", "children"),
    Output("ad-download-area", "children"),
    Input("ad-btn-run", "n_clicks"),
    State("ad-biom-path", "data"),
    State("ad-meta-store", "data"),
    State("ad-sample-id-col", "data"),
    State("ad-metrics", "value"),
    State("ad-group-col", "value"),
    prevent_initial_call=True,
)
def on_run(n_clicks, biom_path, meta_json, sid_col, metrics, group_col):
    if not all([biom_path, meta_json, metrics, group_col]):
        return no_update, no_update, no_update, dbc.Alert("Please fill all inputs.", color="warning"), no_update

    try:
        meta_df = pd.read_json(io.StringIO(meta_json), orient="split")
        table = load_table(biom_path)
        biom_ids = list(table.ids(axis="sample"))
        match_info = validate_metadata_vs_biom(meta_df, sid_col, biom_ids)

        if not match_info["matched"]:
            return no_update, no_update, no_update, dbc.Alert("No matching sample IDs between BIOM and metadata.", color="danger"), no_update

        count_df = biom_to_count_df(biom_path)
        matched = match_info["matched"]
        count_sub = count_df[[s for s in count_df.columns if s in matched]]

        # Compute alpha diversity
        diversity_df = compute_alpha(count_sub, metrics)

        # Build boxplots
        n_metrics = len(metrics)
        fig = make_subplots(
            rows=n_metrics, cols=1,
            subplot_titles=[m.replace("_", " ").title() for m in metrics],
        )

        # Merge with groups
        merged = diversity_df.copy()
        merged.index = merged.index.astype(str)
        meta_map = meta_df.set_index(meta_df[sid_col].astype(str))[group_col]
        merged["group"] = merged.index.map(meta_map)
        merged = merged.dropna(subset=["group"])

        groups = sorted(merged["group"].unique())
        colors = ["#3498db", "#e74c3c", "#2ecc71", "#f39c12", "#9b59b6",
                  "#1abc9c", "#e67e22", "#34495e"]

        for row_idx, metric in enumerate(metrics):
            for g_idx, group in enumerate(groups):
                vals = merged[merged["group"] == group][metric]
                fig.add_trace(
                    go.Box(
                        y=vals, name=group,
                        marker_color=colors[g_idx % len(colors)],
                        showlegend=(row_idx == 0),
                        legendgroup=group,
                    ),
                    row=row_idx + 1, col=1,
                )

        fig.update_layout(
            template="plotly_dark",
            height=350 * n_metrics,
            boxmode="group",
            legend=dict(orientation="h", yanchor="bottom", y=1.02),
        )

        # Run stats for each metric
        stats_rows = []
        pairwise_rows = []
        for metric in metrics:
            stats = run_alpha_stats(diversity_df, meta_df, sid_col, group_col, metric)
            label = metric.replace("_", " ").title()
            stats_rows.append({
                "Metric": label,
                "Kruskal-Wallis H": f"{stats['kruskal_H']:.3f}" if stats["kruskal_H"] is not None else "N/A",
                "p-value": f"{stats['kruskal_p']:.4f}" if stats["kruskal_p"] is not None else "N/A",
                "Significant": "Yes" if stats["kruskal_p"] and stats["kruskal_p"] < 0.05 else "No",
            })
            for pw in stats.get("pairwise", []):
                pairwise_rows.append({
                    "Metric": label,
                    "Group 1": pw["group1"],
                    "Group 2": pw["group2"],
                    "U": f"{pw['U']:.1f}",
                    "p-value": f"{pw['pvalue']:.4f}",
                    "q-value": f"{pw.get('qvalue', 'N/A'):.4f}" if isinstance(pw.get("qvalue"), float) else "N/A",
                })

        # Build stats tables
        tables = []
        if stats_rows:
            tables.append(html.H5("Kruskal-Wallis Test", className="mt-3"))
            tables.append(dbc.Table.from_dataframe(
                pd.DataFrame(stats_rows), striped=True, bordered=True, hover=True,
                color="dark", size="sm",
            ))
        if pairwise_rows:
            tables.append(html.H5("Pairwise Mann-Whitney U (BH-corrected)", className="mt-3"))
            tables.append(dbc.Table.from_dataframe(
                pd.DataFrame(pairwise_rows), striped=True, bordered=True, hover=True,
                color="dark", size="sm",
            ))

        dl = html.Div([
            dcc.Download(id="ad-download"),
            dbc.Button("Download Diversity Table (CSV)", id="ad-btn-download",
                       color="secondary", size="sm", className="mt-2"),
            dcc.Store(id="ad-diversity-csv", data=diversity_df.to_csv()),
        ])

        return fig, {"display": "block"}, html.Div(tables), "", dl

    except Exception as e:
        return no_update, no_update, no_update, dbc.Alert(f"Error: {e}\n{traceback.format_exc()}", color="danger"), no_update


@dash_app.callback(
    Output("ad-btn-run", "disabled", allow_duplicate=True),
    Input("ad-metrics", "value"),
    Input("ad-group-col", "value"),
    State("ad-biom-path", "data"),
    State("ad-meta-store", "data"),
    prevent_initial_call=True,
)
def on_settings_change(_metrics, _group, biom_path, meta_json):
    """Re-enable Run button when the user changes analysis options."""
    if biom_path and meta_json:
        return False
    return True


@dash_app.callback(
    Output("ad-download", "data"),
    Input("ad-btn-download", "n_clicks"),
    State("ad-diversity-csv", "data"),
    prevent_initial_call=True,
)
def on_download(n_clicks, csv_data):
    if csv_data:
        return dcc.send_string(csv_data, "alpha_diversity.csv")
    return no_update


