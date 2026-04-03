"""
MicrobiomeDash — Taxonomy composition page.
"""
import io
import traceback

import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from biom import load_table
from dash import Input, Output, State, callback_context, dcc, html, no_update
from plotly.subplots import make_subplots
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

from app.analysis.shared import (
    find_metadata_for_samples,
    get_dataset_metadata_df,
    get_group_columns,
    get_pipeline_biom_options,
    parse_uploaded_biom,
    parse_uploaded_metadata,
    validate_metadata_vs_biom,
)
from app.analysis.taxonomy import LEVEL_MAP, aggregate_taxonomy, build_heatmap_data
from app.dashboard.app import app as dash_app

LEVEL_OPTIONS = [{"label": k, "value": k} for k in LEVEL_MAP.keys() if k != "ASV"]


def get_layout():
    pipeline_opts = get_pipeline_biom_options()
    return dbc.Container([
        html.H3("Taxonomy Composition", className="mb-2"),
        html.P(
            "Visualize taxonomic composition across samples as a heatmap.",
            className="text-muted mb-4",
        ),
        dcc.Store(id="tx-biom-path"),
        dcc.Store(id="tx-meta-store"),
        dcc.Store(id="tx-sample-id-col"),

        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("Input Data"),
                    dbc.CardBody([
                        dbc.Label("BIOM Table", className="fw-bold"),
                        dbc.Select(
                            id="tx-select-pipeline",
                            options=[{"label": "— none —", "value": ""}]
                            + pipeline_opts,
                            value="",
                            className="mb-2",
                        ),
                        html.Div("— or upload —", className="text-center text-muted small mb-2"),
                        dcc.Upload(
                            id="tx-upload-biom",
                            children=html.Div(["Drag & drop or ", html.A("select .biom")]),
                            style={"borderWidth": "2px", "borderStyle": "dashed",
                                   "borderRadius": "5px", "borderColor": "#555",
                                   "textAlign": "center", "padding": "8px"},
                            multiple=False,
                        ),
                        html.Div(id="tx-biom-status", className="mt-1 small mb-3"),

                        dbc.Label("Metadata (CSV/TSV)", className="fw-bold"),
                        html.Div(
                            "Auto-loaded from pipeline dataset if available. "
                            "Upload to override.",
                            className="text-muted small mb-1",
                        ),
                        dcc.Upload(
                            id="tx-upload-meta",
                            children=html.Div(["Drag & drop or ", html.A("select metadata file")]),
                            style={"borderWidth": "2px", "borderStyle": "dashed",
                                   "borderRadius": "5px", "borderColor": "#555",
                                   "textAlign": "center", "padding": "8px"},
                            multiple=False,
                        ),
                        html.Div(id="tx-meta-status", className="mt-1 small mb-3"),

                        dbc.Label("Taxonomy Level", className="fw-bold"),
                        dbc.Select(
                            id="tx-level",
                            options=LEVEL_OPTIONS,
                            value="Genus",
                            className="mb-2",
                        ),

                        dbc.Label("Top N Taxa", className="fw-bold"),
                        dbc.Input(id="tx-top-n", type="number", value=20, min=5, max=100,
                                  className="mb-2"),

                        dbc.Label("Group Column", className="fw-bold"),
                        dbc.Select(id="tx-group-col", placeholder="Select or upload metadata...",
                                   className="mb-3"),

                        dbc.Switch(
                            id="tx-exclude-unassigned",
                            label="Exclude unassigned",
                            value=True,
                            className="mb-3",
                        ),

                        dbc.Button("Generate Heatmap", id="tx-btn-run", color="primary",
                                   className="w-100", disabled=True),
                    ]),
                ], className="mb-3"),
            ], md=4),

            dbc.Col([
                dbc.Spinner([
                    html.Div(id="tx-error", className="mb-2"),
                    dcc.Graph(id="tx-heatmap", style={"display": "none"},
                             config={"toImageButtonOptions": {"format": "svg", "scale": 2}}),
                    html.Div(id="tx-download-area"),
                ], color="primary"),
            ], md=8),
        ]),
    ], fluid=True)


# ── Callbacks ────────────────────────────────────────────────────────────────


@dash_app.callback(
    Output("tx-biom-path", "data"),
    Output("tx-biom-status", "children"),
    Output("tx-meta-store", "data"),
    Output("tx-sample-id-col", "data"),
    Output("tx-meta-status", "children"),
    Output("tx-group-col", "options"),
    Output("tx-btn-run", "disabled"),
    Input("tx-upload-biom", "contents"),
    Input("tx-select-pipeline", "value"),
    Input("tx-upload-meta", "contents"),
    State("tx-upload-biom", "filename"),
    State("tx-upload-meta", "filename"),
    State("tx-biom-path", "data"),
    State("tx-meta-store", "data"),
    State("tx-sample-id-col", "data"),
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

    if trigger == "tx-select-pipeline" and not pipeline_value:
        biom_path = None
        biom_status = ""
        meta_json = None
        sid_col = None
        meta_status = ""
        group_opts = []
        btn_disabled = True
        return biom_path, biom_status, meta_json, sid_col, meta_status, group_opts, btn_disabled

    if trigger == "tx-select-pipeline" and pipeline_value:
        try:
            table = load_table(pipeline_value)
            n = len(table.ids(axis="sample"))
            obs = table.ids(axis="observation")[0]
            md = table.metadata(obs, axis="observation")
            has_tax = md and "taxonomy" in md
            tax_msg = "" if has_tax else " (no taxonomy metadata found!)"
            color = "success" if has_tax else "warning"
            biom_path = pipeline_value
            biom_status = dbc.Alert(f"Pipeline dataset: {n} samples{tax_msg}", color=color)

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
            else:
                meta_json = None
                sid_col = None
                group_opts = []
                meta_status = dbc.Alert(
                    "No metadata found for this dataset. Upload a metadata file.",
                    color="info",
                )
            # Taxonomy page can run without metadata (no grouping)
            btn_disabled = False
        except Exception as e:
            biom_path = None
            biom_status = dbc.Alert(f"Error: {e}", color="danger")
            btn_disabled = True

    elif trigger == "tx-upload-biom" and biom_contents:
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
            btn_disabled = False

    elif trigger == "tx-upload-meta" and meta_contents:
        df, s_col, error = parse_uploaded_metadata(meta_contents, meta_filename or "")
        if error:
            meta_json = None
            sid_col = None
            meta_status = dbc.Alert(error, color="danger")
            group_opts = []
        else:
            meta_json = df.to_json(date_format="iso", orient="split")
            sid_col = s_col
            gcols = get_group_columns(df, s_col)
            group_opts = [{"label": c, "value": c} for c in gcols]
            meta_status = dbc.Alert(f"Metadata loaded: {len(df)} samples", color="success")
        btn_disabled = biom_path is None

    return biom_path, biom_status, meta_json, sid_col, meta_status, group_opts, btn_disabled


@dash_app.callback(
    Output("tx-heatmap", "figure"),
    Output("tx-heatmap", "style"),
    Output("tx-error", "children"),
    Output("tx-download-area", "children"),
    Input("tx-btn-run", "n_clicks"),
    State("tx-biom-path", "data"),
    State("tx-meta-store", "data"),
    State("tx-sample-id-col", "data"),
    State("tx-level", "value"),
    State("tx-top-n", "value"),
    State("tx-group-col", "value"),
    State("tx-exclude-unassigned", "value"),
    prevent_initial_call=True,
)
def on_run(n_clicks, biom_path, meta_json, sid_col, level, top_n, group_col, exclude_unassigned):
    if not biom_path:
        return no_update, no_update, dbc.Alert("Please select a BIOM table.", color="warning"), no_update

    try:
        top_n = int(top_n) if top_n else 20
        tax_df = aggregate_taxonomy(biom_path, level, top_n)

        if tax_df.empty:
            return no_update, no_update, dbc.Alert(
                "No taxonomy metadata found in this BIOM file. "
                "Taxonomy is embedded during the pipeline's taxonomy assignment step.",
                color="warning",
            ), no_update

        if exclude_unassigned:
            mask = ~tax_df.index.str.lower().str.contains("unassigned")
            tax_df = tax_df.loc[mask]

        if tax_df.empty:
            return no_update, no_update, dbc.Alert(
                "All taxa are unassigned at this level.", color="warning",
            ), no_update

        group_labels = []
        if meta_json and sid_col and group_col:
            meta_df = pd.read_json(io.StringIO(meta_json), orient="split")
            tax_df, group_labels = build_heatmap_data(tax_df, meta_df, sid_col, group_col)

        # ── Hierarchical clustering on samples ──
        n_samples = len(tax_df.columns)
        if n_samples >= 2:
            sample_data = tax_df.values.T  # samples x taxa
            dist = pdist(sample_data, metric="braycurtis")
            dist = np.nan_to_num(dist, nan=0.0)
            Z = linkage(dist, method="average")
            dn = dendrogram(Z, no_plot=True, labels=tax_df.columns.tolist())
            leaf_order = dn["leaves"]

            # Reorder everything by dendrogram leaf order
            ordered_cols = [tax_df.columns[i] for i in leaf_order]
            tax_df = tax_df[ordered_cols]
            if group_labels:
                group_labels = [group_labels[i] for i in leaf_order]

            # Build dendrogram traces — remap x from scipy's 5,15,25,...
            # to 0,1,2,... to align with heatmap integer positions
            dendro_traces = []
            for xs, ys in zip(dn["icoord"], dn["dcoord"]):
                mapped_xs = [(x - 5) / 10 for x in xs]
                dendro_traces.append(go.Scatter(
                    x=mapped_xs, y=ys,
                    mode="lines",
                    line=dict(color="#888", width=1),
                    hoverinfo="skip",
                    showlegend=False,
                ))
        else:
            dendro_traces = []

        # ── Build subplot figure ──
        has_groups = bool(group_labels)
        if has_groups and dendro_traces:
            row_heights = [0.12, 0.03, 0.85]
            n_rows = 3
        elif dendro_traces:
            row_heights = [0.15, 0.85]
            n_rows = 2
        else:
            row_heights = [1.0]
            n_rows = 1

        fig = make_subplots(
            rows=n_rows, cols=1,
            shared_xaxes=True,
            vertical_spacing=0.005,
            row_heights=row_heights,
        )

        current_row = 1

        # Row 1: Dendrogram
        if dendro_traces:
            for trace in dendro_traces:
                fig.add_trace(trace, row=current_row, col=1)
            fig.update_yaxes(
                showticklabels=False, showgrid=False, zeroline=False,
                row=current_row, col=1,
            )
            fig.update_xaxes(
                showticklabels=False, showgrid=False, zeroline=False,
                row=current_row, col=1,
            )
            current_row += 1

        # Row 2: Group color bar
        if has_groups and dendro_traces:
            GROUP_COLORS = [
                "#3498db", "#e74c3c", "#2ecc71", "#f39c12", "#9b59b6",
                "#1abc9c", "#e67e22", "#34495e", "#e84393", "#00cec9",
            ]
            unique_groups = sorted(set(group_labels))
            group_color_map = {
                g: GROUP_COLORS[i % len(GROUP_COLORS)]
                for i, g in enumerate(unique_groups)
            }
            # One-row heatmap encoding groups as integers for coloring
            group_z = [[unique_groups.index(g) for g in group_labels]]
            colorscale = []
            n_groups = len(unique_groups)
            for i, g in enumerate(unique_groups):
                lo = i / max(n_groups, 1)
                hi = (i + 1) / max(n_groups, 1)
                colorscale.append([lo, group_color_map[g]])
                colorscale.append([hi, group_color_map[g]])

            fig.add_trace(go.Heatmap(
                z=group_z,
                x=list(range(len(group_labels))),
                customdata=[[g for g in group_labels]],
                colorscale=colorscale,
                showscale=False,
                hovertemplate="%{customdata}<extra></extra>",
            ), row=current_row, col=1)
            fig.update_yaxes(
                showticklabels=False, showgrid=False, zeroline=False,
                row=current_row, col=1,
            )
            fig.update_xaxes(
                showticklabels=False, showgrid=False, zeroline=False,
                row=current_row, col=1,
            )
            current_row += 1

        # Row 3 (or 2 or 1): Heatmap
        heatmap_x = list(range(n_samples)) if dendro_traces else tax_df.columns.tolist()
        fig.add_trace(go.Heatmap(
            z=tax_df.values,
            x=heatmap_x,
            y=tax_df.index.tolist(),
            colorscale="Viridis",
            hovertemplate="Sample: %{text}<br>Taxon: %{y}<br>Rel. Abundance: %{z:.3f}<extra></extra>",
            text=[tax_df.columns.tolist()] * len(tax_df),
        ), row=current_row, col=1)

        # X-axis labels on the heatmap row
        if dendro_traces:
            fig.update_xaxes(
                tickvals=list(range(n_samples)),
                ticktext=tax_df.columns.tolist(),
                tickangle=45,
                tickfont=dict(size=8 if n_samples > 30 else 10),
                row=current_row, col=1,
            )
        else:
            fig.update_xaxes(
                tickangle=45,
                tickfont=dict(size=8 if n_samples > 30 else 10),
                row=current_row, col=1,
            )

        # ── Group legend annotations ──
        if has_groups and dendro_traces:
            for g in unique_groups:
                fig.add_trace(go.Scatter(
                    x=[None], y=[None],
                    mode="markers",
                    marker=dict(size=10, color=group_color_map[g]),
                    name=g,
                    showlegend=True,
                ))

        heatmap_height = max(400, 50 + len(tax_df) * 25)
        n_unique = len(set(group_labels)) if group_labels else 0
        legend_rows = max(1, (n_unique + 5) // 6) if n_unique else 0
        bot_margin = 120 + legend_rows * 25 if n_unique else 80
        total_height = heatmap_height + (120 if dendro_traces else 0) + 50 + bot_margin
        fig.update_layout(
            title=dict(
                text=f"{level}-level Composition (top {top_n})",
                y=0.99, yanchor="top",
            ),
            template="plotly_dark",
            height=total_height,
            margin=dict(t=50, b=bot_margin),
            legend=dict(
                orientation="h", yanchor="top", y=-0.15, xanchor="left", x=0,
                title=group_col if group_col else None,
            ),
        )

        dl = html.Div([
            dcc.Download(id="tx-download"),
            dbc.Button("Download Table (CSV)", id="tx-btn-download",
                       color="secondary", size="sm", className="mt-2"),
            dcc.Store(id="tx-csv-data", data=tax_df.to_csv()),
        ])

        return fig, {"display": "block"}, "", dl

    except Exception as e:
        return no_update, no_update, dbc.Alert(f"Error: {e}\n{traceback.format_exc()}", color="danger"), no_update


@dash_app.callback(
    Output("tx-download", "data"),
    Input("tx-btn-download", "n_clicks"),
    State("tx-csv-data", "data"),
    prevent_initial_call=True,
)
def on_download(n_clicks, csv_data):
    if csv_data:
        return dcc.send_string(csv_data, "taxonomy_composition.csv")
    return no_update


