"""
Core Microbiome Analysis — Identify taxa present in a given percentage of samples.
"""
import io
import traceback

import dash_bootstrap_components as dbc
import pandas as pd
import plotly.graph_objects as go
from dash import Input, Output, State, dcc, html, no_update

from app.analysis.core_microbiome import compute_core_microbiome
from app.analysis.shared import (
    get_dataset_metadata_df,
    get_group_columns,
    get_pipeline_biom_options,
    parse_uploaded_biom,
    parse_uploaded_metadata,
    validate_metadata_vs_biom,
)
from app.dashboard.app import app as dash_app

LEVEL_OPTIONS = [
    {"label": "Phylum", "value": "Phylum"},
    {"label": "Class", "value": "Class"},
    {"label": "Order", "value": "Order"},
    {"label": "Family", "value": "Family"},
    {"label": "Genus", "value": "Genus"},
    {"label": "Species", "value": "Species"},
    {"label": "ASV", "value": "ASV"},
]


def get_layout():
    pipeline_opts = get_pipeline_biom_options()
    return dbc.Container([
        html.H3("Core Microbiome", className="mb-2"),
        html.P(
            "Identify taxa present in a given percentage of samples within each group.",
            className="text-muted mb-4",
        ),
        dcc.Store(id="cm-biom-path"),
        dcc.Store(id="cm-meta-store"),
        dcc.Store(id="cm-sample-id-col"),
        dcc.Store(id="cm-csv-store"),

        dbc.Row([
            # Left column: inputs
            dbc.Col([
                dbc.Label("Pipeline Dataset"),
                dcc.Dropdown(id="cm-pipeline-select", options=pipeline_opts,
                             placeholder="Select a completed dataset..."),
                html.Hr(),
                dbc.Label("Or upload BIOM file"),
                dcc.Upload(id="cm-upload-biom",
                           children=html.Div(["Drag & drop .biom or ", html.A("browse")]),
                           style={"borderWidth": "1px", "borderStyle": "dashed",
                                  "borderRadius": "5px", "textAlign": "center",
                                  "padding": "10px", "cursor": "pointer"},
                           ),
                html.Hr(),
                dbc.Label("Metadata"),
                dcc.Upload(id="cm-upload-meta",
                           children=html.Div(["Drag & drop .csv/.tsv or ", html.A("browse")]),
                           style={"borderWidth": "1px", "borderStyle": "dashed",
                                  "borderRadius": "5px", "textAlign": "center",
                                  "padding": "10px", "cursor": "pointer"},
                           ),
                html.Div(id="cm-meta-status", className="mt-2"),
                html.Hr(),
                dbc.Label("Group Column"),
                dcc.Dropdown(id="cm-group-col", placeholder="Select group..."),
                dbc.Label("Taxonomy Level", className="mt-3"),
                dcc.Dropdown(id="cm-tax-level", options=LEVEL_OPTIONS, value="Genus"),
                dbc.Label("Prevalence Threshold", className="mt-3"),
                dcc.Slider(
                    id="cm-threshold",
                    min=0, max=100, step=5, value=80,
                    marks={i: f"{i}%" for i in range(0, 101, 20)},
                    tooltip={"placement": "bottom", "always_visible": True},
                ),
                dbc.Button("Run", id="cm-btn-run", color="primary",
                           className="mt-3 w-100", disabled=True),
            ], md=4),

            # Right column: results
            dbc.Col([
                dcc.Loading(html.Div(id="cm-results")),
                html.Div([
                    dbc.Button("Download CSV", id="cm-btn-download", color="outline-secondary",
                               size="sm", disabled=True, className="mt-2"),
                    dcc.Download(id="cm-download"),
                ]),
            ], md=8),
        ]),
    ], fluid=True)


# ── Input callbacks ──────────────────────────────────────────────────────────


@dash_app.callback(
    Output("cm-biom-path", "data"),
    Output("cm-meta-store", "data"),
    Output("cm-sample-id-col", "data"),
    Output("cm-meta-status", "children"),
    Output("cm-group-col", "options"),
    Output("cm-btn-run", "disabled"),
    Input("cm-pipeline-select", "value"),
    Input("cm-upload-biom", "contents"),
    Input("cm-upload-meta", "contents"),
    State("cm-upload-biom", "filename"),
    State("cm-upload-meta", "filename"),
    prevent_initial_call=True,
)
def on_input_change(pipeline_val, biom_contents, meta_contents,
                    biom_filename, meta_filename):
    biom_path = None
    meta_json = None
    sid_col = None
    meta_status = ""
    group_opts = []
    can_run = False

    # Pipeline dataset selected
    if pipeline_val:
        from app.db.database import SessionLocal
        from app.db.models import Dataset
        db = SessionLocal()
        try:
            ds = db.query(Dataset).filter(Dataset.id == pipeline_val).first()
            if ds and ds.asv_table_path:
                biom_path = ds.asv_table_path.replace(".tsv", ".biom")
                if not biom_path.endswith(".biom"):
                    biom_path = ds.asv_table_path

                # Try to auto-load metadata
                db_meta, db_sid = get_dataset_metadata_df(pipeline_val)
                if db_meta is not None:
                    meta_json = db_meta.to_json(date_format="iso", orient="split")
                    sid_col = db_sid
                    group_opts = [{"label": c, "value": c} for c in get_group_columns(db_meta, db_sid)]
                    meta_status = dbc.Badge("Metadata loaded from database", color="success")
                    can_run = True
                else:
                    meta_status = dbc.Badge("No metadata found — please upload", color="warning")
        finally:
            db.close()

    # Uploaded BIOM
    if biom_contents and not pipeline_val:
        result = parse_uploaded_biom(biom_contents, biom_filename)
        if result.get("path"):
            biom_path = result["path"]

    # Uploaded metadata
    if meta_contents:
        result = parse_uploaded_metadata(meta_contents, meta_filename)
        if result.get("df") is not None:
            meta_df = result["df"]
            sid_col = result.get("sample_id_column", "sample-id")
            meta_json = meta_df.to_json(date_format="iso", orient="split")
            group_opts = [{"label": c, "value": c} for c in get_group_columns(meta_df, sid_col)]
            meta_status = dbc.Badge(f"Metadata: {len(meta_df)} samples", color="success")
            can_run = biom_path is not None

    if biom_path and meta_json:
        can_run = True

    return biom_path, meta_json, sid_col, meta_status, group_opts, not can_run


# ── Run callback ─────────────────────────────────────────────────────────────


@dash_app.callback(
    Output("cm-results", "children"),
    Output("cm-csv-store", "data"),
    Output("cm-btn-download", "disabled"),
    Input("cm-btn-run", "n_clicks"),
    State("cm-biom-path", "data"),
    State("cm-meta-store", "data"),
    State("cm-sample-id-col", "data"),
    State("cm-group-col", "value"),
    State("cm-tax-level", "value"),
    State("cm-threshold", "value"),
    prevent_initial_call=True,
)
def on_run(n_clicks, biom_path, meta_json, sid_col, group_col, level, threshold_pct):
    if not n_clicks or not biom_path or not meta_json:
        return no_update, no_update, no_update

    try:
        meta_df = pd.read_json(io.StringIO(meta_json), orient="split")
        threshold = (threshold_pct or 80) / 100.0

        result = compute_core_microbiome(
            biom_path=biom_path,
            meta_df=meta_df,
            sample_id_col=sid_col,
            group_col=group_col,
            level=level,
            threshold=threshold,
        )

        children = []

        # Overall core taxa
        overall = result["overall"]
        if overall.empty:
            children.append(
                dbc.Alert(
                    f"No taxa found at ≥{threshold_pct}% prevalence.",
                    color="warning",
                )
            )
        else:
            children.append(html.H5(f"Core Taxa (≥{threshold_pct}% prevalence): {len(overall)} taxa"))

            # Bar chart
            fig = go.Figure()
            all_prev = result["all_prevalences"]

            if group_col and result["per_group"]:
                # Show per-group bars for core taxa
                core_taxa = overall["taxon"].tolist()[:30]  # Top 30
                for group_name in sorted(result["per_group"].keys()):
                    group_data = all_prev[all_prev["taxon"].isin(core_taxa)]
                    if group_name in group_data.columns:
                        fig.add_trace(go.Bar(
                            y=group_data["taxon"],
                            x=group_data[group_name] * 100,
                            name=group_name,
                            orientation="h",
                        ))
            else:
                core_taxa = overall["taxon"].tolist()[:30]
                fig.add_trace(go.Bar(
                    y=core_taxa,
                    x=(overall["prevalence"].head(30) * 100).tolist(),
                    orientation="h",
                    name="All samples",
                ))

            fig.update_layout(
                barmode="group",
                xaxis_title="Prevalence (%)",
                yaxis_title="",
                yaxis=dict(autorange="reversed"),
                height=max(400, len(overall.head(30)) * 25 + 100),
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)",
                font_color="#fff",
                margin=dict(l=200),
            )

            # Add threshold line
            fig.add_vline(x=threshold_pct, line_dash="dash", line_color="red",
                          annotation_text=f"Threshold ({threshold_pct}%)")

            children.append(dcc.Graph(figure=fig))

            # Summary table
            table_df = overall.copy()
            table_df["prevalence"] = (table_df["prevalence"] * 100).round(1).astype(str) + "%"
            children.append(html.H6("Core Taxa Table", className="mt-3"))
            children.append(
                dbc.Table.from_dataframe(
                    table_df, striped=True, bordered=True, hover=True, size="sm",
                )
            )

            # Per-group summaries
            if result["per_group"]:
                children.append(html.H6("Core Taxa per Group", className="mt-3"))
                for group_name, group_df in sorted(result["per_group"].items()):
                    children.append(
                        dbc.Badge(f"{group_name}: {len(group_df)} core taxa", color="info", className="me-2 mb-2")
                    )

        # Prepare CSV for download
        csv_data = overall.to_csv(index=False) if not overall.empty else ""

        return children, csv_data, overall.empty

    except Exception as e:
        return dbc.Alert(f"Error: {e}\n{traceback.format_exc()}", color="danger"), no_update, True


# ── Download callback ────────────────────────────────────────────────────────


@dash_app.callback(
    Output("cm-download", "data"),
    Input("cm-btn-download", "n_clicks"),
    State("cm-csv-store", "data"),
    prevent_initial_call=True,
)
def on_download(n_clicks, csv_data):
    if not n_clicks or not csv_data:
        return no_update
    return dcc.send_string(csv_data, "core_microbiome.csv")
