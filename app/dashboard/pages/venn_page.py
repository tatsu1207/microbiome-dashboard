"""
Venn Diagram — Show overlap of taxa between groups.
"""
import io
import traceback

import dash_bootstrap_components as dbc
import pandas as pd
from dash import Input, Output, State, dcc, html, no_update

from app.analysis.shared import (
    get_dataset_metadata_df,
    get_group_columns,
    get_pipeline_biom_options,
    parse_uploaded_biom,
    parse_uploaded_metadata,
)
from app.analysis.venn import compute_intersection_table, compute_venn_sets, render_venn_figure
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
        html.H3("Venn Diagram", className="mb-2"),
        html.P(
            "Compare shared and unique taxa between 2–4 groups.",
            className="text-muted mb-4",
        ),
        dcc.Store(id="vn-biom-path"),
        dcc.Store(id="vn-meta-store"),
        dcc.Store(id="vn-sample-id-col"),
        dcc.Store(id="vn-csv-store"),

        dbc.Row([
            # Left column: inputs
            dbc.Col([
                dbc.Label("Pipeline Dataset"),
                dcc.Dropdown(id="vn-pipeline-select", options=pipeline_opts,
                             placeholder="Select a completed dataset..."),
                html.Hr(),
                dbc.Label("Or upload BIOM file"),
                dcc.Upload(id="vn-upload-biom",
                           children=html.Div(["Drag & drop .biom or ", html.A("browse")]),
                           style={"borderWidth": "1px", "borderStyle": "dashed",
                                  "borderRadius": "5px", "textAlign": "center",
                                  "padding": "10px", "cursor": "pointer"}),
                html.Hr(),
                dbc.Label("Metadata"),
                dcc.Upload(id="vn-upload-meta",
                           children=html.Div(["Drag & drop .csv/.tsv or ", html.A("browse")]),
                           style={"borderWidth": "1px", "borderStyle": "dashed",
                                  "borderRadius": "5px", "textAlign": "center",
                                  "padding": "10px", "cursor": "pointer"}),
                html.Div(id="vn-meta-status", className="mt-2"),
                html.Hr(),
                dbc.Label("Group Column"),
                dcc.Dropdown(id="vn-group-col", placeholder="Select group..."),
                dbc.Label("Select Groups (2–4)", className="mt-3"),
                dcc.Dropdown(id="vn-group-select", multi=True,
                             placeholder="Select 2–4 groups to compare..."),
                dbc.Label("Taxonomy Level", className="mt-3"),
                dcc.Dropdown(id="vn-tax-level", options=LEVEL_OPTIONS, value="Genus"),
                dbc.Button("Run", id="vn-btn-run", color="primary",
                           className="mt-3 w-100", disabled=True),
            ], md=4),

            # Right column: results
            dbc.Col([
                dcc.Loading(html.Div(id="vn-results")),
                html.Div([
                    dbc.Button("Download CSV", id="vn-btn-download", color="outline-secondary",
                               size="sm", disabled=True, className="mt-2"),
                    dcc.Download(id="vn-download"),
                ]),
            ], md=8),
        ]),
    ], fluid=True)


# ── Input callbacks ──────────────────────────────────────────────────────────


@dash_app.callback(
    Output("vn-biom-path", "data"),
    Output("vn-meta-store", "data"),
    Output("vn-sample-id-col", "data"),
    Output("vn-meta-status", "children"),
    Output("vn-group-col", "options"),
    Input("vn-pipeline-select", "value"),
    Input("vn-upload-biom", "contents"),
    Input("vn-upload-meta", "contents"),
    State("vn-upload-biom", "filename"),
    State("vn-upload-meta", "filename"),
    prevent_initial_call=True,
)
def on_input_change(pipeline_val, biom_contents, meta_contents,
                    biom_filename, meta_filename):
    biom_path = None
    meta_json = None
    sid_col = None
    meta_status = ""
    group_opts = []

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
                db_meta, db_sid = get_dataset_metadata_df(pipeline_val)
                if db_meta is not None:
                    meta_json = db_meta.to_json(date_format="iso", orient="split")
                    sid_col = db_sid
                    group_opts = [{"label": c, "value": c} for c in get_group_columns(db_meta, db_sid)]
                    meta_status = dbc.Badge("Metadata loaded from database", color="success")
        finally:
            db.close()

    if biom_contents and not pipeline_val:
        result = parse_uploaded_biom(biom_contents, biom_filename)
        if result.get("path"):
            biom_path = result["path"]

    if meta_contents:
        result = parse_uploaded_metadata(meta_contents, meta_filename)
        if result.get("df") is not None:
            meta_df = result["df"]
            sid_col = result.get("sample_id_column", "sample-id")
            meta_json = meta_df.to_json(date_format="iso", orient="split")
            group_opts = [{"label": c, "value": c} for c in get_group_columns(meta_df, sid_col)]
            meta_status = dbc.Badge(f"Metadata: {len(meta_df)} samples", color="success")

    return biom_path, meta_json, sid_col, meta_status, group_opts


@dash_app.callback(
    Output("vn-group-select", "options"),
    Input("vn-group-col", "value"),
    State("vn-meta-store", "data"),
    State("vn-sample-id-col", "data"),
    prevent_initial_call=True,
)
def on_group_col_change(group_col, meta_json, sid_col):
    if not group_col or not meta_json:
        return []
    meta_df = pd.read_json(io.StringIO(meta_json), orient="split")
    if group_col not in meta_df.columns:
        return []
    groups = sorted(meta_df[group_col].dropna().unique().tolist())
    return [{"label": str(g), "value": str(g)} for g in groups]


@dash_app.callback(
    Output("vn-btn-run", "disabled"),
    Input("vn-biom-path", "data"),
    Input("vn-meta-store", "data"),
    Input("vn-group-select", "value"),
)
def enable_run(biom_path, meta_json, selected_groups):
    if not biom_path or not meta_json or not selected_groups:
        return True
    return not (2 <= len(selected_groups) <= 4)


# ── Run callback ─────────────────────────────────────────────────────────────


@dash_app.callback(
    Output("vn-results", "children"),
    Output("vn-csv-store", "data"),
    Output("vn-btn-download", "disabled"),
    Input("vn-btn-run", "n_clicks"),
    State("vn-biom-path", "data"),
    State("vn-meta-store", "data"),
    State("vn-sample-id-col", "data"),
    State("vn-group-col", "value"),
    State("vn-group-select", "value"),
    State("vn-tax-level", "value"),
    prevent_initial_call=True,
)
def on_run(n_clicks, biom_path, meta_json, sid_col, group_col, groups, level):
    if not n_clicks or not biom_path or not meta_json or not groups:
        return no_update, no_update, no_update

    try:
        meta_df = pd.read_json(io.StringIO(meta_json), orient="split")

        sets = compute_venn_sets(biom_path, meta_df, sid_col, group_col, groups, level)

        children = []

        # Summary badges
        for g, s in sets.items():
            children.append(
                dbc.Badge(f"{g}: {len(s)} taxa", color="info", className="me-2 mb-2")
            )

        # Venn diagram image
        img_src = render_venn_figure(sets)
        if img_src:
            children.append(html.Img(src=img_src, style={"maxWidth": "100%"}, className="my-3"))

        # Intersection table
        table_df = compute_intersection_table(sets)
        if not table_df.empty:
            children.append(html.H6("Intersection Details", className="mt-3"))
            display_df = table_df[["region", "count"]].copy()
            children.append(
                dbc.Table.from_dataframe(
                    display_df, striped=True, bordered=True, hover=True, size="sm",
                )
            )

        csv_data = table_df.to_csv(index=False) if not table_df.empty else ""
        return children, csv_data, table_df.empty

    except Exception as e:
        return dbc.Alert(f"Error: {e}\n{traceback.format_exc()}", color="danger"), no_update, True


# ── Download callback ────────────────────────────────────────────────────────


@dash_app.callback(
    Output("vn-download", "data"),
    Input("vn-btn-download", "n_clicks"),
    State("vn-csv-store", "data"),
    prevent_initial_call=True,
)
def on_download(n_clicks, csv_data):
    if not n_clicks or not csv_data:
        return no_update
    return dcc.send_string(csv_data, "venn_intersections.csv")
