"""
MicrobiomeDash — Rarefy & Filter page: filter samples and optionally rarefy.

Allows users to select/deselect individual samples from a BIOM file,
apply a read-depth threshold, and download a new BIOM with only the
selected samples (optionally rarefied to a uniform depth).
"""
import dash_bootstrap_components as dbc
from dash import ALL, Input, Output, State, dcc, html, no_update

from app.analysis.shared import (
    get_dataset_metadata_df,
    get_pipeline_biom_options,
    parse_uploaded_biom,
)
from app.dashboard.app import app as dash_app


def get_layout():
    pipeline_opts = get_pipeline_biom_options()

    return dbc.Container(
        [
            html.H3("Rarefy & Filter", className="mb-2"),
            html.P(
                "Filter samples from a BIOM file and optionally rarefy to a "
                "uniform read depth before downstream analysis.",
                className="text-muted mb-4",
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Label("Pipeline Dataset", className="fw-bold"),
                            dbc.Select(
                                id="ss-pipeline-select",
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
                                id="ss-biom-upload",
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
            html.Div(id="ss-source-status", className="mb-2"),
            dbc.Button(
                "Load BIOM",
                id="ss-btn-load",
                color="primary",
                disabled=True,
                className="mb-4",
            ),
            html.Div(id="ss-main-area"),
            dcc.Download(id="ss-download-biom"),
            dcc.Store(id="ss-biom-path"),
            dcc.Store(id="ss-sample-stats"),
            dcc.Store(id="ss-meta-cols"),
        ],
        fluid=True,
    )


# ── Callback 1: resolve BIOM source ─────────────────────────────────────────


@dash_app.callback(
    Output("ss-biom-path", "data"),
    Output("ss-source-status", "children"),
    Output("ss-btn-load", "disabled"),
    Output("ss-pipeline-select", "value"),
    Input("ss-pipeline-select", "value"),
    Input("ss-biom-upload", "contents"),
    State("ss-biom-upload", "filename"),
    prevent_initial_call=True,
)
def ss_on_biom_input(pipeline_val, upload_contents, upload_filename):
    from dash import ctx

    triggered = ctx.triggered_id

    if triggered == "ss-biom-upload" and upload_contents:
        tmp_path, err = parse_uploaded_biom(upload_contents, upload_filename)
        if err:
            return no_update, dbc.Alert(err, color="danger"), True, ""
        return (
            tmp_path,
            html.Span(f"Uploaded: {upload_filename}", className="text-success"),
            False, "",
        )

    if triggered == "ss-pipeline-select" and pipeline_val:
        return (
            pipeline_val,
            html.Span("Pipeline dataset selected", className="text-success"),
            False, no_update,
        )

    return None, "", True, no_update


# ── Callback 2: load BIOM → build entire controls area ──────────────────────


@dash_app.callback(
    Output("ss-sample-stats", "data"),
    Output("ss-meta-cols", "data"),
    Output("ss-main-area", "children"),
    Input("ss-btn-load", "n_clicks"),
    State("ss-biom-path", "data"),
    prevent_initial_call=True,
)
def ss_on_load(n_clicks, biom_path):
    if not n_clicks or not biom_path:
        return no_update, no_update, no_update

    try:
        from biom import load_table

        table = load_table(biom_path)
        sample_ids = list(table.ids(axis="sample"))
        mat = table.matrix_data.toarray()

        stats = []
        for j, sid in enumerate(sample_ids):
            col = mat[:, j]
            stats.append({
                "Sample Name": str(sid),
                "Reads": int(col.sum()),
                "ASVs": int((col > 0).sum()),
            })

        meta_df, sample_id_col = get_dataset_metadata_df(biom_path)
        meta_cols = []
        if meta_df is not None and sample_id_col:
            meta_cols = [c for c in meta_df.columns if c != sample_id_col]
            meta_lookup = {}
            for _, row in meta_df.iterrows():
                meta_lookup[str(row[sample_id_col])] = {
                    c: str(row[c]) for c in meta_cols
                }
            for s in stats:
                meta = meta_lookup.get(s["Sample Name"], {})
                for c in meta_cols:
                    s[c] = meta.get(c, "")

        min_reads = min(s["Reads"] for s in stats) if stats else ""
        all_ids = [s["Sample Name"] for s in stats]

        tbl = _build_sample_table(stats, meta_cols, all_ids)

        controls = html.Div([
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Label("Rarefaction Depth"),
                            html.Small(
                                "Set to 0 for no rarefaction (filter only)",
                                className="text-muted d-block mb-1",
                            ),
                            dbc.Input(
                                id="ss-rarefy-depth",
                                type="number",
                                min=0, step=1,
                                value=min_reads,
                                placeholder="e.g. 10000",
                            ),
                        ],
                        md=3,
                    ),
                    dbc.Col(
                        dbc.Button("Apply Threshold", id="ss-btn-threshold",
                                   color="secondary", className="mt-4"),
                        md="auto",
                    ),
                    dbc.Col(
                        dbc.Button("Select All", id="ss-btn-select-all",
                                   color="outline-primary", size="sm",
                                   className="mt-4 me-1"),
                        md="auto",
                    ),
                    dbc.Col(
                        dbc.Button("Deselect All", id="ss-btn-deselect-all",
                                   color="outline-secondary", size="sm",
                                   className="mt-4"),
                        md="auto",
                    ),
                ],
                className="mb-3 align-items-end",
            ),
            html.Div(tbl, id="ss-table-wrapper"),
            html.Div(
                html.Span(
                    f"{len(stats)} of {len(stats)} samples selected",
                    className="text-info",
                ),
                id="ss-selection-info",
                className="mt-2 mb-3",
            ),
            dbc.Button(
                "Subsample & Download", id="ss-btn-download",
                color="success", size="lg", className="w-100 mb-3",
            ),
            html.Div(id="ss-result"),
        ])

        return stats, meta_cols, controls

    except Exception as e:
        return (
            no_update, no_update,
            dbc.Alert(f"Error loading BIOM: {e}", color="danger"),
        )


# ── Callback 3: threshold / select-all / deselect-all → rebuild table ───────


@dash_app.callback(
    Output("ss-table-wrapper", "children"),
    Output("ss-selection-info", "children"),
    Input("ss-btn-threshold", "n_clicks"),
    Input("ss-btn-select-all", "n_clicks"),
    Input("ss-btn-deselect-all", "n_clicks"),
    State("ss-rarefy-depth", "value"),
    State("ss-sample-stats", "data"),
    State("ss-meta-cols", "data"),
    prevent_initial_call=True,
)
def ss_on_bulk_selection(n_thresh, n_all, n_none, depth, stats, meta_cols):
    from dash import ctx

    if not stats:
        return no_update, no_update

    triggered = ctx.triggered_id

    if triggered == "ss-btn-threshold":
        if depth is None:
            return no_update, no_update
        d = int(depth)
        checked = [s["Sample Name"] for s in stats if s["Reads"] >= d]
    elif triggered == "ss-btn-select-all":
        checked = [s["Sample Name"] for s in stats]
    elif triggered == "ss-btn-deselect-all":
        checked = []
    else:
        return no_update, no_update

    tbl = _build_sample_table(stats, meta_cols, checked)
    n = len(checked)
    if n == 0:
        summary = html.Span("No samples selected", className="text-warning")
    else:
        summary = html.Span(
            f"{n} of {len(stats)} samples selected", className="text-info",
        )
    return tbl, summary


# ── Callback 4: update summary from checkbox clicks ─────────────────────────


@dash_app.callback(
    Output("ss-selection-info", "children", allow_duplicate=True),
    Input({"type": "ss-chk", "index": ALL}, "value"),
    State("ss-sample-stats", "data"),
    prevent_initial_call=True,
)
def ss_on_checkbox_change(values, stats):
    if not stats or not values:
        return no_update
    n_checked = sum(1 for v in values if v)
    n_total = len(stats)
    if n_checked == 0:
        return html.Span("No samples selected", className="text-warning")
    return html.Span(
        f"{n_checked} of {n_total} samples selected", className="text-info",
    )


# ── Callback 5: subsample & download ────────────────────────────────────────


@dash_app.callback(
    Output("ss-result", "children"),
    Output("ss-download-biom", "data"),
    Input("ss-btn-download", "n_clicks"),
    State({"type": "ss-chk", "index": ALL}, "value"),
    State("ss-sample-stats", "data"),
    State("ss-biom-path", "data"),
    State("ss-rarefy-depth", "value"),
    prevent_initial_call=True,
)
def ss_on_subsample(n_clicks, checkbox_values, stats, biom_path, rarefy_depth):
    if not n_clicks or not stats or not biom_path:
        return no_update, no_update

    checked_ids = [
        s["Sample Name"]
        for s, v in zip(stats, checkbox_values or [])
        if v
    ]
    if not checked_ids:
        return dbc.Alert("No samples selected.", color="warning"), no_update

    try:
        if rarefy_depth and int(rarefy_depth) > 0:
            from app.data_manager.subsample import rarefy_samples

            depth = int(rarefy_depth)
            biom_bytes = rarefy_samples(biom_path, checked_ids, depth)
            filename = f"subsampled_rarefied_{depth}.biom"
            msg = f"Rarefied {len(checked_ids)} samples to {depth:,} reads."
        else:
            from app.data_manager.subsample import filter_samples

            biom_bytes = filter_samples(biom_path, checked_ids)
            filename = f"subsampled_{len(checked_ids)}samples.biom"
            msg = f"Filtered to {len(checked_ids)} samples."

        return (
            dbc.Alert(msg, color="success"),
            dcc.send_bytes(biom_bytes, filename=filename),
        )
    except Exception as e:
        return dbc.Alert(f"Subsample failed: {e}", color="danger"), no_update


# ── Helper ───────────────────────────────────────────────────────────────────


def _build_sample_table(stats, meta_cols, checked_ids):
    """Build a dbc.Table with checkboxes."""
    meta_cols = meta_cols or []
    checked_set = set(checked_ids) if checked_ids else set()

    header = html.Thead(
        html.Tr(
            [html.Th("", style={"width": "40px"}), html.Th("Sample Name"),
             html.Th("Reads"), html.Th("ASVs")]
            + [html.Th(c) for c in meta_cols]
        ),
        className="table-dark",
    )

    rows = []
    for i, s in enumerate(stats):
        sid = s["Sample Name"]
        rows.append(html.Tr([
            html.Td(dbc.Checkbox(
                id={"type": "ss-chk", "index": i},
                value=sid in checked_set,
            )),
            html.Td(sid),
            html.Td(f"{s['Reads']:,}"),
            html.Td(f"{s['ASVs']:,}"),
        ] + [html.Td(s.get(c, "")) for c in meta_cols]))

    return dbc.Table(
        [header, html.Tbody(rows)],
        bordered=True, hover=True, size="sm", color="dark",
        style={"maxHeight": "500px", "overflowY": "auto", "display": "block"},
    )
