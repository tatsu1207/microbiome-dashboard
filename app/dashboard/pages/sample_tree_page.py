"""
MicrobiomeDash — Outlier Detection page.

Upload a BIOM file + metadata → view a Bray-Curtis UPGMA dendrogram and a
PCoA scatter plot with 95 % confidence ellipses coloured by a metadata group.
Toggle individual samples via checkboxes; both plots update instantly.
Download a filtered BIOM file with only the selected (included) samples.

State is persisted in session storage so plots survive page navigation.
"""
import io
import tempfile

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import dash_bootstrap_components as dbc
from biom import load_table
from biom.util import biom_open
from dash import Input, Output, State, callback_context, dcc, html, no_update
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

from app.analysis.beta import compute_confidence_ellipse, compute_distance, run_pcoa
from app.analysis.shared import (
    biom_to_count_df,
    find_metadata_for_samples,
    get_dataset_metadata_df,
    get_group_columns,
    get_pipeline_biom_options,
    parse_uploaded_biom,
    parse_uploaded_metadata,
)
from app.dashboard.app import app as dash_app

COLORS = [
    "#3498db", "#e74c3c", "#2ecc71", "#f39c12", "#9b59b6",
    "#1abc9c", "#e67e22", "#34495e", "#e84393", "#00cec9",
]

HIDDEN = {"display": "none"}
SESSION = "session"


def get_layout():
    pipeline_opts = get_pipeline_biom_options()

    return dbc.Container(
        [
            html.H3("Outlier Detection", className="mb-2"),
            html.P(
                "Identify outlier samples by inspecting a Bray-Curtis UPGMA "
                "dendrogram and PCoA plot.  Uncheck samples to exclude them — "
                "both plots update instantly.  Download a cleaned BIOM file "
                "containing only the included samples.",
                className="text-muted mb-4",
            ),

            # ── Input section ───────────────────────────────────────────
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Label("Pipeline Dataset", className="fw-bold"),
                            dbc.Select(
                                id="st-pipeline-select",
                                options=[{"label": "— none —", "value": ""}]
                                + pipeline_opts,
                                value="",
                            ),
                        ],
                        md=4,
                    ),
                    dbc.Col(
                        [
                            dbc.Label("Or Upload BIOM", className="fw-bold"),
                            dcc.Upload(
                                id="st-biom-upload",
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
                        md=4,
                    ),
                    dbc.Col(
                        [
                            dbc.Label("Metadata (CSV/TSV)", className="fw-bold"),
                            dcc.Upload(
                                id="st-meta-upload",
                                children=html.Div(
                                    [
                                        "Drag & drop or ",
                                        html.A("select file", className="text-info"),
                                    ],
                                    className="text-center py-2",
                                ),
                                style={
                                    "borderWidth": "2px",
                                    "borderStyle": "dashed",
                                    "borderRadius": "5px",
                                    "borderColor": "#555",
                                },
                            ),
                        ],
                        md=4,
                    ),
                ],
                className="mb-2",
            ),
            html.Div(id="st-source-status", className="mb-2"),

            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Label("Color / Group by", className="fw-bold"),
                            dbc.Select(
                                id="st-group-col",
                                placeholder="Select column...",
                                disabled=True,
                            ),
                        ],
                        md=4,
                    ),
                    dbc.Col(
                        [
                            dbc.Button(
                                "Load & Plot",
                                id="st-btn-load",
                                color="primary",
                                disabled=True,
                                className="mt-4 me-2",
                            ),
                            dbc.Button(
                                "Clear",
                                id="st-btn-clear",
                                color="outline-danger",
                                className="mt-4",
                            ),
                        ],
                        md=4,
                    ),
                ],
                className="mb-4",
            ),

            # ── Main area (populated after load or restored from session) ─
            html.Div(id="st-main-area"),

            # ── Session-persisted stores ─────────────────────────────────
            dcc.Store(id="st-biom-path", storage_type=SESSION),
            dcc.Store(id="st-meta-store", storage_type=SESSION),
            dcc.Store(id="st-sample-id-col", storage_type=SESSION),
            dcc.Store(id="st-group-col-val", storage_type=SESSION),
            dcc.Store(id="st-all-samples", storage_type=SESSION),
            dcc.Store(id="st-selected-samples", storage_type=SESSION),
            dcc.Download(id="st-download-biom"),
            # One-shot interval to trigger auto-restore on page load
            dcc.Interval(id="st-restore-tick", interval=300, max_intervals=1),
        ],
        fluid=True,
    )


# ── Callback 1: resolve BIOM + metadata sources ─────────────────────────────


@dash_app.callback(
    Output("st-biom-path", "data"),
    Output("st-meta-store", "data"),
    Output("st-sample-id-col", "data"),
    Output("st-source-status", "children"),
    Output("st-group-col", "options"),
    Output("st-group-col", "disabled"),
    Output("st-btn-load", "disabled"),
    Output("st-pipeline-select", "value"),
    Input("st-pipeline-select", "value"),
    Input("st-biom-upload", "contents"),
    Input("st-meta-upload", "contents"),
    State("st-biom-upload", "filename"),
    State("st-meta-upload", "filename"),
    State("st-biom-path", "data"),
    State("st-meta-store", "data"),
    State("st-sample-id-col", "data"),
    prevent_initial_call=True,
)
def st_on_input(pipeline_val, biom_contents, meta_contents,
                biom_filename, meta_filename,
                prev_biom, prev_meta, prev_sid):
    trigger = callback_context.triggered_id

    biom_path = prev_biom
    meta_json = prev_meta
    sid_col = prev_sid
    status_parts = []
    group_opts = no_update
    group_disabled = no_update
    btn_disabled = True
    pipeline_reset = no_update

    # ── Pipeline select ──────────────────────────────────────────
    if trigger == "st-pipeline-select":
        if not pipeline_val:
            return (None, None, None, "", [], True, True, no_update)
        try:
            table = load_table(pipeline_val)
            n = len(table.ids(axis="sample"))
            biom_path = pipeline_val
            status_parts.append(f"Pipeline dataset: {n} samples")

            db_meta, db_sid = get_dataset_metadata_df(pipeline_val)
            if db_meta is not None:
                meta_json = db_meta.to_json(date_format="iso", orient="split")
                sid_col = db_sid
                gcols = get_group_columns(db_meta, db_sid)
                group_opts = [{"label": c, "value": c} for c in gcols]
                group_disabled = len(gcols) == 0
                status_parts.append(f"Metadata auto-loaded ({len(gcols)} columns)")
                btn_disabled = len(gcols) == 0
            else:
                meta_json = None
                sid_col = None
                group_opts = []
                group_disabled = True
                status_parts.append("No metadata found — upload one")
        except Exception as e:
            return (None, None, None, dbc.Alert(str(e), color="danger"),
                    [], True, True, no_update)

    # ── BIOM upload ──────────────────────────────────────────────
    elif trigger == "st-biom-upload" and biom_contents:
        path, err = parse_uploaded_biom(biom_contents, biom_filename or "")
        if err:
            return (no_update, no_update, no_update,
                    dbc.Alert(err, color="danger"),
                    no_update, no_update, True, "")
        table = load_table(path)
        sample_ids = list(table.ids(axis="sample"))
        biom_path = path
        pipeline_reset = ""
        status_parts.append(f"Uploaded BIOM: {len(sample_ids)} samples")

        # Auto-match metadata from pipeline DB
        match_df, match_sid, match_name = find_metadata_for_samples(sample_ids)
        if match_df is not None:
            meta_json = match_df.to_json(date_format="iso", orient="split")
            sid_col = match_sid
            gcols = get_group_columns(match_df, match_sid)
            group_opts = [{"label": c, "value": c} for c in gcols]
            group_disabled = len(gcols) == 0
            status_parts.append(
                f'Metadata auto-matched from "{match_name}" ({len(gcols)} columns)'
            )
            btn_disabled = len(gcols) == 0
        elif meta_json:
            btn_disabled = False
        else:
            group_opts = []
            group_disabled = True
            status_parts.append("Upload metadata for group coloring")

    # ── Metadata upload ──────────────────────────────────────────
    elif trigger == "st-meta-upload" and meta_contents:
        df, s_col, err = parse_uploaded_metadata(meta_contents, meta_filename or "")
        if err:
            return (no_update, no_update, no_update,
                    dbc.Alert(err, color="danger"),
                    no_update, no_update, no_update, no_update)
        meta_json = df.to_json(date_format="iso", orient="split")
        sid_col = s_col
        gcols = get_group_columns(df, s_col)
        group_opts = [{"label": c, "value": c} for c in gcols]
        group_disabled = len(gcols) == 0
        status_parts.append(f"Metadata loaded: {len(df)} samples, {len(gcols)} columns")
        btn_disabled = biom_path is None or len(gcols) == 0

    status_msg = ""
    if status_parts:
        status_msg = dbc.Alert(" | ".join(status_parts), color="success")

    return (biom_path, meta_json, sid_col, status_msg,
            group_opts, group_disabled, btn_disabled, pipeline_reset)


# ── Callback 2: Load & Plot — build main area + persist session state ────────


@dash_app.callback(
    Output("st-main-area", "children", allow_duplicate=True),
    Output("st-group-col-val", "data"),
    Output("st-all-samples", "data", allow_duplicate=True),
    Output("st-selected-samples", "data", allow_duplicate=True),
    Input("st-btn-load", "n_clicks"),
    State("st-biom-path", "data"),
    State("st-meta-store", "data"),
    State("st-sample-id-col", "data"),
    State("st-group-col", "value"),
    prevent_initial_call=True,
)
def st_on_load(n_clicks, biom_path, meta_json, sid_col, group_col):
    if not n_clicks or not biom_path:
        return no_update, no_update, no_update, no_update

    try:
        count_df = biom_to_count_df(biom_path)
        samples = sorted(count_df.columns.tolist())

        if len(samples) < 3:
            return (
                dbc.Alert("Need at least 3 samples.", color="warning"),
                no_update, no_update, no_update,
            )

        # Resolve group_col
        if not group_col and meta_json and sid_col:
            meta_df = pd.read_json(io.StringIO(meta_json), orient="split")
            gcols = get_group_columns(meta_df, sid_col)
            group_col = gcols[0] if gcols else None

        main_area = _build_main_area(
            biom_path, meta_json, sid_col, group_col, samples, samples,
        )
        return main_area, group_col, samples, samples

    except Exception as e:
        return dbc.Alert(f"Error: {e}", color="danger"), no_update, no_update, no_update


# ── Callback 3: auto-restore from session on page load ──────────────────────


@dash_app.callback(
    Output("st-main-area", "children"),
    Output("st-all-samples", "data"),
    Output("st-selected-samples", "data"),
    Input("st-restore-tick", "n_intervals"),
    State("st-biom-path", "data"),
    State("st-meta-store", "data"),
    State("st-sample-id-col", "data"),
    State("st-group-col-val", "data"),
    State("st-all-samples", "data"),
    State("st-selected-samples", "data"),
)
def st_auto_restore(_n, biom_path, meta_json, sid_col, group_col,
                    all_samples, selected_samples):
    if not biom_path or not all_samples:
        return "", no_update, no_update

    try:
        # Validate the BIOM file still exists
        from pathlib import Path
        if not Path(biom_path).exists():
            # File gone (e.g. temp file from upload) — clear state
            return "", None, None

        selected = selected_samples or all_samples
        main_area = _build_main_area(
            biom_path, meta_json, sid_col, group_col, all_samples, selected,
        )
        return main_area, no_update, no_update

    except Exception:
        return "", None, None


# ── Callback 4: select / deselect all ───────────────────────────────────────


@dash_app.callback(
    Output("st-sample-checklist", "value"),
    Input("st-btn-select-all", "n_clicks"),
    Input("st-btn-deselect-all", "n_clicks"),
    State("st-all-samples", "data"),
    prevent_initial_call=True,
)
def st_toggle_all(sel_clicks, desel_clicks, all_samples):
    triggered = callback_context.triggered_id
    if triggered == "st-btn-select-all":
        return all_samples or []
    if triggered == "st-btn-deselect-all":
        return []
    return no_update


# ── Callback 5: update plots + persist selection when checkboxes change ──────


@dash_app.callback(
    Output("st-tree-graph", "figure"),
    Output("st-pcoa-graph", "figure"),
    Output("st-plot-status", "children"),
    Output("st-selection-info", "children"),
    Output("st-selected-samples", "data", allow_duplicate=True),
    Input("st-sample-checklist", "value"),
    State("st-biom-path", "data"),
    State("st-meta-store", "data"),
    State("st-sample-id-col", "data"),
    State("st-group-col-val", "data"),
    State("st-all-samples", "data"),
    prevent_initial_call=True,
)
def st_update_plots(selected, biom_path, meta_json, sid_col, group_col, all_samples):
    empty_tree = go.Figure()
    empty_pcoa = go.Figure()
    empty_tree.update_layout(template="plotly_dark")
    empty_pcoa.update_layout(template="plotly_dark")

    n_total = len(all_samples) if all_samples else 0
    n_sel = len(selected) if selected else 0
    n_excluded = n_total - n_sel
    info = html.Small(
        f"{n_sel} included / {n_excluded} excluded",
        className="text-muted",
    )

    if not biom_path or not selected or len(selected) < 3:
        for fig in (empty_tree, empty_pcoa):
            fig.add_annotation(
                text="Select at least 3 samples",
                xref="paper", yref="paper", x=0.5, y=0.5,
                showarrow=False, font=dict(size=16),
            )
        return (
            empty_tree, empty_pcoa,
            dbc.Alert("Select at least 3 samples.", color="warning"),
            info, selected,
        )

    try:
        count_df = biom_to_count_df(biom_path)
        keep = [s for s in selected if s in count_df.columns]
        if len(keep) < 3:
            raise ValueError("Fewer than 3 valid samples after filtering.")

        tree_fig, pcoa_fig = _compute_plots(
            count_df[keep], meta_json, sid_col, group_col,
        )
        status = html.Span(
            f"Showing {len(keep)} of {n_total} samples",
            className="text-success",
        )
        return tree_fig, pcoa_fig, status, info, selected

    except Exception as e:
        return (
            empty_tree, empty_pcoa,
            dbc.Alert(f"Error: {e}", color="danger"),
            info, selected,
        )


# ── Callback 6: download filtered BIOM ──────────────────────────────────────


@dash_app.callback(
    Output("st-download-biom", "data"),
    Input("st-btn-download", "n_clicks"),
    State("st-biom-path", "data"),
    State("st-sample-checklist", "value"),
    prevent_initial_call=True,
)
def st_download(n_clicks, biom_path, selected):
    if not n_clicks or not biom_path or not selected:
        return no_update

    table = load_table(biom_path)
    keep = [s for s in selected if s in table.ids(axis="sample")]
    if not keep:
        return no_update

    filtered = table.filter(keep, axis="sample", inplace=False)
    filtered.remove_empty(axis="observation", inplace=True)

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".biom")
    with biom_open(tmp.name, "w") as f:
        filtered.to_hdf5(f, "Filtered by Outlier Detection")
    tmp.close()

    with open(tmp.name, "rb") as f:
        biom_bytes = f.read()

    n_removed = len(table.ids(axis="sample")) - len(keep)
    filename = f"filtered_{len(keep)}samples_excl{n_removed}.biom"
    return dcc.send_bytes(biom_bytes, filename=filename)


# ── Callback 7: clear all session state ──────────────────────────────────────


@dash_app.callback(
    Output("st-main-area", "children", allow_duplicate=True),
    Output("st-biom-path", "data", allow_duplicate=True),
    Output("st-meta-store", "data", allow_duplicate=True),
    Output("st-sample-id-col", "data", allow_duplicate=True),
    Output("st-group-col-val", "data", allow_duplicate=True),
    Output("st-all-samples", "data", allow_duplicate=True),
    Output("st-selected-samples", "data", allow_duplicate=True),
    Output("st-source-status", "children", allow_duplicate=True),
    Output("st-group-col", "options", allow_duplicate=True),
    Output("st-group-col", "disabled", allow_duplicate=True),
    Output("st-group-col", "value", allow_duplicate=True),
    Output("st-btn-load", "disabled", allow_duplicate=True),
    Output("st-pipeline-select", "value", allow_duplicate=True),
    Input("st-btn-clear", "n_clicks"),
    prevent_initial_call=True,
)
def st_clear(n_clicks):
    if not n_clicks:
        return tuple([no_update] * 13)
    return ("", None, None, None, None, None, None,
            "", [], True, None, True, "")


# ── Shared builders ─────────────────────────────────────────────────────────


def _build_main_area(biom_path, meta_json, sid_col, group_col,
                     all_samples, selected_samples):
    """Build the full main area: checklist + plots."""
    count_df = biom_to_count_df(biom_path)
    keep = [s for s in selected_samples if s in count_df.columns]
    n_total = len(all_samples)
    n_sel = len(keep)

    if n_sel >= 3:
        tree_fig, pcoa_fig = _compute_plots(
            count_df[keep], meta_json, sid_col, group_col,
        )
        status_child = html.Span(
            f"Showing {n_sel} of {n_total} samples",
            className="text-success",
        )
    else:
        tree_fig = go.Figure()
        pcoa_fig = go.Figure()
        tree_fig.update_layout(template="plotly_dark")
        pcoa_fig.update_layout(template="plotly_dark")
        status_child = dbc.Alert("Select at least 3 samples.", color="warning")

    checklist_options = [{"label": s, "value": s} for s in all_samples]
    n_excluded = n_total - n_sel

    return html.Div(
        [
            dbc.Row(
                [
                    # Left panel: sample checklist
                    dbc.Col(
                        [
                            html.H5("Samples", className="mb-2"),
                            dbc.Button(
                                "Select All",
                                id="st-btn-select-all",
                                color="outline-secondary",
                                size="sm",
                                className="me-2 mb-2",
                            ),
                            dbc.Button(
                                "Deselect All",
                                id="st-btn-deselect-all",
                                color="outline-secondary",
                                size="sm",
                                className="mb-2",
                            ),
                            html.Div(
                                dbc.Checklist(
                                    id="st-sample-checklist",
                                    options=checklist_options,
                                    value=list(keep),
                                    className="ms-1",
                                ),
                                style={
                                    "maxHeight": "75vh",
                                    "overflowY": "auto",
                                },
                            ),
                            html.Hr(),
                            html.Div(
                                html.Small(
                                    f"{n_sel} included / {n_excluded} excluded",
                                    className="text-muted",
                                ),
                                id="st-selection-info",
                                className="mb-2",
                            ),
                            dbc.Button(
                                "Download Filtered BIOM",
                                id="st-btn-download",
                                color="success",
                                className="w-100",
                            ),
                        ],
                        md=3,
                    ),
                    # Right panel: plots
                    dbc.Col(
                        [
                            html.Div(
                                status_child,
                                id="st-plot-status",
                                className="mb-2",
                            ),
                            dbc.Row(
                                [
                                    dbc.Col(
                                        dcc.Graph(
                                            id="st-tree-graph",
                                            figure=tree_fig,
                                            config={
                                                "displayModeBar": True,
                                                "toImageButtonOptions": {"format": "svg"},
                                            },
                                            style={"height": "75vh"},
                                        ),
                                        md=6,
                                    ),
                                    dbc.Col(
                                        dcc.Graph(
                                            id="st-pcoa-graph",
                                            figure=pcoa_fig,
                                            config={
                                                "displayModeBar": True,
                                                "toImageButtonOptions": {"format": "svg"},
                                            },
                                            style={"height": "75vh"},
                                        ),
                                        md=6,
                                    ),
                                ]
                            ),
                        ],
                        md=9,
                    ),
                ]
            ),
        ]
    )


def _compute_plots(count_df, meta_json, sid_col, group_col):
    """Compute distance, dendrogram, and PCoA from a count DataFrame."""
    dm = compute_distance(count_df, "braycurtis")
    condensed = squareform(dm.data)
    condensed = np.nan_to_num(condensed, nan=0.0)
    Z = linkage(condensed, method="average")
    dn = dendrogram(Z, labels=list(dm.ids), no_plot=True, orientation="left")
    tree_fig = _build_dendrogram(dn, group_col, meta_json, sid_col)
    pcoa_fig = _build_pcoa(dm, group_col, meta_json, sid_col)
    return tree_fig, pcoa_fig


# ── Plot helpers ─────────────────────────────────────────────────────────────


def _resolve_group_map(sample_ids, group_col, meta_json, sid_col):
    """Return a Series mapping sample_id → group label (str), or None."""
    if not meta_json or not sid_col or not group_col:
        return None
    meta_df = pd.read_json(io.StringIO(meta_json), orient="split")
    meta_indexed = meta_df.set_index(meta_df[sid_col].astype(str))
    common = [s for s in sample_ids if s in meta_indexed.index]
    if not common:
        return None
    return meta_indexed.loc[common, group_col].fillna("Unknown").astype(str)


def _build_dendrogram(dn, group_col, meta_json, sid_col):
    """Build a Plotly dendrogram figure, colouring leaves by group."""
    fig = go.Figure()

    ivl = dn["ivl"]
    group_map = _resolve_group_map(ivl, group_col, meta_json, sid_col)
    groups = sorted(group_map.unique().tolist()) if group_map is not None else []
    group_color = {g: COLORS[i % len(COLORS)] for i, g in enumerate(groups)}

    for xs, ys in zip(dn["dcoord"], dn["icoord"]):
        fig.add_trace(
            go.Scatter(
                x=xs, y=ys,
                mode="lines",
                line=dict(color="#888", width=1.2),
                hoverinfo="skip",
                showlegend=False,
            )
        )

    leaf_y = [5 + 10 * i for i in range(len(ivl))]

    if group_map is not None:
        shown_legend = set()
        for i, label in enumerate(ivl):
            grp = group_map.get(label, "Unknown")
            clr = group_color.get(grp, "#aaa")
            show = grp not in shown_legend
            shown_legend.add(grp)
            fig.add_trace(
                go.Scatter(
                    x=[0], y=[leaf_y[i]],
                    mode="markers",
                    marker=dict(color=clr, size=8),
                    name=str(grp),
                    showlegend=show,
                    hovertemplate=f"<b>{label}</b><br>{group_col}: {grp}<extra></extra>",
                )
            )

    fig.update_layout(
        template="plotly_dark",
        title=dict(text="Bray-Curtis UPGMA Dendrogram", y=0.98, yanchor="top"),
        xaxis=dict(title="Bray-Curtis Distance", side="top", zeroline=False),
        yaxis=dict(
            tickvals=leaf_y,
            ticktext=ivl,
            zeroline=False,
            automargin=True,
        ),
        margin=dict(l=10, r=20, t=60, b=30),
        hovermode="closest",
        legend=dict(title=group_col or ""),
    )
    return fig


def _build_pcoa(dm, group_col, meta_json, sid_col):
    """Build a PCoA scatter plot with 95% confidence ellipses."""
    coords, prop_exp = run_pcoa(dm)
    coords.columns = ["Axis1", "Axis2"]
    coords.index = coords.index.astype(str)

    ax1_label = f"PC1 ({prop_exp.get('PC1', 0):.1%})"
    ax2_label = f"PC2 ({prop_exp.get('PC2', 0):.1%})"

    group_map = _resolve_group_map(
        coords.index.tolist(), group_col, meta_json, sid_col
    )

    fig = go.Figure()

    if group_map is not None:
        coords["group"] = coords.index.map(group_map).fillna("Unknown").astype(str)
        groups = sorted(coords["group"].unique())

        for i, grp in enumerate(groups):
            sub = coords[coords["group"] == grp]
            clr = COLORS[i % len(COLORS)]
            fig.add_trace(
                go.Scatter(
                    x=sub["Axis1"], y=sub["Axis2"],
                    mode="markers",
                    name=str(grp),
                    marker=dict(color=clr, size=10, opacity=0.85),
                    text=sub.index,
                    hovertemplate="<b>%{text}</b><br>%{x:.3f}, %{y:.3f}<extra></extra>",
                )
            )

        group_series = coords.index.to_series().map(group_map).fillna("Unknown").astype(str)
        ellipses = compute_confidence_ellipse(coords, group_series)
        for ell in ellipses:
            grp = ell["group"]
            i = groups.index(grp) if grp in groups else 0
            clr = COLORS[i % len(COLORS)]
            fig.add_trace(
                go.Scatter(
                    x=ell["x"], y=ell["y"],
                    mode="lines",
                    fill="toself",
                    fillcolor=f"rgba({_hex_to_rgb(clr)}, 0.10)",
                    line=dict(color=clr, width=1.5, dash="dot"),
                    name=f"{grp} (95% CI)",
                    showlegend=False,
                    hoverinfo="skip",
                )
            )
    else:
        fig.add_trace(
            go.Scatter(
                x=coords["Axis1"], y=coords["Axis2"],
                mode="markers",
                marker=dict(color=COLORS[0], size=10, opacity=0.85),
                text=coords.index,
                hovertemplate="<b>%{text}</b><br>%{x:.3f}, %{y:.3f}<extra></extra>",
                showlegend=False,
            )
        )

    fig.update_layout(
        template="plotly_dark",
        title="PCoA (Bray-Curtis)",
        xaxis_title=ax1_label,
        yaxis_title=ax2_label,
        height=None,
        legend=dict(title=group_col or ""),
    )
    return fig


def _hex_to_rgb(hex_color: str) -> str:
    """Convert '#3498db' to '52, 152, 219'."""
    h = hex_color.lstrip("#")
    return ", ".join(str(int(h[i : i + 2], 16)) for i in (0, 2, 4))
