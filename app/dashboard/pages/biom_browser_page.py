"""
MicrobiomeDash — BIOM Browser page.

Read-only browsing of BIOM files: sample stats, region detection, metadata
lookup, and taxonomy summary. Works with both pipeline datasets (with
UploadMetadata) and standalone uploaded .biom files.
"""
import io

import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
from biom import load_table
from dash import Input, Output, State, dcc, html, no_update

from app.analysis.shared import (
    biom_to_count_df,
    get_dataset_metadata_df,
    get_pipeline_biom_options,
    parse_uploaded_biom,
)
from app.dashboard.app import app as dash_app


# ── Helpers ──────────────────────────────────────────────────────────────────


def _summary_card(title: str, value: str, color: str = "light") -> dbc.Col:
    return dbc.Col(
        dbc.Card(
            dbc.CardBody(
                [
                    html.H6(title, className="card-subtitle text-muted mb-1"),
                    html.H4(value, className="card-title mb-0"),
                ],
                className="py-2 px-3",
            ),
            color=color,
            outline=True,
        ),
        md=2,
        className="mb-2",
    )


# ── Layout ───────────────────────────────────────────────────────────────────


def get_layout():
    pipeline_opts = get_pipeline_biom_options()

    return dbc.Container(
        [
            html.H3("BIOM Browser", className="mb-2"),
            html.P(
                "Inspect a BIOM file: sample names, read counts, ASV counts, "
                "variable region, primer status, and associated upload metadata.",
                className="text-muted mb-4",
            ),
            # ── Section 1: BIOM Input ────────────────────────────────────
            dbc.Card(
                dbc.CardBody(
                    [
                        dbc.Label("Pipeline Dataset", className="fw-bold"),
                        dbc.Select(
                            id="bb-select-pipeline",
                            options=[{"label": "— none —", "value": ""}]
                            + [{"label": o["label"], "value": o["value"]} for o in pipeline_opts],
                            value="",
                        ),
                        html.Hr(className="my-2"),
                        dbc.Label("Or Upload a BIOM File", className="fw-bold"),
                        dcc.Upload(
                            id="bb-upload-biom",
                            children=html.Div(
                                [
                                    "Drag & drop or ",
                                    html.A("select .biom file", className="text-info"),
                                ],
                                className="text-center py-3",
                            ),
                            style={
                                "borderWidth": "2px",
                                "borderStyle": "dashed",
                                "borderRadius": "5px",
                                "borderColor": "#555",
                            },
                            multiple=False,
                            accept=".biom,application/octet-stream,application/x-hdf5",
                        ),
                        html.Div(id="bb-biom-status", className="mt-1 small"),
                    ]
                ),
                className="mb-3",
            ),
            # ── Section 2: Summary Cards ─────────────────────────────────
            html.Div(id="bb-summary-section", style={"display": "none"}),
            # ── Section 3: Sample Table ──────────────────────────────────
            html.Div(
                [
                    dbc.Button(
                        "Download CSV",
                        id="bb-btn-download",
                        color="secondary",
                        size="sm",
                        className="mb-2",
                    ),
                ],
                id="bb-download-btn-wrapper",
                style={"display": "none"},
            ),
            html.Div(id="bb-sample-table-section", style={"display": "none"}),
            # ── Section 4: Taxonomy Summary ──────────────────────────────
            html.Div(id="bb-taxonomy-section", style={"display": "none"}),
            # ── Hidden stores ────────────────────────────────────────────
            dcc.Store(id="bb-biom-path"),
            dcc.Store(id="bb-dataset-id"),
            dcc.Store(id="bb-csv-data"),
            dcc.Store(id="bb-error"),
            dcc.Download(id="bb-download"),
        ],
        fluid=True,
    )


# ── Callback 1: Process BIOM Input ──────────────────────────────────────────


@dash_app.callback(
    Output("bb-biom-path", "data"),
    Output("bb-dataset-id", "data"),
    Output("bb-biom-status", "children"),
    Output("bb-summary-section", "children"),
    Output("bb-summary-section", "style"),
    Output("bb-download-btn-wrapper", "style"),
    Output("bb-sample-table-section", "children"),
    Output("bb-sample-table-section", "style"),
    Output("bb-taxonomy-section", "children"),
    Output("bb-taxonomy-section", "style"),
    Output("bb-csv-data", "data"),
    Output("bb-error", "data"),
    Input("bb-upload-biom", "contents"),
    Input("bb-select-pipeline", "value"),
    State("bb-upload-biom", "filename"),
    prevent_initial_call=True,
)
def on_biom_input(upload_contents, pipeline_value, upload_filename):
    """Resolve BIOM source, compute stats, build all display sections."""
    from dash import ctx

    HIDE = {"display": "none"}
    SHOW = {"display": "block"}
    #        biom_path, dataset_id, status, summary_children, summary_style,
    #        btn_style, table_children, table_style, tax_children, tax_style,
    #        csv_data, error
    EMPTY = (None, None, "", [], HIDE, HIDE, [], HIDE, [], HIDE, None, None)

    triggered = ctx.triggered_id
    biom_path = None
    is_pipeline = False

    def _err(msg):
        return (*EMPTY[:2], dbc.Alert(msg, color="danger", className="py-1"), *EMPTY[3:])

    # ── Resolve BIOM path ────────────────────────────────────────────────
    if triggered == "bb-upload-biom" and upload_contents:
        if not upload_filename or not upload_filename.endswith(".biom"):
            return _err("Please upload a .biom file")
        biom_path, err = parse_uploaded_biom(upload_contents, upload_filename)
        if err:
            return _err(err)
    elif triggered == "bb-select-pipeline" and pipeline_value:
        biom_path = pipeline_value
        is_pipeline = True
    else:
        return EMPTY

    if not biom_path:
        return EMPTY

    # ── Load BIOM and compute stats ──────────────────────────────────────
    try:
        count_df = biom_to_count_df(biom_path)
    except Exception as e:
        return _err(f"Error loading BIOM: {e}")

    sample_ids = list(count_df.columns)
    n_samples = len(sample_ids)
    n_asvs = len(count_df.index)
    sample_reads = count_df.sum(axis=0)
    total_reads = int(sample_reads.sum())
    sample_asv_counts = (count_df > 0).sum(axis=0)

    # ── Region detection ─────────────────────────────────────────────────
    from app.data_manager.biom_ops import detect_region_from_biom

    region_info = detect_region_from_biom(biom_path)
    region = region_info.get("region") or "Unknown"
    confidence = region_info.get("confidence", 0)
    primers_present = region_info.get("primers_present", False)
    median_len = region_info.get("median_len")

    # ── Metadata lookup (pipeline datasets) ─────────────────────────────
    meta_df = None
    meta_cols: list[str] = []
    sample_id_col = None
    if is_pipeline:
        meta_df, sample_id_col = get_dataset_metadata_df(biom_path)
        if meta_df is not None and sample_id_col:
            meta_cols = [c for c in meta_df.columns if c != sample_id_col]
            meta_df = meta_df.set_index(sample_id_col)

    # ── Status text ──────────────────────────────────────────────────────
    source_label = upload_filename if not is_pipeline else "Pipeline dataset"
    status = html.Span(f"Loaded: {source_label}", className="text-success")

    # ── Section 2: Summary cards ─────────────────────────────────────────
    summary_cards = dbc.Row(
        [
            _summary_card("Total Samples", f"{n_samples:,}"),
            _summary_card("Total ASVs", f"{n_asvs:,}"),
            _summary_card("Total Reads", f"{total_reads:,}"),
            _summary_card(
                "Variable Region",
                f"{region} ({confidence:.0%})" if confidence else region,
            ),
            _summary_card("Primers Present", "Yes" if primers_present else "No",
                          color="warning" if primers_present else "light"),
            _summary_card("Median ASV Length", f"{median_len} bp" if median_len else "N/A"),
        ],
        className="mb-3",
    )
    source_line = html.Small(
        f"Source: {'pipeline dataset' if is_pipeline else 'uploaded file'}",
        className="text-muted d-block mb-2",
    )
    summary_section = [summary_cards, source_line]

    # ── Section 3: Sample table ──────────────────────────────────────────
    table_header_cols = ["Sample Name", "Read Count", "ASV Count"] + meta_cols
    header = html.Thead(html.Tr([html.Th(c) for c in table_header_cols]))

    rows = []
    csv_rows = [table_header_cols]
    for sid in sorted(sample_ids):
        reads = int(sample_reads[sid])
        asvs = int(sample_asv_counts[sid])
        cells = [html.Td(sid), html.Td(f"{reads:,}"), html.Td(f"{asvs:,}")]
        csv_row = [sid, str(reads), str(asvs)]

        if meta_df is not None and meta_cols:
            for col in meta_cols:
                val = ""
                if sid in meta_df.index:
                    val = str(meta_df.loc[sid, col])
                cells.append(html.Td(val))
                csv_row.append(val)

        rows.append(html.Tr(cells))
        csv_rows.append(csv_row)

    body = html.Tbody(rows)
    sample_table = dbc.Table(
        [header, body],
        color="dark",
        striped=True,
        hover=True,
        responsive=True,
        style={"maxHeight": "500px", "overflowY": "auto", "display": "block"},
    )

    # CSV as string for download
    csv_str = "\n".join(",".join(r) for r in csv_rows)

    sample_section = [
        html.H5("Sample Table", className="mt-2 mb-2"),
        sample_table,
    ]

    # ── Section 4: Taxonomy summary ──────────────────────────────────────
    taxonomy_section: list = []
    taxonomy_style = HIDE

    def _tax_table(level, top_n):
        from app.analysis.taxonomy import aggregate_taxonomy

        df = aggregate_taxonomy(biom_path, level, top_n=top_n)
        if df.empty:
            return None
        mean_abund = df.mean(axis=1)
        prevalence = (df > 0).mean(axis=1)
        hdr = html.Thead(
            html.Tr([html.Th(level), html.Th("Mean Rel. Abundance"), html.Th("Prevalence")])
        )
        tbody_rows = [
            html.Tr([
                html.Td(taxon),
                html.Td(f"{mean_abund[taxon]:.3f}"),
                html.Td(f"{prevalence[taxon]:.1%}"),
            ])
            for taxon in mean_abund.index
        ]
        return dbc.Table(
            [hdr, html.Tbody(tbody_rows)],
            color="dark",
            striped=True,
            hover=True,
            responsive=True,
        )

    try:
        phylum_table = _tax_table("Phylum", 10)
        if phylum_table:
            taxonomy_section.append(html.H5("Top 10 Phyla", className="mt-3 mb-2"))
            taxonomy_section.append(phylum_table)
            taxonomy_style = SHOW

        genus_table = _tax_table("Genus", 20)
        if genus_table:
            taxonomy_section.append(html.H5("Top 20 Genera", className="mt-3 mb-2"))
            taxonomy_section.append(genus_table)
            taxonomy_style = SHOW
    except Exception:
        pass

    return (
        biom_path,
        None,
        status,
        summary_section,
        SHOW,
        SHOW,
        sample_section,
        SHOW,
        taxonomy_section,
        taxonomy_style,
        csv_str,
        None,
    )


# ── Callback 2: CSV Download ────────────────────────────────────────────────


@dash_app.callback(
    Output("bb-download", "data"),
    Input("bb-btn-download", "n_clicks"),
    State("bb-csv-data", "data"),
    prevent_initial_call=True,
)
def on_download(n_clicks, csv_data):
    if not n_clicks or not csv_data:
        return no_update
    return dcc.send_string(csv_data, "biom_browser_samples.csv")
