"""
MicrobiomeDash — BIOM & MOTHUR bidirectional conversion page.

Tab 1 (BIOM to MOTHUR): Select pipeline dataset or upload BIOM → download
    MOTHUR-compatible FASTA + count_table as ZIP.
Tab 2 (MOTHUR to BIOM): Upload FASTA + count_table → download BIOM file.
"""
import base64
import tempfile
from pathlib import Path

import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, no_update

from app.analysis.shared import get_pipeline_biom_options
from app.dashboard.app import app as dash_app


# ═════════════════════════════════════════════════════════════════════════════
# Layout
# ═════════════════════════════════════════════════════════════════════════════


def get_layout():
    pipeline_opts = get_pipeline_biom_options()
    return dbc.Container(
        [
            html.H3("BIOM & MOTHUR", className="mb-2"),
            html.P(
                "Convert between BIOM and MOTHUR formats.",
                className="text-muted mb-4",
            ),
            dbc.Tabs(
                [
                    dbc.Tab(
                        _tab_biom_to_mothur(pipeline_opts),
                        label="BIOM to MOTHUR",
                        tab_id="tab-b2m",
                    ),
                    dbc.Tab(
                        _tab_mothur_to_biom(),
                        label="MOTHUR to BIOM",
                        tab_id="tab-m2b",
                    ),
                ],
                id="mt-tabs",
                active_tab="tab-b2m",
                className="mb-3",
            ),
        ],
        fluid=True,
    )


def _tab_biom_to_mothur(pipeline_opts):
    """Layout for the BIOM → MOTHUR tab."""
    return html.Div(
        [
            html.P(
                "Select a pipeline dataset or upload a BIOM file. "
                "Downloads a ZIP containing a FASTA (representative sequences) "
                "and a count_table (sample abundances).",
                className="text-muted mb-3 mt-3",
            ),
            # ── Pipeline dataset selector ────────────────────────────────
            dbc.Label("BIOM Table", className="fw-bold"),
            dbc.Select(
                id="mt-select-pipeline",
                options=[{"label": "— none —", "value": ""}]
                + pipeline_opts,
                value="",
                className="mb-2",
            ),
            html.Div(
                "— or upload —",
                className="text-center text-muted small mb-2",
            ),
            # ── Upload BIOM ──────────────────────────────────────────────
            dcc.Upload(
                id="mt-upload-biom",
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
            html.Div(id="mt-upload-status", className="mt-1 small mb-3"),
            # ── Convert button ────────────────────────────────────────────
            dbc.Button(
                "Convert to MOTHUR",
                id="mt-btn-convert",
                color="primary",
                size="lg",
                className="w-100 mb-3",
                disabled=True,
            ),
            html.Div(id="mt-convert-error"),
            html.Div(id="mt-biom-details", style={"display": "none"}),
            # ── Download ──────────────────────────────────────────────────
            html.Div(
                id="mt-download-section",
                children=[
                    dbc.Button(
                        "Download MOTHUR Files (.zip)",
                        id="mt-btn-download",
                        color="success",
                        size="lg",
                        className="w-100 mt-2",
                    ),
                    dcc.Download(id="mt-download-zip"),
                ],
                style={"display": "none"},
            ),
            # ── Stores ───────────────────────────────────────────────────
            dcc.Store(id="mt-biom-temppath"),
            dcc.Store(id="mt-biom-filename"),
            dcc.Store(id="mt-zip-bytes"),
        ]
    )


def _tab_mothur_to_biom():
    """Layout for the MOTHUR → BIOM tab."""
    return html.Div(
        [
            html.P(
                "Upload a FASTA file (representative sequences) and a MOTHUR "
                "count_table to produce a BIOM file.",
                className="text-muted mb-3 mt-3",
            ),
            # ── Upload FASTA ─────────────────────────────────────────────
            dbc.Label("FASTA File", className="fw-bold"),
            dcc.Upload(
                id="m2b-upload-fasta",
                children=html.Div(
                    [
                        "Drag & drop or ",
                        html.A(
                            "select .fasta / .fa / .fna file",
                            className="text-info",
                        ),
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
            ),
            html.Div(id="m2b-fasta-status", className="mt-1 small mb-3"),
            # ── Upload count_table ───────────────────────────────────────
            dbc.Label("Count Table", className="fw-bold"),
            dcc.Upload(
                id="m2b-upload-count",
                children=html.Div(
                    [
                        "Drag & drop or ",
                        html.A(
                            "select .count_table file",
                            className="text-info",
                        ),
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
            ),
            html.Div(id="m2b-count-status", className="mt-1 small mb-3"),
            # ── Convert button ────────────────────────────────────────────
            dbc.Button(
                "Convert to BIOM",
                id="m2b-btn-convert",
                color="primary",
                size="lg",
                className="w-100 mb-3",
                disabled=True,
            ),
            html.Div(id="m2b-convert-error"),
            html.Div(id="m2b-result-info", style={"display": "none"}),
            # ── Download ──────────────────────────────────────────────────
            html.Div(
                id="m2b-download-section",
                children=[
                    dbc.Button(
                        "Download BIOM File",
                        id="m2b-btn-download",
                        color="success",
                        size="lg",
                        className="w-100 mt-2",
                    ),
                    dcc.Download(id="m2b-download-biom"),
                ],
                style={"display": "none"},
            ),
            # ── Stores ───────────────────────────────────────────────────
            dcc.Store(id="m2b-fasta-data"),
            dcc.Store(id="m2b-count-data"),
            dcc.Store(id="m2b-biom-bytes"),
            dcc.Store(id="m2b-output-name"),
        ]
    )


# ═════════════════════════════════════════════════════════════════════════════
# Callbacks — BIOM to MOTHUR (unchanged logic, "mt-" prefix)
# ═════════════════════════════════════════════════════════════════════════════


@dash_app.callback(
    Output("mt-biom-temppath", "data"),
    Output("mt-biom-filename", "data"),
    Output("mt-upload-status", "children"),
    Output("mt-btn-convert", "disabled"),
    Output("mt-select-pipeline", "value"),
    Output("mt-biom-details", "style"),
    Output("mt-download-section", "style"),
    Output("mt-convert-error", "children"),
    Input("mt-upload-biom", "contents"),
    State("mt-upload-biom", "filename"),
    prevent_initial_call=True,
)
def on_upload_biom(contents, filename):
    """Decode uploaded BIOM, save to temp, enable convert button."""
    hide = {"display": "none"}
    if not contents or not filename:
        return None, None, "", True, no_update, hide, hide, ""

    if not filename.endswith(".biom"):
        return (
            None, None,
            dbc.Alert("Please upload a .biom file", color="danger", className="py-1"),
            True, no_update, hide, hide, "",
        )

    try:
        _, content_string = contents.split(",", 1)
        data = base64.b64decode(content_string)

        tmp_dir = tempfile.mkdtemp()
        tmp_path = str(Path(tmp_dir) / filename)
        Path(tmp_path).write_bytes(data)
    except Exception as e:
        return (
            None, None,
            dbc.Alert(f"Error reading BIOM: {e}", color="danger", className="py-1"),
            True, no_update, hide, hide, "",
        )

    return (
        tmp_path,
        filename,
        html.Span(f"  {filename}", className="text-success"),
        False,   # enable convert button
        "",      # clear pipeline selector
        hide,    # hide old details
        hide,    # hide old download
        "",      # clear error
    )


@dash_app.callback(
    Output("mt-biom-temppath", "data", allow_duplicate=True),
    Output("mt-biom-filename", "data", allow_duplicate=True),
    Output("mt-upload-status", "children", allow_duplicate=True),
    Output("mt-btn-convert", "disabled", allow_duplicate=True),
    Output("mt-biom-details", "style", allow_duplicate=True),
    Output("mt-download-section", "style", allow_duplicate=True),
    Output("mt-convert-error", "children", allow_duplicate=True),
    Input("mt-select-pipeline", "value"),
    prevent_initial_call=True,
)
def on_select_pipeline(biom_path):
    """Set BIOM path from pipeline dataset and enable convert button."""
    hide = {"display": "none"}
    if not biom_path:
        return None, None, "", True, hide, hide, ""

    p = Path(biom_path)
    if not p.exists():
        return (
            None, None,
            dbc.Alert("BIOM file not found on disk.", color="danger", className="py-1"),
            True, hide, hide, "",
        )

    return (
        str(p),
        p.name,
        html.Span(f"  {p.name} (pipeline)", className="text-success"),
        False,   # enable convert button
        hide,    # hide old details
        hide,    # hide old download
        "",      # clear error
    )


@dash_app.callback(
    Output("mt-biom-details", "children"),
    Output("mt-biom-details", "style", allow_duplicate=True),
    Output("mt-download-section", "style", allow_duplicate=True),
    Output("mt-convert-error", "children", allow_duplicate=True),
    Output("mt-zip-bytes", "data"),
    Input("mt-btn-convert", "n_clicks"),
    State("mt-biom-temppath", "data"),
    State("mt-biom-filename", "data"),
    prevent_initial_call=True,
)
def on_convert(n_clicks, biom_path, original_filename):
    """Run the conversion: show BIOM details and prepare download."""
    hide = {"display": "none"}
    show = {"display": "block"}
    if not n_clicks or not biom_path:
        return "", hide, hide, "", no_update

    try:
        details = _build_biom_details(biom_path)

        from app.data_manager.mothur_convert import biom_to_mothur_zip

        stem = Path(original_filename).stem if original_filename else "mothur"
        zip_bytes = biom_to_mothur_zip(biom_path, name=stem)
    except Exception as e:
        return (
            "",
            hide,
            hide,
            dbc.Alert(f"Conversion failed: {e}", color="danger"),
            no_update,
        )

    return (
        details,
        show,
        show,
        "",
        base64.b64encode(zip_bytes).decode("ascii"),
    )


@dash_app.callback(
    Output("mt-download-zip", "data"),
    Input("mt-btn-download", "n_clicks"),
    State("mt-zip-bytes", "data"),
    State("mt-biom-filename", "data"),
    prevent_initial_call=True,
)
def on_download(n_clicks, zip_b64, original_filename):
    """Send the already-converted ZIP to the browser."""
    if not n_clicks or not zip_b64:
        return no_update

    zip_bytes = base64.b64decode(zip_b64)
    stem = Path(original_filename).stem if original_filename else "mothur"
    return dcc.send_bytes(zip_bytes, filename=f"{stem}_mothur.zip")


# ═════════════════════════════════════════════════════════════════════════════
# Callbacks — MOTHUR to BIOM ("m2b-" prefix)
# ═════════════════════════════════════════════════════════════════════════════


@dash_app.callback(
    Output("m2b-fasta-data", "data"),
    Output("m2b-fasta-status", "children"),
    Input("m2b-upload-fasta", "contents"),
    State("m2b-upload-fasta", "filename"),
    prevent_initial_call=True,
)
def on_upload_fasta(contents, filename):
    """Store uploaded FASTA contents."""
    if not contents or not filename:
        return None, ""
    return (
        {"contents": contents, "filename": filename},
        html.Span(f"  {filename}", className="text-success"),
    )


@dash_app.callback(
    Output("m2b-count-data", "data"),
    Output("m2b-count-status", "children"),
    Input("m2b-upload-count", "contents"),
    State("m2b-upload-count", "filename"),
    prevent_initial_call=True,
)
def on_upload_count(contents, filename):
    """Store uploaded count_table contents."""
    if not contents or not filename:
        return None, ""
    return (
        {"contents": contents, "filename": filename},
        html.Span(f"  {filename}", className="text-success"),
    )


@dash_app.callback(
    Output("m2b-btn-convert", "disabled"),
    Input("m2b-fasta-data", "data"),
    Input("m2b-count-data", "data"),
)
def toggle_m2b_convert(fasta_data, count_data):
    """Enable convert button when both files are uploaded."""
    return not (fasta_data and count_data)


@dash_app.callback(
    Output("m2b-result-info", "children"),
    Output("m2b-result-info", "style"),
    Output("m2b-download-section", "style"),
    Output("m2b-convert-error", "children"),
    Output("m2b-biom-bytes", "data"),
    Output("m2b-output-name", "data"),
    Input("m2b-btn-convert", "n_clicks"),
    State("m2b-fasta-data", "data"),
    State("m2b-count-data", "data"),
    prevent_initial_call=True,
)
def on_m2b_convert(n_clicks, fasta_data, count_data):
    """Convert MOTHUR files to BIOM."""
    hide = {"display": "none"}
    show = {"display": "block"}
    if not n_clicks or not fasta_data or not count_data:
        return "", hide, hide, "", no_update, no_update

    try:
        tmp_dir = tempfile.mkdtemp()

        # Save FASTA
        fasta_path = Path(tmp_dir) / fasta_data["filename"]
        _save_upload(fasta_data["contents"], fasta_path)

        # Save count_table
        count_path = Path(tmp_dir) / count_data["filename"]
        _save_upload(count_data["contents"], count_path)

        from app.data_manager.mothur_convert import mothur_to_biom

        biom_bytes = mothur_to_biom(str(fasta_path), str(count_path))

        # Build summary from the generated BIOM
        from biom import load_table as _load_biom

        tmp_biom = Path(tmp_dir) / "result.biom"
        tmp_biom.write_bytes(biom_bytes)
        table = _load_biom(str(tmp_biom))
        n_obs = len(table.ids(axis="observation"))
        n_samples = len(table.ids(axis="sample"))

        # Output name from the FASTA filename stem
        output_name = fasta_path.stem

        info_card = dbc.Card(
            dbc.CardBody(
                html.Span(f"{n_obs} ASVs  |  {n_samples} samples")
            ),
            className="mb-3",
        )
    except Exception as e:
        return (
            "",
            hide,
            hide,
            dbc.Alert(f"Conversion failed: {e}", color="danger"),
            no_update,
            no_update,
        )

    return (
        info_card,
        show,
        show,
        "",
        base64.b64encode(biom_bytes).decode("ascii"),
        output_name,
    )


@dash_app.callback(
    Output("m2b-download-biom", "data"),
    Input("m2b-btn-download", "n_clicks"),
    State("m2b-biom-bytes", "data"),
    State("m2b-output-name", "data"),
    prevent_initial_call=True,
)
def on_m2b_download(n_clicks, biom_b64, output_name):
    """Send the converted BIOM file to the browser."""
    if not n_clicks or not biom_b64:
        return no_update

    biom_bytes = base64.b64decode(biom_b64)
    name = output_name or "converted"
    return dcc.send_bytes(biom_bytes, filename=f"{name}.biom")


# ═════════════════════════════════════════════════════════════════════════════
# Helpers
# ═════════════════════════════════════════════════════════════════════════════


def _save_upload(contents_str: str, dest: Path):
    """Decode a Dash Upload contents string and save to disk."""
    _, content_string = contents_str.split(",", 1)
    data = base64.b64decode(content_string)
    dest.write_bytes(data)


def _build_biom_details(biom_path: str):
    """Build the BIOM details card from a file path."""
    from app.data_manager.biom_ops import detect_region_from_biom

    info = detect_region_from_biom(biom_path)
    region = info.get("region")

    card_children = [
        html.Div(
            [
                dbc.Badge(
                    region or "Unknown region",
                    color="primary" if region else "secondary",
                    className="me-2 fs-6",
                ),
                html.Span(
                    f"{info['n_asvs']} ASVs  |  {info['n_samples']} samples"
                    + (
                        f"  |  median length: {info.get('median_len', '?')} bp"
                        if info.get("median_len")
                        else ""
                    ),
                ),
            ]
        ),
        html.Small(
            f"Confidence: {info.get('confidence', 0):.0%}",
            className="text-muted d-block mt-1",
        ),
    ]

    if info.get("primers_present"):
        card_children.append(
            dbc.Alert(
                "Primer sequences detected in ASVs. This BIOM file may not have "
                "been primer-trimmed. Run Cutadapt or trim primers before analysis.",
                color="warning",
                className="py-2 mt-2 mb-0",
            )
        )

    return dbc.Card(dbc.CardBody(card_children), className="mb-3")
