"""
MicrobiomeDash — Datasets page: upload BIOM, detect region, extract sub-regions.
"""
import base64
import tempfile
from pathlib import Path

import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, no_update, ALL, ctx

from app.analysis.shared import get_pipeline_biom_options
from app.config import DATASET_DIR
from app.dashboard.app import app as dash_app


def get_layout():
    return dbc.Container(
        [
            html.H3("V-Region Extraction", className="mb-2"),
            html.P(
                "Upload a BIOM file to inspect its variable region and extract "
                "sub-regions (e.g., V4 from V3-V4) for cross-study comparison.",
                className="text-muted mb-4",
            ),
            # ── Pipeline dataset selector ──────────────────────────────
            dbc.Label("Pipeline Dataset", className="fw-bold"),
            dbc.Select(
                id="ds-select-pipeline",
                options=[{"label": "— none —", "value": ""}]
                + get_pipeline_biom_options(),
                value="",
                className="mb-2",
            ),
            html.Div(
                "— or upload a file —",
                className="text-center text-muted small my-2",
            ),
            # ── Upload BIOM ──────────────────────────────────────────────
            dbc.Label("Upload BIOM", className="fw-bold"),
            dcc.Upload(
                id="ds-upload-biom",
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
            html.Div(id="ds-upload-status", className="mt-1 small mb-3"),
            # ── BIOM Details (hidden initially) ──────────────────────────
            html.Div(
                id="ds-biom-details",
                style={"display": "none"},
            ),
            # ── Extract Region ───────────────────────────────────────────
            html.Div(
                id="ds-extract-section",
                children=[
                    html.H5("Extract Region", className="mt-3"),
                    html.P(
                        "Extract a sub-region from the uploaded BIOM file.",
                        className="text-muted small",
                    ),
                    dbc.InputGroup(
                        [
                            dbc.Select(
                                id="ds-extract-target",
                                options=[],
                                placeholder="Select target region...",
                                style={"maxWidth": "200px"},
                            ),
                            dbc.Button(
                                "Extract",
                                id="ds-btn-extract",
                                color="primary",
                            ),
                        ],
                        className="mb-3",
                    ),
                    html.Div(id="ds-extract-result"),
                    dcc.Download(id="ds-download-extracted"),
                ],
                style={"display": "none"},
            ),
            # ── Pipeline Datasets ────────────────────────────────────────
            html.Hr(className="mt-4"),
            html.H5("Pipeline Datasets"),
            html.P(
                "Completed pipeline runs with BIOM files. "
                "Select one to inspect or extract a sub-region.",
                className="text-muted small",
            ),
            html.Div(
                id="ds-pipeline-table",
                style={"maxHeight": "400px", "overflowY": "auto"},
            ),
            # ── Pipeline extract section ─────────────────────────────────
            html.Div(
                id="ds-pipeline-extract-section",
                children=[
                    html.H5("Extract from Pipeline Dataset", className="mt-3"),
                    html.Div(id="ds-pipeline-biom-info"),
                    dbc.InputGroup(
                        [
                            dbc.Select(
                                id="ds-pipeline-extract-target",
                                options=[],
                                placeholder="Select target region...",
                                style={"maxWidth": "200px"},
                            ),
                            dbc.Button(
                                "Extract",
                                id="ds-btn-pipeline-extract",
                                color="primary",
                            ),
                        ],
                        className="mb-3",
                    ),
                    html.Div(id="ds-pipeline-extract-result"),
                    dcc.Download(id="ds-download-pipeline-extracted"),
                ],
                style={"display": "none"},
            ),
            # ── Stores ───────────────────────────────────────────────────
            dcc.Store(id="ds-biom-info"),           # uploaded BIOM detection info
            dcc.Store(id="ds-biom-temppath"),        # path to saved temp BIOM
            dcc.Store(id="ds-pipeline-dataset-id"),  # selected pipeline dataset ID
            # One-shot interval to build table on page load
            dcc.Interval(
                id="ds-init", interval=200, max_intervals=1, n_intervals=0
            ),
        ],
        fluid=True,
    )


# ── Callbacks ────────────────────────────────────────────────────────────────


@dash_app.callback(
    Output("ds-biom-info", "data"),
    Output("ds-biom-temppath", "data"),
    Output("ds-upload-status", "children"),
    Output("ds-biom-details", "children"),
    Output("ds-biom-details", "style"),
    Output("ds-extract-section", "style"),
    Input("ds-upload-biom", "contents"),
    Input("ds-select-pipeline", "value"),
    State("ds-upload-biom", "filename"),
    prevent_initial_call=True,
)
def on_upload_or_select_biom(contents, pipeline_path, filename):
    """Handle BIOM from upload or pipeline dataset selection."""
    triggered = ctx.triggered_id
    empty = (None, None, "", "", {"display": "none"}, {"display": "none"})

    if triggered == "ds-select-pipeline":
        # ── Pipeline branch ───────────────────────────────────────
        if not pipeline_path:
            return empty
        try:
            from app.data_manager.biom_ops import detect_region_from_biom

            info = detect_region_from_biom(pipeline_path)
        except Exception as e:
            return (
                None, None,
                dbc.Alert(f"Error reading BIOM: {e}", color="danger", className="py-1"),
                "", {"display": "none"}, {"display": "none"},
            )

        details, extract_style = _build_biom_details(info)
        label = Path(pipeline_path).parent.name  # dataset id
        return (
            info,
            pipeline_path,
            html.Span(f"  Pipeline dataset #{label}", className="text-success"),
            details,
            {"display": "block"},
            extract_style,
        )

    # ── Upload branch ─────────────────────────────────────────────
    if not contents or not filename:
        return empty

    if not filename.endswith(".biom"):
        return (
            None, None,
            dbc.Alert("Please upload a .biom file", color="danger", className="py-1"),
            "", {"display": "none"}, {"display": "none"},
        )

    try:
        _, content_string = contents.split(",", 1)
        data = base64.b64decode(content_string)

        tmp_dir = tempfile.mkdtemp()
        tmp_path = str(Path(tmp_dir) / filename)
        Path(tmp_path).write_bytes(data)

        from app.data_manager.biom_ops import detect_region_from_biom

        info = detect_region_from_biom(tmp_path)
    except Exception as e:
        return (
            None, None,
            dbc.Alert(f"Error reading BIOM: {e}", color="danger", className="py-1"),
            "", {"display": "none"}, {"display": "none"},
        )

    details, extract_style = _build_biom_details(info)

    return (
        info,
        tmp_path,
        html.Span(f"  {filename}", className="text-success"),
        details,
        {"display": "block"},
        extract_style,
    )


@dash_app.callback(
    Output("ds-extract-target", "options"),
    Input("ds-biom-info", "data"),
    prevent_initial_call=True,
)
def populate_extract_options(info):
    """Populate target region dropdown from detected region."""
    if not info or not info.get("region"):
        return []

    from app.data_manager.biom_ops import get_valid_extractions

    targets = get_valid_extractions(info["region"])
    return [{"label": t, "value": t} for t in targets]


@dash_app.callback(
    Output("ds-extract-result", "children"),
    Output("ds-download-extracted", "data"),
    Input("ds-btn-extract", "n_clicks"),
    State("ds-biom-temppath", "data"),
    State("ds-biom-info", "data"),
    State("ds-extract-target", "value"),
    prevent_initial_call=True,
)
def on_extract(n_clicks, biom_path, info, target_region):
    """Run region extraction on uploaded BIOM."""
    if not n_clicks or not biom_path or not info or not target_region:
        return no_update, no_update

    source_region = info.get("region")
    if not source_region:
        return dbc.Alert("No source region detected.", color="danger"), no_update

    try:
        from app.data_manager.biom_ops import extract_region

        result = extract_region(biom_path, source_region, target_region)
    except Exception as e:
        return dbc.Alert(f"Extraction failed: {e}", color="danger"), no_update

    stats = dbc.Alert(
        [
            html.Strong(f"Extracted {target_region} from {source_region}"),
            html.Br(),
            f"Input: {result['n_input']} ASVs  ->  Output: {result['n_output']} ASVs",
            html.Br(),
            f"Failed: {result['n_failed']}  |  Collapsed: {result['n_collapsed']} "
            "(identical sequences merged)",
        ],
        color="success",
    )

    # Keep original filename stem, append target region
    original_stem = Path(biom_path).stem  # e.g. "my_samples" from "my_samples.biom"
    dl = dcc.send_bytes(
        result["biom_bytes"],
        filename=f"{original_stem}_{target_region}.biom",
    )

    return stats, dl


# ── Pipeline Datasets table ──────────────────────────────────────────────────


@dash_app.callback(
    Output("ds-pipeline-table", "children"),
    Input("ds-init", "n_intervals"),
)
def load_pipeline_datasets(_):
    return _build_pipeline_table()


@dash_app.callback(
    Output("ds-pipeline-dataset-id", "data"),
    Output("ds-pipeline-extract-section", "style"),
    Output("ds-pipeline-biom-info", "children"),
    Output("ds-pipeline-extract-target", "options"),
    Input({"type": "btn-ds-inspect", "index": ALL}, "n_clicks"),
    prevent_initial_call=True,
)
def on_inspect_pipeline(n_clicks_list):
    """Inspect a pipeline dataset's BIOM file."""
    if not any(n_clicks_list):
        return no_update, no_update, no_update, no_update

    triggered = ctx.triggered_id
    if not triggered:
        return no_update, no_update, no_update, no_update
    dataset_id = triggered["index"]

    biom_path = DATASET_DIR / str(dataset_id) / "asv_table.biom"
    if not biom_path.exists():
        return (
            no_update,
            no_update,
            dbc.Alert("BIOM file not found.", color="danger"),
            [],
        )

    try:
        from app.data_manager.biom_ops import detect_region_from_biom, get_valid_extractions

        info = detect_region_from_biom(str(biom_path))
    except Exception as e:
        return (
            no_update,
            no_update,
            dbc.Alert(f"Error: {e}", color="danger"),
            [],
        )

    region = info.get("region")
    extractions = get_valid_extractions(region) if region else []

    card_body = [
        html.Div(
            [
                html.Strong(f"Dataset #{dataset_id}"),
                html.Span(" — ", className="mx-1"),
                dbc.Badge(
                    region or "Unknown",
                    color="primary" if region else "secondary",
                    className="me-1",
                ),
                html.Span(
                    f"{info['n_asvs']} ASVs  |  {info['n_samples']} samples"
                ),
            ]
        ),
    ]
    if info.get("primers_present"):
        card_body.append(
            dbc.Alert(
                "Primer sequences detected in ASVs. This BIOM file may not have "
                "been primer-trimmed. Run Cutadapt or trim primers before analysis.",
                color="warning",
                className="py-2 mt-2 mb-0",
            )
        )
    info_card = dbc.Card(dbc.CardBody(card_body), className="mb-3")

    options = [{"label": t, "value": t} for t in extractions]
    show = {"display": "block"} if extractions else {"display": "none"}

    return dataset_id, show, info_card, options


@dash_app.callback(
    Output("ds-pipeline-extract-result", "children"),
    Output("ds-download-pipeline-extracted", "data"),
    Input("ds-btn-pipeline-extract", "n_clicks"),
    State("ds-pipeline-dataset-id", "data"),
    State("ds-pipeline-extract-target", "value"),
    prevent_initial_call=True,
)
def on_pipeline_extract(n_clicks, dataset_id, target_region):
    """Run region extraction on a pipeline dataset's BIOM."""
    if not n_clicks or not dataset_id or not target_region:
        return no_update, no_update

    biom_path = DATASET_DIR / str(dataset_id) / "asv_table.biom"
    if not biom_path.exists():
        return dbc.Alert("BIOM file not found.", color="danger"), no_update

    try:
        from app.data_manager.biom_ops import detect_region_from_biom, extract_region

        info = detect_region_from_biom(str(biom_path))
        source_region = info.get("region")
        if not source_region:
            return dbc.Alert("Could not detect source region.", color="danger"), no_update

        result = extract_region(str(biom_path), source_region, target_region)
    except Exception as e:
        return dbc.Alert(f"Extraction failed: {e}", color="danger"), no_update

    stats = dbc.Alert(
        [
            html.Strong(f"Extracted {target_region} from {source_region}"),
            html.Br(),
            f"Input: {result['n_input']} ASVs  ->  Output: {result['n_output']} ASVs",
            html.Br(),
            f"Failed: {result['n_failed']}  |  Collapsed: {result['n_collapsed']} "
            "(identical sequences merged)",
        ],
        color="success",
    )

    from app.db.database import get_session
    from app.db.models import Dataset

    with get_session() as db:
        ds = db.query(Dataset).filter(Dataset.id == dataset_id).first()

    name_part = ds.name if ds else f"dataset{dataset_id}"
    dl = dcc.send_bytes(
        result["biom_bytes"],
        filename=f"{name_part}_{target_region}.biom",
    )

    return stats, dl


# ── Helpers ──────────────────────────────────────────────────────────────────


def _build_biom_details(info):
    """Build details card and extract section style from detection info."""
    from app.data_manager.biom_ops import get_valid_extractions

    region = info.get("region")
    extractions = get_valid_extractions(region) if region else []

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
                    + (f"  |  median length: {info.get('median_len', '?')} bp" if info.get("median_len") else ""),
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

    if not extractions:
        card_children.append(
            html.Small(
                f"No sub-region extraction available for {region or 'this region'}. "
                "Extraction is supported for multi-region amplicons (V3-V4, V4-V5).",
                className="text-muted d-block mt-2",
            )
        )

    details = dbc.Card(dbc.CardBody(card_children), className="mb-3")
    extract_style = {"display": "block"} if extractions else {"display": "none"}
    return details, extract_style


def _build_pipeline_table():
    """Build table of completed pipeline datasets with BIOM files."""
    from app.db.database import get_session
    from app.db.models import Dataset

    with get_session() as db:
        datasets = (
            db.query(Dataset)
            .filter(
                Dataset.source_type == "pipeline",
                Dataset.status == "complete",
            )
            .order_by(Dataset.created_at.desc())
            .all()
        )

    if not datasets:
        return html.P("No completed pipeline datasets yet.", className="text-muted")

    # Filter to those that actually have a BIOM file
    rows = []
    for d in datasets:
        biom_path = DATASET_DIR / str(d.id) / "asv_table.biom"
        if not biom_path.exists():
            continue

        rows.append(
            html.Tr(
                [
                    html.Td(d.name),
                    html.Td(d.variable_region or "—"),
                    html.Td(d.sample_count or "—"),
                    html.Td(d.asv_count or "—"),
                    html.Td(
                        d.created_at.strftime("%Y-%m-%d %H:%M")
                        if d.created_at else ""
                    ),
                    html.Td(
                        dbc.Button(
                            "Inspect",
                            id={"type": "btn-ds-inspect", "index": d.id},
                            color="info",
                            size="sm",
                            outline=True,
                        )
                    ),
                ]
            )
        )

    if not rows:
        return html.P("No pipeline datasets with BIOM files.", className="text-muted")

    return dbc.Table(
        [
            html.Thead(
                html.Tr(
                    [
                        html.Th("Name"),
                        html.Th("Region"),
                        html.Th("Samples"),
                        html.Th("ASVs"),
                        html.Th("Created"),
                        html.Th(""),
                    ]
                )
            ),
            html.Tbody(rows),
        ],
        bordered=True,
        hover=True,
        size="sm",
    )
