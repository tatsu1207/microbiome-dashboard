"""
MicrobiomeDash — Combine Datasets page: merge multiple BIOM files.

Supports two merge strategies:
- Same region: merge by sequence identity (produces BIOM)
- Different regions: merge by taxonomy at a chosen level (produces TSV)
"""
import base64
import tempfile
from pathlib import Path

import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, no_update

from app.config import DATASET_DIR
from app.dashboard.app import app as dash_app


def _get_pipeline_options():
    """Return checklist options for completed pipeline datasets with BIOM files."""
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

    options = []
    for d in datasets:
        biom_path = DATASET_DIR / str(d.id) / "asv_table.biom"
        if not biom_path.exists():
            continue
        region = d.variable_region or "?"
        label = (
            f"#{d.id} {d.name} ({region}, "
            f"{d.sample_count or '?'} samples, {d.asv_count or '?'} ASVs)"
        )
        options.append({"label": label, "value": d.id})
    return options


def get_layout():
    pipeline_options = _get_pipeline_options()

    return dbc.Container(
        [
            html.H3("Combine Datasets", className="mb-2"),
            html.P(
                "Merge multiple BIOM files into one. Same-region datasets are "
                "merged by sequence identity; different-region datasets are "
                "merged by taxonomy.",
                className="text-muted mb-4",
            ),
            # ── Upload BIOMs ─────────────────────────────────────────────
            dbc.Label("Upload BIOM Files", className="fw-bold"),
            dcc.Upload(
                id="cb-upload-bioms",
                children=html.Div(
                    [
                        "Drag & drop or ",
                        html.A("select .biom files", className="text-info"),
                    ],
                    className="text-center py-3",
                ),
                style={
                    "borderWidth": "2px",
                    "borderStyle": "dashed",
                    "borderRadius": "5px",
                    "borderColor": "#555",
                },
                multiple=True,
                accept=".biom,application/octet-stream,application/x-hdf5",
            ),
            html.Div(id="cb-upload-list", className="mt-2 mb-3"),
            # ── Pipeline Datasets ────────────────────────────────────────
            html.Hr(),
            dbc.Label("Pipeline Datasets", className="fw-bold"),
            html.P(
                "Check completed pipeline datasets to include in the merge.",
                className="text-muted small",
            ),
            (
                dbc.Checklist(
                    id="cb-pipeline-checklist",
                    options=pipeline_options,
                    value=[],
                    className="mb-3",
                )
                if pipeline_options
                else html.Div(
                    [
                        html.P("No pipeline datasets with BIOM files.", className="text-muted"),
                        # Hidden checklist so callbacks don't break
                        dbc.Checklist(
                            id="cb-pipeline-checklist",
                            options=[],
                            value=[],
                            style={"display": "none"},
                        ),
                    ]
                )
            ),
            # ── Selection Summary ────────────────────────────────────────
            html.Hr(),
            html.Div(id="cb-summary", className="mb-3"),
            # ── Merge Options ────────────────────────────────────────────
            html.Div(
                id="cb-merge-options",
                children=[
                    dbc.Label("Merge Strategy", className="fw-bold"),
                    dbc.RadioItems(
                        id="cb-merge-mode",
                        options=[
                            {"label": "By sequence (same region)", "value": "sequence"},
                            {"label": "By taxonomy (cross-region)", "value": "taxonomy"},
                        ],
                        value="sequence",
                        className="mb-3",
                    ),
                    html.Div(
                        id="cb-tax-level-section",
                        children=[
                            dbc.Label("Taxonomy Level", className="fw-bold"),
                            dbc.Select(
                                id="cb-tax-level",
                                options=[
                                    {"label": level, "value": level}
                                    for level in [
                                        "Kingdom", "Phylum", "Class",
                                        "Order", "Family", "Genus", "Species",
                                    ]
                                ],
                                value="Genus",
                                style={"maxWidth": "200px"},
                            ),
                        ],
                        style={"display": "none"},
                        className="mb-3",
                    ),
                    dbc.Button(
                        "Combine",
                        id="cb-btn-combine",
                        color="success",
                        size="lg",
                        disabled=True,
                        className="w-100 mb-3",
                    ),
                ],
                style={"display": "none"},
            ),
            html.Div(id="cb-result"),
            dcc.Download(id="cb-download-combined"),
            # ── Stores ───────────────────────────────────────────────────
            dcc.Store(id="cb-uploaded-files"),  # [{path, name, region, n_asvs, n_samples}]
        ],
        fluid=True,
    )


# ── Callbacks ────────────────────────────────────────────────────────────────


@dash_app.callback(
    Output("cb-uploaded-files", "data"),
    Output("cb-upload-list", "children"),
    Input("cb-upload-bioms", "contents"),
    State("cb-upload-bioms", "filename"),
    State("cb-uploaded-files", "data"),
    prevent_initial_call=True,
)
def on_upload_bioms(contents_list, filenames, existing_files):
    """Handle multi-file BIOM upload, accumulating across uploads."""
    if not contents_list or not filenames:
        return no_update, no_update

    from app.data_manager.biom_ops import detect_region_from_biom

    # Start from existing uploads
    uploaded = list(existing_files) if existing_files else []
    existing_names = {f["name"] for f in uploaded}

    # Rebuild display items for existing files
    items = []
    for f in uploaded:
        row_children = [
            html.Span(f"{f['name']}.biom", className="text-success"),
            html.Span(" — "),
            dbc.Badge(
                f.get("region") or "Unknown",
                color="primary" if f.get("region") else "secondary",
                className="ms-1 me-1",
            ),
            html.Span(
                f"{f.get('n_asvs', '?')} ASVs, "
                f"{f.get('n_samples', '?')} samples",
                className="small text-muted",
            ),
        ]
        if f.get("primers_present"):
            row_children.append(
                dbc.Badge("primers not trimmed", color="warning", className="ms-2")
            )
        items.append(html.Div(row_children, className="mb-1"))

    for contents, filename in zip(contents_list, filenames):
        if not filename.endswith(".biom"):
            items.append(
                html.Div(
                    f"Skipped {filename} (not .biom)",
                    className="text-warning small",
                )
            )
            continue

        name = filename.rsplit(".", 1)[0]
        if name in existing_names:
            items.append(
                html.Div(
                    f"Skipped {filename} (already added)",
                    className="text-warning small",
                )
            )
            continue

        try:
            _, content_string = contents.split(",", 1)
            data = base64.b64decode(content_string)

            tmp_dir = tempfile.mkdtemp()
            tmp_path = str(Path(tmp_dir) / filename)
            Path(tmp_path).write_bytes(data)

            info = detect_region_from_biom(tmp_path)
            region = info.get("region", "Unknown")

            primers_present = info.get("primers_present", False)
            uploaded.append({
                "path": tmp_path,
                "name": name,
                "region": region,
                "n_asvs": info.get("n_asvs", 0),
                "n_samples": info.get("n_samples", 0),
                "primers_present": primers_present,
            })
            existing_names.add(name)

            row_children = [
                html.Span(f"{filename}", className="text-success"),
                html.Span(" — "),
                dbc.Badge(
                    region or "Unknown",
                    color="primary" if region else "secondary",
                    className="ms-1 me-1",
                ),
                html.Span(
                    f"{info.get('n_asvs', '?')} ASVs, "
                    f"{info.get('n_samples', '?')} samples",
                    className="small text-muted",
                ),
            ]
            if primers_present:
                row_children.append(
                    dbc.Badge("primers not trimmed", color="warning", className="ms-2")
                )
            items.append(html.Div(row_children, className="mb-1"))
        except Exception as e:
            items.append(
                html.Div(
                    f"Error reading {filename}: {e}",
                    className="text-danger small",
                )
            )

    return uploaded or None, html.Div(items)


@dash_app.callback(
    Output("cb-summary", "children"),
    Output("cb-merge-options", "style"),
    Output("cb-merge-mode", "value"),
    Output("cb-btn-combine", "disabled"),
    Input("cb-uploaded-files", "data"),
    Input("cb-pipeline-checklist", "value"),
    State("cb-merge-mode", "value"),
)
def update_summary(uploaded_files, pipeline_selected, current_mode):
    """Update the selection summary and auto-detect merge strategy."""
    sources = []

    if uploaded_files:
        for f in uploaded_files:
            sources.append({
                "name": f["name"],
                "region": f.get("region"),
                "n_asvs": f.get("n_asvs", 0),
                "n_samples": f.get("n_samples", 0),
            })

    # Include pipeline datasets
    if pipeline_selected:
        from app.data_manager.biom_ops import detect_region_from_biom
        from app.db.database import get_session
        from app.db.models import Dataset

        with get_session() as db:
            for ds_id in pipeline_selected:
                ds = db.query(Dataset).filter(Dataset.id == ds_id).first()
                if not ds:
                    continue
                biom_path = DATASET_DIR / str(ds_id) / "asv_table.biom"
                if not biom_path.exists():
                    continue
                info = detect_region_from_biom(str(biom_path))
                sources.append({
                    "name": ds.name,
                    "region": info.get("region") or ds.variable_region,
                    "n_asvs": info.get("n_asvs", ds.asv_count or 0),
                    "n_samples": info.get("n_samples", ds.sample_count or 0),
                })

    if len(sources) < 2:
        return (
            html.P(
                "Upload or select at least 2 datasets to combine.",
                className="text-muted",
            ),
            {"display": "none"},
            current_mode or "sequence",
            True,
        )

    regions = set(s["region"] for s in sources if s.get("region"))
    same_region = len(regions) == 1
    total_asvs = sum(s.get("n_asvs", 0) for s in sources)
    total_samples = sum(s.get("n_samples", 0) for s in sources)

    summary_items = [
        html.Div(
            [
                html.Strong(f"{len(sources)} datasets selected"),
                html.Span(
                    f"  |  {total_asvs} total ASVs  |  {total_samples} total samples"
                ),
            ]
        ),
        html.Div(
            [
                html.Span("Regions: "),
                *[
                    dbc.Badge(r, color="primary", className="me-1")
                    for r in sorted(regions) if r
                ],
            ],
            className="mt-1",
        ),
    ]

    if same_region:
        summary_items.append(
            html.Small(
                "All datasets share the same region — sequence-based merge recommended.",
                className="text-success d-block mt-1",
            )
        )
        mode = "sequence"
    else:
        summary_items.append(
            html.Small(
                "Mixed regions detected — taxonomy-based merge recommended.",
                className="text-warning d-block mt-1",
            )
        )
        mode = "taxonomy"

    summary = dbc.Card(dbc.CardBody(summary_items), className="mb-3")
    return summary, {"display": "block"}, mode, False


@dash_app.callback(
    Output("cb-tax-level-section", "style"),
    Input("cb-merge-mode", "value"),
    prevent_initial_call=True,
)
def toggle_tax_level(mode):
    if mode == "taxonomy":
        return {"display": "block"}
    return {"display": "none"}


@dash_app.callback(
    Output("cb-result", "children"),
    Output("cb-download-combined", "data"),
    Input("cb-btn-combine", "n_clicks"),
    State("cb-uploaded-files", "data"),
    State("cb-pipeline-checklist", "value"),
    State("cb-merge-mode", "value"),
    State("cb-tax-level", "value"),
    prevent_initial_call=True,
)
def on_combine(n_clicks, uploaded_files, pipeline_selected, merge_mode, tax_level):
    """Run the combine operation."""
    if not n_clicks:
        return no_update, no_update

    # Collect all BIOM paths and names
    biom_paths = []
    source_names = []

    if uploaded_files:
        for f in uploaded_files:
            biom_paths.append(f["path"])
            source_names.append(f["name"])

    # Include pipeline datasets
    if pipeline_selected:
        from app.db.database import get_session
        from app.db.models import Dataset

        with get_session() as db:
            for ds_id in pipeline_selected:
                ds = db.query(Dataset).filter(Dataset.id == ds_id).first()
                if not ds:
                    continue
                biom_path = DATASET_DIR / str(ds_id) / "asv_table.biom"
                if biom_path.exists():
                    biom_paths.append(str(biom_path))
                    source_names.append(ds.name)

    if len(biom_paths) < 2:
        return dbc.Alert("Need at least 2 datasets.", color="warning"), no_update

    try:
        if merge_mode == "sequence":
            from app.data_manager.biom_ops import combine_biom_same_region

            result = combine_biom_same_region(biom_paths, source_names)

            stats_children = [
                html.Strong("Combined by sequence identity"),
                html.Br(),
                f"Total ASVs: {result['total_asvs']}  ->  "
                f"Unique ASVs: {result['unique_asvs']}",
                html.Br(),
                f"Total samples: {result['total_samples']}",
            ]
            if result["collisions"]:
                stats_children.extend([
                    html.Br(),
                    html.Small(
                        "Sample name collisions detected — names were prefixed with source.",
                        className="text-warning",
                    ),
                ])

            stats = dbc.Alert(stats_children, color="success")
            dl = dcc.send_bytes(
                result["biom_bytes"],
                filename="combined_sequences.biom",
            )

        else:  # taxonomy
            from app.data_manager.biom_ops import combine_biom_by_taxonomy

            result = combine_biom_by_taxonomy(biom_paths, source_names, tax_level)

            stats = dbc.Alert(
                [
                    html.Strong(f"Combined by taxonomy ({result['tax_level']} level)"),
                    html.Br(),
                    f"Total taxa: {result['total_taxa']}  |  "
                    f"Total samples: {result['total_samples']}",
                ],
                color="success",
            )
            dl = dcc.send_bytes(
                result["tsv_bytes"],
                filename=f"combined_{tax_level.lower()}.tsv",
            )

        return stats, dl

    except Exception as e:
        return dbc.Alert(f"Combine failed: {e}", color="danger"), no_update
