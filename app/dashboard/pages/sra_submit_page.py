"""
SRA Submission Helper — Generate NCBI SRA metadata spreadsheets.
"""
import dash_bootstrap_components as dbc
import dash
from dash import Input, Output, State, dcc, html, no_update, ALL

from app.dashboard.app import app as dash_app
from app.db.database import SessionLocal
from app.db.models import Upload
from app.sra.submission import (
    BIOSAMPLE_PACKAGES,
    INSTRUMENT_OPTIONS,
    generate_biosample_metadata,
    generate_sra_metadata,
    get_upload_info,
)


def get_layout():
    # Fetch uploads for dropdown
    db = SessionLocal()
    try:
        uploads = (
            db.query(Upload)
            .order_by(Upload.created_at.desc())
            .all()
        )
        upload_options = [
            {
                "label": (
                    f"Upload #{u.id} — "
                    f"{u.sequencing_type or '?'}, "
                    f"{u.variable_region or '?'}, "
                    f"{u.total_files or 0} files"
                    f"{' — ' + u.study if u.study else ''}"
                ),
                "value": u.id,
            }
            for u in uploads
        ]
    finally:
        db.close()

    return dbc.Container(
        [
            html.H3("SRA Submission Helper", className="mb-3"),
            html.P(
                "Generate NCBI SRA submission metadata spreadsheets from your uploaded files. "
                "Select an upload, review the auto-detected settings, fill in the required fields, "
                "and download the TSV file for submission to NCBI.",
                className="text-muted mb-4",
            ),
            # Upload selector
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Label("Select Upload"),
                            dcc.Dropdown(
                                id="sra-sub-select-upload",
                                options=upload_options,
                                placeholder="Choose an upload...",
                            ),
                        ],
                        md=8,
                    ),
                ],
                className="mb-3",
            ),
            # Detection info (shown after upload selection)
            html.Div(id="sra-sub-detection-info", className="mb-4"),
            # Editable fields
            html.Div(
                id="sra-sub-fields-section",
                children=[
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    dbc.Label("BioProject Accession (optional)"),
                                    dbc.Input(
                                        id="sra-sub-bioproject",
                                        placeholder="PRJNA...",
                                        type="text",
                                    ),
                                ],
                                md=4,
                            ),
                            dbc.Col(
                                [
                                    dbc.Label("Instrument Model"),
                                    dcc.Dropdown(
                                        id="sra-sub-instrument",
                                        placeholder="Select or type instrument...",
                                    ),
                                    dbc.Input(
                                        id="sra-sub-instrument-custom",
                                        placeholder="Or type custom model...",
                                        type="text",
                                        className="mt-1",
                                        size="sm",
                                    ),
                                ],
                                md=4,
                            ),
                        ],
                        className="mb-3",
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    dbc.Label("Title (applied to all samples)"),
                                    dbc.Input(
                                        id="sra-sub-title",
                                        placeholder="e.g., 16S rRNA sequencing of gut microbiome",
                                        type="text",
                                    ),
                                ],
                                md=8,
                            ),
                        ],
                        className="mb-3",
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    dbc.Label("Design Description"),
                                    dbc.Textarea(
                                        id="sra-sub-design-desc",
                                        placeholder="e.g., 16S rRNA V3-V4 amplicon sequencing using Illumina MiSeq",
                                        style={"height": "80px"},
                                    ),
                                ],
                                md=8,
                            ),
                        ],
                        className="mb-3",
                    ),
                ],
                style={"display": "none"},
            ),
            # Preview table with include/exclude
            html.Div(id="sra-sub-preview", className="mb-4"),
            dcc.Store(id="sra-sub-sample-list"),
            # Download button
            html.Div(
                [
                    dbc.Button(
                        "Download SRA Metadata TSV",
                        id="sra-sub-btn-download",
                        color="success",
                        disabled=True,
                        className="me-2",
                    ),
                ],
                className="mb-4",
            ),
            dcc.Download(id="sra-sub-download-file"),

            # ── BioSample Metadata Section ───────────────────────────────
            html.Hr(className="my-4"),
            html.H4("BioSample Metadata", className="mb-3"),
            html.P(
                "Generate a BioSample submission TSV for NCBI. "
                "Select a MIMS package matching your sample type, fill in shared values, "
                "and download the template with sample names pre-filled. "
                "All required and common optional columns are included in the output — "
                "fill in any remaining blank fields directly in the TSV before submitting to NCBI.",
                className="text-muted mb-3",
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Label("MIMS Package"),
                            dcc.Dropdown(
                                id="sra-bio-package",
                                options=[
                                    {"label": v["label"], "value": k}
                                    for k, v in BIOSAMPLE_PACKAGES.items()
                                ],
                                placeholder="Select sample type...",
                            ),
                        ],
                        md=4,
                    ),
                    dbc.Col(
                        [
                            dbc.Label("Organism"),
                            dbc.Input(
                                id="sra-bio-organism",
                                placeholder="Auto-filled from package...",
                                type="text",
                            ),
                        ],
                        md=4,
                    ),
                ],
                className="mb-3",
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Label("Collection Date"),
                            dbc.Input(
                                id="sra-bio-collection-date",
                                placeholder="e.g., 2024-06-15",
                                type="text",
                            ),
                        ],
                        md=3,
                    ),
                    dbc.Col(
                        [
                            dbc.Label("Geographic Location"),
                            dbc.Input(
                                id="sra-bio-geo-loc",
                                placeholder="e.g., South Korea: Seoul",
                                type="text",
                            ),
                        ],
                        md=3,
                    ),
                    dbc.Col(
                        [
                            dbc.Label("Latitude and Longitude"),
                            dbc.Input(
                                id="sra-bio-lat-lon",
                                placeholder="e.g., 37.57 N 126.98 E",
                                type="text",
                            ),
                        ],
                        md=3,
                    ),
                ],
                className="mb-3",
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Label("Broad-scale Environmental Context"),
                            dbc.Input(
                                id="sra-bio-env-broad",
                                placeholder="e.g., forest biome [ENVO:01000174]",
                                type="text",
                            ),
                        ],
                        md=4,
                    ),
                    dbc.Col(
                        [
                            dbc.Label("Local-scale Environmental Context"),
                            dbc.Input(
                                id="sra-bio-env-local",
                                placeholder="e.g., agricultural field [ENVO:00000114]",
                                type="text",
                            ),
                        ],
                        md=4,
                    ),
                    dbc.Col(
                        [
                            dbc.Label("Environmental Medium"),
                            dbc.Input(
                                id="sra-bio-env-medium",
                                placeholder="e.g., soil [ENVO:00001998]",
                                type="text",
                            ),
                        ],
                        md=4,
                    ),
                ],
                className="mb-3",
            ),
            # Host field (shown for host-associated/human-associated)
            html.Div(
                id="sra-bio-host-section",
                children=dbc.Row(
                    dbc.Col(
                        [
                            dbc.Label("Host"),
                            dbc.Input(
                                id="sra-bio-host",
                                placeholder="e.g., Homo sapiens",
                                type="text",
                            ),
                        ],
                        md=4,
                    ),
                    className="mb-3",
                ),
                style={"display": "none"},
            ),
            # Depth/Elevation fields (shown for soil/water)
            html.Div(
                id="sra-bio-depth-section",
                children=dbc.Row(
                    [
                        dbc.Col(
                            [
                                dbc.Label("Depth"),
                                dbc.Input(
                                    id="sra-bio-depth",
                                    placeholder="e.g., 0.1 m",
                                    type="text",
                                ),
                            ],
                            md=3,
                        ),
                        dbc.Col(
                            [
                                dbc.Label("Elevation"),
                                dbc.Input(
                                    id="sra-bio-elevation",
                                    placeholder="e.g., 150 m",
                                    type="text",
                                ),
                            ],
                            md=3,
                        ),
                    ],
                    className="mb-3",
                ),
                style={"display": "none"},
            ),
            # Preview + Download
            html.Div(id="sra-bio-preview", className="mb-3"),
            dbc.Button(
                "Download BioSample TSV",
                id="sra-bio-btn-download",
                color="success",
                disabled=True,
                className="mb-4",
            ),
            dcc.Download(id="sra-bio-download-file"),
            # Status message
            html.Div(id="sra-sub-status"),
        ],
        fluid=True,
    )


# ── Callback: on upload selection ────────────────────────────────────────────


@dash_app.callback(
    Output("sra-sub-detection-info", "children"),
    Output("sra-sub-instrument", "options"),
    Output("sra-sub-instrument", "value"),
    Output("sra-sub-design-desc", "value"),
    Output("sra-sub-fields-section", "style"),
    Output("sra-sub-preview", "children"),
    Output("sra-sub-btn-download", "disabled"),
    Output("sra-sub-sample-list", "data"),
    Input("sra-sub-select-upload", "value"),
    prevent_initial_call=True,
)
def on_select_upload(upload_id):
    if not upload_id:
        return no_update, no_update, no_update, no_update, {"display": "none"}, "", True, no_update

    info = get_upload_info(upload_id)
    if not info:
        return (
            dbc.Alert("Upload not found.", color="danger"),
            [], None, "", {"display": "none"}, "", True, None,
        )

    # Detection info badges
    detection_info = html.Div(
        [
            dbc.Badge(info["sequencing_type"], color="primary", className="me-2"),
            dbc.Badge(info["variable_region"], color="info", className="me-2"),
            dbc.Badge(info["platform"], color="secondary", className="me-2"),
            dbc.Badge(f"{info['n_samples']} samples", color="success", className="me-2"),
            dbc.Badge(f"{info['n_files']} files", color="light", text_color="dark", className="me-2"),
        ],
        className="mb-2",
    )

    # Instrument options based on platform
    sra_platform = info["sra_platform"]
    instrument_opts = [{"label": m, "value": m} for m in INSTRUMENT_OPTIONS.get(sra_platform, [])]

    # Design description hint
    region = info["variable_region"]
    design_hint = f"16S rRNA {region} amplicon sequencing" if region != "unknown" else "16S rRNA amplicon sequencing"

    # Generate preview table
    sample_ids = []
    try:
        df = generate_sra_metadata(upload_id)
        sample_ids = df["library_ID"].tolist() if "library_ID" in df.columns else []
        preview = _make_preview_table(df)
    except Exception as e:
        preview = dbc.Alert(f"Error generating preview: {e}", color="danger")

    return (
        detection_info,
        instrument_opts,
        None,
        design_hint,
        {"display": "block"},
        preview,
        False,
        sample_ids,
    )


# ── Callback: disable dropdown when custom input is provided ─────────────────


@dash_app.callback(
    Output("sra-sub-instrument", "disabled"),
    Input("sra-sub-instrument-custom", "value"),
)
def toggle_instrument_dropdown(custom_value):
    return bool((custom_value or "").strip())


# ── Callback: download TSV ───────────────────────────────────────────────────


@dash_app.callback(
    Output("sra-sub-download-file", "data"),
    Input("sra-sub-btn-download", "n_clicks"),
    State("sra-sub-select-upload", "value"),
    State("sra-sub-bioproject", "value"),
    State("sra-sub-instrument", "value"),
    State("sra-sub-instrument-custom", "value"),
    State("sra-sub-title", "value"),
    State("sra-sub-design-desc", "value"),
    State({"type": "sra-sub-include", "index": dash.ALL}, "value"),
    State("sra-sub-sample-list", "data"),
    prevent_initial_call=True,
)
def on_download(n_clicks, upload_id, bioproject, instrument, instrument_custom,
                title, design_desc, include_values, sample_list):
    if not n_clicks or not upload_id:
        return no_update

    df = generate_sra_metadata(upload_id)

    # Filter to only included samples
    if include_values and sample_list and "library_ID" in df.columns:
        included = set()
        for vals in include_values:
            if vals:
                included.update(vals)
        if included:
            df = df[df["library_ID"].isin(included)]
        if df.empty:
            return no_update

    # Apply user-provided values
    if bioproject:
        df["biosample_accession"] = bioproject
    # Custom input takes priority over dropdown
    final_instrument = (instrument_custom or "").strip() or instrument
    if final_instrument:
        df["instrument_model"] = final_instrument
    if title:
        df["title"] = title
    if design_desc:
        df["design_description"] = design_desc

    tsv_str = df.to_csv(sep="\t", index=False)
    return dcc.send_string(tsv_str, filename=f"SRA_metadata_upload{upload_id}.tsv")


# ── Helper ───────────────────────────────────────────────────────────────────


def _make_preview_table(df):
    """Create a Dash Bootstrap table with include/exclude switches."""
    if df.empty:
        return html.P("No samples found.", className="text-muted")

    # Show a subset of columns for preview
    preview_cols = [
        "library_ID", "library_layout", "platform",
        "instrument_model", "filename", "filename2",
    ]
    cols = [c for c in preview_cols if c in df.columns]
    subset = df[cols]

    header = html.Thead(html.Tr(
        [html.Th("Include")] + [html.Th(c) for c in cols]
    ))

    sample_ids = df["library_ID"].tolist() if "library_ID" in df.columns else []

    body = html.Tbody(
        [
            html.Tr(
                [
                    html.Td(
                        dbc.Checklist(
                            options=[{"label": "", "value": sample_ids[i]}],
                            value=[sample_ids[i]],
                            id={"type": "sra-sub-include", "index": i},
                            switch=True,
                        ),
                    ),
                ]
                + [html.Td(str(row[c])) for c in cols]
            )
            for i, (_, row) in enumerate(subset.iterrows())
        ]
    )

    return html.Div(
        [
            html.H6(f"Preview ({len(df)} samples)", className="mb-2"),
            dbc.Table(
                [header, body],
                bordered=True,
                striped=True,
                hover=True,
                size="sm",
                className="mb-0",
            ),
        ]
    )


# ── BioSample callbacks ─────────────────────────────────────────────────────


@dash_app.callback(
    Output("sra-bio-organism", "value"),
    Output("sra-bio-host-section", "style"),
    Output("sra-bio-depth-section", "style"),
    Output("sra-bio-btn-download", "disabled"),
    Input("sra-bio-package", "value"),
    State("sra-sub-select-upload", "value"),
)
def on_package_select(package, upload_id):
    if not package:
        return "", {"display": "none"}, {"display": "none"}, True

    pkg = BIOSAMPLE_PACKAGES.get(package, {})
    organism = pkg.get("organism_default", "metagenome")

    show_host = "host" in pkg.get("required", [])
    show_depth = "depth" in pkg.get("required", [])

    host_style = {"display": "block"} if show_host else {"display": "none"}
    depth_style = {"display": "block"} if show_depth else {"display": "none"}

    return organism, host_style, depth_style, not bool(upload_id)


@dash_app.callback(
    Output("sra-bio-download-file", "data"),
    Input("sra-bio-btn-download", "n_clicks"),
    State("sra-sub-select-upload", "value"),
    State("sra-bio-package", "value"),
    State("sra-bio-organism", "value"),
    State("sra-bio-collection-date", "value"),
    State("sra-bio-geo-loc", "value"),
    State("sra-bio-lat-lon", "value"),
    State("sra-bio-env-broad", "value"),
    State("sra-bio-env-local", "value"),
    State("sra-bio-env-medium", "value"),
    State("sra-bio-host", "value"),
    State("sra-bio-depth", "value"),
    State("sra-bio-elevation", "value"),
    prevent_initial_call=True,
)
def on_bio_download(n_clicks, upload_id, package,
                    organism, collection_date, geo_loc, lat_lon,
                    env_broad, env_local, env_medium,
                    host, depth, elevation):
    if not n_clicks or not upload_id or not package:
        return no_update

    shared = {}
    if organism:
        shared["organism"] = organism.strip()
    if collection_date:
        shared["collection_date"] = collection_date.strip()
    if geo_loc:
        shared["geo_loc_name"] = geo_loc.strip()
    if lat_lon:
        shared["lat_lon"] = lat_lon.strip()
    if env_broad:
        shared["env_broad_scale"] = env_broad.strip()
    if env_local:
        shared["env_local_scale"] = env_local.strip()
    if env_medium:
        shared["env_medium"] = env_medium.strip()
    if host:
        shared["host"] = host.strip()
    if depth:
        shared["depth"] = depth.strip()
    if elevation:
        shared["elevation"] = elevation.strip()

    df = generate_biosample_metadata(upload_id, package, shared)
    if df.empty:
        return no_update

    tsv_str = df.to_csv(sep="\t", index=False)
    pkg_short = package.split(".")[-1]
    return dcc.send_string(tsv_str, filename=f"biosample_{pkg_short}_upload{upload_id}.tsv")
