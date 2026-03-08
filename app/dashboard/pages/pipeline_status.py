"""
MicrobiomeDash — Pipeline page: select samples, launch DADA2 pipeline.
"""
import os
import re
from collections import defaultdict
from pathlib import Path

import dash_bootstrap_components as dbc
from dash import ALL, Input, Output, State, ctx, dcc, html, no_update

from app.dashboard.app import app as dash_app

MAX_CPUS = os.cpu_count() or 1

STEPS = ["fastqc", "cutadapt", "dada2", "taxonomy"]
STEPS_LONGREAD = ["fastqc", "dada2_longread", "taxonomy"]
STEP_LABELS = {
    "fastqc": "FastQC",
    "cutadapt": "Cutadapt",
    "dada2": "DADA2",
    "dada2_longread": "DADA2 (Long-Read)",
    "taxonomy": "Taxonomy",
}


# ── Page Layout ───────────────────────────────────────────────────────────────

layout = dbc.Container(
    [
        html.H3("DADA2", className="mb-4"),
        html.P(
            "Select samples to run through the DADA2 pipeline. "
            "Use column filters to narrow the list, then check individual rows "
            "or use Select All.",
            className="text-muted mb-3",
        ),
        # ── Hidden stores for sample selection ─────────────────────────────
        dcc.Store(id="ps-checked-samples", data=[]),
        dcc.Store(id="ps-sort-field", data="sample_name"),
        dcc.Store(id="ps-sort-asc", data=True),
        # ── Sample Selection Card ─────────────────────────────────────────
        dbc.Card(
            [
                dbc.CardHeader(
                    dbc.Row(
                        [
                            dbc.Col(
                                html.Span("Select Samples", className="fw-bold"),
                                width="auto",
                            ),
                            dbc.Col(
                                dbc.ButtonGroup(
                                    [
                                        dbc.Button("Sample", id="ps-sort-sample", color="primary", size="sm", outline=True),
                                        dbc.Button("Type", id="ps-sort-type", color="primary", size="sm", outline=True),
                                        dbc.Button("Reads", id="ps-sort-reads", color="primary", size="sm", outline=True),
                                        dbc.Button("Region", id="ps-sort-region", color="primary", size="sm", outline=True),
                                        dbc.Button("Read Len", id="ps-sort-readlen", color="primary", size="sm", outline=True),
                                        dbc.Button("Study", id="ps-sort-study", color="primary", size="sm", outline=True),
                                        dbc.Button("Source", id="ps-sort-source", color="primary", size="sm", outline=True),
                                        dbc.Button("Metadata", id="ps-sort-metadata", color="primary", size="sm", outline=True),
                                        dbc.Button("Date", id="ps-sort-date", color="primary", size="sm", outline=True),
                                    ],
                                    size="sm",
                                ),
                                width="auto",
                                className="ms-auto",
                            ),
                        ],
                        align="center",
                        className="g-0",
                    ),
                ),
                dbc.CardBody(
                    [
                        dbc.Row(
                            [
                                dbc.Col(width={"size": 1, "order": 0}, style={"maxWidth": "30px"}),  # checkbox col
                                dbc.Col(dbc.Input(id="ps-f-sample", placeholder="Sample", debounce=True, size="sm"), width=2),
                                dbc.Col(dbc.Input(id="ps-f-type", placeholder="Type", debounce=True, size="sm"), width=1),
                                dbc.Col(width=1),  # Reads (no filter)
                                dbc.Col(dbc.Input(id="ps-f-region", placeholder="Region", debounce=True, size="sm"), width=1),
                                dbc.Col(width=1),  # Primers (no filter)
                                dbc.Col(width=1),  # Avg Read Len (no filter)
                                dbc.Col(dbc.Input(id="ps-f-study", placeholder="Study", debounce=True, size="sm"), width=1),
                                dbc.Col(dbc.Input(id="ps-f-source", placeholder="Source", debounce=True, size="sm"), width=1),
                                dbc.Col(dbc.Input(id="ps-f-metadata", placeholder="Metadata", debounce=True, size="sm"), width=2),
                                dbc.Col(width=1),  # Date (no filter)
                            ],
                            className="mb-2 g-1",
                        ),
                        html.Div(id="ps-sample-table", children="No samples found."),
                    ],
                ),
            ],
            className="mb-3",
        ),
        # ── Selection Summary ─────────────────────────────────────────────
        html.Div(id="sample-selection-summary", className="mb-3"),
        # ── Controls Row ──────────────────────────────────────────────────
        dbc.Row(
            [
                dbc.Col(
                    [
                        dbc.Label("Run Name", className="fw-bold"),
                        dbc.Input(
                            id="input-run-name",
                            type="text",
                            placeholder="e.g. Stool_V3V4_run1",
                            debounce=True,
                        ),
                    ],
                    width=3,
                ),
                dbc.Col(
                    [
                        dbc.Label("CPU Threads", className="fw-bold"),
                        dbc.InputGroup(
                            [
                                dbc.Input(
                                    id="input-threads",
                                    type="text",
                                    value="default",
                                    placeholder="default",
                                    style={"maxWidth": "120px"},
                                ),
                                dbc.InputGroupText(
                                    f"default = samples x 2 (max {MAX_CPUS} cores)",
                                    className="text-muted small",
                                ),
                            ],
                        ),
                    ],
                    width=5,
                ),
                dbc.Col(
                    dbc.Button(
                        "Start Pipeline",
                        id="btn-start-pipeline",
                        color="success",
                        disabled=True,
                        className="mt-4",
                    ),
                    width="auto",
                    className="ms-auto",
                ),
            ],
            align="end",
            className="mb-3",
        ),
        # ── DADA2 Parameters (optional) ────────────────────────────────
        dbc.Card(
            [
                dbc.CardHeader(
                    html.Span(
                        "DADA2 Parameters (optional — auto-detected by default)",
                        className="fw-bold small",
                    ),
                ),
                dbc.CardBody(
                    [
                        html.P(
                            "Leave blank for auto-detection. Truncation positions "
                            "are determined from quality profiles; trim-left removes "
                            "leading bases (e.g. N-containing positions).",
                            className="text-muted small mb-2",
                        ),
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        dbc.Label("Trunc Fwd", className="small"),
                                        dbc.Input(
                                            id="input-trunc-f",
                                            type="number",
                                            placeholder="auto",
                                            min=0,
                                        ),
                                    ],
                                    width=2,
                                ),
                                dbc.Col(
                                    [
                                        dbc.Label("Trunc Rev", className="small"),
                                        dbc.Input(
                                            id="input-trunc-r",
                                            type="number",
                                            placeholder="auto",
                                            min=0,
                                        ),
                                    ],
                                    width=2,
                                ),
                                dbc.Col(
                                    [
                                        dbc.Label("Trim Left Fwd", className="small"),
                                        dbc.Input(
                                            id="input-trimleft-f",
                                            type="number",
                                            placeholder="auto",
                                            min=0,
                                        ),
                                    ],
                                    width=2,
                                ),
                                dbc.Col(
                                    [
                                        dbc.Label("Trim Left Rev", className="small"),
                                        dbc.Input(
                                            id="input-trimleft-r",
                                            type="number",
                                            placeholder="auto",
                                            min=0,
                                        ),
                                    ],
                                    width=2,
                                ),
                            ],
                            className="mb-3",
                        ),
                        html.Hr(className="my-2"),
                        html.P(
                            "Provide custom primer sequences to override auto-detection. "
                            "When set, Cutadapt will always run with these primers "
                            "regardless of the auto-detected region.",
                            className="text-muted small mb-2",
                        ),
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        dbc.Label("Forward Primer (5'→3')", className="small"),
                                        dbc.Input(
                                            id="input-fwd-primer",
                                            type="text",
                                            placeholder="e.g. CCTACGGGNGGCWGCAG",
                                            debounce=True,
                                            className="font-monospace",
                                        ),
                                    ],
                                    width=6,
                                ),
                                dbc.Col(
                                    [
                                        dbc.Label("Reverse Primer (5'→3')", className="small"),
                                        dbc.Input(
                                            id="input-rev-primer",
                                            type="text",
                                            placeholder="e.g. GACTACHVGGGTATCTAATCC",
                                            debounce=True,
                                            className="font-monospace",
                                        ),
                                    ],
                                    width=6,
                                ),
                            ],
                        ),
                    ],
                ),
            ],
            className="mb-3",
        ),
        html.Div(id="pipeline-launch-error"),
        # Detection summary (shown after scan+launch)
        html.Div(id="pipeline-detection-summary"),
        # ── Progress Section ─────────────────────────────────────────────
        html.Div(
            id="pipeline-progress-section",
            children=[
                html.Div(
                    [
                        html.H5("Progress", className="d-inline"),
                        dbc.Button(
                            "Cancel Pipeline",
                            id="btn-cancel-pipeline",
                            color="danger",
                            size="sm",
                            className="ms-3",
                        ),
                    ],
                    className="d-flex align-items-center mb-2",
                ),
                html.Div(id="pipeline-status-badge", className="mb-3"),
                html.Div(
                    dbc.Progress(
                        id="pipeline-progress-bar",
                        value=0,
                        striped=True,
                        animated=True,
                        className="mb-3",
                    ),
                    id="pipeline-progress-bar-wrapper",
                ),
                html.Div(id="pipeline-steps-display", className="mb-3"),
                html.H6("Pipeline Log"),
                html.Pre(
                    id="pipeline-log-viewer",
                    style={
                        "maxHeight": "400px",
                        "overflowY": "auto",
                        "backgroundColor": "#1a1a2e",
                        "color": "#e0e0e0",
                        "padding": "10px",
                        "borderRadius": "5px",
                        "fontSize": "0.8rem",
                        "whiteSpace": "pre-wrap",
                    },
                ),
                # ── Download Section (visible after completion) ───────
                html.Div(
                    id="pipeline-download-section",
                    children=[
                        dbc.Button(
                            "Download ASV Table (.biom)",
                            id="btn-download-asv",
                            color="success",
                            size="sm",
                            className="mt-3",
                        ),
                        dcc.Download(id="download-asv-table"),
                    ],
                    style={"display": "none"},
                ),
            ],
            style={"display": "none"},
        ),
        html.Hr(className="mt-4"),
        # ── Pipeline History ─────────────────────────────────────────────
        html.H5("Pipeline History"),
        html.Div(
            id="pipeline-history-table",
            style={"maxHeight": "300px", "overflowY": "auto"},
        ),
        dcc.Download(id="download-history-biom"),
        dcc.Download(id="download-history-qc"),
        dcc.Download(id="download-history-log"),
        # Stores and intervals
        dcc.Store(id="pipeline-active-dataset", storage_type="session"),
        dcc.Store(id="file-ids-map-store"),
        dcc.Interval(
            id="pipeline-poll", interval=3000, disabled=True, n_intervals=0
        ),
    ],
    fluid=True,
)


# ── Helpers ───────────────────────────────────────────────────────────────────


def _build_sample_data(sort_by="sample_name", filters=None, ascending=True, checked_samples=None):
    """Query all FastqFile + Upload + UploadMetadata records, group by sample.

    Returns (html_table_div, file_ids_map) where file_ids_map maps
    sample_name → list of FastqFile IDs.
    """
    if filters is None:
        filters = {}
    if checked_samples is None:
        checked_samples = set()

    from app.db.database import get_session
    from app.db.models import FastqFile, Upload, UploadMetadata

    with get_session() as db:
        results = (
            db.query(FastqFile, Upload.created_at, Upload.variable_region, Upload.primers_detected, Upload.study)
            .join(Upload, FastqFile.upload_id == Upload.id)
            .all()
        )
        if not results:
            return "No samples found.", {}

        # Load metadata keys per sample + extract "source" and "study" values
        meta_rows = db.query(UploadMetadata).all()

    meta_by_sample = {}
    source_by_sample = {}
    study_by_sample = {}
    for m in meta_rows:
        if m.sample_name not in meta_by_sample:
            meta_by_sample[m.sample_name] = []
        meta_by_sample[m.sample_name].append(m.key)
        if m.key == "source":
            source_by_sample[m.sample_name] = m.value or ""
        if m.key == "study":
            study_by_sample[m.sample_name] = m.value or ""

    sample_map = defaultdict(lambda: {
        "files": [], "region": None, "date": None, "primers_detected": None,
        "study": None,
    })
    for f, upload_date, region, primers_detected, study in results:
        s = sample_map[f.sample_name]
        s["files"].append(f)
        if region:
            s["region"] = region
        if upload_date:
            s["date"] = upload_date
        if primers_detected is not None:
            s["primers_detected"] = primers_detected
        if study:
            s["study"] = study

    # Build row dicts
    sample_rows = []
    file_ids_map = {}
    for sample_name, info in sample_map.items():
        files = info["files"]
        file_ids_map[sample_name] = [f.id for f in files]

        directions = set(f.read_direction for f in files)
        if "R1" in directions and "R2" in directions:
            seq_type = "PE"
        elif directions == {"single"}:
            seq_type = "SE"
        else:
            seq_type = ", ".join(sorted(d for d in directions if d))

        total_reads = sum(f.read_count or 0 for f in files)
        avg_lengths = [f.avg_read_length for f in files if f.avg_read_length]
        avg_len = round(sum(avg_lengths) / len(avg_lengths)) if avg_lengths else None

        meta_keys = meta_by_sample.get(sample_name, [])
        meta_str = ", ".join(sorted(meta_keys)) if meta_keys else "—"
        source_val = source_by_sample.get(sample_name, "")

        study_val = info["study"] or study_by_sample.get(sample_name, "")

        sample_rows.append({
            "sample_name": sample_name,
            "direction": seq_type,
            "total_reads": total_reads,
            "region": info["region"],
            "avg_len": avg_len,
            "study": study_val,
            "source": source_val,
            "metadata": meta_str,
            "date": info["date"],
            "primers_detected": info["primers_detected"],
        })

    # Per-column filters (AND logic, case-insensitive substring)
    def _matches(r):
        if filters.get("sample") and filters["sample"] not in r["sample_name"].lower():
            return False
        if filters.get("type") and filters["type"] not in r["direction"].lower():
            return False
        if filters.get("region") and filters["region"] not in (r["region"] or "").lower():
            return False
        if filters.get("study") and filters["study"] not in r["study"].lower():
            return False
        if filters.get("source") and filters["source"] not in r["source"].lower():
            return False
        if filters.get("metadata") and filters["metadata"] not in r["metadata"].lower():
            return False
        return True

    sample_rows = [r for r in sample_rows if _matches(r)]

    if not sample_rows:
        return html.P("No matching samples.", className="text-muted"), file_ids_map

    # Sort
    reverse = not ascending
    sort_keys = {
        "sample_name": lambda r: r["sample_name"].lower(),
        "direction": lambda r: r["direction"],
        "total_reads": lambda r: r["total_reads"],
        "region": lambda r: r["region"] or "",
        "avg_read_length": lambda r: r["avg_len"] or 0,
        "study": lambda r: r["study"].lower(),
        "source": lambda r: r["source"].lower(),
        "metadata": lambda r: r["metadata"],
        "created_at": lambda r: r["date"] or "",
    }
    numeric_fields = {"total_reads", "avg_read_length", "created_at"}
    key_fn = sort_keys.get(sort_by, sort_keys["sample_name"])
    effective_reverse = (not reverse) if sort_by in numeric_fields else reverse
    sample_rows.sort(key=key_fn, reverse=effective_reverse)

    # Build HTML rows
    rows = []
    for r in sample_rows:
        date_str = r["date"].strftime("%Y-%m-%d") if r["date"] else ""
        reads_str = f"{r['total_reads']:,}" if r["total_reads"] else "—"
        is_checked = r["sample_name"] in checked_samples
        pd_val = r["primers_detected"]
        if pd_val is True:
            primers_badge = dbc.Badge("Yes", color="success")
        elif pd_val is False:
            primers_badge = dbc.Badge("No", color="info")
        else:
            primers_badge = dbc.Badge("?", color="secondary")
        rows.append(
            html.Tr(
                [
                    html.Td(
                        dbc.Checkbox(
                            id={"type": "ps-check", "index": r["sample_name"]},
                            value=is_checked,
                        ),
                        style={"width": "30px"},
                    ),
                    html.Td(r["sample_name"]),
                    html.Td(
                        dbc.Badge(
                            r["direction"],
                            color="primary" if r["direction"] == "PE" else "info",
                        )
                    ),
                    html.Td(reads_str),
                    html.Td(r["region"] or "—", className="small"),
                    html.Td(primers_badge),
                    html.Td(f"{r['avg_len']} bp" if r["avg_len"] else "—", className="small"),
                    html.Td(r["study"] or "—", className="small"),
                    html.Td(r["source"] or "—", className="small"),
                    html.Td(r["metadata"], className="small"),
                    html.Td(date_str, className="text-muted small"),
                ]
            )
        )

    total_reads = sum(r["total_reads"] for r in sample_rows)

    # Arrow indicator for sorted column
    arrow = " ^" if ascending else " v"
    def _th(label, field):
        suffix = arrow if sort_by == field else ""
        return html.Th(f"{label}{suffix}")

    table_div = html.Div(
        [
            html.P(
                f"{len(sample_rows)} samples, {total_reads:,} reads total",
                className="text-muted mb-2",
            ),
            dbc.Table(
                [
                    html.Thead(
                        html.Tr(
                            [
                                html.Th(
                                    dbc.Checkbox(id="ps-check-all", value=False),
                                    style={"width": "30px"},
                                ),
                                _th("Sample", "sample_name"),
                                _th("Type", "direction"),
                                _th("# Reads", "total_reads"),
                                _th("Region", "region"),
                                html.Th("Primers"),
                                _th("Avg Read Len", "avg_read_length"),
                                _th("Study", "study"),
                                _th("Source", "source"),
                                _th("Metadata", "metadata"),
                                _th("Date", "created_at"),
                            ]
                        )
                    ),
                    html.Tbody(rows),
                ],
                bordered=True,
                color="dark",
                hover=True,
                size="sm",
            ),
        ]
    )

    return table_div, file_ids_map


def _download_filename(ds, prefix: str, ext: str) -> str:
    """Build a download filename from dataset name, variable region, and date."""
    if not ds:
        return f"{prefix}.{ext}"
    name = re.sub(r"[^\w\-.]", "_", ds.name or prefix).strip("_") or prefix
    region = ds.variable_region or ""
    date = ds.created_at.strftime("%Y-%m-%d") if ds.created_at else ""
    parts = [p for p in (name, region, date) if p]
    return f"{'_'.join(parts)}.{ext}"


def _dir_size(path: Path) -> int:
    """Return total size in bytes of all files under *path*."""
    if not path.is_dir():
        return 0
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())


def _format_size(nbytes: int) -> str:
    """Human-readable file size."""
    if nbytes == 0:
        return "—"
    for unit in ("B", "KB", "MB", "GB"):
        if nbytes < 1024:
            return f"{nbytes:.1f} {unit}" if unit != "B" else f"{nbytes} B"
        nbytes /= 1024
    return f"{nbytes:.1f} TB"


def _build_history_table():
    """Build the pipeline history table."""
    from app.db.database import get_session
    from app.db.models import Dataset

    with get_session() as db:
        datasets = (
            db.query(Dataset)
            .filter(Dataset.source_type == "pipeline")
            .order_by(Dataset.created_at.desc())
            .all()
        )

    if not datasets:
        return html.P("No pipeline runs yet.", className="text-muted")

    badge_color = {
        "pending": "secondary",
        "processing": "primary",
        "complete": "success",
        "failed": "danger",
        "cancelled": "warning",
    }

    from app.config import DATASET_DIR

    rows = []
    for d in datasets:
        ds_dir = DATASET_DIR / str(d.id)
        disk_bytes = _dir_size(ds_dir)
        disk_label = _format_size(disk_bytes)
        biom_exists = (
            d.status == "complete" and (ds_dir / "asv_table.biom").exists()
        )
        qc_exists = (
            d.status in ("complete", "failed")
            and (ds_dir / "qc_report.pdf").exists()
        )
        biom_cell = html.Td(
            dbc.Button(
                ".biom",
                id={"type": "btn-history-dl", "index": d.id},
                color="success",
                size="sm",
                outline=True,
            )
            if biom_exists else ""
        )
        qc_cell = html.Td(
            dbc.Button(
                "QC",
                id={"type": "btn-history-qc", "index": d.id},
                color="info",
                size="sm",
                outline=True,
            )
            if qc_exists else ""
        )

        is_active = d.status in ("pending", "processing")
        del_cell = html.Td(
            dbc.Button(
                "x",
                id={"type": "btn-history-del", "index": d.id},
                color="danger",
                size="sm",
                outline=True,
                disabled=is_active,
                title="Cannot delete a running pipeline" if is_active else "Delete",
            )
        )

        log_exists = (
            d.status in ("failed", "cancelled")
            and (ds_dir / "pipeline.log").exists()
        )
        if log_exists:
            status_cell = html.Td(
                dbc.Button(
                    d.status.upper(),
                    id={"type": "btn-history-log", "index": d.id},
                    color=badge_color.get(d.status, "secondary"),
                    size="sm",
                    title="Download pipeline log",
                )
            )
        else:
            status_cell = html.Td(
                dbc.Badge(
                    d.status.upper(),
                    color=badge_color.get(d.status, "secondary"),
                )
            )

        rows.append(
            html.Tr(
                [
                    html.Td(d.name),
                    status_cell,
                    html.Td(
                        f"{d.sequencing_type or ''} ({d.platform.title()})"
                        if d.platform and d.platform != "illumina"
                        else (d.sequencing_type or "")
                    ),
                    html.Td(d.variable_region or ""),
                    html.Td(d.sample_count or "—"),
                    html.Td(d.asv_count or "—"),
                    html.Td(
                        d.created_at.strftime("%Y-%m-%d %H:%M")
                        if d.created_at
                        else ""
                    ),
                    html.Td(disk_label, className="text-muted small"),
                    biom_cell,
                    qc_cell,
                    del_cell,
                ]
            )
        )

    return dbc.Table(
        [
            html.Thead(
                html.Tr(
                    [
                        html.Th("Name"),
                        html.Th("Status"),
                        html.Th("Type"),
                        html.Th("Region"),
                        html.Th("Samples"),
                        html.Th("ASVs"),
                        html.Th("Created"),
                        html.Th("Disk"),
                        html.Th("ASV Table"),
                        html.Th("QC Report"),
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


# ── Callbacks ─────────────────────────────────────────────────────────────────


# Sort button definitions (id → sort field)
_PS_SORT_BUTTONS = [
    ("ps-sort-sample", "sample_name"),
    ("ps-sort-type", "direction"),
    ("ps-sort-reads", "total_reads"),
    ("ps-sort-region", "region"),
    ("ps-sort-readlen", "avg_read_length"),
    ("ps-sort-study", "study"),
    ("ps-sort-source", "source"),
    ("ps-sort-metadata", "metadata"),
    ("ps-sort-date", "created_at"),
]


@dash_app.callback(
    Output("ps-sort-field", "data"),
    Output("ps-sort-asc", "data"),
    *[Output(btn_id, "outline") for btn_id, _ in _PS_SORT_BUTTONS],
    *[Input(btn_id, "n_clicks") for btn_id, _ in _PS_SORT_BUTTONS],
    State("ps-sort-field", "data"),
    State("ps-sort-asc", "data"),
    prevent_initial_call=True,
)
def on_ps_sort_click(*args):
    """Update sort field and toggle direction. Highlight active button."""
    current_field = args[-2]
    current_asc = args[-1]
    triggered = ctx.triggered_id
    sort_map = {btn_id: field for btn_id, field in _PS_SORT_BUTTONS}
    new_field = sort_map.get(triggered, "sample_name")

    if new_field == current_field:
        new_asc = not current_asc
    else:
        new_asc = True

    outlines = [new_field != field for _, field in _PS_SORT_BUTTONS]
    return new_field, new_asc, *outlines


@dash_app.callback(
    Output("ps-sample-table", "children"),
    Output("file-ids-map-store", "data"),
    Input("ps-sort-field", "data"),
    Input("ps-sort-asc", "data"),
    Input("ps-f-sample", "value"),
    Input("ps-f-type", "value"),
    Input("ps-f-region", "value"),
    Input("ps-f-study", "value"),
    Input("ps-f-source", "value"),
    Input("ps-f-metadata", "value"),
    Input("url", "pathname"),
    State("ps-checked-samples", "data"),
)
def refresh_ps_table(sort_by, sort_asc, f_sample, f_type, f_region, f_study, f_source,
                     f_metadata, pathname, checked):
    """Rebuild the sample table whenever sort, filter, or page changes."""
    if pathname != "/pipeline":
        return no_update, no_update
    filters = {
        "sample": (f_sample or "").strip().lower(),
        "type": (f_type or "").strip().lower(),
        "region": (f_region or "").strip().lower(),
        "study": (f_study or "").strip().lower(),
        "source": (f_source or "").strip().lower(),
        "metadata": (f_metadata or "").strip().lower(),
    }
    checked_set = set(checked) if checked else set()
    table_div, file_ids_map = _build_sample_data(
        sort_by or "sample_name",
        filters,
        sort_asc if sort_asc is not None else True,
        checked_samples=checked_set,
    )
    return table_div, file_ids_map


@dash_app.callback(
    Output("ps-checked-samples", "data"),
    Input({"type": "ps-check", "index": ALL}, "value"),
    State({"type": "ps-check", "index": ALL}, "id"),
    prevent_initial_call=True,
)
def on_ps_row_check(values, ids):
    """Sync checkbox states to checked-samples store."""
    checked = [
        id_dict["index"] for id_dict, val in zip(ids, values) if val
    ]
    return checked


@dash_app.callback(
    Output({"type": "ps-check", "index": ALL}, "value"),
    Input("ps-check-all", "value"),
    State({"type": "ps-check", "index": ALL}, "id"),
    prevent_initial_call=True,
)
def on_ps_check_all(select_all, ids):
    """Toggle all visible row checkboxes."""
    return [bool(select_all)] * len(ids)


@dash_app.callback(
    Output("sample-selection-summary", "children"),
    Output("btn-start-pipeline", "disabled"),
    Input("ps-checked-samples", "data"),
    Input("input-run-name", "value"),
    State("file-ids-map-store", "data"),
)
def update_selection_summary(checked_samples, run_name, file_ids_map):
    """Show 'N samples selected' + type/region badges + mixed warnings."""
    if not checked_samples or not file_ids_map:
        return "", True

    # Re-query DB to get type/region info for checked samples
    from app.db.database import get_session
    from app.db.models import FastqFile, Upload

    with get_session() as db:
        results = (
            db.query(FastqFile, Upload.variable_region)
            .join(Upload, FastqFile.upload_id == Upload.id)
            .filter(FastqFile.sample_name.in_(checked_samples))
            .all()
        )

    sample_info = defaultdict(lambda: {"directions": set(), "region": None})
    for f, region in results:
        s = sample_info[f.sample_name]
        s["directions"].add(f.read_direction)
        if region:
            s["region"] = region

    n = len(sample_info)
    if n == 0:
        return "", True

    types = set()
    regions = set()
    for info in sample_info.values():
        dirs = info["directions"]
        if "R1" in dirs and "R2" in dirs:
            types.add("PE")
        elif dirs == {"single"}:
            types.add("SE")
        else:
            types.add(", ".join(sorted(dirs)))
        if info["region"]:
            regions.add(info["region"])

    badges = []
    for t in sorted(types):
        badges.append(dbc.Badge(
            t, color="success" if t == "PE" else "info", className="me-1",
        ))
    for r in sorted(regions):
        badges.append(dbc.Badge(r, color="primary", className="me-1"))

    warnings = []
    if len(types) > 1:
        warnings.append(dbc.Badge(
            "Mixed SE/PE", color="warning", className="me-1",
        ))
    if len(regions) > 1:
        warnings.append(dbc.Badge(
            "Mixed regions", color="warning", className="me-1",
        ))

    summary = html.Div([
        html.Span(f"{n} sample{'s' if n != 1 else ''} selected  ", className="me-2"),
        *badges,
        *warnings,
    ])

    # Enable button only when samples are selected AND run name is provided
    enabled = bool(run_name and run_name.strip())
    return summary, not enabled


@dash_app.callback(
    Output("pipeline-active-dataset", "data"),
    Output("pipeline-progress-section", "style"),
    Output("pipeline-poll", "disabled"),
    Output("pipeline-launch-error", "children"),
    Output("pipeline-detection-summary", "children"),
    Output("pipeline-history-table", "children", allow_duplicate=True),
    Input("btn-start-pipeline", "n_clicks"),
    State("ps-checked-samples", "data"),
    State("file-ids-map-store", "data"),
    State("input-threads", "value"),
    State("input-run-name", "value"),
    State("input-fwd-primer", "value"),
    State("input-rev-primer", "value"),
    State("input-trunc-f", "value"),
    State("input-trunc-r", "value"),
    State("input-trimleft-f", "value"),
    State("input-trimleft-r", "value"),
    prevent_initial_call=True,
)
def on_start_pipeline(
    n_clicks, checked_samples, file_ids_map, threads_val, run_name,
    fwd_primer, rev_primer, trunc_f_val, trunc_r_val, trimleft_f_val, trimleft_r_val,
):
    """Launch the DADA2 pipeline on the selected samples."""
    if not n_clicks or not checked_samples or not file_ids_map:
        return no_update, no_update, no_update, no_update, no_update, no_update

    no_change = (no_update, no_update, no_update)

    # Parse threads
    threads_override = None
    use_default_threads = True
    if threads_val and str(threads_val).strip().lower() != "default":
        use_default_threads = False
        try:
            threads_override = max(1, int(threads_val))
        except ValueError:
            return (
                *no_change,
                dbc.Alert(
                    f"Invalid CPU threads value: '{threads_val}'. "
                    "Enter a number or 'default'.",
                    color="danger",
                ),
                no_update, no_update,
            )

    # Parse optional DADA2 truncation/trim parameters
    def _int_or_none(val):
        if val is None or val == "":
            return None
        try:
            return int(val)
        except (ValueError, TypeError):
            return None

    trunc_f = _int_or_none(trunc_f_val)
    trunc_r = _int_or_none(trunc_r_val)
    trimleft_f = _int_or_none(trimleft_f_val)
    trimleft_r = _int_or_none(trimleft_r_val)

    # ── 0. Check for duplicate run name ─────────────────────────────
    dataset_name = run_name.strip() if run_name and run_name.strip() else f"Samples_{len(checked_samples)}"

    from app.db.database import get_session
    from app.db.models import Dataset

    with get_session() as db:
        existing = (
            db.query(Dataset)
            .filter(Dataset.name == dataset_name, Dataset.source_type == "pipeline")
            .first()
        )
    if existing:
        return (
            *no_change,
            dbc.Alert(
                f"A pipeline run named '{dataset_name}' already exists. "
                "Please choose a different name.",
                color="danger",
            ),
            no_update, no_update,
        )

    # ── 1. Gather file_ids from checked samples ──────────────────────
    all_file_ids = []
    for sample_name in checked_samples:
        ids = file_ids_map.get(sample_name, [])
        all_file_ids.extend(ids)

    if not all_file_ids:
        return (
            *no_change,
            dbc.Alert("No files found for selected samples.", color="danger"),
            no_update, no_update,
        )

    n_samples = len(checked_samples)

    # ── 2. Build detection summary ────────────────────────────────────
    from app.db.database import get_session
    from app.db.models import FastqFile, Upload

    with get_session() as db:
        sel_results = (
            db.query(FastqFile, Upload.variable_region)
            .join(Upload, FastqFile.upload_id == Upload.id)
            .filter(FastqFile.sample_name.in_(checked_samples))
            .all()
        )

    sel_info = defaultdict(lambda: {"directions": set(), "region": None})
    for f, region in sel_results:
        s = sel_info[f.sample_name]
        s["directions"].add(f.read_direction)
        if region:
            s["region"] = region

    types = set()
    regions = set()
    for info in sel_info.values():
        dirs = info["directions"]
        if "R1" in dirs and "R2" in dirs:
            types.add("PE")
        elif dirs == {"single"}:
            types.add("SE")
        else:
            types.add(", ".join(sorted(dirs)))
        if info["region"]:
            regions.add(info["region"])

    # Determine display values
    if len(types) == 1:
        type_str = next(iter(types))
        type_color = "success" if type_str == "PE" else "info"
    else:
        type_str = "MIXED"
        type_color = "warning"

    if len(regions) == 1:
        region_str = next(iter(regions))
        region_color = "primary"
    elif len(regions) > 1:
        region_str = "Mixed regions"
        region_color = "warning"
    else:
        region_str = "Region unknown"
        region_color = "secondary"

    card_items = [
        html.Div([
            dbc.Badge(type_str.upper(), color=type_color, className="me-2 fs-6"),
            dbc.Badge(region_str, color=region_color, className="me-2 fs-6"),
            html.Span(
                f"{n_samples} sample(s)  |  {len(all_file_ids)} file(s)"
            ),
        ]),
        html.Small(
            "Truncation: auto-detected" if trunc_f is None and trunc_r is None
            else f"Truncation: F={trunc_f or 'auto'}, R={trunc_r or 'auto'}",
            className="text-muted d-block mt-2",
        ),
    ]
    summary_card = dbc.Card(dbc.CardBody(card_items), className="mb-3")

    # ── 3. Compute thread count ───────────────────────────────────────
    if use_default_threads:
        threads = min(n_samples * 2, max(1, MAX_CPUS - 1))
    else:
        threads = threads_override

    # ── 4. Launch pipeline ────────────────────────────────────────────
    from app.pipeline.runner import launch_pipeline

    # Normalize custom primers (strip whitespace, uppercase, None if empty)
    custom_fwd = fwd_primer.strip().upper() if fwd_primer and fwd_primer.strip() else None
    custom_rev = rev_primer.strip().upper() if rev_primer and rev_primer.strip() else None

    try:
        dataset_id = launch_pipeline(
            file_ids=all_file_ids,
            dataset_name=dataset_name,
            threads=threads,
            trunc_len_f=trunc_f,
            trunc_len_r=trunc_r,
            trim_left_f=trimleft_f,
            trim_left_r=trimleft_r,
            custom_fwd_primer=custom_fwd,
            custom_rev_primer=custom_rev,
        )
    except Exception as e:
        return (
            *no_change,
            dbc.Alert(f"Pipeline launch failed: {e}", color="danger"),
            summary_card,
            no_update,
        )

    return (
        dataset_id,
        {"display": "block"},
        False,
        "",
        summary_card,
        _build_history_table(),
    )


@dash_app.callback(
    Output("pipeline-launch-error", "children", allow_duplicate=True),
    Output("pipeline-progress-section", "style", allow_duplicate=True),
    Output("pipeline-status-badge", "children", allow_duplicate=True),
    Output("pipeline-progress-bar", "value", allow_duplicate=True),
    Output("pipeline-progress-bar", "label", allow_duplicate=True),
    Output("pipeline-progress-bar", "animated", allow_duplicate=True),
    Output("pipeline-log-viewer", "children", allow_duplicate=True),
    Output("pipeline-poll", "disabled", allow_duplicate=True),
    Output("pipeline-history-table", "children", allow_duplicate=True),
    Input("btn-cancel-pipeline", "n_clicks"),
    State("pipeline-active-dataset", "data"),
    prevent_initial_call=True,
)
def on_cancel_pipeline(n_clicks, dataset_id):
    """Cancel the running pipeline and immediately update the UI."""
    if not n_clicks or not dataset_id:
        return (no_update,) * 9

    import time

    from app.pipeline.runner import cancel_pipeline, get_pipeline_status

    ok = cancel_pipeline(dataset_id)
    if not ok:
        return (
            dbc.Alert("No running pipeline found to cancel.", color="secondary"),
            {"display": "none"},
            "", 0, "", False, "",
            True,
            _build_history_table(),
        )

    time.sleep(0.5)

    status = get_pipeline_status(dataset_id)
    badge = dbc.Badge("CANCELLED", color="warning", className="fs-6 p-2")
    log_tail = status.get("log_tail", "Pipeline cancelled by user.")

    return (
        dbc.Alert("Pipeline cancelled.", color="warning"),
        {"display": "block"},
        badge,
        0,
        "Cancelled",
        False,
        log_tail,
        True,
        _build_history_table(),
    )


@dash_app.callback(
    Output("pipeline-status-badge", "children"),
    Output("pipeline-progress-bar", "value"),
    Output("pipeline-progress-bar", "label"),
    Output("pipeline-steps-display", "children"),
    Output("pipeline-log-viewer", "children"),
    Output("pipeline-poll", "disabled", allow_duplicate=True),
    Output("pipeline-history-table", "children", allow_duplicate=True),
    Output("pipeline-download-section", "style"),
    Output("pipeline-progress-bar-wrapper", "style"),
    Input("pipeline-poll", "n_intervals"),
    State("pipeline-active-dataset", "data"),
    prevent_initial_call=True,
)
def on_poll(n_intervals, dataset_id):
    """Poll pipeline status and update the UI."""
    if not dataset_id:
        return (no_update, no_update, no_update, no_update, no_update,
                True, no_update, no_update, no_update)

    from app.pipeline.runner import get_pipeline_status

    status = get_pipeline_status(dataset_id)
    db_status = status["status"]
    pct = status["progress_pct"]

    badge_color = {
        "pending": "secondary",
        "processing": "primary",
        "complete": "success",
        "failed": "danger",
        "cancelled": "warning",
    }.get(db_status, "secondary")
    badge = dbc.Badge(
        db_status.upper(), color=badge_color, className="fs-6 p-2"
    )

    completed = set(status.get("steps_completed", []))
    current = status.get("current_step")

    # Use long-read step list if any long-read step is in progress or completed
    is_longread = bool({"dada2_longread"} & (completed | {current}))
    if not is_longread:
        from app.db.database import SessionLocal
        from app.db.models import Dataset
        _db = SessionLocal()
        try:
            _ds = _db.query(Dataset).filter(Dataset.id == dataset_id).first()
            if _ds and _ds.platform in ("pacbio", "nanopore"):
                is_longread = True
        finally:
            _db.close()
    active_steps = STEPS_LONGREAD if is_longread else STEPS

    step_items = []
    for step in active_steps:
        if step in completed and step != current:
            icon = "  "
            color = "text-success"
        elif step == current and db_status == "processing":
            icon = "  "
            color = "text-primary"
        elif step == current and db_status == "failed":
            icon = "  "
            color = "text-danger"
        else:
            icon = "  "
            color = "text-muted"
        step_items.append(
            html.Div(
                f"{icon} {STEP_LABELS.get(step, step)}",
                className=f"{color} mb-1",
            )
        )

    auto_info = None
    if status.get("auto_trunc_details"):
        auto_info = dbc.Alert(
            [
                html.Strong("Auto-detected: "),
                html.Span(status["auto_trunc_details"]),
            ],
            color="info",
            className="mt-2 py-2 small",
        )

    stop_poll = db_status in ("complete", "failed", "cancelled")
    history = _build_history_table()

    steps_display = html.Div(step_items + ([auto_info] if auto_info else []))

    download_style = (
        {"display": "block"} if db_status == "complete" else {"display": "none"}
    )

    finished = db_status in ("complete", "failed", "cancelled")
    bar_style = {"display": "none"} if finished else {"display": "block"}

    return (
        badge,
        pct,
        f"{pct}%",
        steps_display,
        status.get("log_tail", ""),
        stop_poll,
        history,
        download_style,
        bar_style,
    )


@dash_app.callback(
    Output("pipeline-history-table", "children"),
    Input("pipeline-active-dataset", "id"),  # fires on page load
)
def load_history(_):
    return _build_history_table()


@dash_app.callback(
    Output("pipeline-active-dataset", "data", allow_duplicate=True),
    Output("pipeline-progress-section", "style", allow_duplicate=True),
    Output("pipeline-poll", "disabled", allow_duplicate=True),
    Output("pipeline-status-badge", "children", allow_duplicate=True),
    Output("pipeline-progress-bar", "value", allow_duplicate=True),
    Output("pipeline-progress-bar", "label", allow_duplicate=True),
    Output("pipeline-steps-display", "children", allow_duplicate=True),
    Output("pipeline-log-viewer", "children", allow_duplicate=True),
    Output("pipeline-download-section", "style", allow_duplicate=True),
    Output("pipeline-progress-bar-wrapper", "style", allow_duplicate=True),
    Input("pipeline-active-dataset", "id"),  # fires on page load
    State("pipeline-active-dataset", "data"),
    prevent_initial_call="initial_duplicate",
)
def restore_progress_on_load(_, dataset_id):
    """Re-attach to a running pipeline on page load.

    Populates all progress fields immediately so the user sees current
    state without waiting for the first poll cycle.
    """
    from app.db.database import get_session
    from app.db.models import Dataset
    from app.pipeline.runner import get_pipeline_status

    no_all = (no_update,) * 10

    if dataset_id:
        with get_session() as db:
            ds = db.query(Dataset).filter(Dataset.id == dataset_id).first()
            if not ds or ds.status not in ("pending", "processing"):
                return (None, {"display": "none"}, True,
                        no_update, no_update, no_update, no_update,
                        no_update, no_update, no_update)
    else:
        with get_session() as db:
            running = (
                db.query(Dataset)
                .filter(
                    Dataset.source_type == "pipeline",
                    Dataset.status.in_(["pending", "processing"]),
                )
                .order_by(Dataset.created_at.desc())
                .first()
            )
            if running:
                dataset_id = running.id

    if not dataset_id:
        return no_all

    # Populate all progress fields immediately
    status = get_pipeline_status(dataset_id)
    db_status = status["status"]
    pct = status["progress_pct"]

    badge_color = {
        "pending": "secondary",
        "processing": "primary",
        "complete": "success",
        "failed": "danger",
        "cancelled": "warning",
    }.get(db_status, "secondary")
    badge = dbc.Badge(db_status.upper(), color=badge_color, className="fs-6 p-2")

    completed = set(status.get("steps_completed", []))
    current = status.get("current_step")

    is_longread = bool({"dada2_longread"} & (completed | {current}))
    if not is_longread:
        from app.db.database import SessionLocal
        from app.db.models import Dataset
        _db = SessionLocal()
        try:
            _ds = _db.query(Dataset).filter(Dataset.id == dataset_id).first()
            if _ds and _ds.platform in ("pacbio", "nanopore"):
                is_longread = True
        finally:
            _db.close()
    active_steps = STEPS_LONGREAD if is_longread else STEPS

    step_items = []
    for step in active_steps:
        if step in completed and step != current:
            icon = "  "
            color = "text-success"
        elif step == current and db_status == "processing":
            icon = "  "
            color = "text-primary"
        elif step == current and db_status == "failed":
            icon = "  "
            color = "text-danger"
        else:
            icon = "  "
            color = "text-muted"
        step_items.append(
            html.Div(
                f"{icon} {STEP_LABELS.get(step, step)}",
                className=f"{color} mb-1",
            )
        )
    steps_display = html.Div(step_items)
    log_tail = status.get("log_tail", "")
    bar_style = {"display": "block"}

    # Enable polling for active jobs
    keep_polling = db_status in ("pending", "processing")

    return (
        dataset_id,
        {"display": "block"},
        not keep_polling,   # poll disabled = True when finished
        badge,
        pct,
        f"{pct}%",
        steps_display,
        log_tail,
        {"display": "block"} if db_status == "complete" else {"display": "none"},
        {"display": "none"} if not keep_polling else bar_style,
    )


# ── Download callbacks ────────────────────────────────────────────────────


@dash_app.callback(
    Output("download-asv-table", "data"),
    Input("btn-download-asv", "n_clicks"),
    State("pipeline-active-dataset", "data"),
    prevent_initial_call=True,
)
def on_download_asv(n_clicks, dataset_id):
    """Serve the ASV table as a BIOM file."""
    if not n_clicks or not dataset_id:
        return no_update

    from app.config import DATASET_DIR
    from app.db.database import get_session
    from app.db.models import Dataset

    with get_session() as db:
        ds = db.query(Dataset).filter(Dataset.id == dataset_id).first()

    biom_path = DATASET_DIR / str(dataset_id) / "asv_table.biom"
    if not biom_path.exists():
        return no_update

    return dcc.send_file(
        str(biom_path),
        filename=_download_filename(ds, "asv_table", "biom"),
    )


@dash_app.callback(
    Output("download-history-biom", "data"),
    Input({"type": "btn-history-dl", "index": ALL}, "n_clicks"),
    prevent_initial_call=True,
)
def on_history_download(n_clicks_list):
    """Download BIOM file from the pipeline history table."""
    if not any(n_clicks_list):
        return no_update

    triggered = ctx.triggered_id
    if not triggered:
        return no_update
    dataset_id = triggered["index"]

    from app.config import DATASET_DIR
    from app.db.database import get_session
    from app.db.models import Dataset

    with get_session() as db:
        ds = db.query(Dataset).filter(Dataset.id == dataset_id).first()

    biom_path = DATASET_DIR / str(dataset_id) / "asv_table.biom"
    if not biom_path.exists():
        return no_update

    return dcc.send_file(
        str(biom_path),
        filename=_download_filename(ds, "asv_table", "biom"),
    )


@dash_app.callback(
    Output("download-history-qc", "data"),
    Input({"type": "btn-history-qc", "index": ALL}, "n_clicks"),
    prevent_initial_call=True,
)
def on_history_qc_download(n_clicks_list):
    """Download QC report PDF from the pipeline history table."""
    if not any(n_clicks_list):
        return no_update

    triggered = ctx.triggered_id
    if not triggered:
        return no_update
    dataset_id = triggered["index"]

    from app.config import DATASET_DIR
    from app.db.database import get_session
    from app.db.models import Dataset

    with get_session() as db:
        ds = db.query(Dataset).filter(Dataset.id == dataset_id).first()

    pdf_path = DATASET_DIR / str(dataset_id) / "qc_report.pdf"
    if not pdf_path.exists():
        return no_update

    return dcc.send_file(
        str(pdf_path),
        filename=_download_filename(ds, "qc_report", "pdf"),
    )


@dash_app.callback(
    Output("download-history-log", "data"),
    Input({"type": "btn-history-log", "index": ALL}, "n_clicks"),
    prevent_initial_call=True,
)
def on_history_log_download(n_clicks_list):
    """Download pipeline log from the history table (for failed/cancelled runs)."""
    if not any(n_clicks_list):
        return no_update

    triggered = ctx.triggered_id
    if not triggered:
        return no_update
    dataset_id = triggered["index"]

    from app.config import DATASET_DIR
    from app.db.database import get_session
    from app.db.models import Dataset

    with get_session() as db:
        ds = db.query(Dataset).filter(Dataset.id == dataset_id).first()

    log_path = DATASET_DIR / str(dataset_id) / "pipeline.log"
    if not log_path.exists():
        return no_update

    return dcc.send_file(
        str(log_path),
        filename=_download_filename(ds, "pipeline_log", "txt"),
    )


# ── Delete callback ────────────────────────────────────────────────────────


@dash_app.callback(
    Output("pipeline-history-table", "children", allow_duplicate=True),
    Input({"type": "btn-history-del", "index": ALL}, "n_clicks"),
    prevent_initial_call=True,
)
def on_history_delete(n_clicks_list):
    """Delete a pipeline run: remove dataset DB record and files on disk.

    Uploads are NOT deleted here — they are managed in File Manager.
    """
    if not any(n_clicks_list):
        return no_update

    triggered = ctx.triggered_id
    if not triggered:
        return no_update
    dataset_id = triggered["index"]

    import shutil

    from app.config import DATASET_DIR
    from app.db.database import get_session
    from app.db.models import Dataset

    with get_session() as db:
        ds = db.query(Dataset).filter(Dataset.id == dataset_id).first()
        if not ds:
            return _build_history_table()

        if ds.status in ("pending", "processing"):
            return no_update

        # Delete the dataset (cascades to samples, analysis_results,
        # DatasetFastqFile, etc.)
        db.delete(ds)
        db.commit()

    # Remove dataset files from disk
    ds_path = DATASET_DIR / str(dataset_id)
    if ds_path.exists():
        shutil.rmtree(ds_path, ignore_errors=True)

    return _build_history_table()
