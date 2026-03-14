"""
MicrobiomeDash — File Manager page: upload FASTQ files, manage registered files, attach metadata.
"""
import base64
import io
import shutil
from pathlib import Path

import dash_bootstrap_components as dbc
import dash_uploader as du
import pandas as pd
from dash import ALL, Input, Output, State, ctx, dcc, html, no_update

from app.config import UPLOAD_DIR
from app.dashboard.app import app as dash_app
from app.db.database import SessionLocal
from app.db.models import FastqFile, Upload, UploadMetadata
from app.pipeline.detect import (
    PRIMERS,
    _primer_matches,
    _read_fastq_sequences,
    detect_platform,
    detect_sequencing_type,
    detect_variable_region,
    extract_sample_name,
)

# ── Layout ───────────────────────────────────────────────────────────────────

layout = dbc.Container(
    [
        html.H3("File Manager", className="mb-4"),
        # ── FASTQ Upload ─────────────────────────────────────────────────
        dbc.Card(
            [
                dbc.CardHeader(html.H5("Upload FASTQ Files", className="mb-0")),
                dbc.CardBody(
                    [
                        du.Upload(
                            id="fm-du-upload",
                            text="Drag & drop FASTQ files here or click to browse",
                            text_completed="Uploaded: ",
                            max_file_size=5120,  # 5 GB per file
                            chunk_size=10,  # 10 MB chunks
                            max_files=100,
                            filetypes=["gz", "fastq", "fq"],
                            cancel_button=True,
                            pause_button=True,
                        ),
                        html.Div(id="fm-upload-file-status", className="mt-2"),
                        dbc.InputGroup(
                            [
                                dbc.InputGroupText("Study"),
                                dbc.Input(
                                    id="fm-study-name",
                                    placeholder="e.g. CRC_cohort_2024",
                                ),
                            ],
                            className="mt-2",
                            style={"maxWidth": "400px"},
                        ),
                        # Hidden button kept for callback wiring
                        dbc.Button(
                            id="fm-btn-register",
                            style={"display": "none"},
                        ),
                        html.Div(id="fm-upload-status", className="mt-3"),
                    ]
                ),
            ],
            className="mb-4",
        ),
        # Hidden stores
        dcc.Store(id="fm-upload-trigger", data=0),
        dcc.Store(id="fm-sort-field", data="sample_name"),
        dcc.Store(id="fm-sort-asc", data=True),
        dcc.Store(id="fm-checked-samples", data=[]),
        dcc.Store(id="fm-pending-file-count", data=0),
        dcc.Interval(id="fm-register-interval", interval=2000, disabled=True),
        dcc.Download(id="fm-download-meta-selected"),
        dcc.ConfirmDialog(id="fm-confirm-delete", message=""),
        # ── Metadata Upload ──────────────────────────────────────────────
        dbc.Card(
            [
                dbc.CardHeader(
                    dbc.Row(
                        [
                            dbc.Col(html.H5("Upload Sample Metadata", className="mb-0"), width="auto"),
                            dbc.Col(
                                [
                                    dbc.Button(
                                        "Download Template",
                                        id="fm-btn-dl-template",
                                        color="outline-secondary",
                                        size="sm",
                                        className="me-2",
                                    ),
                                    dbc.Button(
                                        "Download Full Metadata (.tsv)",
                                        id="fm-btn-dl-meta",
                                        color="info",
                                        size="sm",
                                        disabled=True,
                                    ),
                                ],
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
                        html.P(
                            "Upload a TSV file with sample names in the first column. "
                            "Metadata will be matched to registered FASTQ files by sample name.",
                            className="text-muted small mb-3",
                        ),
                        dcc.Upload(
                            id="fm-upload-metadata",
                            children=html.Div(
                                [
                                    "Drag & drop a metadata file (.csv or .tsv), or ",
                                    html.A("click to browse"),
                                ]
                            ),
                            style={
                                "borderWidth": "2px",
                                "borderStyle": "dashed",
                                "borderRadius": "10px",
                                "textAlign": "center",
                                "padding": "20px",
                                "cursor": "pointer",
                            },
                            multiple=False,
                            accept=".tsv,.csv,.txt",
                        ),
                        html.Div(id="fm-meta-status", className="mt-3"),
                        html.Div(id="fm-meta-preview", className="mt-3"),
                    ]
                ),
            ],
            className="mb-4",
        ),
        # ── Registered Files Table ─────────────────────────────────────────
        dbc.Card(
            [
                dbc.CardHeader(
                    dbc.Row(
                        [
                            dbc.Col(html.H5("Registered Files", className="mb-0"), width="auto"),
                            dbc.Col(
                                [
                                    dbc.Button(
                                        "Delete Selected",
                                        id="fm-btn-delete-selected",
                                        color="danger",
                                        size="sm",
                                        disabled=True,
                                        className="me-2",
                                    ),
                                    dbc.Button(
                                        "Metadata",
                                        id="fm-btn-dl-meta-selected",
                                        color="info",
                                        size="sm",
                                        disabled=True,
                                        className="me-3",
                                    ),
                                ],
                                width="auto",
                                className="ms-auto",
                            ),
                            dbc.Col(
                                dbc.ButtonGroup(
                                    [
                                        dbc.Button("Sample", id="fm-sort-sample", color="primary", size="sm", outline=True),
                                        dbc.Button("Type", id="fm-sort-type", color="primary", size="sm", outline=True),
                                        dbc.Button("Reads", id="fm-sort-reads", color="primary", size="sm", outline=True),
                                        dbc.Button("Region", id="fm-sort-region", color="primary", size="sm", outline=True),
                                        dbc.Button("Read Len", id="fm-sort-readlen", color="primary", size="sm", outline=True),
                                        dbc.Button("Study", id="fm-sort-study", color="primary", size="sm", outline=True),
                                        dbc.Button("Source", id="fm-sort-source", color="primary", size="sm", outline=True),
                                        dbc.Button("Metadata", id="fm-sort-metadata", color="primary", size="sm", outline=True),
                                        dbc.Button("Date", id="fm-sort-date", color="primary", size="sm", outline=True),
                                    ],
                                    size="sm",
                                ),
                                width="auto",
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
                                dbc.Col(dbc.Input(id="fm-f-sample", placeholder="Sample", debounce=True, size="sm"), width=2),
                                dbc.Col(dbc.Input(id="fm-f-type", placeholder="Type", debounce=True, size="sm"), width=1),
                                dbc.Col(width=1),  # # Reads (no filter)
                                dbc.Col(dbc.Input(id="fm-f-region", placeholder="Region", debounce=True, size="sm"), width=1),
                                dbc.Col(width=1),  # Primers (no filter)
                                dbc.Col(width=1),  # Avg Read Len (no filter)
                                dbc.Col(dbc.Input(id="fm-f-study", placeholder="Study", debounce=True, size="sm"), width=1),
                                dbc.Col(dbc.Input(id="fm-f-source", placeholder="Source", debounce=True, size="sm"), width=1),
                                dbc.Col(dbc.Input(id="fm-f-metadata", placeholder="Metadata", debounce=True, size="sm"), width=2),
                                dbc.Col(width=1),  # Date (no filter)
                            ],
                            className="mb-2 g-1",
                        ),
                        html.Div(id="fm-files-table", children="No files registered yet."),
                    ]
                ),
            ],
            className="mb-4",
        ),
        # Hidden download components
        dcc.Download(id="fm-download-metadata"),
        dcc.Download(id="fm-download-template"),
    ],
    fluid=True,
)


# ── Helpers ──────────────────────────────────────────────────────────────────


def _human_size(size_mb):
    """Format size in MB to human-readable string."""
    if size_mb is None:
        return "?"
    if size_mb < 1:
        return f"{size_mb * 1024:.0f} KB"
    if size_mb >= 1024:
        return f"{size_mb / 1024:.1f} GB"
    return f"{size_mb:.1f} MB"


def _count_reads(fastq_path: Path) -> int | None:
    """Count total reads in a FASTQ file (each 4 lines = 1 read)."""
    import gzip
    try:
        opener = gzip.open if str(fastq_path).endswith(".gz") else open
        count = 0
        with opener(fastq_path, "rt") as f:
            for _ in f:
                count += 1
        return count // 4
    except Exception:
        return None


def _avg_read_length(fastq_path: Path, n_reads: int = 200) -> int | None:
    """Compute average read length from the first n reads of a FASTQ file."""
    import gzip
    try:
        opener = gzip.open if str(fastq_path).endswith(".gz") else open
        lengths = []
        with opener(fastq_path, "rt") as f:
            line_num = 0
            for line in f:
                line_num += 1
                if line_num % 4 == 2:  # sequence line
                    lengths.append(len(line.strip()))
                    if len(lengths) >= n_reads:
                        break
        return round(sum(lengths) / len(lengths)) if lengths else None
    except Exception:
        return None


def _build_files_table(sort_by="sample_name", filters=None, ascending=True, checked_samples=None):
    """Build a sample-level table of all registered FASTQ files with metadata."""
    if filters is None:
        filters = {}
    if checked_samples is None:
        checked_samples = set()
    db = SessionLocal()
    try:
        query = (
            db.query(FastqFile, Upload.created_at, Upload.variable_region, Upload.primers_detected, Upload.study)
            .join(Upload, FastqFile.upload_id == Upload.id)
        )

        results = query.all()
        if not results:
            return "No files registered yet."

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

        # Group files by sample_name
        from collections import defaultdict
        sample_map = defaultdict(lambda: {
            "files": [], "total_size": 0.0, "date": None,
            "region": None, "upload_id": None, "primers_detected": None,
            "study": None,
        })
        for f, upload_date, region, primers_detected, study in results:
            s = sample_map[f.sample_name]
            s["files"].append(f)
            s["total_size"] += f.file_size_mb or 0
            s["upload_id"] = f.upload_id
            if upload_date:
                s["date"] = upload_date
            if region:
                s["region"] = region
            if primers_detected is not None:
                s["primers_detected"] = primers_detected
            if study:
                s["study"] = study

        # Build rows
        sample_rows = []
        for sample_name, info in sample_map.items():
            files = sorted(info["files"], key=lambda x: x.filename)
            directions = set(f.read_direction for f in files)
            if "R1" in directions and "R2" in directions:
                direction = "PE"
            elif directions == {"single"}:
                direction = "SE"
            else:
                direction = ", ".join(sorted(directions))

            file_links = [
                html.A(
                    f.filename,
                    href=f"/api/uploads/{f.upload_id}/files/{f.filename}",
                    target="_blank",
                    className="text-info me-2",
                )
                for f in files
            ]

            total_reads = sum(f.read_count or 0 for f in files)
            avg_lengths = [f.avg_read_length for f in files if f.avg_read_length]
            avg_len = round(sum(avg_lengths) / len(avg_lengths)) if avg_lengths else None

            meta_keys = meta_by_sample.get(sample_name, [])
            meta_str = ", ".join(sorted(meta_keys)) if meta_keys else "—"
            source_val = source_by_sample.get(sample_name, "")

            study_val = info["study"] or study_by_sample.get(sample_name, "")

            sample_rows.append({
                "sample_name": sample_name,
                "direction": direction,
                "file_links": file_links,
                "total_reads": total_reads,
                "region": info["region"],
                "avg_len": avg_len,
                "study": study_val,
                "source": source_val,
                "date": info["date"],
                "n_files": len(files),
                "metadata": meta_str,
                "primers_detected": info["primers_detected"],
            })

        # Per-column filters (AND logic)
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
            return html.P("No matching samples.", className="text-muted")

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
        # Numeric fields: higher = "first" in ascending, so invert reverse
        numeric_fields = {"total_reads", "avg_read_length", "created_at"}
        key_fn = sort_keys.get(sort_by, sort_keys["sample_name"])
        effective_reverse = (not reverse) if sort_by in numeric_fields else reverse
        sample_rows.sort(key=key_fn, reverse=effective_reverse)

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
                                id={"type": "fm-check", "index": r["sample_name"]},
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

        total_files = sum(r["n_files"] for r in sample_rows)
        total_reads = sum(r["total_reads"] for r in sample_rows)

        # Arrow indicator for sorted column
        arrow = " ^" if ascending else " v"
        def _th(label, field):
            suffix = arrow if sort_by == field else ""
            return html.Th(f"{label}{suffix}")

        return html.Div(
            [
                html.P(
                    f"{len(sample_rows)} samples, {total_files} files, {total_reads:,} reads total",
                    className="text-muted mb-2",
                ),
                dbc.Table(
                    [
                        html.Thead(
                            html.Tr(
                                [
                                    html.Th(
                                        dbc.Checkbox(id="fm-check-all", value=False),
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
    finally:
        db.close()


# ── Callback 1a: Track file uploads (fires per-file) ────────────────────────


@du.callback(
    output=Output("fm-upload-file-status", "children"),
    id="fm-du-upload",
)
def on_file_uploaded(file_paths):
    """Called by dash-uploader each time a file finishes uploading."""
    if not file_paths:
        return no_update
    # file_paths only contains the latest file; scan the upload directory
    # to get the real total count of all uploaded files so far.
    latest = Path(file_paths[0])
    upload_dir = latest.parent  # UPLOAD_DIR/{upload_id}/
    all_files = [
        f.name for f in sorted(upload_dir.iterdir())
        if f.is_file() and (f.name.endswith(".fastq.gz") or f.name.endswith(".fq.gz"))
    ]
    n = len(all_files)
    return html.Small(
        f"{n} file(s) transferred: {', '.join(all_files)}",
        className="text-success",
    )


# ── Callback 1b: Start debounce timer when files arrive ─────────────────────


@dash_app.callback(
    Output("fm-register-interval", "disabled"),
    Output("fm-pending-file-count", "data"),
    Output("fm-btn-register", "disabled"),
    Input("fm-upload-file-status", "children"),
    State("fm-du-upload", "upload_id"),
    prevent_initial_call=True,
)
def on_file_arrived(file_status, du_upload_id):
    """When a new file finishes uploading, (re)start the debounce interval."""
    if not file_status or not du_upload_id:
        return no_update, no_update, no_update
    staging_dir = UPLOAD_DIR / du_upload_id
    if not staging_dir.exists():
        return no_update, no_update, no_update
    n = sum(
        1 for f in staging_dir.iterdir()
        if f.is_file() and (f.name.endswith(".fastq.gz") or f.name.endswith(".fq.gz"))
    )
    # Store current count and enable the 2s interval
    return False, n, no_update


# ── Callback 1c: Debounced auto-registration ────────────────────────────────


@dash_app.callback(
    Output("fm-upload-status", "children"),
    Output("fm-upload-trigger", "data"),
    Output("fm-register-interval", "disabled", allow_duplicate=True),
    Output("fm-pending-file-count", "data", allow_duplicate=True),
    Input("fm-register-interval", "n_intervals"),
    State("fm-du-upload", "upload_id"),
    State("fm-study-name", "value"),
    State("fm-upload-trigger", "data"),
    State("fm-pending-file-count", "data"),
    prevent_initial_call=True,
)
def on_register_upload(n_intervals, du_upload_id, study_name, trigger, prev_count):
    """Register once the file count has stabilised (no new files for 2s)."""
    if not du_upload_id:
        return no_update, no_update, True, no_update

    staging_dir = UPLOAD_DIR / du_upload_id
    if not staging_dir.exists():
        return no_update, no_update, True, no_update

    # Count current files
    staged_files = [
        f for f in sorted(staging_dir.iterdir())
        if f.is_file() and (f.name.endswith(".fastq.gz") or f.name.endswith(".fq.gz"))
    ]
    current_count = len(staged_files)

    if not staged_files:
        shutil.rmtree(staging_dir, ignore_errors=True)
        return no_update, no_update, True, no_update

    # If count changed since last tick, wait another cycle
    if current_count != prev_count:
        return no_update, no_update, no_update, current_count

    saved_filenames = [f.name for f in staged_files]

    # Check for duplicate filenames already registered in the DB
    db = SessionLocal()
    try:
        existing = (
            db.query(FastqFile.filename)
            .filter(FastqFile.filename.in_(saved_filenames))
            .all()
        )
        dupes = [row.filename for row in existing]
    finally:
        db.close()

    if dupes:
        # Clean up staging directory
        shutil.rmtree(staging_dir, ignore_errors=True)
        dupe_list = html.Ul([html.Li(f, className="text-warning") for f in sorted(dupes)])
        return (
            dbc.Alert(
                [
                    html.P("These filenames are already registered:", className="mb-1 fw-bold"),
                    dupe_list,
                    html.P("Please rename the files and upload again.", className="mb-0 mt-2"),
                ],
                color="danger",
            ),
            no_update,
            True,
            0,
        )

    # Move files from staging to flat UPLOAD_DIR
    total_size = 0.0
    for fpath in staged_files:
        dest = UPLOAD_DIR / fpath.name
        fpath.rename(dest)
        total_size += dest.stat().st_size / (1024 * 1024)
    # Remove empty staging directory
    shutil.rmtree(staging_dir, ignore_errors=True)

    # Detect sequencing type
    detection = detect_sequencing_type(saved_filenames)

    # Detect variable region from first R1 / single-end file
    variable_region = None
    first_sample = next(iter(detection["samples"].values()), {})
    r1_name = first_sample.get("R1")
    if r1_name:
        try:
            region_result = detect_variable_region(UPLOAD_DIR / r1_name)
            variable_region = region_result["region"]
        except Exception:
            pass

    # Detect platform (Illumina / PacBio / Nanopore)
    platform = None
    if r1_name:
        try:
            platform_result = detect_platform(UPLOAD_DIR / r1_name)
            platform = platform_result["platform"]
        except Exception:
            pass

    # For long-read data, force region to V1-V9
    if platform in ("pacbio", "nanopore"):
        variable_region = "V1-V9"
        # Override sequencing type to single-end for long reads
        detection["type"] = "single-end"

    # Detect primer presence (check if reads start with expected forward primer)
    primers_detected = None
    if variable_region and r1_name:
        try:
            fwd_primer = PRIMERS[variable_region]["forward"]
            seqs = _read_fastq_sequences(UPLOAD_DIR / r1_name, n_reads=100)
            if seqs:
                matches = sum(1 for s in seqs if _primer_matches(s, fwd_primer))
                primers_detected = (matches / len(seqs)) >= 0.30
        except Exception:
            pass

    # Create DB records
    db = SessionLocal()
    try:
        upload = Upload(
            upload_dir=str(UPLOAD_DIR),
            sequencing_type=detection["type"],
            variable_region=variable_region,
            platform=platform,
            primers_detected=primers_detected,
            study=(study_name or "").strip() or None,
            total_files=len(saved_filenames),
            total_size_mb=round(total_size, 2),
            status="uploaded",
        )
        db.add(upload)
        db.flush()

        for filename in saved_filenames:
            sample_name = extract_sample_name(filename)
            sample_info = detection["samples"].get(sample_name, {})
            if sample_info.get("R1") == filename:
                read_direction = "R1"
            elif sample_info.get("R2") == filename:
                read_direction = "R2"
            else:
                read_direction = "single"

            fpath = UPLOAD_DIR / filename
            avg_len = _avg_read_length(fpath)
            n_reads = _count_reads(fpath)
            db.add(
                FastqFile(
                    upload_id=upload.id,
                    sample_name=sample_name,
                    filename=filename,
                    file_path=str(fpath),
                    read_direction=read_direction,
                    file_size_mb=round(fpath.stat().st_size / (1024 * 1024), 2),
                    read_count=n_reads,
                    avg_read_length=avg_len,
                )
            )

        db.commit()
    except Exception as e:
        db.rollback()
        return dbc.Alert(f"Database error: {e}", color="danger"), no_update, True, 0
    finally:
        db.close()

    # Build success message
    region_str = f", region: {variable_region}" if variable_region else ""
    platform_str = f", platform: {platform}" if platform and platform != "illumina" else ""
    if primers_detected is True:
        primer_str = ", primers: detected"
    elif primers_detected is False:
        primer_str = ", primers: not detected"
    else:
        primer_str = ""
    msg = (
        f"Registered {len(saved_filenames)} files "
        f"({_human_size(total_size)}), {detection['type']}{region_str}{platform_str}{primer_str}"
    )
    warnings = detection.get("errors", [])
    alert_children = [html.P(msg, className="mb-0")]
    if warnings:
        alert_children.append(
            html.Small(
                " | ".join(warnings), className="text-warning d-block mt-1"
            )
        )

    return (
        dbc.Alert(alert_children, color="success"),
        (trigger or 0) + 1,
        True,
        0,
    )


# ── Callback 2: Sort buttons ─────────────────────────────────────────────────


_SORT_BUTTONS = [
    ("fm-sort-sample", "sample_name"),
    ("fm-sort-type", "direction"),
    ("fm-sort-reads", "total_reads"),
    ("fm-sort-region", "region"),
    ("fm-sort-readlen", "avg_read_length"),
    ("fm-sort-study", "study"),
    ("fm-sort-source", "source"),
    ("fm-sort-metadata", "metadata"),
    ("fm-sort-date", "created_at"),
]


@dash_app.callback(
    Output("fm-sort-field", "data"),
    Output("fm-sort-asc", "data"),
    *[Output(btn_id, "outline") for btn_id, _ in _SORT_BUTTONS],
    *[Input(btn_id, "n_clicks") for btn_id, _ in _SORT_BUTTONS],
    State("fm-sort-field", "data"),
    State("fm-sort-asc", "data"),
    prevent_initial_call=True,
)
def on_sort_click(*args):
    """Update sort field and toggle direction. Highlight active button."""
    current_field = args[-2]
    current_asc = args[-1]
    triggered = ctx.triggered_id
    sort_map = {btn_id: field for btn_id, field in _SORT_BUTTONS}
    new_field = sort_map.get(triggered, "sample_name")

    # Toggle direction if same button clicked again
    if new_field == current_field:
        new_asc = not current_asc
    else:
        new_asc = True

    # Highlight: active button is solid (outline=False), others outlined
    outlines = [new_field != field for _, field in _SORT_BUTTONS]
    return new_field, new_asc, *outlines


# ── Callback 3: Refresh files table ──────────────────────────────────────────


@dash_app.callback(
    Output("fm-files-table", "children"),
    Input("fm-upload-trigger", "data"),
    Input("fm-sort-field", "data"),
    Input("fm-sort-asc", "data"),
    Input("fm-f-sample", "value"),
    Input("fm-f-type", "value"),
    Input("fm-f-region", "value"),
    Input("fm-f-study", "value"),
    Input("fm-f-source", "value"),
    Input("fm-f-metadata", "value"),
    Input("url", "pathname"),
    State("fm-checked-samples", "data"),
)
def refresh_files_table(trigger, sort_by, sort_asc, f_sample, f_type, f_region, f_study, f_source, f_metadata, pathname, checked):
    """Rebuild the flat files table whenever data changes, sort, or filter changes."""
    if pathname != "/files":
        return no_update
    filters = {
        "sample": (f_sample or "").strip().lower(),
        "type": (f_type or "").strip().lower(),
        "region": (f_region or "").strip().lower(),
        "study": (f_study or "").strip().lower(),
        "source": (f_source or "").strip().lower(),
        "metadata": (f_metadata or "").strip().lower(),
    }
    checked_set = set(checked) if checked else set()
    return _build_files_table(sort_by or "sample_name", filters, sort_asc if sort_asc is not None else True, checked_samples=checked_set)


# ── Callback 4: Metadata Upload ─────────────────────────────────────────────


@dash_app.callback(
    Output("fm-meta-status", "children"),
    Output("fm-meta-preview", "children"),
    Output("fm-upload-trigger", "data", allow_duplicate=True),
    Input("fm-upload-metadata", "contents"),
    State("fm-upload-metadata", "filename"),
    State("fm-upload-trigger", "data"),
    prevent_initial_call=True,
)
def on_metadata_upload(content, filename, trigger):
    if not content:
        return no_update, no_update, no_update

    # Decode the metadata file (auto-detect CSV or TSV)
    try:
        _, content_encoded = content.split(",", 1)
        decoded = base64.b64decode(content_encoded)
        text = decoded.decode("utf-8")
        first_line = text.split("\n")[0]
        sep = "\t" if ("\t" in first_line or (filename and filename.endswith(".tsv"))) else ","
        df = pd.read_csv(io.StringIO(text), sep=sep)
    except Exception as e:
        return dbc.Alert(f"Error parsing metadata file: {e}", color="danger"), no_update, no_update

    if df.empty or len(df.columns) < 2:
        return (
            dbc.Alert(
                "Metadata file must have at least 2 columns (sample_name + metadata).",
                color="warning",
            ),
            no_update,
            no_update,
        )

    # Lowercase column headers and string values
    df.columns = [c.lower() for c in df.columns]
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].str.lower()

    # First column = sample name
    sample_col = df.columns[0]
    meta_cols = list(df.columns[1:])

    # Cross-reference with all registered FASTQ sample names
    db = SessionLocal()
    try:
        all_fastq = db.query(FastqFile).all()
        fastq_samples = set(f.sample_name for f in all_fastq)
        meta_samples = set(df[sample_col].astype(str))

        matched = fastq_samples & meta_samples
        meta_only = meta_samples - fastq_samples
        fastq_only = fastq_samples - meta_samples

        # Build validation messages
        messages = []
        if matched:
            messages.append(
                dbc.Alert(
                    f"Matched {len(matched)} sample(s).",
                    color="success",
                    className="py-1 mb-1",
                )
            )
        if meta_only:
            messages.append(
                dbc.Alert(
                    f"Metadata-only (no FASTQ): {', '.join(sorted(meta_only))}",
                    color="warning",
                    className="py-1 mb-1",
                )
            )
        if fastq_only:
            messages.append(
                dbc.Alert(
                    f"FASTQ-only (no metadata): {', '.join(sorted(fastq_only))}",
                    color="info",
                    className="py-1 mb-1",
                )
            )

        # Save metadata to each upload that has matching samples
        # Group fastq files by upload_id
        upload_ids = set(f.upload_id for f in all_fastq)
        for uid in upload_ids:
            upload_samples = set(
                f.sample_name for f in all_fastq if f.upload_id == uid
            )
            matching = upload_samples & meta_samples
            if not matching:
                continue

            # Delete existing metadata rows for this upload, insert fresh
            db.query(UploadMetadata).filter(
                UploadMetadata.upload_id == uid
            ).delete()

            for _, row in df.iterrows():
                sample_name = str(row[sample_col])
                if sample_name not in upload_samples:
                    continue
                for col in meta_cols:
                    val = row[col]
                    if pd.notna(val):
                        db.add(
                            UploadMetadata(
                                upload_id=uid,
                                sample_name=sample_name,
                                key=col,
                                value=str(val),
                            )
                        )

        db.commit()

        # Preview table
        preview_df = df.head(20)
        table = dbc.Table.from_dataframe(
            preview_df, bordered=True, color="dark", hover=True, size="sm", striped=True
        )
        preview = html.Div(
            [
                html.H6(f"Metadata Preview ({len(df)} rows, {len(meta_cols)} columns)"),
                table,
            ]
        )

        return html.Div(messages), preview, (trigger or 0) + 1

    except Exception as e:
        db.rollback()
        return dbc.Alert(f"Error saving metadata: {e}", color="danger"), no_update, no_update
    finally:
        db.close()


# ── Callback 5: Enable/disable metadata download button ──────────────────────


@dash_app.callback(
    Output("fm-btn-dl-meta", "disabled"),
    Input("fm-upload-trigger", "data"),
    Input("url", "pathname"),
)
def toggle_meta_download(trigger, pathname):
    """Enable download button when metadata exists."""
    if pathname != "/files":
        return no_update
    db = SessionLocal()
    try:
        has_meta = db.query(UploadMetadata).first() is not None
        return not has_meta
    finally:
        db.close()


# ── Callback 6: Download full metadata TSV ───────────────────────────────────


@dash_app.callback(
    Output("fm-download-metadata", "data"),
    Input("fm-btn-dl-meta", "n_clicks"),
    prevent_initial_call=True,
)
def on_download_metadata(n_clicks):
    """Export all stored metadata as a TSV file."""
    if not n_clicks:
        return no_update

    db = SessionLocal()
    try:
        rows = db.query(UploadMetadata).order_by(UploadMetadata.sample_name).all()
        if not rows:
            return no_update

        # Pivot to dataframe
        data = {}
        for meta in rows:
            if meta.sample_name not in data:
                data[meta.sample_name] = {"sample_name": meta.sample_name}
            data[meta.sample_name][meta.key] = meta.value

        df = pd.DataFrame(data.values())
        cols = ["sample_name"] + [c for c in df.columns if c != "sample_name"]
        df = df[cols]

        return dcc.send_data_frame(df.to_csv, "metadata.tsv", sep="\t", index=False)
    finally:
        db.close()


# ── Callback 6b: Download metadata template CSV ──────────────────────────────


@dash_app.callback(
    Output("fm-download-template", "data"),
    Input("fm-btn-dl-template", "n_clicks"),
    prevent_initial_call=True,
)
def on_download_template(n_clicks):
    """Generate a metadata template CSV with sample names pre-filled."""
    if not n_clicks:
        return no_update

    db = SessionLocal()
    try:
        # Get all unique sample names from registered uploads
        sample_names = (
            db.query(FastqFile.sample_name)
            .distinct()
            .order_by(FastqFile.sample_name)
            .all()
        )
        if not sample_names:
            return no_update

        names = [row.sample_name for row in sample_names]
        df = pd.DataFrame({
            "sample-id": names,
            "group": [""] * len(names),
            "treatment": [""] * len(names),
        })

        return dcc.send_data_frame(df.to_csv, "metadata_template.csv", index=False)
    finally:
        db.close()


# ── Callback 7: Row checkbox → update checked list + button state ─────────


@dash_app.callback(
    Output("fm-checked-samples", "data"),
    Output("fm-btn-delete-selected", "disabled"),
    Output("fm-btn-delete-selected", "children"),
    Output("fm-btn-dl-meta-selected", "disabled"),
    Output("fm-btn-dl-meta-selected", "children"),
    Input({"type": "fm-check", "index": ALL}, "value"),
    State({"type": "fm-check", "index": ALL}, "id"),
    prevent_initial_call=True,
)
def on_row_check(values, ids):
    """Sync checkbox states to checked-samples store and update button labels."""
    checked = [
        id_dict["index"] for id_dict, val in zip(ids, values) if val
    ]
    n = len(checked)
    if n:
        return (
            checked,
            False,
            f"Delete Selected ({n})",
            False,
            f"Metadata ({n})",
        )
    return checked, True, "Delete Selected", True, "Metadata"


# ── Callback 8: Select All checkbox ──────────────────────────────────────


@dash_app.callback(
    Output({"type": "fm-check", "index": ALL}, "value"),
    Input("fm-check-all", "value"),
    State({"type": "fm-check", "index": ALL}, "id"),
    prevent_initial_call=True,
)
def on_check_all(select_all, ids):
    """Toggle all visible row checkboxes."""
    return [bool(select_all)] * len(ids)


# ── Callback 9: Delete Selected click → show confirm dialog ─────────────


@dash_app.callback(
    Output("fm-confirm-delete", "displayed"),
    Output("fm-confirm-delete", "message"),
    Input("fm-btn-delete-selected", "n_clicks"),
    State("fm-checked-samples", "data"),
    prevent_initial_call=True,
)
def on_delete_click(n_clicks, checked):
    """Show a confirmation dialog before deleting."""
    if not n_clicks or not checked:
        return False, ""
    n = len(checked)
    names = ", ".join(sorted(checked)[:10])
    suffix = f"... and {n - 10} more" if n > 10 else ""
    return True, (
        f"Permanently delete {n} sample(s) and their FASTQ files?\n\n"
        f"{names}{suffix}\n\n"
        "This cannot be undone."
    )


# ── Callback 10: Delete confirmed → remove files from DB + disk ─────────


@dash_app.callback(
    Output("fm-upload-trigger", "data", allow_duplicate=True),
    Output("fm-checked-samples", "data", allow_duplicate=True),
    Input("fm-confirm-delete", "submit_n_clicks"),
    State("fm-checked-samples", "data"),
    State("fm-upload-trigger", "data"),
    prevent_initial_call=True,
)
def on_delete_confirmed(submit_n_clicks, checked, trigger):
    """Delete selected samples: FASTQ files from disk + DB records."""
    if not submit_n_clicks or not checked:
        return no_update, no_update

    checked_set = set(checked)
    db = SessionLocal()
    try:
        # Find all FastqFile records for checked samples
        fastq_records = (
            db.query(FastqFile)
            .filter(FastqFile.sample_name.in_(checked_set))
            .all()
        )
        affected_upload_ids = set()
        for ff in fastq_records:
            affected_upload_ids.add(ff.upload_id)
            # Delete actual file from disk
            p = Path(ff.file_path)
            p.unlink(missing_ok=True)
            db.delete(ff)

        # Delete UploadMetadata rows for these samples
        db.query(UploadMetadata).filter(
            UploadMetadata.sample_name.in_(checked_set)
        ).delete(synchronize_session="fetch")

        db.flush()

        # Clean up orphaned Uploads (no remaining FastqFiles)
        for uid in affected_upload_ids:
            remaining = (
                db.query(FastqFile)
                .filter(FastqFile.upload_id == uid)
                .count()
            )
            if remaining == 0:
                upload = db.query(Upload).filter(Upload.id == uid).first()
                if upload:
                    db.delete(upload)

        db.commit()
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()

    return (trigger or 0) + 1, []


# ── Callback 11: Download metadata for selected samples ──────────────────


@dash_app.callback(
    Output("fm-download-meta-selected", "data"),
    Input("fm-btn-dl-meta-selected", "n_clicks"),
    State("fm-checked-samples", "data"),
    prevent_initial_call=True,
)
def on_dl_meta_selected(n_clicks, checked):
    """Export metadata for only the checked samples as a TSV file."""
    if not n_clicks or not checked:
        return no_update

    checked_set = set(checked)
    db = SessionLocal()
    try:
        rows = (
            db.query(UploadMetadata)
            .filter(UploadMetadata.sample_name.in_(checked_set))
            .order_by(UploadMetadata.sample_name)
            .all()
        )
        if not rows:
            return no_update

        # Pivot to dataframe
        data = {}
        for meta in rows:
            if meta.sample_name not in data:
                data[meta.sample_name] = {"sample_name": meta.sample_name}
            data[meta.sample_name][meta.key] = meta.value

        df = pd.DataFrame(data.values())
        cols = ["sample_name"] + [c for c in df.columns if c != "sample_name"]
        df = df[cols]

        return dcc.send_data_frame(df.to_csv, "metadata_selected.tsv", sep="\t", index=False)
    finally:
        db.close()
