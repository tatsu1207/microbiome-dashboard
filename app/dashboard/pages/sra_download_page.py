"""
SRA Downloader — Download FASTQ files from NCBI SRA by accession number.
"""
import uuid

import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, no_update

from app.config import SRA_CACHE_DIR
from app.dashboard.app import app as dash_app
from app.sra.downloader import (
    cancel_download,
    get_job_status,
    launch_download,
    parse_accessions,
    validate_accessions,
)


def get_layout():
    return dbc.Container(
        [
            html.H3("SRA Downloader", className="mb-3"),
            html.P(
                "Download 16S rRNA FASTQ files from NCBI SRA. "
                "Enter SRR accession numbers (individual runs), SRP (study), "
                "or PRJNA (BioProject) accessions. "
                "Downloaded files will be registered automatically in the File Manager.",
                className="text-muted mb-4",
            ),
            dbc.Alert(
                [
                    html.Strong("Tip: "),
                    "Download one dataset per job. Do not mix accessions from different "
                    "sequencing platforms (e.g. Illumina + PacBio) or different variable "
                    "regions (e.g. V4 + V3-V4) in the same download — they will be "
                    "registered as a single upload and the pipeline cannot handle mixed types.",
                ],
                color="info",
                className="mb-3",
            ),
            # Input section
            dbc.Row(
                dbc.Col(
                    [
                        dbc.Label("Accession Numbers"),
                        dbc.Textarea(
                            id="sra-dl-input-accessions",
                            placeholder=(
                                "Enter accessions, one per line or comma-separated.\n"
                                "Examples:\n"
                                "  SRR12345678\n"
                                "  SRP123456\n"
                                "  PRJNA123456\n\n"
                                "Note: Use one download per dataset/platform."
                            ),
                            style={"height": "150px"},
                        ),
                    ],
                    md=6,
                ),
                className="mb-3",
            ),
            # Launch / Cancel buttons
            html.Div(
                [
                    dbc.Button(
                        "Download",
                        id="sra-dl-btn-launch",
                        color="primary",
                        className="me-2",
                    ),
                    dbc.Button(
                        "Cancel",
                        id="sra-dl-btn-cancel",
                        color="danger",
                        outline=True,
                        style={"display": "none"},
                    ),
                ],
                className="mb-3",
            ),
            # Status / error messages
            html.Div(id="sra-dl-status", className="mb-3"),
            # Progress section (hidden by default)
            html.Div(
                id="sra-dl-progress-section",
                children=[
                    dbc.Progress(
                        id="sra-dl-progress-bar",
                        value=0,
                        striped=True,
                        animated=True,
                        className="mb-2",
                    ),
                    html.P(id="sra-dl-progress-text", className="text-muted small"),
                    # Live log area
                    html.Pre(
                        id="sra-dl-log",
                        style={
                            "maxHeight": "300px",
                            "overflowY": "auto",
                            "fontSize": "0.8rem",
                            "padding": "10px",
                            "borderRadius": "4px",
                            "whiteSpace": "pre-wrap",
                        },
                        className="bg-dark text-light mt-2",
                    ),
                ],
                style={"display": "none"},
            ),
            # Session store for active job
            dcc.Store(id="sra-dl-active-job", storage_type="session"),
            # Polling interval (disabled by default)
            dcc.Interval(id="sra-dl-poll", interval=3000, disabled=True),
            # One-shot init interval to restore state on page load
            dcc.Interval(id="sra-dl-init", interval=500, max_intervals=1),
        ],
        fluid=True,
    )


# ── Callback: Launch download ────────────────────────────────────────────────


@dash_app.callback(
    Output("sra-dl-active-job", "data"),
    Output("sra-dl-poll", "disabled"),
    Output("sra-dl-status", "children"),
    Output("sra-dl-progress-section", "style"),
    Output("sra-dl-btn-launch", "style"),
    Output("sra-dl-btn-cancel", "style"),
    Input("sra-dl-btn-launch", "n_clicks"),
    State("sra-dl-input-accessions", "value"),
    prevent_initial_call=True,
)
def on_launch(n_clicks, accession_text):
    if not n_clicks or not accession_text:
        return (
            no_update, no_update,
            dbc.Alert("Please enter at least one accession number.", color="warning"),
            no_update, no_update, no_update,
        )

    # Parse and validate
    raw = parse_accessions(accession_text)
    valid, invalid = validate_accessions(raw)

    if invalid:
        return (
            no_update, no_update,
            dbc.Alert(
                f"Invalid accessions: {', '.join(invalid)}. "
                "Use SRR, SRP, or PRJNA format.",
                color="danger",
            ),
            no_update, no_update, no_update,
        )

    if not valid:
        return (
            no_update, no_update,
            dbc.Alert("No valid accessions found.", color="warning"),
            no_update, no_update, no_update,
        )

    # Warn if mixing multiple PRJNA/SRP accessions (likely different platforms)
    import re
    project_accs = [a for a in valid if re.match(r"^(PRJ|[SED]RP)", a, re.IGNORECASE)]
    if len(project_accs) > 1:
        return (
            no_update, no_update,
            dbc.Alert(
                [
                    html.Strong("Multiple project accessions detected: "),
                    f"{', '.join(project_accs)}. ",
                    "Each project may use a different platform or region. "
                    "Please download one project per job to avoid mixing "
                    "incompatible data types.",
                ],
                color="danger",
            ),
            no_update, no_update, no_update,
        )

    # Create job
    job_id = str(uuid.uuid4())[:8]
    output_dir = SRA_CACHE_DIR / f"job_{job_id}"

    launch_download(
        job_id=job_id,
        accessions=valid,
        output_dir=output_dir,
    )

    job_data = {"job_id": job_id, "output_dir": str(output_dir)}

    return (
        job_data,
        False,  # Enable polling
        dbc.Alert(
            f"Download started for {len(valid)} accession(s)...",
            color="info",
        ),
        {"display": "block"},
        {"display": "none"},  # Hide launch button
        {"display": "inline-block"},  # Show cancel button
    )


# ── Callback: Poll progress ──────────────────────────────────────────────────


@dash_app.callback(
    Output("sra-dl-progress-bar", "value"),
    Output("sra-dl-progress-text", "children"),
    Output("sra-dl-log", "children"),
    Output("sra-dl-status", "children", allow_duplicate=True),
    Output("sra-dl-poll", "disabled", allow_duplicate=True),
    Output("sra-dl-progress-section", "style", allow_duplicate=True),
    Output("sra-dl-btn-launch", "style", allow_duplicate=True),
    Output("sra-dl-btn-cancel", "style", allow_duplicate=True),
    Input("sra-dl-poll", "n_intervals"),
    State("sra-dl-active-job", "data"),
    prevent_initial_call=True,
)
def on_poll(n_intervals, job_data):
    if not job_data:
        return (no_update, no_update, no_update, no_update,
                True, no_update, no_update, no_update)

    job_id = job_data["job_id"]
    output_dir = job_data["output_dir"]
    from pathlib import Path

    status = get_job_status(job_id, Path(output_dir))

    total = status.get("total_accessions", 1) or 1
    completed = status.get("completed", 0)
    current = status.get("current", "")
    pct = int((completed / total) * 100) if total > 0 else 0
    job_status = status.get("status", "pending")

    # Detect dead jobs: status says running but thread is gone AND
    # status.json hasn't been updated recently (>60s stale)
    if job_status in ("downloading", "resolving", "registering"):
        from app.sra.downloader import _download_threads
        thread = _download_threads.get(job_id)
        if thread is None or not thread.is_alive():
            import time
            status_path = Path(output_dir) / "status.json"
            try:
                age = time.time() - status_path.stat().st_mtime
            except OSError:
                age = 999
            if age > 60:
                job_status = "failed"
                status["error"] = status.get("error") or "Download process died unexpectedly"
    log_text = "\n".join(status.get("log", []))

    if job_status == "complete":
        upload_id = status.get("upload_id", "?")
        n_files = status.get("n_files", 0)
        return (
            100,
            "Complete!",
            log_text,
            dbc.Alert(
                [
                    f"Download complete! {n_files} file(s) registered as Upload #{upload_id}. ",
                    html.A("Go to File Manager", href="/files", className="alert-link"),
                    " to view the files.",
                ],
                color="success",
            ),
            True,  # Disable polling
            {"display": "block"},  # Keep log visible
            {"display": "inline-block"},  # Show launch button
            {"display": "none"},  # Hide cancel button
        )

    if job_status == "failed":
        error = status.get("error", "Unknown error")
        return (
            0, "",
            log_text,
            dbc.Alert(f"Download failed: {error}", color="danger"),
            True,
            {"display": "block"},
            {"display": "inline-block"},
            {"display": "none"},
        )

    if job_status == "cancelled":
        return (
            0, "",
            log_text,
            dbc.Alert("Download cancelled.", color="warning"),
            True,
            {"display": "block"},
            {"display": "inline-block"},
            {"display": "none"},
        )

    # Still in progress
    return (
        pct, current, log_text,
        no_update, no_update, no_update, no_update, no_update,
    )


# ── Callback: Cancel download ────────────────────────────────────────────────


@dash_app.callback(
    Output("sra-dl-status", "children", allow_duplicate=True),
    Input("sra-dl-btn-cancel", "n_clicks"),
    State("sra-dl-active-job", "data"),
    prevent_initial_call=True,
)
def on_cancel(n_clicks, job_data):
    if not n_clicks or not job_data:
        return no_update

    cancel_download(job_data["job_id"])
    return dbc.Alert("Cancelling download...", color="warning")


# ── Callback: Restore state on page load ─────────────────────────────────────


@dash_app.callback(
    Output("sra-dl-poll", "disabled", allow_duplicate=True),
    Output("sra-dl-progress-section", "style", allow_duplicate=True),
    Output("sra-dl-btn-launch", "style", allow_duplicate=True),
    Output("sra-dl-btn-cancel", "style", allow_duplicate=True),
    Input("sra-dl-init", "n_intervals"),
    State("sra-dl-active-job", "data"),
    prevent_initial_call=True,
)
def restore_on_load(n, job_data):
    if not job_data:
        return True, {"display": "none"}, {"display": "inline-block"}, {"display": "none"}

    from pathlib import Path
    status = get_job_status(job_data["job_id"], Path(job_data["output_dir"]))
    job_status = status.get("status", "pending")

    if job_status in ("complete", "failed", "cancelled"):
        return True, {"display": "none"}, {"display": "inline-block"}, {"display": "none"}

    # Job still running — re-enable polling
    return (
        False,
        {"display": "block"},
        {"display": "none"},
        {"display": "inline-block"},
    )
