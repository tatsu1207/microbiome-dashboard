"""
MicrobiomeDash — Standalone PICRUSt2 page.

Upload a BIOM table to run PICRUSt2 functional prediction independently
of the main pipeline. Representative sequences are extracted automatically
from the BIOM observation metadata.
"""
import base64
import os
import time
from pathlib import Path

import dash_bootstrap_components as dbc
from dash import ALL, Input, Output, State, ctx, dcc, html, no_update

from app.config import DATASET_DIR, PICRUST2_RUNS_DIR
from app.dashboard.app import app as dash_app

MAX_CPUS = os.cpu_count() or 1


_picrust2_available: bool | None = None


def _is_picrust2_available() -> bool:
    """Check whether the picrust2_16S conda env has picrust2 installed.

    Result is cached after first call to avoid repeated subprocess spawns.
    """
    global _picrust2_available
    if _picrust2_available is not None:
        return _picrust2_available
    import subprocess
    from app.config import PICRUST2_ENV_NAME
    try:
        result = subprocess.run(
            ["conda", "run", "-n", PICRUST2_ENV_NAME, "which", "picrust2_pipeline.py"],
            capture_output=True, text=True, timeout=10,
        )
        _picrust2_available = result.returncode == 0
    except Exception:
        _picrust2_available = False
    return _picrust2_available

# Compute RAM-aware default for hint text
def _default_picrust2_threads() -> int:
    """RAM-aware default: each PICRUSt2 HSP worker uses ~2 GB."""
    try:
        total_ram_gb = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / (1024 ** 3)
    except (ValueError, OSError):
        total_ram_gb = 8
    max_by_ram = max(1, int((total_ram_gb - 4) / 2))
    return max(1, min(max(1, MAX_CPUS - 1), max_by_ram))

_DEFAULT_THREADS = _default_picrust2_threads()


def _get_completed_datasets():
    """Return dropdown options for datasets that have a BIOM table."""
    from app.db.database import get_session
    from app.db.models import Dataset

    with get_session() as db:
        datasets = (
            db.query(Dataset)
            .filter(Dataset.status == "complete", Dataset.asv_table_path.isnot(None))
            .order_by(Dataset.created_at.desc())
            .all()
        )
    return [
        {"label": f"#{d.id} — {d.name}", "value": d.id}
        for d in datasets
    ]


# ── Page Layout ───────────────────────────────────────────────────────────────


def get_layout():
    """Build a fresh layout each time the page is visited.

    Pre-renders the history table server-side so it's always visible,
    and includes a one-shot interval to restore active-run progress.
    """
    picrust2_available = _is_picrust2_available()
    dataset_options = _get_completed_datasets()

    # Show unavailability banner when PICRUSt2 is not installed
    unavailable_banner = []
    if not picrust2_available:
        unavailable_banner = [
            dbc.Alert(
                [
                    html.H5("PICRUSt2 is not available", className="alert-heading"),
                    html.P(
                        "The PICRUSt2 conda environment is not installed. "
                        "This typically happens on ARM64/Apple Silicon systems "
                        "where PICRUSt2 has no compatible package available.",
                        className="mb-2",
                    ),
                    html.P(
                        "All other features of 16S Pipeline work normally. "
                        "PICRUSt2 analysis pages (Pathway Comparison, KEGG Map) "
                        "can still be used if you have pre-computed PICRUSt2 output.",
                        className="mb-0",
                    ),
                ],
                color="warning",
                className="mb-4",
            ),
        ]

    return dbc.Container(
        [
            html.H3("PICRUSt2 Functional Prediction", className="mb-2"),
            html.P(
                "Upload a BIOM table or select a completed DADA2 dataset to run "
                "PICRUSt2. Representative sequences are extracted automatically "
                "from BIOM observation metadata. "
                "Produces EC, KO (optional), and MetaCyc pathway predictions.",
                className="text-muted mb-4",
            ),
            *unavailable_banner,
            # ── DADA2 Dataset Selector ───────────────────────────────────
            dbc.Label("Use DADA2 Dataset", className="fw-bold"),
            dcc.Dropdown(
                id="picrust2-dataset-select",
                options=dataset_options,
                placeholder="Select a completed DADA2 run..." if dataset_options
                else "No completed DADA2 runs available",
                disabled=not dataset_options or not picrust2_available,
                className="mb-2",
            ),
            html.Div(
                "— or —",
                className="text-center text-muted my-2",
            ),
            # ── BIOM Upload ──────────────────────────────────────────────
            dbc.Label("Upload BIOM File", className="fw-bold"),
            dcc.Upload(
                id="picrust2-upload-biom",
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
                disabled=not picrust2_available,
                multiple=False,
                accept=".biom,application/octet-stream,application/x-hdf5",
            ),
            html.Div(id="picrust2-biom-status", className="mt-1 small mb-3"),
            # ── CPU Threads ──────────────────────────────────────────────
            dbc.Label("CPU Threads", className="fw-bold"),
            dbc.InputGroup(
                [
                    dbc.Input(
                        id="picrust2-input-threads",
                        type="text",
                        value="default",
                        placeholder="default",
                        style={"maxWidth": "120px"},
                    ),
                    dbc.InputGroupText(
                        f"default = {_DEFAULT_THREADS} (RAM-aware, max {MAX_CPUS} cores)",
                        className="text-muted small",
                    ),
                ],
                className="mb-3",
            ),
            # ── Options + Launch ─────────────────────────────────────────
            dbc.Checkbox(
                id="picrust2-skip-ko",
                label="Skip KO prediction (faster — EC + MetaCyc only)",
                value=False,
                className="mb-3",
            ),
            dbc.Button(
                "Run PICRUSt2",
                id="picrust2-btn-launch",
                color="warning",
                size="lg",
                disabled=True,
                className="w-100 mb-4",
            ),
            html.Div(id="picrust2-launch-error"),
            # ── Stores ───────────────────────────────────────────────────
            dcc.Store(id="picrust2-biom-data"),
            dcc.Store(id="picrust2-active-run", storage_type="session"),
            dcc.Interval(
                id="picrust2-poll", interval=3000, disabled=True, n_intervals=0
            ),
            # One-shot interval: fires once on page load to restore state
            dcc.Interval(
                id="picrust2-init", interval=200, max_intervals=1, n_intervals=0
            ),
            # ── Progress Section ─────────────────────────────────────────
            html.Div(
                id="picrust2-progress-section",
                children=[
                    html.Div(
                        [
                            html.H5("Progress", className="d-inline"),
                            dbc.Button(
                                "Cancel",
                                id="picrust2-btn-cancel",
                                color="danger",
                                size="sm",
                                className="ms-3",
                            ),
                        ],
                    ),
                    html.Div(id="picrust2-cancel-msg"),
                    html.Div(id="picrust2-status-badge", className="mb-3"),
                    html.Div(
                        id="picrust2-substep-label",
                        className="text-muted small mb-1",
                    ),
                    dbc.Progress(
                        id="picrust2-progress-bar",
                        value=0,
                        striped=True,
                        animated=True,
                        color="warning",
                        className="mb-3",
                        style={"height": "20px"},
                    ),
                    html.H6("Log"),
                    html.Pre(
                        id="picrust2-log-viewer",
                        style={
                            "maxHeight": "300px",
                            "overflowY": "auto",
                            "backgroundColor": "#1a1a2e",
                            "color": "#e0e0e0",
                            "padding": "10px",
                            "borderRadius": "5px",
                            "fontSize": "0.8rem",
                            "whiteSpace": "pre-wrap",
                        },
                    ),
                ],
                style={"display": "none"},
            ),
            # ── Results Section ──────────────────────────────────────────
            html.Div(
                id="picrust2-results-section",
                className="mt-4",
                style={"display": "none"},
            ),
            # ── Download targets ───────────────────────────────────────────
            dcc.Download(id="download-p2-results"),
            # ── History ──────────────────────────────────────────────────
            html.Hr(className="mt-4"),
            html.H5("Run History"),
            html.Div(
                id="picrust2-history-table",
                children=_build_history_table(),
                style={"maxHeight": "300px", "overflowY": "auto"},
            ),
        ],
        fluid=True,
    )


# ── Callbacks ─────────────────────────────────────────────────────────────────


@dash_app.callback(
    Output("picrust2-biom-data", "data"),
    Output("picrust2-biom-status", "children"),
    Output("picrust2-dataset-select", "value"),
    Input("picrust2-upload-biom", "contents"),
    State("picrust2-upload-biom", "filename"),
    prevent_initial_call=True,
)
def on_upload_biom(contents, filename):
    """Validate and store uploaded BIOM file, clear dataset selection."""
    if not contents or not filename:
        return None, "", no_update
    if not filename.endswith(".biom"):
        return None, dbc.Alert("Please upload a .biom file", color="danger", className="py-1"), no_update
    return {"contents": contents, "filename": filename}, html.Span(
        f"  {filename}", className="text-success"
    ), None  # clear dataset dropdown


@dash_app.callback(
    Output("picrust2-biom-data", "data", allow_duplicate=True),
    Output("picrust2-biom-status", "children", allow_duplicate=True),
    Input("picrust2-dataset-select", "value"),
    prevent_initial_call=True,
)
def on_select_dataset(dataset_id):
    """Load BIOM path from a completed DADA2 dataset."""
    if not dataset_id:
        # Return no_update so that a programmatic clear (e.g. from on_upload_biom)
        # does not wipe out biom_data that was just set by the upload callback.
        return no_update, no_update

    from app.db.database import get_session
    from app.db.models import Dataset

    with get_session() as db:
        ds = db.query(Dataset).filter(Dataset.id == dataset_id).first()
        if not ds or not ds.asv_table_path:
            return None, dbc.Alert("Dataset not found or missing BIOM.", color="danger", className="py-1")

        biom_path = Path(ds.asv_table_path)
        if not biom_path.is_absolute():
            biom_path = DATASET_DIR / str(ds.id) / biom_path
        if not biom_path.exists():
            return None, dbc.Alert(
                f"BIOM file not found: {biom_path}", color="danger", className="py-1"
            )

        data = {
            "dataset_id": ds.id,
            "dataset_biom_path": str(biom_path),
            "filename": biom_path.name,
        }
        # Include rep_seqs path if available
        if ds.rep_seqs_path:
            fasta = Path(ds.rep_seqs_path)
            if not fasta.is_absolute():
                fasta = DATASET_DIR / str(ds.id) / fasta
            if fasta.exists():
                data["dataset_fasta_path"] = str(fasta)

        return data, html.Span(
            f"  DADA2 #{ds.id}: {biom_path.name}", className="text-success"
        )


@dash_app.callback(
    Output("picrust2-btn-launch", "disabled"),
    Input("picrust2-biom-data", "data"),
)
def toggle_launch_button(biom_data):
    """Enable launch button when BIOM file is uploaded or dataset selected."""
    return not biom_data


@dash_app.callback(
    Output("picrust2-active-run", "data"),
    Output("picrust2-progress-section", "style"),
    Output("picrust2-poll", "disabled"),
    Output("picrust2-launch-error", "children"),
    Output("picrust2-history-table", "children", allow_duplicate=True),
    Input("picrust2-btn-launch", "n_clicks"),
    State("picrust2-biom-data", "data"),
    State("picrust2-skip-ko", "value"),
    State("picrust2-input-threads", "value"),
    prevent_initial_call=True,
)
def on_launch(n_clicks, biom_data, skip_ko, threads_val):
    """Save uploaded BIOM file and launch PICRUSt2."""
    if not n_clicks or not biom_data:
        return no_update, no_update, no_update, no_update, no_update

    from app.db.database import get_session
    from app.db.models import Dataset, Picrust2Run
    from app.pipeline.runner import launch_picrust2_standalone

    # Parse threads
    threads = None  # None → default in runner
    if threads_val and str(threads_val).strip().lower() != "default":
        try:
            threads = max(1, int(threads_val))
        except ValueError:
            return (
                no_update,
                no_update,
                no_update,
                dbc.Alert(
                    f"Invalid CPU threads value: '{threads_val}'. "
                    "Enter a number or 'default'.",
                    color="danger",
                ),
                no_update,
            )

    try:
        import shutil

        # Create DB record first to get the run ID
        with get_session() as db:
            # Use dataset name if selected from DADA2, otherwise use filename
            if biom_data.get("dataset_id"):
                _ds = db.query(Dataset).filter(Dataset.id == biom_data["dataset_id"]).first()
                run_label = _ds.name if _ds else biom_data["filename"].rsplit(".", 1)[0]
            else:
                run_label = biom_data["filename"].rsplit(".", 1)[0]
            run = Picrust2Run(
                name=run_label,
                status="pending",
                biom_path="",  # will be set below
                skip_ko=bool(skip_ko),
            )
            db.add(run)
            db.flush()
            run_id = run.id

            # Save BIOM file
            run_dir = PICRUST2_RUNS_DIR / str(run_id) / "input"
            run_dir.mkdir(parents=True, exist_ok=True)

            if "dataset_biom_path" in biom_data:
                # Source from completed DADA2 dataset — copy files
                biom_path = run_dir / biom_data["filename"]
                shutil.copy2(biom_data["dataset_biom_path"], biom_path)
                # Copy FASTA if available (avoids re-extracting from BIOM)
                if biom_data.get("dataset_fasta_path"):
                    fasta_dst = run_dir / "rep_seqs.fasta"
                    shutil.copy2(biom_data["dataset_fasta_path"], fasta_dst)
                    run.fasta_path = str(fasta_dst)
            else:
                # Uploaded file — decode base64
                biom_path = run_dir / biom_data["filename"]
                _save_upload(biom_data["contents"], biom_path)

            run.biom_path = str(biom_path)
            db.commit()

        # Launch (FASTA is extracted from BIOM in the runner)
        success = launch_picrust2_standalone(run_id, threads=threads)
        if not success:
            return (
                no_update,
                no_update,
                no_update,
                dbc.Alert("Failed to launch PICRUSt2", color="danger"),
                no_update,
            )

        return (
            run_id,
            {"display": "block"},
            False,
            "",
            _build_history_table(),
        )

    except Exception as e:
        return (
            no_update,
            no_update,
            no_update,
            dbc.Alert(f"Error: {e}", color="danger"),
            no_update,
        )


# ── Cancel callback ──────────────────────────────────────────────────────────


@dash_app.callback(
    Output("picrust2-cancel-msg", "children"),
    Output("picrust2-status-badge", "children", allow_duplicate=True),
    Output("picrust2-substep-label", "children", allow_duplicate=True),
    Output("picrust2-progress-bar", "value", allow_duplicate=True),
    Output("picrust2-progress-bar", "label", allow_duplicate=True),
    Output("picrust2-progress-bar", "animated", allow_duplicate=True),
    Output("picrust2-log-viewer", "children", allow_duplicate=True),
    Output("picrust2-poll", "disabled", allow_duplicate=True),
    Output("picrust2-history-table", "children", allow_duplicate=True),
    Output("picrust2-progress-section", "style", allow_duplicate=True),
    Input("picrust2-btn-cancel", "n_clicks"),
    State("picrust2-active-run", "data"),
    prevent_initial_call=True,
)
def on_cancel(n_clicks, run_id):
    """Cancel the running PICRUSt2 job and immediately update the UI."""
    if not n_clicks or not run_id:
        return (no_update,) * 10

    from app.pipeline.runner import cancel_picrust2, get_picrust2_status

    ok = cancel_picrust2(run_id)
    if not ok:
        return (
            dbc.Alert("No running PICRUSt2 job found to cancel.", color="secondary"),
            *((no_update,) * 8),
            _build_history_table(),
        )

    # Brief wait for the thread to register cancellation
    time.sleep(0.5)

    status = get_picrust2_status(run_id)
    badge = dbc.Badge("FAILED", color="danger", className="fs-6 p-2")
    log_tail = status.get("log_tail", "PICRUSt2 cancelled by user.")

    return (
        dbc.Alert("PICRUSt2 cancelled.", color="warning"),
        badge,
        "Cancelled",
        0,
        "Cancelled",
        False,
        log_tail,
        True,  # stop polling
        _build_history_table(),
        {"display": "none"},  # hide progress section
    )


# ── Poll callback ────────────────────────────────────────────────────────────


@dash_app.callback(
    Output("picrust2-status-badge", "children"),
    Output("picrust2-substep-label", "children"),
    Output("picrust2-progress-bar", "value"),
    Output("picrust2-progress-bar", "label"),
    Output("picrust2-log-viewer", "children"),
    Output("picrust2-poll", "disabled", allow_duplicate=True),
    Output("picrust2-results-section", "style"),
    Output("picrust2-results-section", "children"),
    Output("picrust2-history-table", "children", allow_duplicate=True),
    Output("picrust2-btn-cancel", "style"),
    Output("picrust2-progress-section", "style", allow_duplicate=True),
    Input("picrust2-poll", "n_intervals"),
    State("picrust2-active-run", "data"),
    prevent_initial_call=True,
)
def on_poll(n_intervals, run_id):
    """Poll PICRUSt2 status and update the UI."""
    no_results = ({"display": "none"}, "")
    cancel_visible = {"display": "inline-block"}
    cancel_hidden = {"display": "none"}
    if not run_id:
        return ("", "", 0, "", "", True, *no_results, no_update, cancel_hidden, {"display": "none"})

    from app.pipeline.runner import get_picrust2_status

    status = get_picrust2_status(run_id)
    db_status = status["status"]

    badge_color = {
        "pending": "secondary",
        "processing": "primary",
        "complete": "success",
        "failed": "danger",
    }.get(db_status, "secondary")
    badge = dbc.Badge(db_status.upper(), color=badge_color, className="fs-6 p-2")

    substep = status.get("picrust2_substep", "")
    pct = status.get("picrust2_pct", 0)
    log_tail = status.get("log_tail", "")

    finished = db_status in ("complete", "failed")
    stop_poll = finished

    # Show/hide cancel button
    cancel_style = cancel_hidden if finished else cancel_visible

    # Hide progress section when finished
    progress_style = {"display": "none"} if finished else {"display": "block"}

    # Build results section if complete
    if db_status == "complete":
        results_style = {"display": "block"}
        results_children = _build_results(run_id)
    else:
        results_style, results_children = no_results

    return (
        badge,
        substep,
        pct,
        f"{pct}%",
        log_tail,
        stop_poll,
        results_style,
        results_children,
        _build_history_table(),
        cancel_style,
        progress_style,
    )


# ── Restore on page load ────────────────────────────────────────────────────


@dash_app.callback(
    Output("picrust2-progress-section", "style", allow_duplicate=True),
    Output("picrust2-poll", "disabled", allow_duplicate=True),
    Output("picrust2-history-table", "children", allow_duplicate=True),
    Output("picrust2-status-badge", "children", allow_duplicate=True),
    Output("picrust2-substep-label", "children", allow_duplicate=True),
    Output("picrust2-progress-bar", "value", allow_duplicate=True),
    Output("picrust2-progress-bar", "label", allow_duplicate=True),
    Output("picrust2-progress-bar", "animated", allow_duplicate=True),
    Output("picrust2-log-viewer", "children", allow_duplicate=True),
    Output("picrust2-btn-cancel", "style", allow_duplicate=True),
    Output("picrust2-results-section", "style", allow_duplicate=True),
    Output("picrust2-results-section", "children", allow_duplicate=True),
    Input("picrust2-init", "n_intervals"),
    State("picrust2-active-run", "data"),
    prevent_initial_call=True,
)
def restore_on_load(n_intervals, run_id):
    """Restore progress section with full data on page load.

    Triggered by a one-shot interval (fires once 200ms after mount).
    Populates badge, progress bar, substep, and log immediately so the
    user sees current state without waiting for the first poll cycle.
    """
    history = _build_history_table()
    no_restore = (no_update, no_update, history,
                  no_update, no_update, no_update, no_update, no_update,
                  no_update, no_update, no_update, no_update)
    if not run_id:
        return no_restore

    from app.pipeline.runner import get_picrust2_status

    status = get_picrust2_status(run_id)
    db_status = status["status"]

    badge_color = {
        "pending": "secondary",
        "processing": "primary",
        "complete": "success",
        "failed": "danger",
    }.get(db_status, "secondary")
    badge = dbc.Badge(db_status.upper(), color=badge_color, className="fs-6 p-2")
    substep = status.get("picrust2_substep", "")
    pct = status.get("picrust2_pct", 0)
    log_tail = status.get("log_tail", "")

    if db_status in ("pending", "processing"):
        # Still running — show progress and enable polling
        return (
            {"display": "block"}, False, history,
            badge, substep, pct, f"{pct}%", True, log_tail,
            {"display": "inline-block"},  # cancel visible
            {"display": "none"}, "",      # results hidden
        )

    if db_status == "complete":
        # Finished — show results, hide progress
        return (
            {"display": "none"}, True, history,
            badge, substep, 100, "100%", False, log_tail,
            {"display": "none"},                        # cancel hidden
            {"display": "block"}, _build_results(run_id),  # results shown
        )

    # Failed or other terminal — show progress section with final state
    return (
        {"display": "block"}, True, history,
        badge, substep, pct, f"{pct}%", False, log_tail,
        {"display": "none"},          # cancel hidden
        {"display": "none"}, "",      # results hidden
    )


# ── Helpers ───────────────────────────────────────────────────────────────────


def _save_upload(contents_str: str, dest: Path):
    """Decode a Dash Upload contents string and save to disk."""
    # Format: "data:<mime>;base64,<data>"
    _, content_string = contents_str.split(",", 1)
    data = base64.b64decode(content_string)
    dest.write_bytes(data)


def _build_history_table():
    """Build the PICRUSt2 run history table."""
    from app.db.database import get_session
    from app.db.models import Picrust2Run

    with get_session() as db:
        runs = (
            db.query(Picrust2Run)
            .order_by(Picrust2Run.created_at.desc())
            .limit(20)
            .all()
        )

    if not runs:
        return html.P("No PICRUSt2 runs yet.", className="text-muted")

    badge_color = {
        "pending": "secondary",
        "processing": "primary",
        "complete": "success",
        "failed": "danger",
    }

    rows = []
    for r in runs:
        run_dir = PICRUST2_RUNS_DIR / str(r.id) / "picrust2"
        has_output = r.status == "complete" and run_dir.exists()

        # Download button — only for completed runs with output
        dl_cell = html.Td(
            dbc.Button(
                "Download",
                id={"type": "btn-p2-history-dl", "index": r.id},
                color="success",
                size="sm",
                outline=True,
            )
            if has_output else ""
        )

        # Delete button — disabled for active runs
        is_active = r.status in ("pending", "processing")
        del_cell = html.Td(
            dbc.Button(
                "x",
                id={"type": "btn-p2-history-del", "index": r.id},
                color="danger",
                size="sm",
                outline=True,
                disabled=is_active,
                title="Cannot delete a running job" if is_active else "Delete",
            )
        )

        rows.append(
            html.Tr(
                [
                    html.Td(r.id),
                    html.Td(r.name),
                    html.Td(
                        dbc.Badge(
                            r.status.upper(),
                            color=badge_color.get(r.status, "secondary"),
                        )
                    ),
                    html.Td("No KO" if r.skip_ko else "EC + KO"),
                    html.Td(
                        r.created_at.strftime("%Y-%m-%d %H:%M")
                        if r.created_at
                        else ""
                    ),
                    dl_cell,
                    del_cell,
                ]
            )
        )

    return dbc.Table(
        [
            html.Thead(
                html.Tr(
                    [
                        html.Th("ID"),
                        html.Th("Name"),
                        html.Th("Status"),
                        html.Th("Mode"),
                        html.Th("Created"),
                        html.Th("Results"),
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


def _build_results(run_id: int):
    """Build download buttons for completed PICRUSt2 outputs."""
    run_dir = PICRUST2_RUNS_DIR / str(run_id) / "picrust2"
    if not run_dir.exists():
        return html.P("Output directory not found.", className="text-muted")

    # Map output directories to user-friendly labels and icons
    sections = [
        ("EC_metagenome_out", "EC Metagenome"),
        ("KO_metagenome_out", "KO Metagenome"),
        ("pathways_out", "MetaCyc Pathways"),
    ]

    buttons = []
    for dir_name, label in sections:
        path = run_dir / dir_name
        if not path.exists() or not path.is_dir():
            continue
        files = sorted(path.glob("*.tsv.gz")) + sorted(path.glob("*.tsv"))
        for f in files:
            rel = f"/static-data/picrust2_runs/{run_id}/picrust2/{dir_name}/{f.name}"
            buttons.append(
                html.A(
                    dbc.Button(
                        f"{label}: {f.name}",
                        color="info",
                        size="sm",
                        className="me-2 mb-2",
                    ),
                    href=rel,
                    download=f.name,
                )
            )

    if not buttons:
        return html.P("No output files found.", className="text-muted")

    return html.Div(
        [
            html.H5("Download Results"),
            html.Div(buttons),
        ]
    )


# ── History download callback ────────────────────────────────────────────


@dash_app.callback(
    Output("download-p2-results", "data"),
    Input({"type": "btn-p2-history-dl", "index": ALL}, "n_clicks"),
    prevent_initial_call=True,
)
def on_p2_history_download(n_clicks_list):
    """Download key PICRUSt2 prediction files as a zip."""
    if not any(n_clicks_list):
        return no_update

    triggered = ctx.triggered_id
    if not triggered:
        return no_update
    run_id = triggered["index"]

    import tempfile
    import zipfile

    from app.db.database import get_session
    from app.db.models import Picrust2Run

    run_dir = PICRUST2_RUNS_DIR / str(run_id) / "picrust2"
    if not run_dir.exists():
        return no_update

    # Look up the run name for the download filename
    with get_session() as db:
        run = db.query(Picrust2Run).filter(Picrust2Run.id == run_id).first()
        run_name = run.name if run else f"run{run_id}"

    # Collect the key prediction files (prefer described versions)
    # (subdir, described_name, fallback_name, zip_prefix)
    targets = [
        ("pathways_out", "path_abun_unstrat_described.tsv.gz", "path_abun_unstrat.tsv.gz", "metacyc"),
        ("EC_metagenome_out", "pred_metagenome_unstrat_described.tsv.gz", "pred_metagenome_unstrat.tsv.gz", "EC"),
        ("KO_metagenome_out", "pred_metagenome_unstrat_described.tsv.gz", "pred_metagenome_unstrat.tsv.gz", "KO"),
    ]

    tmp = tempfile.mkdtemp()
    zip_path = Path(tmp) / f"picrust2_run{run_id}.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for subdir, described, fallback, prefix in targets:
            fpath = run_dir / subdir / described
            if not fpath.exists():
                fpath = run_dir / subdir / fallback
            if fpath.exists():
                zf.write(fpath, arcname=f"{prefix}_{fpath.name}")

    if zip_path.stat().st_size <= 22:  # empty zip
        return no_update

    return dcc.send_file(str(zip_path), filename=f"{run_name}_picrust2.zip")


# ── History delete callback ──────────────────────────────────────────────


@dash_app.callback(
    Output("picrust2-history-table", "children", allow_duplicate=True),
    Input({"type": "btn-p2-history-del", "index": ALL}, "n_clicks"),
    prevent_initial_call=True,
)
def on_p2_history_delete(n_clicks_list):
    """Delete a PICRUSt2 run: remove DB record and files on disk."""
    if not any(n_clicks_list):
        return no_update

    triggered = ctx.triggered_id
    if not triggered:
        return no_update
    run_id = triggered["index"]

    import shutil

    from app.db.database import get_session
    from app.db.models import Picrust2Run

    with get_session() as db:
        run = db.query(Picrust2Run).filter(Picrust2Run.id == run_id).first()
        if not run:
            return _build_history_table()

        # Don't delete active jobs
        if run.status in ("pending", "processing"):
            return no_update

        db.delete(run)
        db.commit()

    # Remove run directory from disk
    run_path = PICRUST2_RUNS_DIR / str(run_id)
    if run_path.exists():
        shutil.rmtree(run_path, ignore_errors=True)

    return _build_history_table()
