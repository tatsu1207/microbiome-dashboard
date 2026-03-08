"""
MicrobiomeDash — Pipeline orchestrator.

Launches the 5-step pipeline in a background thread and tracks progress
via a status.json file in the dataset output directory.
"""
import json
import logging
import os
import signal
import subprocess
import threading
from datetime import datetime
from pathlib import Path

from app.config import DADA2_DEFAULTS, DATASET_DIR, LONGREAD_DADA2_DEFAULTS, PICRUST2_RUNS_DIR

# Active pipeline threads keyed by dataset_id (or "picrust2-{run_id}")
_running_pipelines: dict[int | str, threading.Thread] = {}

# Cancellation infrastructure — DADA2 pipeline
_cancel_events: dict[int, threading.Event] = {}
_active_procs: dict[int, subprocess.Popen] = {}

# Cancellation infrastructure — PICRUSt2
_picrust2_cancel_events: dict[int, threading.Event] = {}
_picrust2_active_procs: dict[int, subprocess.Popen] = {}


class PipelineCancelled(Exception):
    """Raised when a pipeline is cancelled by the user."""
    pass


def cancel_pipeline(dataset_id: int) -> bool:
    """Cancel a running pipeline.

    Sets the cancel event and kills the active subprocess (if any).
    Returns True if a running pipeline was found and signalled.
    """
    event = _cancel_events.get(dataset_id)
    if not event:
        return False

    event.set()

    # Kill active subprocess and its entire process group.
    # Do this in a short-lived thread so the caller (UI) isn't blocked.
    proc = _active_procs.get(dataset_id)
    if proc and proc.poll() is None:
        def _kill():
            try:
                os.killpg(proc.pid, signal.SIGTERM)
            except (OSError, ProcessLookupError):
                pass
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(proc.pid, signal.SIGKILL)
                except (OSError, ProcessLookupError):
                    pass

        threading.Thread(target=_kill, daemon=True).start()

    return True


def cancel_picrust2(run_id: int) -> bool:
    """Cancel a running PICRUSt2 run.

    Sets the cancel event and kills the active subprocess.
    Returns True if a running job was found and signalled.
    """
    event = _picrust2_cancel_events.get(run_id)
    if not event:
        return False

    event.set()

    proc = _picrust2_active_procs.get(run_id)
    if proc and proc.poll() is None:
        def _kill():
            try:
                os.killpg(proc.pid, signal.SIGTERM)
            except (OSError, ProcessLookupError):
                pass
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(proc.pid, signal.SIGKILL)
                except (OSError, ProcessLookupError):
                    pass

        threading.Thread(target=_kill, daemon=True).start()

    return True


def _check_cancel(dataset_id: int):
    """Raise PipelineCancelled if cancellation was requested."""
    event = _cancel_events.get(dataset_id)
    if event and event.is_set():
        raise PipelineCancelled()


def _register_proc(dataset_id: int, proc: subprocess.Popen):
    """Track the active subprocess for a pipeline so it can be killed on cancel."""
    _active_procs[dataset_id] = proc


def _is_pid_alive(pid: int) -> bool:
    """Check if a process with the given PID is still running."""
    if not pid:
        return False
    try:
        os.kill(pid, 0)
        return True
    except (OSError, ProcessLookupError):
        return False


def _save_pid_to_status(status_file: Path, pid: int):
    """Save a subprocess PID into status.json for recovery after server restart."""
    try:
        data = json.loads(status_file.read_text()) if status_file.exists() else {}
        data["pid"] = pid
        status_file.write_text(json.dumps(data, indent=2))
    except Exception:
        pass


def _reattach_picrust2(run_id: int, pid: int):
    """Spawn a lightweight thread that waits for an orphaned PICRUSt2 subprocess
    to finish and updates the DB status accordingly."""
    from app.pipeline.picrust2 import _monitor_progress

    run_dir = PICRUST2_RUNS_DIR / str(run_id)

    def _waiter():
        from app.db.database import SessionLocal
        from app.db.models import Picrust2Run

        logger = logging.getLogger(f"picrust2-{run_id}-reattach")
        logger.setLevel(logging.INFO)
        logger.handlers.clear()
        fh = logging.FileHandler(run_dir / "pipeline.log")
        fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
        logger.addHandler(fh)
        logger.info(f"Re-attached to PICRUSt2 subprocess (PID {pid})")

        # Read skip_ko from DB
        db = SessionLocal()
        try:
            run = db.query(Picrust2Run).filter(Picrust2Run.id == run_id).first()
            skip_ko = run.skip_ko if run else False
        finally:
            db.close()

        def _progress(substep, pct):
            _update_status(run_dir, "picrust2", pct, [],
                           picrust2_substep=substep, picrust2_pct=pct, pid=pid)

        # Re-start the progress monitor
        stop_event = threading.Event()
        p2_out = run_dir / "picrust2"
        monitor = threading.Thread(
            target=_monitor_progress,
            args=(p2_out, skip_ko, _progress, logger, stop_event),
            daemon=True,
        )
        monitor.start()

        # Wait for the process to finish
        try:
            _, exit_code = os.waitpid(pid, 0)
            returncode = os.waitstatus_to_exitcode(exit_code)
        except ChildProcessError:
            # Not our child — poll until it's gone
            import time
            while _is_pid_alive(pid):
                time.sleep(3)
            returncode = 0  # Assume success if we can't get exit code

        stop_event.set()
        monitor.join(timeout=5)

        # Update DB
        db = SessionLocal()
        try:
            run = db.query(Picrust2Run).filter(Picrust2Run.id == run_id).first()
            if run:
                p2_dir = run_dir / "picrust2"
                # Check for key output files to determine success
                has_pathways = (p2_dir / "pathways_out").exists()
                has_ec = (p2_dir / "EC_metagenome_out").exists()
                if returncode == 0 and (has_pathways or has_ec):
                    run.status = "complete"
                    run.output_dir = str(p2_dir)
                    logger.info("PICRUSt2 completed successfully (re-attached).")
                    _update_status(run_dir, "complete", 100, ["picrust2"],
                                   picrust2_substep="Complete", picrust2_pct=100)
                else:
                    run.status = "failed"
                    logger.error(f"PICRUSt2 failed (exit code {returncode}).")
                    _update_status(run_dir, "failed", 0, [],
                                   picrust2_substep="Failed", picrust2_pct=0)
                db.commit()
        finally:
            db.close()

        _running_pipelines.pop(f"picrust2-{run_id}", None)
        logger.handlers.clear()

    key = f"picrust2-{run_id}"
    t = threading.Thread(target=_waiter, daemon=True, name=f"reattach-{key}")
    _running_pipelines[key] = t
    t.start()


def _reattach_pipeline(dataset_id: int, pid: int):
    """Spawn a lightweight thread that waits for an orphaned pipeline subprocess
    to finish and updates the DB status accordingly."""
    output_dir = DATASET_DIR / str(dataset_id)

    def _waiter():
        from app.db.database import SessionLocal
        from app.db.models import Dataset

        logger = logging.getLogger(f"pipeline-{dataset_id}-reattach")
        logger.setLevel(logging.INFO)
        logger.handlers.clear()
        fh = logging.FileHandler(output_dir / "pipeline.log")
        fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
        logger.addHandler(fh)
        logger.info(f"Re-attached to pipeline subprocess (PID {pid})")

        # Wait for the process to finish
        try:
            _, exit_code = os.waitpid(pid, 0)
        except ChildProcessError:
            import time
            while _is_pid_alive(pid):
                time.sleep(3)

        # Update DB based on output files
        db = SessionLocal()
        try:
            dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
            if dataset and dataset.status == "processing":
                if dataset.asv_table_path and Path(dataset.asv_table_path).exists():
                    dataset.status = "complete"
                    logger.info("Pipeline completed successfully (re-attached).")
                else:
                    dataset.status = "failed"
                    logger.error("Pipeline failed (re-attached, no ASV table found).")
                db.commit()
        finally:
            db.close()

        _running_pipelines.pop(dataset_id, None)
        logger.handlers.clear()

    t = threading.Thread(target=_waiter, daemon=True, name=f"reattach-pipeline-{dataset_id}")
    _running_pipelines[dataset_id] = t
    t.start()


STEPS = ["fastqc", "cutadapt", "dada2", "taxonomy"]
STEPS_LONGREAD = ["fastqc", "dada2_longread", "taxonomy"]


def launch_pipeline(
    file_ids: list[int],
    dataset_name: str,
    description: str | None = None,
    project_id: int | None = None,
    trim_left_f: int | None = None,
    trim_left_r: int | None = None,
    trunc_len_f: int | None = None,
    trunc_len_r: int | None = None,
    min_overlap: int | None = None,
    threads: int | None = None,
    custom_fwd_primer: str | None = None,
    custom_rev_primer: str | None = None,
) -> int:
    """Create a Dataset record and launch the pipeline in a background thread.

    Accepts a list of FastqFile IDs (potentially spanning multiple uploads).
    Returns the new dataset_id.
    """
    from app.db.database import SessionLocal
    from app.db.models import Dataset, DatasetFastqFile, FastqFile, Upload

    if not file_ids:
        raise ValueError("No files selected")

    db = SessionLocal()
    try:
        fastq_files = db.query(FastqFile).filter(FastqFile.id.in_(file_ids)).all()
        if not fastq_files:
            raise ValueError("No matching FASTQ files found")

        # Determine upload_id: set it if all files come from one upload, else None
        upload_ids = list({f.upload_id for f in fastq_files})
        single_upload_id = upload_ids[0] if len(upload_ids) == 1 else None

        # Determine sequencing_type and variable_region from the associated uploads
        uploads = db.query(Upload).filter(Upload.id.in_(upload_ids)).all()
        upload_map = {u.id: u for u in uploads}

        seq_types = {u.sequencing_type for u in uploads if u.sequencing_type}
        regions = {u.variable_region for u in uploads if u.variable_region}
        platforms = {u.platform for u in uploads if u.platform}

        sequencing_type = seq_types.pop() if len(seq_types) == 1 else (
            "mixed" if len(seq_types) > 1 else None
        )
        variable_region = regions.pop() if len(regions) == 1 else (
            "mixed" if len(regions) > 1 else None
        )
        platform = platforms.pop() if len(platforms) == 1 else None

        # Derive project_id from uploads if not provided
        if not project_id:
            pids = {u.project_id for u in uploads if u.project_id}
            project_id = pids.pop() if len(pids) == 1 else None

        # Count unique samples
        sample_names = list({f.sample_name for f in fastq_files})
        n_samples = len(sample_names)

        # Apply defaults for unset parameters
        tlf = trim_left_f if trim_left_f is not None else DADA2_DEFAULTS["trim_left_f"]
        tlr = trim_left_r if trim_left_r is not None else DADA2_DEFAULTS["trim_left_r"]
        trf = trunc_len_f if trunc_len_f is not None else DADA2_DEFAULTS["trunc_len_f"]
        trr = trunc_len_r if trunc_len_r is not None else DADA2_DEFAULTS["trunc_len_r"]
        mo = min_overlap if min_overlap is not None else DADA2_DEFAULTS["min_overlap"]

        dataset = Dataset(
            project_id=project_id,
            upload_id=single_upload_id,
            name=dataset_name,
            description=description or "",
            source_type="pipeline",
            status="pending",
            variable_region=variable_region,
            sequencing_type=sequencing_type,
            platform=platform,
            trim_left_f=tlf,
            trim_left_r=tlr,
            trunc_len_f=trf,
            trunc_len_r=trr,
            min_overlap=mo,
            custom_fwd_primer=custom_fwd_primer,
            custom_rev_primer=custom_rev_primer,
        )
        db.add(dataset)
        db.flush()

        dataset_id = dataset.id

        # Create DatasetFastqFile association records
        for fid in file_ids:
            db.add(DatasetFastqFile(dataset_id=dataset_id, fastq_file_id=fid))

        # Create output directory
        output_dir = DATASET_DIR / str(dataset_id)
        output_dir.mkdir(parents=True, exist_ok=True)

        dataset.pipeline_log_path = str(output_dir / "pipeline.log")
        db.commit()
    finally:
        db.close()

    # Set up cancellation event
    _cancel_events[dataset_id] = threading.Event()

    # Resolve thread count: None → samples × 2 (capped at CPUs - 1)
    if threads is not None:
        resolved_threads = threads
    else:
        max_cpus = max(1, (os.cpu_count() or 1) - 1)
        resolved_threads = min(max(1, n_samples * 2), max_cpus)

    # Launch background thread
    t = threading.Thread(
        target=_run_pipeline,
        args=(dataset_id, resolved_threads),
        daemon=True,
        name=f"pipeline-{dataset_id}",
    )
    _running_pipelines[dataset_id] = t
    t.start()

    return dataset_id


# ── Standalone PICRUSt2 ──────────────────────────────────────────────────────


def launch_picrust2_standalone(run_id: int, threads: int | None = None) -> bool:
    """Launch a standalone PICRUSt2 run in a background thread.

    Returns True if launched successfully.
    """
    from app.db.database import SessionLocal
    from app.db.models import Picrust2Run

    db = SessionLocal()
    try:
        run = db.query(Picrust2Run).filter(Picrust2Run.id == run_id).first()
        if not run or run.status != "pending":
            return False

        run.status = "processing"
        db.commit()
    finally:
        db.close()

    # Set up cancellation event
    _picrust2_cancel_events[run_id] = threading.Event()

    t = threading.Thread(
        target=_run_picrust2_standalone,
        args=(run_id, threads),
        daemon=True,
        name=f"picrust2-{run_id}",
    )
    key = f"picrust2-{run_id}"
    _running_pipelines[key] = t
    t.start()
    return True


def _run_picrust2_standalone(run_id: int, threads: int | None = None):
    """Execute standalone PICRUSt2 in a background thread."""
    from app.db.database import SessionLocal
    from app.db.models import Picrust2Run
    from app.pipeline.biom_convert import extract_fasta_from_biom
    from app.pipeline.picrust2 import run_picrust2

    if threads is None:
        threads = DADA2_DEFAULTS["threads"]

    run_dir = PICRUST2_RUNS_DIR / str(run_id)
    run_dir.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger(f"picrust2-{run_id}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fh = logging.FileHandler(run_dir / "pipeline.log")
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logger.addHandler(fh)

    was_cancelled = False

    db = SessionLocal()
    try:
        run = db.query(Picrust2Run).filter(Picrust2Run.id == run_id).first()

        def _progress(substep, pct):
            _update_status(run_dir, "picrust2", pct, [],
                           picrust2_substep=substep, picrust2_pct=pct)

        _update_status(run_dir, "picrust2", 0, [],
                       picrust2_substep="Extracting FASTA", picrust2_pct=0)

        # Extract FASTA from BIOM observation metadata
        if not run.fasta_path:
            input_dir = run_dir / "input"
            input_dir.mkdir(parents=True, exist_ok=True)
            fasta_path = input_dir / "rep_seqs.fasta"
            logger.info("Extracting representative sequences from BIOM metadata...")
            extract_fasta_from_biom(Path(run.biom_path), fasta_path, logger)
            run.fasta_path = str(fasta_path)
            db.commit()

        _update_status(run_dir, "picrust2", 0, [],
                       picrust2_substep="Starting", picrust2_pct=0)

        if run.skip_ko:
            logger.info("Running PICRUSt2 (EC + MetaCyc only, skipping KO)...")
        else:
            logger.info("Running PICRUSt2 (full: EC + KO + MetaCyc)...")

        logger.info(f"Using {threads} threads")

        result = run_picrust2(
            asv_table_path=Path(run.biom_path),  # already BIOM
            rep_seqs_path=Path(run.fasta_path),
            output_dir=run_dir,
            threads=threads,
            logger=logger,
            skip_ko=run.skip_ko,
            progress_callback=_progress,
            biom_ready=True,  # skip TSV→BIOM conversion
            proc_callback=lambda proc: (
                _picrust2_active_procs.__setitem__(run_id, proc),
                _save_pid_to_status(run_dir / "status.json", proc.pid),
            ),
        )

        # Check if cancelled while subprocess was running
        cancel_event = _picrust2_cancel_events.get(run_id)
        if cancel_event and cancel_event.is_set():
            raise PipelineCancelled()

        if result["success"]:
            run.output_dir = result["picrust_dir"]
            run.status = "complete"
            logger.info("PICRUSt2 completed successfully!")
            _update_status(run_dir, "complete", 100, ["picrust2"],
                           picrust2_substep="Complete", picrust2_pct=100)
        else:
            run.status = "failed"
            logger.error(f"PICRUSt2 failed: {result['error']}")
            _update_status(run_dir, "failed", 0, [],
                           picrust2_substep="Failed", picrust2_pct=0)
        db.commit()

    except PipelineCancelled:
        was_cancelled = True
        logger.info("PICRUSt2 cancelled by user.")
        try:
            run = db.query(Picrust2Run).filter(Picrust2Run.id == run_id).first()
            if run:
                run.status = "failed"
                db.commit()
        except Exception:
            pass
        _update_status(run_dir, "failed", 0, [],
                       picrust2_substep="Cancelled", picrust2_pct=0)

    except Exception as e:
        # If cancel event is set, treat subprocess errors as cancellation
        cancel_event = _picrust2_cancel_events.get(run_id)
        if cancel_event and cancel_event.is_set():
            was_cancelled = True
            logger.info("PICRUSt2 cancelled by user.")
            try:
                run = db.query(Picrust2Run).filter(Picrust2Run.id == run_id).first()
                if run:
                    run.status = "failed"
                    db.commit()
            except Exception:
                pass
            _update_status(run_dir, "failed", 0, [],
                           picrust2_substep="Cancelled", picrust2_pct=0)
        else:
            logger.error(f"PICRUSt2 failed: {e}", exc_info=True)
            try:
                run = db.query(Picrust2Run).filter(Picrust2Run.id == run_id).first()
                if run:
                    run.status = "failed"
                    db.commit()
            except Exception:
                pass
            _update_status(run_dir, "failed", 0, [],
                           picrust2_substep="Failed", picrust2_pct=0)
    finally:
        db.close()
        _running_pipelines.pop(f"picrust2-{run_id}", None)
        _picrust2_cancel_events.pop(run_id, None)
        _picrust2_active_procs.pop(run_id, None)
        logger.handlers.clear()


def get_picrust2_status(run_id: int) -> dict:
    """Get the current status of a standalone PICRUSt2 run."""
    from app.db.database import SessionLocal
    from app.db.models import Picrust2Run

    run_dir = PICRUST2_RUNS_DIR / str(run_id)
    status_file = run_dir / "status.json"

    step_info = {}
    if status_file.exists():
        try:
            step_info = json.loads(status_file.read_text())
        except Exception:
            pass

    log_tail = ""
    log_file = run_dir / "pipeline.log"
    if log_file.exists():
        try:
            lines = log_file.read_text().splitlines()
            log_tail = "\n".join(lines[-50:])
        except Exception:
            pass

    db = SessionLocal()
    try:
        run = db.query(Picrust2Run).filter(Picrust2Run.id == run_id).first()
        db_status = run.status if run else "unknown"
    finally:
        db.close()

    # Check for dead thread — but first check if the subprocess is still alive
    key = f"picrust2-{run_id}"
    thread = _running_pipelines.get(key)
    thread_alive = thread is not None and thread.is_alive()

    if db_status == "processing" and not thread_alive:
        # Check if the subprocess PID is still alive (survives server restart)
        pid = step_info.get("pid")
        if pid and _is_pid_alive(pid):
            # Subprocess still running — re-attach a monitor thread
            _reattach_picrust2(run_id, pid)
        else:
            # Both thread and subprocess are dead — mark as failed
            db = SessionLocal()
            try:
                run = db.query(Picrust2Run).filter(Picrust2Run.id == run_id).first()
                if run and run.status == "processing":
                    run.status = "failed"
                    db.commit()
                    db_status = "failed"
            finally:
                db.close()

    return {
        "run_id": run_id,
        "status": db_status,
        "picrust2_substep": step_info.get("picrust2_substep", ""),
        "picrust2_pct": step_info.get("picrust2_pct", 0),
        "log_tail": log_tail,
    }


# ── Pipeline status helpers ──────────────────────────────────────────────────


def get_pipeline_status(dataset_id: int) -> dict:
    """Get the current pipeline status for a dataset."""
    from app.db.database import SessionLocal
    from app.db.models import Dataset

    output_dir = DATASET_DIR / str(dataset_id)
    status_file = output_dir / "status.json"

    # Read status.json
    step_info = {
        "current_step": None,
        "progress_pct": 0,
        "steps_completed": [],
        "started_at": None,
    }
    if status_file.exists():
        try:
            step_info = json.loads(status_file.read_text())
        except Exception:
            pass

    # Read log tail
    log_tail = ""
    log_file = output_dir / "pipeline.log"
    if log_file.exists():
        try:
            lines = log_file.read_text().splitlines()
            log_tail = "\n".join(lines[-50:])
        except Exception:
            pass

    # Get DB status
    db = SessionLocal()
    try:
        dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
        db_status = dataset.status if dataset else "unknown"
    finally:
        db.close()

    # Check for dead thread — but first check if the subprocess is still alive
    thread = _running_pipelines.get(dataset_id)
    thread_alive = thread is not None and thread.is_alive()

    if db_status == "processing" and not thread_alive:
        # Check if the subprocess PID is still alive (survives server restart)
        pid = step_info.get("pid")
        if pid and _is_pid_alive(pid):
            # Subprocess still running — re-attach a monitor thread
            _reattach_pipeline(dataset_id, pid)
        else:
            # Both thread and subprocess are dead — mark as failed or complete
            db = SessionLocal()
            try:
                dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
                if dataset and dataset.status == "processing":
                    if dataset.asv_table_path and Path(dataset.asv_table_path).exists():
                        dataset.status = "complete"
                        db_status = "complete"
                    else:
                        dataset.status = "failed"
                        db_status = "failed"
                    db.commit()
            finally:
                db.close()

    return {
        "dataset_id": dataset_id,
        "status": db_status,
        "current_step": step_info.get("current_step"),
        "progress_pct": step_info.get("progress_pct", 0),
        "steps_completed": step_info.get("steps_completed", []),
        "started_at": step_info.get("started_at"),
        "log_tail": log_tail,
        "thread_alive": thread_alive,
        "auto_trunc_len_f": step_info.get("auto_trunc_len_f"),
        "auto_trunc_len_r": step_info.get("auto_trunc_len_r"),
        "auto_trunc_details": step_info.get("auto_trunc_details"),
    }


def _update_status(output_dir: Path, step: str, pct: int, completed: list[str], **extra):
    """Write the status.json file."""
    status_file = output_dir / "status.json"
    data = {
        "current_step": step,
        "progress_pct": pct,
        "steps_completed": completed,
        "updated_at": datetime.utcnow().isoformat(),
    }
    if status_file.exists():
        try:
            existing = json.loads(status_file.read_text())
            data["started_at"] = existing.get("started_at")
            # Preserve extra fields from previous writes
            for key in ("auto_trunc_len_f", "auto_trunc_len_r", "auto_trunc_details"):
                if key in existing:
                    data[key] = existing[key]
        except Exception:
            pass
    else:
        data["started_at"] = datetime.utcnow().isoformat()

    data.update(extra)
    status_file.write_text(json.dumps(data, indent=2))


def _run_pipeline(dataset_id: int, threads: int | None = None):
    """Execute the full pipeline in a background thread.

    Creates its own DB session for thread safety.
    """
    from app.db.database import SessionLocal
    from app.db.models import Dataset, FastqFile, QcMetric, Sample, Upload

    output_dir = DATASET_DIR / str(dataset_id)
    completed_steps: list[str] = []
    was_cancelled = False

    # Set up logger
    log_file = output_dir / "pipeline.log"
    logger = logging.getLogger(f"pipeline-{dataset_id}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fh = logging.FileHandler(log_file)
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logger.addHandler(fh)

    db = SessionLocal()
    try:
        dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()

        dataset.status = "processing"
        db.commit()

        # Get selected files via DatasetFastqFile association
        from app.db.models import DatasetFastqFile

        dsff_rows = (
            db.query(DatasetFastqFile)
            .filter(DatasetFastqFile.dataset_id == dataset_id)
            .all()
        )
        selected_file_ids = [r.fastq_file_id for r in dsff_rows]
        fastq_files = (
            db.query(FastqFile).filter(FastqFile.id.in_(selected_file_ids)).all()
        )

        # Create a working directory with symlinks to the selected files
        fastq_dir = output_dir / "fastq_input"
        fastq_dir.mkdir(parents=True, exist_ok=True)
        for ff in fastq_files:
            src = Path(ff.file_path)
            link = fastq_dir / ff.filename
            if not link.exists() and src.exists():
                link.symlink_to(src)

        seq_type = dataset.sequencing_type or "single-end"
        var_region = dataset.variable_region
        platform = dataset.platform
        custom_fwd = dataset.custom_fwd_primer
        custom_rev = dataset.custom_rev_primer

        # threads should already be resolved by launch_pipeline; fallback just in case
        if threads is None:
            sample_names = list({ff.sample_name for ff in fastq_files})
            n_samples = len(sample_names)
            max_cpus = max(1, (os.cpu_count() or 1) - 1)
            threads = min(max(1, n_samples * 2), max_cpus)

        # Collect all upload IDs linked to the selected files
        upload_ids = list({ff.upload_id for ff in fastq_files})
        uploads = db.query(Upload).filter(Upload.id.in_(upload_ids)).all()

        logger.info(f"Pipeline started for dataset {dataset_id}")
        logger.info(f"  Files: {len(fastq_files)} from {len(uploads)} upload(s)")
        logger.info(f"  Type: {seq_type}, Region: {var_region}, Platform: {platform or 'illumina'}, Threads: {threads}")
        _update_status(output_dir, "fastqc", 0, completed_steps)

        # ── Step 1: FastQC ─────────────────────────────────────────────

        from app.pipeline.qc import run_fastqc

        _check_cancel(dataset_id)
        _update_status(output_dir, "fastqc", 5, completed_steps)
        logger.info("Running FastQC...")
        qc_result = run_fastqc(fastq_dir, output_dir, logger)
        completed_steps.append("fastqc")

        # ── Branch: Long-read (PacBio) vs Short-read (Illumina) ──────

        if platform in ("pacbio", "nanopore"):
            # Long-read branch: DADA2 with learned error model (no Cutadapt, no truncation)
            from app.pipeline.dada2 import run_dada2_longread
            from app.pipeline.detect import PRIMERS

            _check_cancel(dataset_id)
            _update_status(output_dir, "dada2_longread", 25, completed_steps)

            # Get primers for full-length 16S
            primer_region = var_region or "V1-V9"
            primer_info = PRIMERS.get(primer_region, PRIMERS["V1-V9"])
            fwd_primer = custom_fwd or primer_info["forward"]
            rev_primer = custom_rev or primer_info["reverse"]

            dada2_result = run_dada2_longread(
                input_dir=fastq_dir,
                output_dir=output_dir,
                fwd_primer=fwd_primer,
                rev_primer=rev_primer,
                platform=platform,
                threads=threads,
                logger=logger,
                proc_callback=lambda proc: (
                    _register_proc(dataset_id, proc),
                    _save_pid_to_status(DATASET_DIR / str(dataset_id) / "status.json", proc.pid),
                ),
            )

            _check_cancel(dataset_id)

        else:
            # Short-read branch: Cutadapt → auto-truncation → DADA2

            # ── Step 2: Cutadapt ───────────────────────────────────────

            from app.pipeline.trim import run_cutadapt

            _check_cancel(dataset_id)
            _update_status(output_dir, "cutadapt", 20, completed_steps)
            trim_result = run_cutadapt(
                fastq_dir, output_dir, seq_type, var_region, logger, threads,
                custom_fwd_primer=custom_fwd, custom_rev_primer=custom_rev,
            )
            trimmed_dir = Path(trim_result["trimmed_dir"])
            completed_steps.append("cutadapt")

            # ── Auto-detect truncation params if needed ──────────────

            _check_cancel(dataset_id)

            trunc_f = dataset.trunc_len_f
            trunc_r = dataset.trunc_len_r

            from app.pipeline.quality import detect_truncation_params

            auto = detect_truncation_params(
                trimmed_dir=trimmed_dir,
                sequencing_type=seq_type,
                variable_region=var_region,
                logger=logger,
            )

            # Always apply auto-detected trim_left for leading N bases
            if dataset.trim_left_f == 0 and auto.get("trim_left_f"):
                dataset.trim_left_f = auto["trim_left_f"]
                logger.info(f"Auto trim_left_f={auto['trim_left_f']} (leading N)")
            if dataset.trim_left_r == 0 and auto.get("trim_left_r"):
                dataset.trim_left_r = auto["trim_left_r"]
                logger.info(f"Auto trim_left_r={auto['trim_left_r']} (leading N)")

            # Auto-detect truncation if not manually set
            if trunc_f == 0 and auto["trunc_len_f"]:
                trunc_f = auto["trunc_len_f"]
                dataset.trunc_len_f = trunc_f
            if trunc_r == 0 and auto["trunc_len_r"]:
                trunc_r = auto["trunc_len_r"]
                dataset.trunc_len_r = trunc_r

            db.commit()
            logger.info(f"Auto-detected: {auto['details']}")
            _update_status(
                output_dir, "dada2", 35, completed_steps,
                auto_trunc_len_f=trunc_f,
                auto_trunc_len_r=trunc_r,
                auto_trunc_details=auto["details"],
            )

            # ── Pre-DADA2: trim leading bases and fix N bases ──────────
            # DADA2's maxN=0 rejects any read containing N. Two fixes:
            # 1. Trim leading bases (maxN applies before trimLeft).
            # 2. Replace remaining N bases with A at Q0 so maxEE
            #    penalizes them without outright rejection.
            tlf = dataset.trim_left_f
            tlr = dataset.trim_left_r
            from app.pipeline.trim import trim_leading_bases
            logger.info(f"Pre-DADA2 processing: trim_left F={tlf}/R={tlr}, fix N bases")
            trimmed_dir = trim_leading_bases(
                trimmed_dir, tlf, tlr, seq_type, logger, threads=threads,
            )
            # Reset trim_left since we've already trimmed
            tlf = 0
            tlr = 0

            # ── Step 3: DADA2 ──────────────────────────────────────────

            from app.pipeline.dada2 import run_dada2

            _check_cancel(dataset_id)
            _update_status(output_dir, "dada2", 40, completed_steps)
            dada2_result = run_dada2(
                input_dir=trimmed_dir,
                output_dir=output_dir,
                sequencing_type=seq_type,
                trim_left_f=tlf,
                trim_left_r=tlr,
                trunc_len_f=trunc_f,
                trunc_len_r=trunc_r,
                min_overlap=dataset.min_overlap,
                threads=threads,
                logger=logger,
                proc_callback=lambda proc: (
                    _register_proc(dataset_id, proc),
                    _save_pid_to_status(DATASET_DIR / str(dataset_id) / "status.json", proc.pid),
                ),
            )

            _check_cancel(dataset_id)

        # Generate rep_seqs.fasta from the ASV table (sequences are in the TSV)
        asv_tsv = Path(dada2_result["asv_table_path"])
        rep_seqs_path = output_dir / "rep_seqs.fasta"
        _extract_rep_seqs(asv_tsv, rep_seqs_path, logger)

        # Update dataset paths and counts
        dataset.rep_seqs_path = str(rep_seqs_path)
        dataset.asv_count = dada2_result["asv_count"]
        dataset.sample_count = dada2_result["sample_count"]

        # Create Sample rows from track_reads
        for sample_name, counts in dada2_result["track_reads"].items():
            sample = Sample(
                dataset_id=dataset_id,
                sample_name=sample_name,
                read_count_raw=counts.get("input"),
                read_count_filtered=counts.get("filtered"),
                read_count_denoised=counts.get("denoised"),
                read_count_nonchimeric=counts.get("nonchim"),
            )
            db.add(sample)

        db.flush()

        # Propagate UploadMetadata → SampleMetadata (from all linked uploads)
        for upload in uploads:
            _propagate_upload_metadata(db, upload, dataset_id, logger)

        # Store FastQC metrics linked to samples
        if qc_result.get("sample_metrics"):
            _store_qc_metrics(db, dataset_id, qc_result["sample_metrics"])

        db.commit()
        completed_steps.append("dada2_longread" if platform in ("pacbio", "nanopore") else "dada2")

        # ── Step 4: Taxonomy Assignment ───────────────────────────────

        from app.pipeline.taxonomy import run_taxonomy

        _check_cancel(dataset_id)
        _update_status(output_dir, "taxonomy", 70, completed_steps)
        logger.info("Running taxonomy assignment...")

        tax_result = run_taxonomy(
            rep_seqs_path=rep_seqs_path,
            output_dir=output_dir,
            threads=threads,
            logger=logger,
            proc_callback=lambda proc: (
                _register_proc(dataset_id, proc),
                _save_pid_to_status(DATASET_DIR / str(dataset_id) / "status.json", proc.pid),
            ),
            skip_species=platform not in ("pacbio", "nanopore"),
        )

        _check_cancel(dataset_id)

        taxonomy_path = Path(tax_result["taxonomy_path"])
        dataset.taxonomy_path = str(taxonomy_path)
        db.commit()
        completed_steps.append("taxonomy")
        _update_status(output_dir, "taxonomy", 90, completed_steps)

        # ── Convert ASV table to BIOM (with taxonomy) ─────────────────

        from app.pipeline.biom_convert import tsv_to_biom

        biom_path = tsv_to_biom(
            tsv_path=Path(dada2_result["asv_table_path"]),
            output_dir=output_dir,
            logger=logger,
            taxonomy_path=taxonomy_path,
        )
        dataset.asv_table_path = str(biom_path)
        db.commit()

        # ── Generate QC report PDF ────────────────────────────────────

        try:
            from app.pipeline.qc_pdf import generate_qc_pdf
            generate_qc_pdf(dataset_id, output_dir, logger)
        except Exception as e:
            logger.warning(f"QC PDF generation failed (non-fatal): {e}")

        # ── Done ───────────────────────────────────────────────────────

        dataset.status = "complete"
        db.commit()
        _update_status(output_dir, "complete", 100, completed_steps)
        logger.info("Pipeline completed successfully!")

    except PipelineCancelled:
        was_cancelled = True
        logger.info("Pipeline cancelled by user.")
        try:
            dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
            if dataset:
                dataset.status = "cancelled"
                db.commit()
        except Exception:
            pass
        _update_status(output_dir, "cancelled", 0, completed_steps)

    except Exception as e:
        # If the cancel event is set, treat subprocess errors as cancellation
        cancel_event = _cancel_events.get(dataset_id)
        if cancel_event and cancel_event.is_set():
            was_cancelled = True
            logger.info("Pipeline cancelled by user.")
            try:
                dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
                if dataset:
                    dataset.status = "cancelled"
                    db.commit()
            except Exception:
                pass
            _update_status(output_dir, "cancelled", 0, completed_steps)
        else:
            logger.error(f"Pipeline failed: {e}", exc_info=True)
            try:
                dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
                if dataset:
                    dataset.status = "failed"
                    db.commit()
            except Exception:
                pass
            _update_status(output_dir, "failed", 0, completed_steps)
    finally:
        db.close()
        _running_pipelines.pop(dataset_id, None)
        _cancel_events.pop(dataset_id, None)
        _active_procs.pop(dataset_id, None)
        logger.handlers.clear()
        if was_cancelled and log_file.exists():
            log_file.unlink()


def _extract_rep_seqs(asv_tsv: Path, fasta_path: Path, logger: logging.Logger):
    """Write representative sequences from the ASV table TSV to a FASTA file."""
    import pandas as pd

    df = pd.read_csv(asv_tsv, sep="\t")
    lines = []
    for _, row in df.iterrows():
        lines.append(f">{row['ASV_ID']}")
        lines.append(row["sequence"])
    fasta_path.write_text("\n".join(lines) + "\n")
    logger.info(f"Extracted {len(df)} rep seqs to {fasta_path}")


def _propagate_upload_metadata(db, upload, dataset_id: int, logger):
    """Copy UploadMetadata rows to SampleMetadata for newly created Samples."""
    from app.db.models import Sample, SampleMetadata, UploadMetadata

    meta_rows = (
        db.query(UploadMetadata)
        .filter(UploadMetadata.upload_id == upload.id)
        .all()
    )
    if not meta_rows:
        return

    samples = db.query(Sample).filter(Sample.dataset_id == dataset_id).all()
    sample_map = {s.sample_name: s.id for s in samples}

    count = 0
    for row in meta_rows:
        sample_id = sample_map.get(row.sample_name)
        if sample_id:
            db.add(
                SampleMetadata(
                    sample_id=sample_id,
                    key=row.key,
                    value=row.value,
                )
            )
            count += 1

    if count:
        logger.info(f"Propagated {count} metadata entries from upload to samples")


def _store_qc_metrics(db, dataset_id: int, sample_metrics: dict):
    """Store FastQC metrics as QcMetric rows linked to the dataset's samples."""
    from app.db.models import QcMetric, Sample

    samples = db.query(Sample).filter(Sample.dataset_id == dataset_id).all()
    sample_map = {s.sample_name: s.id for s in samples}

    for qc_name, metrics in sample_metrics.items():
        from app.pipeline.detect import extract_sample_name
        clean_name = extract_sample_name(qc_name)
        sample_id = sample_map.get(clean_name)
        if not sample_id:
            continue

        if "_R1" in qc_name:
            direction = "forward"
        elif "_R2" in qc_name:
            direction = "reverse"
        else:
            direction = "single"

        for metric_name, value in metrics.items():
            if metric_name.startswith("module_"):
                continue
            try:
                db.add(
                    QcMetric(
                        sample_id=sample_id,
                        metric_name=metric_name,
                        metric_value=float(value),
                        read_direction=direction,
                    )
                )
            except (ValueError, TypeError):
                pass
