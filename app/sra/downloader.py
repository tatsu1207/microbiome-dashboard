"""
SRA Downloader — Download FASTQ files from NCBI SRA by accession number.
"""
import json
import logging
import os
import re
import shutil
import signal
import subprocess
import threading
import urllib.request
import urllib.parse
from pathlib import Path

from app.config import UPLOAD_DIR, SRA_CACHE_DIR, conda_cmd
from app.db.database import SessionLocal
from app.db.models import FastqFile, Upload
from app.pipeline.detect import (
    PRIMERS,
    _primer_matches,
    _read_fastq_sequences,
    detect_platform,
    detect_sequencing_type,
    detect_variable_region,
    extract_sample_name,
)

logger = logging.getLogger(__name__)

# ── Global state for cancellation ────────────────────────────────────────────

_cancel_events: dict[str, threading.Event] = {}
_active_procs: dict[str, subprocess.Popen] = {}
_download_threads: dict[str, threading.Thread] = {}

# ── Accession validation ─────────────────────────────────────────────────────

_SRR_RE = re.compile(r"^[SED]RR\d+$", re.IGNORECASE)
_SRP_RE = re.compile(r"^[SED]RP\d+$", re.IGNORECASE)
_PRJNA_RE = re.compile(r"^PRJ[A-Z]{2}\d+$", re.IGNORECASE)


def parse_accessions(text: str) -> list[str]:
    """Parse accession input text into a list of clean accession IDs."""
    # Split on newlines, commas, spaces, tabs
    tokens = re.split(r"[,\s]+", text.strip())
    return [t.strip() for t in tokens if t.strip()]


def validate_accessions(accessions: list[str]) -> tuple[list[str], list[str]]:
    """Validate accessions. Returns (valid, invalid) lists."""
    valid, invalid = [], []
    for acc in accessions:
        if _SRR_RE.match(acc) or _SRP_RE.match(acc) or _PRJNA_RE.match(acc):
            valid.append(acc.upper())
        else:
            invalid.append(acc)
    return valid, invalid


# ── Resolve SRP/PRJNA to SRR list ────────────────────────────────────────────


def resolve_to_srr(accessions: list[str]) -> list[str]:
    """
    Expand any SRP/PRJNA accessions to individual SRR IDs using NCBI E-utilities.
    SRR accessions pass through unchanged.
    """
    srr_list = []
    for acc in accessions:
        if _SRR_RE.match(acc):
            srr_list.append(acc)
        elif _SRP_RE.match(acc) or _PRJNA_RE.match(acc):
            expanded = _eutils_search_srr(acc)
            srr_list.extend(expanded)
    return srr_list


def _eutils_search_srr(accession: str) -> list[str]:
    """Use NCBI E-utilities to find SRR IDs for an SRP or PRJNA accession."""
    try:
        # Step 1: esearch to get WebEnv/QueryKey
        params = urllib.parse.urlencode({
            "db": "sra",
            "term": accession,
            "usehistory": "y",
            "retmax": 0,
        })
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{params}&retmode=json"
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = json.loads(resp.read())

        result = data.get("esearchresult", {})
        count = int(result.get("count", 0))
        webenv = result.get("webenv", "")
        query_key = result.get("querykey", "")

        if count == 0 or not webenv:
            logger.warning(f"No SRA results for {accession}")
            return []

        # Step 2: efetch to get run accessions
        params2 = urllib.parse.urlencode({
            "db": "sra",
            "query_key": query_key,
            "WebEnv": webenv,
            "rettype": "runinfo",
            "retmode": "text",
            "retmax": count,
        })
        url2 = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{params2}"
        with urllib.request.urlopen(url2, timeout=60) as resp2:
            text = resp2.read().decode("utf-8")

        # Parse CSV-like runinfo output — "Run" column has SRR IDs
        srr_ids = []
        lines = text.strip().split("\n")
        if len(lines) < 2:
            return []
        header = lines[0].split(",")
        try:
            run_idx = header.index("Run")
        except ValueError:
            return []
        for line in lines[1:]:
            fields = line.split(",")
            if len(fields) > run_idx:
                srr_id = fields[run_idx].strip()
                if _SRR_RE.match(srr_id):
                    srr_ids.append(srr_id)
        return srr_ids

    except Exception as e:
        logger.error(f"E-utilities lookup failed for {accession}: {e}")
        return []


# ── Download a single SRR ────────────────────────────────────────────────────


def download_srr(
    srr_id: str,
    output_dir: Path,
    job_id: str,
    cancel_event: threading.Event,
) -> list[Path]:
    """
    Download and convert a single SRR accession to FASTQ files.
    Tries prefetch + fasterq-dump first; if prefetch fails, falls back to
    fasterq-dump directly (which downloads on its own).
    Returns list of gzipped FASTQ file paths.
    """
    if cancel_event.is_set():
        return []

    cache_dir = SRA_CACHE_DIR / srr_id
    cache_dir.mkdir(parents=True, exist_ok=True)
    threads = min(8, max(1, os.cpu_count() - 1))

    # Try prefetch first, fall back to direct fasterq-dump on failure
    prefetch_ok = False
    try:
        logger.info(f"Prefetching {srr_id}...")
        _run_cmd(
            conda_cmd(["prefetch", srr_id, "-O", str(SRA_CACHE_DIR)]),
            job_id, cancel_event,
        )
        prefetch_ok = True
    except RuntimeError as e:
        logger.warning(f"prefetch failed for {srr_id}, falling back to direct download: {e}")

    if cancel_event.is_set():
        return []

    # Determine SRA source path
    if prefetch_ok:
        sra_file = SRA_CACHE_DIR / srr_id / f"{srr_id}.sra"
        if not sra_file.exists():
            sra_file = SRA_CACHE_DIR / srr_id / srr_id
        sra_source = str(sra_file)
    else:
        sra_source = srr_id  # fasterq-dump/fastq-dump will download directly

    # Try fasterq-dump first; if it fails (e.g. PacBio), fall back to fastq-dump
    logger.info(f"Converting {srr_id} to FASTQ...")
    try:
        _run_cmd(
            conda_cmd([
                "fasterq-dump", sra_source,
                "-O", str(output_dir),
                "--split-3",
                "-e", str(threads),
            ]),
            job_id, cancel_event,
        )
    except RuntimeError as e:
        if "PACBIO" in str(e) or "fastq-dump instead" in str(e):
            logger.warning(f"fasterq-dump does not support this format, using fastq-dump for {srr_id}")
            _run_cmd(
                conda_cmd([
                    "fastq-dump", sra_source,
                    "--outdir", str(output_dir),
                    "--split-3",
                    "--gzip",
                ]),
                job_id, cancel_event,
            )
        else:
            raise

    if cancel_event.is_set():
        return []

    # Gzip the output FASTQ files (fastq-dump --gzip already produces .gz)
    gzipped = list(output_dir.glob(f"{srr_id}*.fastq.gz"))
    fastq_files = list(output_dir.glob(f"{srr_id}*.fastq"))
    for fq in fastq_files:
        if cancel_event.is_set():
            return []
        logger.info(f"Compressing {fq.name}...")
        _run_cmd(["gzip", str(fq)], job_id, cancel_event)
        gz_path = Path(str(fq) + ".gz")
        if gz_path.exists():
            gzipped.append(gz_path)

    # Clean up SRA cache for this accession
    shutil.rmtree(cache_dir, ignore_errors=True)

    return gzipped


def _run_cmd(cmd: list[str], job_id: str, cancel_event: threading.Event):
    """Run a subprocess with cancellation support."""
    if cancel_event.is_set():
        return

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        start_new_session=True,
    )
    _active_procs[job_id] = proc

    try:
        proc.wait()
        if proc.returncode != 0:
            output = proc.stdout.read().decode("utf-8", errors="replace") if proc.stdout else ""
            raise RuntimeError(f"Command failed (exit {proc.returncode}): {' '.join(cmd[:3])}...\n{output[-500:]}")
    finally:
        _active_procs.pop(job_id, None)


# ── Register downloaded files ────────────────────────────────────────────────


def _count_reads(fastq_path: Path) -> int | None:
    """Count reads in a gzipped FASTQ file."""
    import gzip
    try:
        count = 0
        with gzip.open(fastq_path, "rt") as f:
            for _ in f:
                count += 1
        return count // 4
    except Exception:
        return None


def _avg_read_length(fastq_path: Path, n_reads: int = 200) -> int | None:
    """Compute average read length from first n reads."""
    import gzip
    try:
        lengths = []
        with gzip.open(fastq_path, "rt") as f:
            line_num = 0
            for line in f:
                line_num += 1
                if line_num % 4 == 2:
                    lengths.append(len(line.strip()))
                    if len(lengths) >= n_reads:
                        break
        return int(sum(lengths) / len(lengths)) if lengths else None
    except Exception:
        return None


def register_downloaded_files(
    fastq_dir: Path,
    project_id: int | None = None,
    study: str | None = None,
) -> int:
    """
    Register downloaded FASTQ files as a new Upload.
    Moves files to UPLOAD_DIR and creates DB records.
    Returns the upload ID.
    """
    gz_files = sorted(fastq_dir.glob("*.fastq.gz"))
    if not gz_files:
        raise ValueError("No .fastq.gz files found in download directory")

    filenames = [f.name for f in gz_files]

    # Move to flat UPLOAD_DIR
    total_size = 0.0
    for fpath in gz_files:
        dest = UPLOAD_DIR / fpath.name
        # Handle name conflicts
        if dest.exists():
            stem = fpath.stem.replace(".fastq", "")
            dest = UPLOAD_DIR / f"{stem}_sra.fastq.gz"
        shutil.move(str(fpath), str(dest))
        total_size += dest.stat().st_size / (1024 * 1024)

    # Re-collect the actual filenames after move
    saved_filenames = [f.name for f in sorted(UPLOAD_DIR.glob("*.fastq.gz")) if f.name in filenames or f.name.replace("_sra.fastq.gz", ".fastq.gz") in filenames]
    # Simpler: re-scan based on what we actually moved
    saved_filenames = []
    for fpath in gz_files:
        dest = UPLOAD_DIR / fpath.name
        alt_dest = UPLOAD_DIR / f"{fpath.stem.replace('.fastq', '')}_sra.fastq.gz"
        if dest.exists():
            saved_filenames.append(dest.name)
        elif alt_dest.exists():
            saved_filenames.append(alt_dest.name)

    # Detect sequencing type
    detection = detect_sequencing_type(saved_filenames)

    # Detect variable region from first R1 file
    variable_region = None
    first_sample = next(iter(detection["samples"].values()), {})
    r1_name = first_sample.get("R1") or first_sample.get("single")
    if r1_name:
        try:
            region_result = detect_variable_region(UPLOAD_DIR / r1_name)
            variable_region = region_result["region"]
        except Exception:
            pass

    # Detect platform
    platform = None
    if r1_name:
        try:
            platform_result = detect_platform(UPLOAD_DIR / r1_name)
            platform = platform_result["platform"]
        except Exception:
            pass

    # Long-read adjustments
    if platform in ("pacbio", "nanopore"):
        variable_region = "V1-V9"
        detection["type"] = "single-end"

    # Detect primers
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
            project_id=project_id,
            upload_dir=str(UPLOAD_DIR),
            sequencing_type=detection["type"],
            variable_region=variable_region,
            platform=platform,
            primers_detected=primers_detected,
            study=(study or "").strip() or None,
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
        upload_id = upload.id
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()

    return upload_id


# ── Main download job (runs in background thread) ────────────────────────────


def run_sra_download(
    job_id: str,
    accessions: list[str],
    output_dir: Path,
    project_id: int | None = None,
    study: str | None = None,
):
    """
    Main download job function. Runs in a background thread.
    Resolves accessions, downloads each SRR, registers files.
    Progress tracked via status.json in output_dir.
    """
    cancel_event = _cancel_events.setdefault(job_id, threading.Event())
    status_path = output_dir / "status.json"

    def _update_status(**kwargs):
        current = {}
        if status_path.exists():
            try:
                current = json.loads(status_path.read_text())
            except Exception:
                pass
        current.update(kwargs)
        status_path.write_text(json.dumps(current))

    log_lines = []

    def _log(msg: str):
        """Append a log line and persist it to status.json."""
        log_lines.append(msg)
        logger.info(msg)
        _update_status(log=log_lines)

    try:
        _update_status(
            status="resolving",
            total_accessions=len(accessions),
            completed=0,
            current="Resolving accessions...",
            error=None,
            log=[],
        )
        _log(f"Resolving {len(accessions)} accession(s)...")

        # Resolve SRP/PRJNA to SRR
        srr_list = resolve_to_srr(accessions)
        if not srr_list:
            _log("No SRR accessions found")
            _update_status(status="failed", error="No SRR accessions found")
            return

        _log(f"Resolved to {len(srr_list)} SRR run(s)")
        _update_status(
            status="downloading",
            total_accessions=len(srr_list),
            srr_list=srr_list,
        )

        # Download SRRs in parallel (up to 3 concurrent)
        from concurrent.futures import ThreadPoolExecutor, as_completed

        max_workers = min(3, len(srr_list))
        all_files = []
        completed_count = 0
        lock = threading.Lock()

        def _download_one(idx_srr):
            idx, srr_id = idx_srr
            if cancel_event.is_set():
                return srr_id, [], None
            try:
                files = download_srr(srr_id, output_dir, job_id, cancel_event)
                return srr_id, files, None
            except Exception as e:
                return srr_id, [], e

        _log(f"Downloading {len(srr_list)} run(s) with {max_workers} parallel workers...")

        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = {
                pool.submit(_download_one, (i, srr_id)): (i, srr_id)
                for i, srr_id in enumerate(srr_list)
            }

            for future in as_completed(futures):
                if cancel_event.is_set():
                    _update_status(status="cancelled")
                    pool.shutdown(wait=False, cancel_futures=True)
                    return

                srr_id, files, error = future.result()
                with lock:
                    completed_count += 1
                    if error:
                        _log(f"[{completed_count}/{len(srr_list)}] {srr_id} FAILED: {error}")
                        raise error
                    all_files.extend(files)
                    file_names = ", ".join(f.name for f in files) if files else "no files"
                    _log(f"[{completed_count}/{len(srr_list)}] {srr_id} done ({len(files)} file(s): {file_names})")
                    _update_status(
                        current=f"Downloaded {completed_count}/{len(srr_list)}",
                        completed=completed_count,
                    )

        if cancel_event.is_set():
            _update_status(status="cancelled")
            return

        # Register downloaded files
        _log("Registering downloaded files...")
        _update_status(
            status="registering",
            current="Registering files...",
            completed=len(srr_list),
        )

        upload_id = register_downloaded_files(output_dir, project_id, study)
        _log(f"Registered as Upload #{upload_id}")

        # Fetch BioSample metadata (non-blocking — failure won't affect download)
        try:
            _log("Fetching BioSample metadata...")
            _update_status(current="Fetching BioSample metadata...")
            from app.sra.metadata_fetcher import fetch_biosample_metadata, store_biosample_metadata
            bio_meta = fetch_biosample_metadata(srr_list)
            if bio_meta:
                store_biosample_metadata(upload_id, bio_meta)
                _log(f"Stored BioSample metadata for {len(bio_meta)} sample(s)")
            else:
                _log("No BioSample metadata found")
        except Exception as e:
            _log(f"BioSample metadata fetch failed (non-critical): {e}")

        _log(f"Complete! {len(all_files)} file(s) downloaded")
        _update_status(
            status="complete",
            current="Done",
            upload_id=upload_id,
            n_files=len(all_files),
        )

    except Exception as e:
        logger.exception(f"SRA download job {job_id} failed")
        _log(f"FAILED: {e}")
        _update_status(status="failed", error=str(e))
    finally:
        _cancel_events.pop(job_id, None)
        _download_threads.pop(job_id, None)


def launch_download(
    job_id: str,
    accessions: list[str],
    output_dir: Path,
    project_id: int | None = None,
    study: str | None = None,
):
    """Launch a download job in a background thread."""
    output_dir.mkdir(parents=True, exist_ok=True)
    _cancel_events[job_id] = threading.Event()

    t = threading.Thread(
        target=run_sra_download,
        args=(job_id, accessions, output_dir, project_id, study),
        daemon=True,
    )
    _download_threads[job_id] = t
    t.start()


def cancel_download(job_id: str):
    """Cancel a running download job."""
    event = _cancel_events.get(job_id)
    if event:
        event.set()
    proc = _active_procs.get(job_id)
    if proc and proc.poll() is None:
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
        except (ProcessLookupError, PermissionError):
            pass


def get_job_status(job_id: str, output_dir: Path) -> dict:
    """Read status.json for a download job."""
    status_path = output_dir / "status.json"
    if status_path.exists():
        try:
            return json.loads(status_path.read_text())
        except Exception:
            return {"status": "unknown"}
    return {"status": "pending"}
