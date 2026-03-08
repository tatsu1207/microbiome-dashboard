"""
MicrobiomeDash — Register local FASTQ files into the database.

Files are copied into data/uploads/ so pipeline I/O runs on the fast WSL filesystem.
"""
import shutil
import uuid
from pathlib import Path

from sqlalchemy.orm import Session

from app.config import UPLOAD_DIR
from app.db.models import FastqFile, Upload
from app.pipeline.detect import (
    detect_sequencing_type,
    detect_variable_region,
    extract_sample_name,
)

FASTQ_EXTENSIONS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def _unique_dest(dest: Path) -> Path:
    """Return a path that doesn't collide with existing files.

    For 'sample_R1.fastq.gz', tries 'sample_R1_1.fastq.gz', '_2', etc.
    """
    if not dest.exists():
        return dest

    # Handle compound extensions like .fastq.gz
    name = dest.name
    for ext in FASTQ_EXTENSIONS:
        if name.endswith(ext):
            stem = name[: -len(ext)]
            suffix = ext
            break
    else:
        stem, suffix = dest.stem, dest.suffix

    counter = 1
    while True:
        candidate = dest.parent / f"{stem}_{counter}{suffix}"
        if not candidate.exists():
            return candidate
        counter += 1


def scan_fastq_directory(directory: str | Path) -> list[Path]:
    """Find all FASTQ files in a directory (non-recursive)."""
    d = Path(directory)
    if not d.is_dir():
        raise FileNotFoundError(f"Directory not found: {d}")
    files = []
    for ext in FASTQ_EXTENSIONS:
        files.extend(d.glob(f"*{ext}"))
    return sorted(files)


def register_upload(
    fastq_dir: str | Path,
    db: Session,
    project_id: int | None = None,
    *,
    symlink: bool = False,
) -> Upload:
    """Register local FASTQ files in the database.

    Files are copied (or symlinked) into data/uploads/{id}/fastq/ so that
    pipeline I/O uses a consistent path layout.

    Args:
        fastq_dir: Path to directory containing FASTQ files.
        db: SQLAlchemy session.
        project_id: Optional project to associate with.
        symlink: If True, create symlinks instead of copying files.
                 Much faster for large datasets already on a local filesystem.

    Returns:
        The created Upload ORM object.
    """
    fastq_dir = Path(fastq_dir)
    fastq_files = scan_fastq_directory(fastq_dir)

    if not fastq_files:
        raise ValueError(f"No FASTQ files found in {fastq_dir}")

    filenames = [f.name for f in fastq_files]

    # Detect SE/PE
    detection = detect_sequencing_type(filenames)

    # Create destination directory
    upload_id = uuid.uuid4().hex[:8]
    dest_dir = UPLOAD_DIR / upload_id
    dest_fastq_dir = dest_dir / "fastq"
    dest_fastq_dir.mkdir(parents=True, exist_ok=True)

    # Copy or symlink FASTQ files into the project
    total_size = 0.0
    copied_names: dict[Path, Path] = {}  # original path -> destination path
    for fpath in fastq_files:
        dest_path = _unique_dest(dest_fastq_dir / fpath.name)
        if symlink:
            dest_path.symlink_to(fpath.resolve())
        else:
            shutil.copy2(fpath, dest_path)
        copied_names[fpath] = dest_path
        total_size += fpath.stat().st_size / (1024 * 1024)

    # Detect variable region from the first R1 / single-end file
    r1_file = next(
        (copied_names[f] for f in fastq_files
         if detection["samples"].get(extract_sample_name(f.name), {}).get("R1") == f.name),
        next(iter(copied_names.values())),
    )
    region_result = detect_variable_region(r1_file)
    variable_region = region_result["region"]

    # Create Upload record pointing to the local copy
    upload = Upload(
        project_id=project_id,
        source_dir=str(fastq_dir),
        upload_dir=str(dest_fastq_dir),
        sequencing_type=detection["type"],
        variable_region=variable_region,
        total_files=len(fastq_files),
        total_size_mb=round(total_size, 2),
        status="uploaded",
    )
    db.add(upload)
    db.flush()

    # Create FastqFile records with paths to destination files
    for fpath in fastq_files:
        sample_name = extract_sample_name(fpath.name)
        sample_info = detection["samples"].get(sample_name, {})

        if sample_info.get("R1") == fpath.name:
            read_direction = "R1"
        elif sample_info.get("R2") == fpath.name:
            read_direction = "R2"
        else:
            read_direction = "single"

        copied_path = copied_names[fpath]
        db.add(
            FastqFile(
                upload_id=upload.id,
                sample_name=sample_name,
                filename=fpath.name,
                file_path=str(copied_path),
                read_direction=read_direction,
                file_size_mb=round(fpath.stat().st_size / (1024 * 1024), 2),
            )
        )

    db.commit()
    return upload
