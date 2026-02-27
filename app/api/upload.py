"""
MicrobiomeDash — FastAPI upload endpoint for external/CLI access.
"""
import base64
import shutil
import uuid
from pathlib import Path

from fastapi import APIRouter, Depends, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse
from sqlalchemy.orm import Session

from app.config import UPLOAD_DIR
from app.db.database import get_db
from app.db.models import FastqFile, Upload
from app.pipeline.detect import detect_sequencing_type, extract_sample_name

router = APIRouter(prefix="/api", tags=["upload"])


@router.post("/upload")
async def upload_files(
    fastq_files: list[UploadFile] = File(...),
    metadata: UploadFile | None = File(None),
    project_id: int | None = Form(None),
    db: Session = Depends(get_db),
):
    """Upload FASTQ files and optional metadata via multipart form."""
    upload_id = uuid.uuid4().hex[:8]
    upload_dir = UPLOAD_DIR / upload_id
    fastq_dir = upload_dir / "fastq"
    fastq_dir.mkdir(parents=True, exist_ok=True)

    # Save FASTQ files
    filenames = []
    total_size = 0.0
    for f in fastq_files:
        filepath = fastq_dir / f.filename
        with open(filepath, "wb") as out:
            shutil.copyfileobj(f.file, out)
        filenames.append(f.filename)
        total_size += filepath.stat().st_size / (1024 * 1024)

    # Detect SE/PE
    detection = detect_sequencing_type(filenames)

    # Create DB records
    upload = Upload(
        project_id=project_id,
        upload_dir=str(upload_dir),
        sequencing_type=detection["type"],
        total_files=len(filenames),
        total_size_mb=round(total_size, 2),
        status="uploaded",
    )
    db.add(upload)
    db.flush()

    for filename in filenames:
        sample_name = extract_sample_name(filename)
        sample_info = detection["samples"].get(sample_name, {})
        if sample_info.get("R1") == filename:
            read_direction = "R1"
        elif sample_info.get("R2") == filename:
            read_direction = "R2"
        else:
            read_direction = "single"

        db.add(
            FastqFile(
                upload_id=upload.id,
                sample_name=sample_name,
                filename=filename,
                file_path=str(fastq_dir / filename),
                read_direction=read_direction,
                file_size_mb=round(
                    (fastq_dir / filename).stat().st_size / (1024 * 1024), 2
                ),
            )
        )

    db.commit()

    return {
        "upload_id": upload.id,
        "upload_dir": str(upload_dir),
        "sequencing_type": detection["type"],
        "total_files": len(filenames),
        "total_size_mb": round(total_size, 2),
        "samples": list(detection["samples"].keys()),
        "errors": detection["errors"],
    }


@router.get("/uploads")
async def list_uploads(db: Session = Depends(get_db)):
    """List all uploads."""
    uploads = db.query(Upload).order_by(Upload.created_at.desc()).all()
    return [
        {
            "id": u.id,
            "sequencing_type": u.sequencing_type,
            "variable_region": u.variable_region,
            "total_files": u.total_files,
            "total_size_mb": u.total_size_mb,
            "status": u.status,
            "created_at": u.created_at.isoformat() if u.created_at else None,
        }
        for u in uploads
    ]


@router.get("/uploads/{upload_id}/files/{filename}")
async def download_fastq(upload_id: int, filename: str, db: Session = Depends(get_db)):
    """Download a FASTQ file from an upload."""
    upload = db.query(Upload).filter(Upload.id == upload_id).first()
    if not upload:
        raise HTTPException(status_code=404, detail="Upload not found")

    # Build path and protect against traversal
    upload_dir = Path(upload.upload_dir)

    file_path = (upload_dir / filename).resolve()
    try:
        file_path.relative_to(upload_dir.resolve())
    except ValueError:
        raise HTTPException(status_code=400, detail="Invalid filename")

    if not file_path.is_file():
        raise HTTPException(status_code=404, detail="File not found")

    return FileResponse(
        path=str(file_path),
        filename=filename,
        media_type="application/gzip" if filename.endswith(".gz") else "application/octet-stream",
    )
