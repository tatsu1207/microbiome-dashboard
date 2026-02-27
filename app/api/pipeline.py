"""
MicrobiomeDash — FastAPI pipeline endpoints.
"""
from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel
from sqlalchemy.orm import Session

from app.db.database import get_db
from app.db.models import Dataset, FastqFile, Upload
from app.pipeline.runner import cancel_pipeline, get_pipeline_status, launch_pipeline

router = APIRouter(prefix="/api", tags=["pipeline"])


class PipelineLaunchRequest(BaseModel):
    file_ids: list[int] | None = None
    upload_id: int | None = None  # convenience: selects all files from this upload
    dataset_name: str
    description: str | None = None
    project_id: int | None = None
    trim_left_f: int | None = None
    trim_left_r: int | None = None
    trunc_len_f: int | None = None
    trunc_len_r: int | None = None
    min_overlap: int | None = None


@router.post("/pipeline/launch")
async def launch(req: PipelineLaunchRequest, db: Session = Depends(get_db)):
    """Launch a new pipeline run for the given files or upload."""
    # Resolve file_ids: either provided directly, or from upload_id
    file_ids = req.file_ids
    if not file_ids:
        if not req.upload_id:
            raise HTTPException(400, "Provide file_ids or upload_id")
        upload = db.query(Upload).filter(Upload.id == req.upload_id).first()
        if not upload:
            raise HTTPException(404, f"Upload {req.upload_id} not found")
        file_ids = [
            f.id for f in
            db.query(FastqFile).filter(FastqFile.upload_id == req.upload_id).all()
        ]
        if not file_ids:
            raise HTTPException(400, f"No files in upload {req.upload_id}")

    dataset_id = launch_pipeline(
        file_ids=file_ids,
        dataset_name=req.dataset_name,
        description=req.description,
        project_id=req.project_id,
        trim_left_f=req.trim_left_f,
        trim_left_r=req.trim_left_r,
        trunc_len_f=req.trunc_len_f,
        trunc_len_r=req.trunc_len_r,
        min_overlap=req.min_overlap,
    )
    return {"dataset_id": dataset_id, "status": "launched"}


@router.post("/pipeline/cancel/{dataset_id}")
async def cancel(dataset_id: int):
    """Cancel a running pipeline."""
    ok = cancel_pipeline(dataset_id)
    if not ok:
        raise HTTPException(404, "No running pipeline found for this dataset")
    return {"dataset_id": dataset_id, "status": "cancelling"}


@router.get("/pipeline/status/{dataset_id}")
async def status(dataset_id: int):
    """Get pipeline status for a dataset."""
    return get_pipeline_status(dataset_id)


@router.get("/pipeline/datasets")
async def list_pipeline_datasets(db: Session = Depends(get_db)):
    """List all pipeline datasets with their status."""
    datasets = (
        db.query(Dataset)
        .filter(Dataset.source_type == "pipeline")
        .order_by(Dataset.created_at.desc())
        .all()
    )
    return [
        {
            "id": d.id,
            "name": d.name,
            "upload_id": d.upload_id,
            "status": d.status,
            "sequencing_type": d.sequencing_type,
            "variable_region": d.variable_region,
            "sample_count": d.sample_count,
            "asv_count": d.asv_count,
            "created_at": d.created_at.isoformat() if d.created_at else None,
        }
        for d in datasets
    ]
