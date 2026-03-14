"""
SRA Submission Helper — Generate NCBI SRA metadata spreadsheets from upload data.
"""
import pandas as pd

from app.db.database import SessionLocal
from app.db.models import FastqFile, Upload

# Map internal platform names to NCBI SRA platform names
_PLATFORM_MAP = {
    "illumina": "ILLUMINA",
    "pacbio": "PACBIO_SMRT",
    "nanopore": "OXFORD_NANOPORE",
}

# Instrument model options per platform
INSTRUMENT_OPTIONS = {
    "ILLUMINA": [
        "Illumina MiSeq",
        "Illumina HiSeq 2500",
        "Illumina HiSeq 4000",
        "Illumina NextSeq 500",
        "Illumina NextSeq 2000",
        "Illumina NovaSeq 6000",
        "Illumina iSeq 100",
    ],
    "PACBIO_SMRT": [
        "PacBio RS II",
        "Sequel",
        "Sequel II",
        "Sequel IIe",
        "Revio",
    ],
    "OXFORD_NANOPORE": [
        "MinION",
        "GridION",
        "PromethION",
    ],
}


def generate_sra_metadata(upload_id: int) -> pd.DataFrame:
    """
    Generate an NCBI SRA submission metadata DataFrame from an Upload record.

    Each row = one sample. For paired-end, filename and filename2 are populated.
    """
    db = SessionLocal()
    try:
        upload = db.query(Upload).filter(Upload.id == upload_id).first()
        if not upload:
            raise ValueError(f"Upload {upload_id} not found")

        files = (
            db.query(FastqFile)
            .filter(FastqFile.upload_id == upload_id)
            .order_by(FastqFile.sample_name, FastqFile.read_direction)
            .all()
        )
    finally:
        db.close()

    # Group files by sample
    samples: dict[str, dict] = {}
    for f in files:
        if f.sample_name not in samples:
            samples[f.sample_name] = {}
        if f.read_direction == "R1":
            samples[f.sample_name]["R1"] = f.filename
        elif f.read_direction == "R2":
            samples[f.sample_name]["R2"] = f.filename
        else:
            samples[f.sample_name]["single"] = f.filename

    is_paired = upload.sequencing_type == "paired-end"
    sra_platform = _PLATFORM_MAP.get(upload.platform or "illumina", "ILLUMINA")

    # Build design description hint
    region = upload.variable_region or ""
    design_hint = f"16S rRNA {region} amplicon sequencing".strip() if region else "16S rRNA amplicon sequencing"

    rows = []
    for sample_name, file_info in sorted(samples.items()):
        row = {
            "biosample_accession": "",
            "library_ID": sample_name,
            "title": "",
            "library_strategy": "AMPLICON",
            "library_source": "METAGENOMIC",
            "library_selection": "PCR",
            "library_layout": "paired" if is_paired else "single",
            "platform": sra_platform,
            "instrument_model": "",
            "design_description": design_hint,
            "filetype": "fastq",
            "filename": file_info.get("R1") or file_info.get("single", ""),
            "filename2": file_info.get("R2", "") if is_paired else "",
        }
        rows.append(row)

    return pd.DataFrame(rows)


# ── BioSample metadata generation ────────────────────────────────────────────

# MIMS package definitions: required fields per package type
BIOSAMPLE_PACKAGES = {
    "MIMS.me.host-associated": {
        "label": "Host-associated (gut, skin, etc.)",
        "organism_default": "host-associated metagenome",
        "required": [
            "sample_name", "organism", "collection_date", "env_broad_scale",
            "env_local_scale", "env_medium", "geo_loc_name", "host", "lat_lon",
        ],
        "common_optional": [
            "host_tissue_sampled", "host_body_product", "host_age",
            "host_sex", "host_disease", "host_common_name",
            "isolation_source", "sample_collection_device_or_method",
            "description",
        ],
    },
    "MIMS.me.human-associated": {
        "label": "Human-associated",
        "organism_default": "human gut metagenome",
        "required": [
            "sample_name", "organism", "collection_date", "env_broad_scale",
            "env_local_scale", "env_medium", "geo_loc_name", "host", "lat_lon",
        ],
        "common_optional": [
            "host_tissue_sampled", "host_body_product", "host_age",
            "host_sex", "host_disease", "host_subject_id",
            "isolation_source", "description",
        ],
    },
    "MIMS.me.soil": {
        "label": "Soil",
        "organism_default": "soil metagenome",
        "required": [
            "sample_name", "organism", "collection_date", "depth",
            "elevation", "env_broad_scale", "env_local_scale", "env_medium",
            "geo_loc_name", "lat_lon",
        ],
        "common_optional": [
            "isolation_source", "pH", "soil_type", "soil_horizon",
            "current_land_use", "crop_rotation", "tillage",
            "total_organic_carbon", "total_nitrogen_content",
            "collection_method", "description",
        ],
    },
    "MIMS.me.water": {
        "label": "Water",
        "organism_default": "water metagenome",
        "required": [
            "sample_name", "organism", "collection_date", "depth",
            "env_broad_scale", "env_local_scale", "env_medium",
            "geo_loc_name", "lat_lon",
        ],
        "common_optional": [
            "isolation_source", "temperature", "salinity",
            "collection_method", "size_fraction_selected",
            "description",
        ],
    },
    "MIMS.me.microbial": {
        "label": "Microbial (general environmental)",
        "organism_default": "metagenome",
        "required": [
            "sample_name", "organism", "collection_date",
            "env_broad_scale", "env_local_scale", "env_medium",
            "geo_loc_name", "lat_lon",
        ],
        "common_optional": [
            "isolation_source", "depth", "elevation", "altitude",
            "temperature", "collection_method", "description",
        ],
    },
}


def generate_biosample_metadata(
    upload_id: int,
    package: str,
    shared_values: dict | None = None,
) -> pd.DataFrame:
    """
    Generate a BioSample submission metadata DataFrame.

    Args:
        upload_id: Upload ID to pull sample names from.
        package: MIMS package key (e.g. "MIMS.me.soil").
        shared_values: Dict of {attribute: value} applied to all samples.

    Returns:
        DataFrame with one row per sample, columns = required + optional attributes.
    """
    pkg = BIOSAMPLE_PACKAGES.get(package)
    if not pkg:
        raise ValueError(f"Unknown package: {package}")

    db = SessionLocal()
    try:
        sample_names = (
            db.query(FastqFile.sample_name)
            .filter(FastqFile.upload_id == upload_id)
            .distinct()
            .order_by(FastqFile.sample_name)
            .all()
        )
        names = sorted(set(n[0] for n in sample_names))
    finally:
        db.close()

    if not names:
        return pd.DataFrame()

    shared = shared_values or {}
    all_cols = pkg["required"] + pkg["common_optional"]

    rows = []
    for name in names:
        row = {col: "" for col in all_cols}
        row["sample_name"] = name
        row["organism"] = shared.get("organism") or pkg["organism_default"]
        # Apply shared values
        for k, v in shared.items():
            if k in row and v:
                row[k] = v
        rows.append(row)

    return pd.DataFrame(rows, columns=all_cols)


def get_upload_info(upload_id: int) -> dict:
    """Return summary info about an upload for display."""
    db = SessionLocal()
    try:
        upload = db.query(Upload).filter(Upload.id == upload_id).first()
        if not upload:
            return {}
        n_files = db.query(FastqFile).filter(FastqFile.upload_id == upload_id).count()
        # Count unique samples
        sample_names = (
            db.query(FastqFile.sample_name)
            .filter(FastqFile.upload_id == upload_id)
            .distinct()
            .all()
        )
        return {
            "id": upload.id,
            "sequencing_type": upload.sequencing_type or "unknown",
            "variable_region": upload.variable_region or "unknown",
            "platform": upload.platform or "illumina",
            "sra_platform": _PLATFORM_MAP.get(upload.platform or "illumina", "ILLUMINA"),
            "n_files": n_files,
            "n_samples": len(sample_names),
            "study": upload.study or "",
            "created_at": str(upload.created_at),
        }
    finally:
        db.close()
