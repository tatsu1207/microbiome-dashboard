#!/usr/bin/env python3
"""Register MST FASTQ files + metadata into the 16S Analyzer database.

Usage:
    conda run -n microbiome_16S python register_mst.py
"""
import csv
import sys
from pathlib import Path

# Ensure project root is on the path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from app.db.database import init_db, get_session
from app.db.models import UploadMetadata
from app.utils.file_handler import register_upload

FASTQ_DIR = Path.home() / "github/MST-Pipeline/DB_dev/data"
METADATA_FILE = Path.home() / "github/MST-Pipeline/DB/MST.design"


def main():
    # Ensure DB tables exist
    init_db()

    with get_session() as db:
        # 1. Register FASTQ files (symlinked to avoid copying 38GB)
        print(f"Registering FASTQ files from {FASTQ_DIR} ...")
        upload = register_upload(FASTQ_DIR, db, symlink=True)
        print(f"  Upload ID: {upload.id}")
        print(f"  Type: {upload.sequencing_type}")
        print(f"  Region: {upload.variable_region}")
        print(f"  Files: {upload.total_files}")
        print(f"  Size: {upload.total_size_mb:.0f} MB")

        # 2. Load metadata (Sample → Group as "source")
        print(f"\nLoading metadata from {METADATA_FILE} ...")
        meta_count = 0
        with open(METADATA_FILE) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                sample = row["Sample"]
                group = row["Group"]
                db.add(UploadMetadata(
                    upload_id=upload.id,
                    sample_name=sample,
                    key="source",
                    value=group,
                ))
                meta_count += 1

        db.flush()
        print(f"  Loaded metadata for {meta_count} samples")

    print("\nDone! You can now launch the DADA2 pipeline from the web UI.")


if __name__ == "__main__":
    main()
