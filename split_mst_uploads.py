#!/usr/bin/env python3
"""Split the MST upload into V4 and V3-V4 uploads based on per-sample region detection.

Usage:
    conda run -n microbiome_16S python split_mst_uploads.py
"""
import csv
import sys
import uuid
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from app.config import UPLOAD_DIR
from app.db.database import init_db, SessionLocal
from app.db.models import FastqFile, Upload, UploadMetadata
from app.pipeline.detect import (
    PRIMERS, _primer_matches, _read_fastq_sequences,
    detect_variable_region, detect_sequencing_type, extract_sample_name,
)

FASTQ_DIR = Path.home() / "github/MST-Pipeline/DB_dev/data"
METADATA_FILE = Path.home() / "github/MST-Pipeline/DB/MST.design"


def main():
    init_db()
    db = SessionLocal()

    try:
        # Find current MST upload
        old_upload = (
            db.query(Upload)
            .filter(Upload.source_dir.contains("MST-Pipeline"))
            .order_by(Upload.id.desc())
            .first()
        )
        if not old_upload:
            print("No MST upload found.")
            return

        old_id = old_upload.id
        print(f"Found upload {old_id} with {old_upload.total_files} files")

        # Get all files grouped by sample
        all_files = db.query(FastqFile).filter(FastqFile.upload_id == old_id).all()
        samples = defaultdict(list)
        for f in all_files:
            samples[f.sample_name].append(f)

        # Detect region per sample using R1 file
        print(f"\nDetecting region for {len(samples)} samples...")
        region_groups = defaultdict(list)  # region -> [sample_name, ...]
        done = 0
        for sample_name, files in sorted(samples.items()):
            r1 = next((f for f in files if f.read_direction == "R1"), files[0])
            fpath = Path(r1.file_path)
            # Resolve symlink to actual file
            if fpath.is_symlink():
                fpath = fpath.resolve()
            result = detect_variable_region(fpath)
            region = result["region"] or "unknown"
            region_groups[region].append(sample_name)
            done += 1
            if done % 200 == 0:
                print(f"  {done}/{len(samples)} samples detected")

        print(f"  {done}/{len(samples)} done")
        print("\nRegion breakdown:")
        for region, snames in sorted(region_groups.items()):
            print(f"  {region}: {len(snames)} samples")

        # Load metadata
        meta_map = {}  # sample_name -> group
        with open(METADATA_FILE) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                meta_map[row["Sample"]] = row["Group"]

        old_study = old_upload.study

        # Delete old upload (cascades to FastqFile and UploadMetadata)
        print(f"\nDeleting old upload {old_id}...")
        # Also clean up the symlink directory
        old_dir = Path(old_upload.upload_dir)
        db.delete(old_upload)
        db.flush()

        # Remove old symlink directory
        if old_dir.exists():
            import shutil
            shutil.rmtree(old_dir.parent)  # remove {upload_id}/ dir
            print(f"  Removed {old_dir.parent}")

        # Create new uploads per region
        for region, sample_names in sorted(region_groups.items()):
            if region == "unknown":
                print(f"\nSkipping {len(sample_names)} samples with unknown region")
                continue

            sample_set = set(sample_names)
            upload_id = uuid.uuid4().hex[:8]
            dest_dir = UPLOAD_DIR / upload_id / "fastq"
            dest_dir.mkdir(parents=True, exist_ok=True)

            # Symlink files and build FastqFile records
            total_size = 0.0
            file_count = 0
            fastq_records = []

            for sname in sorted(sample_names):
                for old_f in samples[sname]:
                    # Resolve original path from old file_path (was a symlink)
                    orig = FASTQ_DIR / old_f.filename
                    dest = dest_dir / old_f.filename
                    dest.symlink_to(orig.resolve())
                    total_size += orig.stat().st_size / (1024 * 1024)
                    file_count += 1
                    fastq_records.append({
                        "sample_name": sname,
                        "filename": old_f.filename,
                        "file_path": str(dest),
                        "read_direction": old_f.read_direction,
                        "file_size_mb": old_f.file_size_mb,
                        "read_count": old_f.read_count,
                        "avg_read_length": old_f.avg_read_length,
                    })

            # Detect sequencing type from filenames
            filenames = [r["filename"] for r in fastq_records]
            detection = detect_sequencing_type(filenames)

            # Detect primers from first R1 file in this group
            primers_detected = None
            if region in PRIMERS:
                fwd_primer = PRIMERS[region]["forward"]
                first_r1 = next(
                    (r for r in fastq_records if r["read_direction"] == "R1"), None
                )
                if first_r1:
                    seqs = _read_fastq_sequences(Path(first_r1["file_path"]), n_reads=100)
                    if seqs:
                        matches = sum(1 for s in seqs if _primer_matches(s, fwd_primer))
                        primers_detected = (matches / len(seqs)) >= 0.30
                        print(f"  Primers: {matches}/{len(seqs)} matched {region} fwd → {primers_detected}")

            new_upload = Upload(
                source_dir=str(FASTQ_DIR),
                upload_dir=str(dest_dir),
                sequencing_type=detection["type"],
                variable_region=region,
                total_files=file_count,
                total_size_mb=round(total_size, 2),
                primers_detected=primers_detected,
                study=old_study,
                status="uploaded",
            )
            db.add(new_upload)
            db.flush()

            for rec in fastq_records:
                db.add(FastqFile(upload_id=new_upload.id, **rec))

            # Add metadata
            meta_count = 0
            for sname in sample_names:
                if sname in meta_map:
                    db.add(UploadMetadata(
                        upload_id=new_upload.id,
                        sample_name=sname,
                        key="source",
                        value=meta_map[sname],
                    ))
                    meta_count += 1

            db.flush()
            print(f"\nCreated upload {new_upload.id}:")
            print(f"  Region: {region}")
            print(f"  Samples: {len(sample_names)}")
            print(f"  Files: {file_count}")
            print(f"  Size: {total_size:.0f} MB")
            print(f"  Metadata: {meta_count} samples")

        db.commit()
        print("\nDone!")

    except Exception:
        db.rollback()
        raise
    finally:
        db.close()


if __name__ == "__main__":
    main()
