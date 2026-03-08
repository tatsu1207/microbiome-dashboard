#!/usr/bin/env python3
"""Backfill read_count, avg_read_length, primers_detected, and study for
the MST upload that was registered without these fields.

Usage:
    conda run -n microbiome_16S python backfill_mst.py
"""
import gzip
import sys
from collections import defaultdict
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, str(Path(__file__).resolve().parent))

from app.db.database import init_db, SessionLocal
from app.db.models import FastqFile, Upload
from app.pipeline.detect import PRIMERS, _primer_matches, _read_fastq_sequences


def _count_reads_and_avg_len(file_path: str) -> tuple[int | None, int | None]:
    """Count reads and compute avg read length for one FASTQ file."""
    p = Path(file_path)
    if not p.exists():
        return None, None
    try:
        opener = gzip.open if str(p).endswith(".gz") else open
        total_lines = 0
        lengths = []
        sample_done = False
        with opener(p, "rt") as f:
            for line in f:
                total_lines += 1
                if not sample_done and total_lines % 4 == 2:
                    lengths.append(len(line.strip()))
                    if len(lengths) >= 200:
                        sample_done = True
        read_count = total_lines // 4
        avg_len = round(sum(lengths) / len(lengths)) if lengths else None
        return read_count, avg_len
    except Exception as e:
        print(f"  Error reading {p.name}: {e}")
        return None, None


def _process_file(args: tuple[int, str, str]) -> tuple[int, int | None, int | None]:
    """Worker: returns (file_id, read_count, avg_read_length)."""
    file_id, file_path, _filename = args
    rc, al = _count_reads_and_avg_len(file_path)
    return file_id, rc, al


def main():
    init_db()
    db = SessionLocal()

    try:
        # Find the MST upload (most recent, or the one with source_dir matching)
        upload = (
            db.query(Upload)
            .filter(Upload.source_dir.contains("MST-Pipeline"))
            .order_by(Upload.id.desc())
            .first()
        )
        if not upload:
            print("No MST upload found in database.")
            return

        print(f"Upload ID: {upload.id}, files: {upload.total_files}")

        # 1. Set study = "MST"
        upload.study = "MST"
        db.flush()
        print("Set study = 'MST'")

        # 2. Detect primers
        print("\nDetecting primers...")
        region = upload.variable_region
        if region and region in PRIMERS:
            fwd_primer = PRIMERS[region]["forward"]
            r1 = (
                db.query(FastqFile)
                .filter(FastqFile.upload_id == upload.id, FastqFile.read_direction == "R1")
                .first()
            )
            if r1 and Path(r1.file_path).exists():
                seqs = _read_fastq_sequences(Path(r1.file_path), n_reads=100)
                if seqs:
                    matches = sum(1 for s in seqs if _primer_matches(s, fwd_primer))
                    upload.primers_detected = (matches / len(seqs)) >= 0.30
                    print(f"  Primers detected: {upload.primers_detected} ({matches}/{len(seqs)} matched {region} fwd)")
        db.flush()

        # 3. Compute read_count + avg_read_length for all files (parallel)
        files = (
            db.query(FastqFile.id, FastqFile.file_path, FastqFile.filename)
            .filter(FastqFile.upload_id == upload.id)
            .all()
        )
        total = len(files)
        print(f"\nCounting reads & avg length for {total} files (parallel)...")

        done = 0
        workers = min(8, total)
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_process_file, (f.id, f.file_path, f.filename)): f.id for f in files}
            for future in as_completed(futures):
                file_id, rc, al = future.result()
                db.query(FastqFile).filter(FastqFile.id == file_id).update(
                    {"read_count": rc, "avg_read_length": al}
                )
                done += 1
                if done % 200 == 0:
                    db.flush()
                    print(f"  {done}/{total} files processed")

        db.commit()
        print(f"  {done}/{total} files processed — done!")

        # Summary
        sample_reads = defaultdict(int)
        for f in db.query(FastqFile).filter(FastqFile.upload_id == upload.id).all():
            sample_reads[f.sample_name] += f.read_count or 0
        total_reads = sum(sample_reads.values())
        print(f"\nSummary:")
        print(f"  Total reads across all files: {total_reads:,}")
        print(f"  Study: {upload.study}")
        print(f"  Primers detected: {upload.primers_detected}")
        print(f"  Region: {upload.variable_region}")

    except Exception:
        db.rollback()
        raise
    finally:
        db.close()


if __name__ == "__main__":
    main()
