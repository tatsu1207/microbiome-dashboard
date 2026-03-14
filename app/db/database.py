"""
MicrobiomeDash — Database engine, session management, and initialization.
"""
from contextlib import contextmanager

from sqlalchemy import create_engine
from sqlalchemy.orm import Session, sessionmaker

from app.config import (
    COMBINED_DIR,
    DATABASE_URL,
    DATASET_DIR,
    EXPORT_DIR,
    SRA_CACHE_DIR,
    UPLOAD_DIR,
)
from app.db.models import Base

engine = create_engine(
    DATABASE_URL,
    connect_args={"check_same_thread": False, "timeout": 30},
)

SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine, expire_on_commit=False)


def get_db():
    """FastAPI dependency: yields a DB session, auto-closes."""
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


@contextmanager
def get_session() -> Session:
    """Context manager for Dash callbacks (which bypass FastAPI DI)."""
    db = SessionLocal()
    try:
        yield db
        db.commit()
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()


def init_db():
    """Create all tables and ensure data directories exist."""
    Base.metadata.create_all(bind=engine)
    # Lightweight migrations for columns added after initial schema
    _migrate_add_columns()
    _backfill_primers_detected()
    for d in [UPLOAD_DIR, DATASET_DIR, COMBINED_DIR, EXPORT_DIR, SRA_CACHE_DIR]:
        d.mkdir(parents=True, exist_ok=True)


def _migrate_add_columns():
    """Add columns that were introduced after the initial schema."""
    import sqlite3
    url = str(engine.url)
    if not url.startswith("sqlite"):
        return
    db_path = url.replace("sqlite:///", "")
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    # (table, column, type)
    migrations = [
        ("fastq_files", "avg_read_length", "INTEGER"),
        ("fastq_files", "read_count", "INTEGER"),
        ("uploads", "primers_detected", "BOOLEAN"),
        ("datasets", "custom_fwd_primer", "TEXT"),
        ("datasets", "custom_rev_primer", "TEXT"),
        ("uploads", "platform", "TEXT"),
        ("datasets", "platform", "TEXT"),
    ]
    for table, col, col_type in migrations:
        try:
            cursor.execute(f"ALTER TABLE {table} ADD COLUMN {col} {col_type}")
        except sqlite3.OperationalError:
            pass  # column already exists
    # Drop columns removed from the schema
    drops = [
        ("uploads", "metadata_path"),
        ("datasets", "metadata_path"),
    ]
    for table, col in drops:
        try:
            cursor.execute(f"ALTER TABLE {table} DROP COLUMN {col}")
        except sqlite3.OperationalError:
            pass  # column already gone
    conn.commit()
    conn.close()


def _backfill_primers_detected():
    """Populate primers_detected for existing uploads that have NULL."""
    from pathlib import Path
    from app.db.models import FastqFile, Upload
    from app.pipeline.detect import PRIMERS, _primer_matches, _read_fastq_sequences

    db = SessionLocal()
    try:
        uploads = db.query(Upload).filter(
            Upload.primers_detected.is_(None),
            Upload.variable_region.isnot(None),
        ).all()
        for upload in uploads:
            region = upload.variable_region
            if region not in PRIMERS:
                continue
            # Find the first R1 file for this upload
            r1 = (
                db.query(FastqFile)
                .filter(FastqFile.upload_id == upload.id, FastqFile.read_direction == "R1")
                .first()
            )
            if not r1:
                r1 = (
                    db.query(FastqFile)
                    .filter(FastqFile.upload_id == upload.id)
                    .first()
                )
            if not r1 or not Path(r1.file_path).exists():
                continue
            try:
                fwd_primer = PRIMERS[region]["forward"]
                seqs = _read_fastq_sequences(Path(r1.file_path), n_reads=100)
                if seqs:
                    matches = sum(1 for s in seqs if _primer_matches(s, fwd_primer))
                    upload.primers_detected = (matches / len(seqs)) >= 0.30
            except Exception:
                pass
        db.commit()
    except Exception:
        db.rollback()
    finally:
        db.close()
