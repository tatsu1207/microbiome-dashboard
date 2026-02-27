"""
MicrobiomeDash — SQLAlchemy 2.0 ORM Models
"""
from datetime import datetime

from sqlalchemy import (
    Boolean,
    DateTime,
    Float,
    ForeignKey,
    Index,
    Integer,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship


class Base(DeclarativeBase):
    pass


# ── Projects ─────────────────────────────────────────────────────────────────


class Project(Base):
    __tablename__ = "projects"

    id: Mapped[int] = mapped_column(primary_key=True)
    name: Mapped[str] = mapped_column(String, nullable=False)
    description: Mapped[str | None] = mapped_column(Text)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow)
    updated_at: Mapped[datetime] = mapped_column(
        DateTime, default=datetime.utcnow, onupdate=datetime.utcnow
    )

    uploads: Mapped[list["Upload"]] = relationship(back_populates="project")
    datasets: Mapped[list["Dataset"]] = relationship(back_populates="project")


# ── Uploads ───────────────────────────────────────────────────────────────────


class Upload(Base):
    __tablename__ = "uploads"

    id: Mapped[int] = mapped_column(primary_key=True)
    project_id: Mapped[int | None] = mapped_column(ForeignKey("projects.id"))
    source_dir: Mapped[str | None] = mapped_column(String)
    upload_dir: Mapped[str] = mapped_column(String, nullable=False)
    sequencing_type: Mapped[str | None] = mapped_column(String)
    variable_region: Mapped[str | None] = mapped_column(String)
    total_files: Mapped[int | None] = mapped_column(Integer)
    total_size_mb: Mapped[float | None] = mapped_column(Float)
    primers_detected: Mapped[bool | None] = mapped_column(Boolean, nullable=True)
    study: Mapped[str | None] = mapped_column(String)
    status: Mapped[str] = mapped_column(String, default="uploaded")
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow)

    project: Mapped["Project | None"] = relationship(back_populates="uploads")
    fastq_files: Mapped[list["FastqFile"]] = relationship(
        back_populates="upload", cascade="all, delete-orphan"
    )
    upload_metadata: Mapped[list["UploadMetadata"]] = relationship(
        back_populates="upload", cascade="all, delete-orphan"
    )


# ── FASTQ Files ───────────────────────────────────────────────────────────────


class FastqFile(Base):
    __tablename__ = "fastq_files"
    __table_args__ = (Index("idx_fastq_upload", "upload_id"),)

    id: Mapped[int] = mapped_column(primary_key=True)
    upload_id: Mapped[int] = mapped_column(
        ForeignKey("uploads.id", ondelete="CASCADE")
    )
    sample_name: Mapped[str] = mapped_column(String, nullable=False)
    filename: Mapped[str] = mapped_column(String, nullable=False)
    file_path: Mapped[str] = mapped_column(String, nullable=False)
    read_direction: Mapped[str | None] = mapped_column(String)
    file_size_mb: Mapped[float | None] = mapped_column(Float)
    read_count: Mapped[int | None] = mapped_column(Integer)
    avg_read_length: Mapped[int | None] = mapped_column(Integer)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow)

    upload: Mapped["Upload"] = relationship(back_populates="fastq_files")


# ── Upload Metadata (key-value, pre-pipeline) ────────────────────────────


class UploadMetadata(Base):
    __tablename__ = "upload_metadata"
    __table_args__ = (
        UniqueConstraint("upload_id", "sample_name", "key"),
        Index("idx_upload_metadata_upload", "upload_id"),
    )

    id: Mapped[int] = mapped_column(primary_key=True)
    upload_id: Mapped[int] = mapped_column(
        ForeignKey("uploads.id", ondelete="CASCADE")
    )
    sample_name: Mapped[str] = mapped_column(String, nullable=False)
    key: Mapped[str] = mapped_column(String, nullable=False)
    value: Mapped[str | None] = mapped_column(Text)

    upload: Mapped["Upload"] = relationship(back_populates="upload_metadata")


# ── Dataset–FastqFile association ─────────────────────────────────────────────


class DatasetFastqFile(Base):
    __tablename__ = "dataset_fastq_files"
    __table_args__ = (
        UniqueConstraint("dataset_id", "fastq_file_id"),
        Index("idx_dsff_dataset", "dataset_id"),
    )

    id: Mapped[int] = mapped_column(primary_key=True)
    dataset_id: Mapped[int] = mapped_column(
        ForeignKey("datasets.id", ondelete="CASCADE")
    )
    fastq_file_id: Mapped[int] = mapped_column(
        ForeignKey("fastq_files.id", ondelete="CASCADE")
    )

    dataset: Mapped["Dataset"] = relationship(back_populates="fastq_file_links")
    fastq_file: Mapped["FastqFile"] = relationship()


# ── Datasets ──────────────────────────────────────────────────────────────────


class Dataset(Base):
    __tablename__ = "datasets"
    __table_args__ = (Index("idx_dataset_project", "project_id"),)

    id: Mapped[int] = mapped_column(primary_key=True)
    project_id: Mapped[int | None] = mapped_column(ForeignKey("projects.id"))
    upload_id: Mapped[int | None] = mapped_column(ForeignKey("uploads.id"))
    name: Mapped[str] = mapped_column(String, nullable=False)
    description: Mapped[str | None] = mapped_column(Text)
    source_type: Mapped[str] = mapped_column(String, default="pipeline")
    status: Mapped[str] = mapped_column(String, default="pending")
    variable_region: Mapped[str | None] = mapped_column(String)
    sequencing_type: Mapped[str | None] = mapped_column(String)
    sample_count: Mapped[int | None] = mapped_column(Integer)
    asv_count: Mapped[int | None] = mapped_column(Integer)
    asv_table_path: Mapped[str | None] = mapped_column(String)
    taxonomy_path: Mapped[str | None] = mapped_column(String)
    tree_path: Mapped[str | None] = mapped_column(String)
    rep_seqs_path: Mapped[str | None] = mapped_column(String)
    picrust_dir_path: Mapped[str | None] = mapped_column(String)
    pipeline_log_path: Mapped[str | None] = mapped_column(String)
    trim_left_f: Mapped[int | None] = mapped_column(Integer)
    trim_left_r: Mapped[int | None] = mapped_column(Integer)
    trunc_len_f: Mapped[int | None] = mapped_column(Integer)
    trunc_len_r: Mapped[int | None] = mapped_column(Integer)
    min_overlap: Mapped[int | None] = mapped_column(Integer)
    custom_fwd_primer: Mapped[str | None] = mapped_column(String)
    custom_rev_primer: Mapped[str | None] = mapped_column(String)
    parent_dataset_ids: Mapped[str | None] = mapped_column(Text)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow)
    updated_at: Mapped[datetime] = mapped_column(
        DateTime, default=datetime.utcnow, onupdate=datetime.utcnow
    )

    project: Mapped["Project | None"] = relationship(back_populates="datasets")
    upload: Mapped["Upload | None"] = relationship()
    samples: Mapped[list["Sample"]] = relationship(
        back_populates="dataset", cascade="all, delete-orphan"
    )
    analysis_results: Mapped[list["AnalysisResult"]] = relationship(
        back_populates="dataset", cascade="all, delete-orphan"
    )
    fastq_file_links: Mapped[list["DatasetFastqFile"]] = relationship(
        back_populates="dataset", cascade="all, delete-orphan"
    )
    combinations: Mapped[list["DatasetCombination"]] = relationship(
        back_populates="combined_dataset",
        cascade="all, delete-orphan",
        foreign_keys="DatasetCombination.combined_dataset_id",
    )


# ── Samples ───────────────────────────────────────────────────────────────────


class Sample(Base):
    __tablename__ = "samples"
    __table_args__ = (Index("idx_samples_dataset", "dataset_id"),)

    id: Mapped[int] = mapped_column(primary_key=True)
    dataset_id: Mapped[int] = mapped_column(
        ForeignKey("datasets.id", ondelete="CASCADE")
    )
    sample_name: Mapped[str] = mapped_column(String, nullable=False)
    source_study: Mapped[str | None] = mapped_column(String)
    read_count_raw: Mapped[int | None] = mapped_column(Integer)
    read_count_filtered: Mapped[int | None] = mapped_column(Integer)
    read_count_denoised: Mapped[int | None] = mapped_column(Integer)
    read_count_nonchimeric: Mapped[int | None] = mapped_column(Integer)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow)

    dataset: Mapped["Dataset"] = relationship(back_populates="samples")
    metadata_entries: Mapped[list["SampleMetadata"]] = relationship(
        back_populates="sample", cascade="all, delete-orphan"
    )
    qc_metrics: Mapped[list["QcMetric"]] = relationship(
        back_populates="sample", cascade="all, delete-orphan"
    )


# ── Sample Metadata (key-value) ──────────────────────────────────────────────


class SampleMetadata(Base):
    __tablename__ = "sample_metadata"
    __table_args__ = (
        UniqueConstraint("sample_id", "key"),
        Index("idx_metadata_sample", "sample_id"),
    )

    id: Mapped[int] = mapped_column(primary_key=True)
    sample_id: Mapped[int] = mapped_column(
        ForeignKey("samples.id", ondelete="CASCADE")
    )
    key: Mapped[str] = mapped_column(String, nullable=False)
    value: Mapped[str | None] = mapped_column(Text)

    sample: Mapped["Sample"] = relationship(back_populates="metadata_entries")


# ── QC Metrics ────────────────────────────────────────────────────────────────


class QcMetric(Base):
    __tablename__ = "qc_metrics"
    __table_args__ = (Index("idx_qc_sample", "sample_id"),)

    id: Mapped[int] = mapped_column(primary_key=True)
    sample_id: Mapped[int] = mapped_column(
        ForeignKey("samples.id", ondelete="CASCADE")
    )
    metric_name: Mapped[str] = mapped_column(String, nullable=False)
    metric_value: Mapped[float | None] = mapped_column(Float)
    read_direction: Mapped[str | None] = mapped_column(String)

    sample: Mapped["Sample"] = relationship(back_populates="qc_metrics")


# ── Dataset Combinations ─────────────────────────────────────────────────────


class DatasetCombination(Base):
    __tablename__ = "dataset_combinations"

    id: Mapped[int] = mapped_column(primary_key=True)
    combined_dataset_id: Mapped[int] = mapped_column(
        ForeignKey("datasets.id", ondelete="CASCADE")
    )
    source_dataset_id: Mapped[int | None] = mapped_column(ForeignKey("datasets.id"))
    source_dataset_name: Mapped[str | None] = mapped_column(String)
    sample_count: Mapped[int | None] = mapped_column(Integer)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow)

    combined_dataset: Mapped["Dataset"] = relationship(
        back_populates="combinations",
        foreign_keys=[combined_dataset_id],
    )
    source_dataset: Mapped["Dataset | None"] = relationship(
        foreign_keys=[source_dataset_id],
    )


# ── Analysis Results ──────────────────────────────────────────────────────────


class AnalysisResult(Base):
    __tablename__ = "analysis_results"
    __table_args__ = (Index("idx_analysis_dataset", "dataset_id"),)

    id: Mapped[int] = mapped_column(primary_key=True)
    dataset_id: Mapped[int] = mapped_column(
        ForeignKey("datasets.id", ondelete="CASCADE")
    )
    analysis_type: Mapped[str] = mapped_column(String, nullable=False)
    parameters: Mapped[str | None] = mapped_column(Text)
    result_path: Mapped[str | None] = mapped_column(String)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow)

    dataset: Mapped["Dataset"] = relationship(back_populates="analysis_results")


# ── PICRUSt2 Standalone Runs ────────────────────────────────────────────────


class Picrust2Run(Base):
    __tablename__ = "picrust2_runs"

    id: Mapped[int] = mapped_column(primary_key=True)
    name: Mapped[str] = mapped_column(String, nullable=False)
    status: Mapped[str] = mapped_column(String, default="pending")
    biom_path: Mapped[str] = mapped_column(String, nullable=False)
    fasta_path: Mapped[str | None] = mapped_column(String, nullable=True)
    skip_ko: Mapped[bool] = mapped_column(Boolean, default=False)
    output_dir: Mapped[str | None] = mapped_column(String)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow)
