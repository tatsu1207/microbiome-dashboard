# MicrobiomeDash — Claude Code Briefing Document

## What This Is

A **three-tool web application** for 16S rRNA microbiome analysis, built with **Plotly Dash + FastAPI + SQLite**, running on **Windows WSL**. The user is a bioinformatician building a personal/hobby project.

---

## What's Already Done

The project scaffold is complete at `~/microbiome-dashboard/`:

- `setup.sh` — Automated installer (detects/installs mamba, creates conda env, installs all deps, skip-if-installed logic)
- `run.sh` — Starts the app with `uvicorn`
- `README.md` — Full setup guide including WSL installation from scratch
- `requirements.txt` — Python pip dependencies
- `environment.yml` — Conda environment spec (Python 3.11, R 4.3, bioinformatics tools)
- `Makefile` — Common dev commands (run, setup, db-reset, clean)
- `.gitignore`, `LICENSE` (MIT)
- `app/config.py` — Auto-generated config with paths, DB URL, SILVA paths, DADA2 defaults
- Empty directory structure with `__init__.py` files:
  ```
  app/api/
  app/pipeline/
  app/data_manager/
  app/analysis/
  app/dashboard/pages/
  app/dashboard/components/
  app/db/
  app/utils/
  r_scripts/
  data/uploads/
  data/datasets/
  data/combined/
  data/exports/
  data/references/
  ```

**No application code has been written yet.** All three tools need to be built from scratch.

---

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│               TOOL 1: Pipeline Engine                            │
│          "FASTQ → Microbiome Dataset Builder"                   │
│                                                                  │
│  Input: fastq.gz + metadata (CSV/TSV)                           │
│                                                                  │
│  ┌────────────┐   ┌──────────┐   ┌────────┐   ┌─────────────┐ │
│  │ Auto-Detect │──▶│  FastQC   │──▶│Cutadapt│──▶│   DADA2     │ │
│  │ • SE / PE   │   │  (QC)     │   │ (Trim) │   │  (Denoise)  │ │
│  │ • Var Region│   └──────────┘   └────────┘   └─────────────┘ │
│  │   (V1-V9)  │                                      │          │
│  └────────────┘                                      ▼          │
│                    ┌───────────┐   ┌──────────┐   ┌──────────┐ │
│                    │ PICRUSt2  │◀──│Phylogeny │◀──│ Taxonomy │ │
│                    │ (Pathways)│   │MAFFT+FT  │   │ SILVA138 │ │
│                    └───────────┘   └──────────┘   └──────────┘ │
│                                                                  │
│  Output: Microbiome Dataset (ASV table, taxonomy, tree,         │
│          PICRUSt2 output, QC metrics, variable region info)     │
└──────────────────────────────┬──────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────┐
│               TOOL 2: Data Management Hub                        │
│          "Organize, Browse, Combine, Export"                    │
│                                                                  │
│  • File Manager: browse/delete raw FASTQ & metadata files       │
│  • Dataset Manager: browse/inspect/download/delete datasets     │
│    - Tracks variable region per dataset                         │
│    - Download as ZIP (ASV + taxonomy + tree + PICRUSt2 + meta)  │
│  • Dataset Combiner: merge datasets from different studies      │
│    - Requires same variable region                              │
│    - Merges ASV tables, metadata, rebuilds tree                 │
│    - Optionally re-runs taxonomy & PICRUSt2                     │
└──────────────────────────────┬──────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────┐
│               TOOL 3: Analysis & Visualization                   │
│          "Microbiome Dataset → Insights"                        │
│                                                                  │
│  • Alpha diversity (Shannon, Simpson, Chao1, Faith's PD)        │
│  • Beta diversity (PCoA, NMDS, PERMANOVA)                       │
│  • Taxonomic composition (stacked bar, heatmap, sunburst)       │
│  • Differential abundance (ANCOM-BC / ALDEx2, volcano plots)    │
│  • PICRUSt2 pathway visualization                               │
│  • Export: publication-ready figures (SVG/PDF) + tables (CSV)   │
└─────────────────────────────────────────────────────────────────┘
```

---

## Key Technical Decisions (Already Agreed)

| Decision | Choice | Reason |
|----------|--------|--------|
| Backend framework | **FastAPI** | Async, fast, good Python integration |
| Dashboard UI | **Plotly Dash** (mounted inside FastAPI) | Python-native, interactive, science-friendly |
| Database | **SQLite** via SQLAlchemy | Zero setup, file-based, personal project |
| DADA2 integration | **R scripts via subprocess** (NOT rpy2) | Cleaner, easier to debug, independently testable |
| R stat packages | **ANCOM-BC / ALDEx2 via R subprocess** | Same approach as DADA2 |
| Taxonomy reference | **SILVA 138.1** | Gold standard for 16S |
| PICRUSt2 | **Subprocess** (may be in separate conda env) | Standalone CLI, often has dependency conflicts |
| SE/PE detection | **Auto-detect from filenames** (_R1/_R2, _1/_2 patterns) | User requested |
| Variable region | **Auto-detect from primer sequences + amplicon length** | User requested, stored per upload and dataset |
| Dataset combination | **Same variable region required** | Mixing regions produces meaningless results |
| Package manager | **Mamba** (auto-installed if missing, falls back to conda) | Faster dependency resolution |

---

## Database Schema

```sql
-- Projects group related datasets
CREATE TABLE projects (
    id              INTEGER PRIMARY KEY,
    name            TEXT NOT NULL,
    description     TEXT,
    created_at      DATETIME DEFAULT CURRENT_TIMESTAMP,
    updated_at      DATETIME DEFAULT CURRENT_TIMESTAMP
);

-- Each upload batch (a set of FASTQ files + metadata)
CREATE TABLE uploads (
    id              INTEGER PRIMARY KEY,
    project_id      INTEGER REFERENCES projects(id),
    upload_dir      TEXT NOT NULL,
    metadata_path   TEXT,
    sequencing_type TEXT,                   -- 'single-end' | 'paired-end' (auto-detected)
    variable_region TEXT,                   -- 'V1-V2' | 'V3-V4' | 'V4' | etc. (auto-detected)
    total_files     INTEGER,
    total_size_mb   REAL,
    status          TEXT DEFAULT 'uploaded', -- uploaded | processing | complete
    created_at      DATETIME DEFAULT CURRENT_TIMESTAMP
);

-- Individual FASTQ files within an upload
CREATE TABLE fastq_files (
    id              INTEGER PRIMARY KEY,
    upload_id       INTEGER REFERENCES uploads(id) ON DELETE CASCADE,
    sample_name     TEXT NOT NULL,
    filename        TEXT NOT NULL,
    file_path       TEXT NOT NULL,
    read_direction  TEXT,                   -- 'R1' | 'R2' | 'single'
    file_size_mb    REAL,
    read_count      INTEGER,
    created_at      DATETIME DEFAULT CURRENT_TIMESTAMP
);

-- Each pipeline run produces a dataset
CREATE TABLE datasets (
    id                  INTEGER PRIMARY KEY,
    project_id          INTEGER REFERENCES projects(id),
    upload_id           INTEGER REFERENCES uploads(id),
    name                TEXT NOT NULL,
    description         TEXT,
    source_type         TEXT DEFAULT 'pipeline',  -- 'pipeline' | 'combined' | 'imported'
    status              TEXT DEFAULT 'pending',   -- pending | processing | complete | failed
    variable_region     TEXT,
    sequencing_type     TEXT,
    sample_count        INTEGER,
    asv_count           INTEGER,
    asv_table_path      TEXT,
    taxonomy_path       TEXT,
    tree_path           TEXT,
    rep_seqs_path       TEXT,
    picrust_dir_path    TEXT,
    metadata_path       TEXT,
    pipeline_log_path   TEXT,
    trim_left_f         INTEGER,
    trim_left_r         INTEGER,
    trunc_len_f         INTEGER,
    trunc_len_r         INTEGER,
    min_overlap         INTEGER,
    parent_dataset_ids  TEXT,               -- JSON array for combined datasets
    created_at          DATETIME DEFAULT CURRENT_TIMESTAMP,
    updated_at          DATETIME DEFAULT CURRENT_TIMESTAMP
);

-- Individual samples within a dataset
CREATE TABLE samples (
    id                      INTEGER PRIMARY KEY,
    dataset_id              INTEGER REFERENCES datasets(id) ON DELETE CASCADE,
    sample_name             TEXT NOT NULL,
    source_study            TEXT,           -- for combined datasets
    read_count_raw          INTEGER,
    read_count_filtered     INTEGER,
    read_count_denoised     INTEGER,
    read_count_nonchimeric  INTEGER,
    created_at              DATETIME DEFAULT CURRENT_TIMESTAMP
);

-- Flexible key-value metadata per sample
CREATE TABLE sample_metadata (
    id          INTEGER PRIMARY KEY,
    sample_id   INTEGER REFERENCES samples(id) ON DELETE CASCADE,
    key         TEXT NOT NULL,
    value       TEXT,
    UNIQUE(sample_id, key)
);

-- QC metrics per sample
CREATE TABLE qc_metrics (
    id              INTEGER PRIMARY KEY,
    sample_id       INTEGER REFERENCES samples(id) ON DELETE CASCADE,
    metric_name     TEXT NOT NULL,
    metric_value    REAL,
    read_direction  TEXT                    -- 'forward' | 'reverse' | 'single'
);

-- Dataset combination history
CREATE TABLE dataset_combinations (
    id                  INTEGER PRIMARY KEY,
    combined_dataset_id INTEGER REFERENCES datasets(id) ON DELETE CASCADE,
    source_dataset_id   INTEGER REFERENCES datasets(id),
    source_dataset_name TEXT,
    sample_count        INTEGER,
    created_at          DATETIME DEFAULT CURRENT_TIMESTAMP
);

-- Cached analysis results
CREATE TABLE analysis_results (
    id              INTEGER PRIMARY KEY,
    dataset_id      INTEGER REFERENCES datasets(id) ON DELETE CASCADE,
    analysis_type   TEXT NOT NULL,           -- 'alpha_div', 'beta_div', 'diff_abundance', 'picrust'
    parameters      TEXT,                    -- JSON
    result_path     TEXT,
    created_at      DATETIME DEFAULT CURRENT_TIMESTAMP
);

-- Indexes
CREATE INDEX idx_fastq_upload ON fastq_files(upload_id);
CREATE INDEX idx_dataset_project ON datasets(project_id);
CREATE INDEX idx_samples_dataset ON samples(dataset_id);
CREATE INDEX idx_metadata_sample ON sample_metadata(sample_id);
CREATE INDEX idx_qc_sample ON qc_metrics(sample_id);
CREATE INDEX idx_analysis_dataset ON analysis_results(dataset_id);
```

---

## Project File Map

```
microbiome-dashboard/
├── app/
│   ├── main.py                      # FastAPI app + Dash mounting
│   ├── config.py                    # ✅ EXISTS — auto-generated paths/settings
│   │
│   ├── api/                         # FastAPI REST endpoints
│   │   ├── upload.py                # File upload
│   │   ├── pipeline.py              # Pipeline trigger/status
│   │   ├── datasets.py              # Dataset CRUD
│   │   ├── data_manager.py          # Browse, download, combine
│   │   └── analysis.py              # Analysis trigger
│   │
│   ├── pipeline/                    # TOOL 1
│   │   ├── detect.py                # Auto-detect SE/PE + variable region
│   │   ├── qc.py                    # FastQC wrapper
│   │   ├── trim.py                  # Cutadapt wrapper
│   │   ├── dada2.py                 # Calls r_scripts/run_dada2.R
│   │   ├── taxonomy.py              # Calls r_scripts/run_taxonomy.R
│   │   ├── phylogeny.py             # MAFFT + FastTree subprocess
│   │   ├── picrust2.py              # PICRUSt2 subprocess
│   │   └── runner.py                # Pipeline orchestrator
│   │
│   ├── data_manager/                # TOOL 2
│   │   ├── file_browser.py          # Browse/delete FASTQ & metadata
│   │   ├── dataset_browser.py       # Browse/inspect/delete datasets
│   │   ├── dataset_download.py      # ZIP packaging for download
│   │   ├── dataset_combiner.py      # Merge datasets across studies
│   │   └── compatibility.py         # Variable region compatibility check
│   │
│   ├── analysis/                    # TOOL 3
│   │   ├── alpha_diversity.py
│   │   ├── beta_diversity.py
│   │   ├── diff_abundance.py        # Calls r_scripts/run_ancombc.R or run_aldex2.R
│   │   ├── taxonomy_comp.py
│   │   ├── picrust_analysis.py
│   │   └── stats.py                 # PERMANOVA, Kruskal-Wallis
│   │
│   ├── dashboard/                   # Plotly Dash UI
│   │   ├── app.py                   # Dash init
│   │   ├── layout.py                # Main layout with sidebar nav
│   │   ├── pages/
│   │   │   ├── upload.py            # Upload & metadata page (Tool 1)
│   │   │   ├── pipeline_status.py   # Pipeline progress (Tool 1)
│   │   │   ├── qc_report.py         # QC visualization (Tool 1)
│   │   │   ├── file_manager.py      # File browser (Tool 2)
│   │   │   ├── dataset_manager.py   # Dataset browser (Tool 2)
│   │   │   ├── dataset_combiner.py  # Combine UI (Tool 2)
│   │   │   ├── alpha_div.py         # Alpha diversity (Tool 3)
│   │   │   ├── beta_div.py          # Beta diversity (Tool 3)
│   │   │   ├── taxonomy.py          # Taxonomic composition (Tool 3)
│   │   │   ├── diff_abundance.py    # Differential abundance (Tool 3)
│   │   │   └── picrust.py           # Pathway analysis (Tool 3)
│   │   └── components/
│   │       ├── filters.py           # Shared filter widgets
│   │       ├── tables.py            # Reusable data tables
│   │       ├── plots.py             # Common plot helpers
│   │       └── dataset_selector.py  # Reusable dataset picker
│   │
│   ├── db/
│   │   ├── models.py                # SQLAlchemy ORM models
│   │   ├── database.py              # Engine, session, init_db()
│   │   └── crud.py                  # Database operations
│   │
│   └── utils/
│       ├── file_handler.py          # FASTQ parsing, file management
│       ├── metadata_parser.py       # Metadata CSV/TSV validation
│       ├── variable_region.py       # Variable region detection logic
│       └── zip_export.py            # ZIP packaging for dataset download
│
├── r_scripts/
│   ├── run_dada2.R                  # DADA2 pipeline (accepts --mode paired|single + CLI args)
│   ├── run_taxonomy.R               # SILVA 138.1 assignTaxonomy + addSpecies
│   ├── run_ancombc.R                # ANCOM-BC differential abundance
│   └── run_aldex2.R                 # ALDEx2 differential abundance
│
├── data/
│   ├── uploads/{upload_id}/fastq/ + metadata.tsv
│   ├── datasets/{dataset_id}/asv_table.tsv, taxonomy.tsv, tree.nwk, etc.
│   ├── combined/{combined_id}/...
│   ├── exports/                     # Temporary ZIP files
│   └── references/silva_*.fa.gz     # SILVA 138.1 (downloaded by setup.sh)
│
├── setup.sh                         # ✅ EXISTS
├── run.sh                           # ✅ EXISTS
├── README.md                        # ✅ EXISTS
├── requirements.txt                 # ✅ EXISTS
├── environment.yml                  # ✅ EXISTS
├── Makefile                         # ✅ EXISTS
├── .gitignore                       # ✅ EXISTS
└── LICENSE                          # ✅ EXISTS
```

---

## Pipeline Detail: Auto-Detection

### SE/PE Detection (from filenames)
```
Paired-end patterns:  *_R1_*.fastq.gz / *_R2_*.fastq.gz
                      *_1.fastq.gz   / *_2.fastq.gz
                      *_R1.fastq.gz  / *_R2.fastq.gz
Single-end:           Everything else (no R2 pair found)

Logic: Group files by sample name → if all have R1+R2 = paired; all have 1 file = single; mixed = error
```

### Variable Region Detection
```
1. PRIMER MATCHING (primary): Read first 100 reads, match against known 16S primers:
   V1-V2: 27F  (AGAGTTTGATCMTGGCTCAG) / 338R
   V3-V4: 341F (CCTACGGGNGGCWGCAG)    / 805R
   V4:    515F (GTGYCAGCMGCCGCGGTAA)  / 806R
   V4-V5: 515F / 926R
   V5-V6: 784F / 1061R
   Allow 2-3 mismatches for degenerate bases.

2. AMPLICON LENGTH (secondary validation):
   V1-V2: ~300-350bp, V3-V4: ~460bp, V4: ~253bp, V4-V5: ~400bp

3. USER OVERRIDE: Show detected region with confidence, allow manual change.
```

### R Script Interface (DADA2)
```
Rscript r_scripts/run_dada2.R \
  --input_dir /path/to/fastq/ \
  --output_dir /path/to/output/ \
  --mode paired|single \
  --trim_left_f 17 \
  --trim_left_r 21 \
  --trunc_len_f 250 \
  --trunc_len_r 200 \
  --silva_train /path/to/silva_nr99_v138.1_train_set.fa.gz \
  --silva_species /path/to/silva_species_assignment_v138.1.fa.gz \
  --threads 4

Outputs: asv_table.tsv, taxonomy.tsv, representative_seqs.fasta, track_reads.tsv
```

### Dataset Combination Rules
- Same variable region required (enforced)
- Merge ASV tables (prefix ASV IDs with dataset ID to avoid collision)
- Merge metadata with `source_study` column added
- Optionally re-run: taxonomy, phylogeny (MAFFT+FastTree), PICRUSt2
- Combined dataset gets `source_type = 'combined'` and `parent_dataset_ids` JSON

---

## Phased Build Plan

### Phase 1 — Foundation (current priority)
- [ ] SQLAlchemy models in `app/db/models.py` (from schema above)
- [ ] Database engine/session in `app/db/database.py` with `init_db()`
- [ ] FastAPI + Dash integration in `app/main.py` (Dash mounted as sub-app)
- [ ] Dash layout shell in `app/dashboard/layout.py` (sidebar with 3 tool sections)
- [ ] Upload page `app/dashboard/pages/upload.py` (drag & drop FASTQ + metadata)
- [ ] Auto-detect logic in `app/pipeline/detect.py` (SE/PE + variable region)
- [ ] Metadata parser in `app/utils/metadata_parser.py` (validate CSV/TSV)
- [ ] File handler in `app/utils/file_handler.py` (store uploads, register in DB)
- [ ] Upload API endpoint in `app/api/upload.py`

### Phase 2 — Pipeline Engine (Tool 1)
- [ ] FastQC wrapper, Cutadapt wrapper
- [ ] DADA2 R script + Python caller
- [ ] Taxonomy R script
- [ ] Phylogeny (MAFFT + FastTree)
- [ ] PICRUSt2 wrapper
- [ ] Pipeline orchestrator with background execution + status tracking
- [ ] Pipeline status page + QC report page

### Phase 3 — Data Management Hub (Tool 2)
- [ ] File manager page (browse/delete uploads)
- [ ] Dataset manager page (browse/inspect/download/delete)
- [ ] ZIP export for dataset download
- [ ] Dataset combiner (compatibility check, merge, re-run)

### Phase 4 — Analysis & Visualization (Tool 3)
- [ ] Alpha/beta diversity computation + plots
- [ ] Taxonomic composition visualizations
- [ ] Differential abundance (R scripts + volcano plots)
- [ ] PICRUSt2 pathway visualization
- [ ] Export (SVG/PDF/CSV)

### Phase 5 — Polish
- [ ] Error handling, theming, responsive layout
- [ ] Help/documentation tooltips
- [ ] Cross-study comparison helpers

---

## How to Start

```bash
cd ~/microbiome-dashboard
conda activate microbiome

# Start building Phase 1:
# 1. Create app/db/models.py with SQLAlchemy models from the schema
# 2. Create app/db/database.py with engine, session, init_db()
# 3. Create app/main.py wiring FastAPI + Dash
# 4. Create app/dashboard/layout.py with sidebar navigation
# 5. Create app/dashboard/pages/upload.py
# 6. Create app/pipeline/detect.py
# 7. Test with: uvicorn app.main:app --reload --host 0.0.0.0 --port 8050
```
