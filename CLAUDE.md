# 16S Analyzer — Claude Reference Guide

Read this file before working on this project. It explains how the web tool works end-to-end.

## What This Tool Does

An end-to-end microbiome analysis platform: upload raw 16S rRNA FASTQ files → denoise with DADA2 → assign taxonomy → build phylogenetic tree → run functional prediction → perform statistical analysis. All through a browser UI.

**Supported inputs**: Illumina paired-end/single-end (V1-V2, V3-V4, V4, V4-V5, V5-V6 regions), full-length 16S long reads (PacBio HiFi, Nanopore).

## Architecture

- **Frontend**: Plotly Dash (Python interactive dashboard)
- **Backend API**: FastAPI (REST endpoints)
- **Database**: SQLite via SQLAlchemy ORM (`microbiome.db` at project root)
- **R Integration**: Subprocess calls to conda environments (`conda run -n <env> Rscript ...`)
- **Pipeline execution**: Background Python threads with subprocess spawning
- **Server**: Uvicorn, port = 7000 + user UID

### 4 Conda Environments

| Environment | Purpose |
|---|---|
| `microbiome_16S` | Web app, FastQC, Cutadapt, MAFFT, FastTree, vsearch |
| `dada2_16S` | R 4.3 + DADA2 for denoising and taxonomy |
| `analysis_16S` | R 4.3 + ALDEx2, DESeq2, ANCOM-BC2, LinDA, MaAsLin2, vegan |
| `picrust2_16S` | PICRUSt2 functional prediction (isolated due to dependency conflicts) |

## Project Structure

```
app/
├── main.py                    # FastAPI + Dash entry point, startup init
├── config.py                  # All paths, conda env names, defaults, port
├── api/
│   ├── upload.py              # POST /api/upload — file upload
│   └── pipeline.py            # POST/GET /api/pipeline/* — launch/cancel/status
├── pipeline/
│   ├── runner.py              # Main orchestrator (~1150 lines), _run_pipeline()
│   ├── detect.py              # SE/PE detection, variable region detection, platform detection
│   ├── quality.py             # Auto-truncation parameter detection (Q20 sliding window)
│   ├── dada2.py               # DADA2 subprocess wrapper
│   ├── taxonomy.py            # SILVA taxonomy wrapper
│   ├── phylogeny.py           # MAFFT alignment + FastTree
│   ├── biom_convert.py        # ASV table → BIOM (HDF5) conversion
│   ├── qc.py                  # FastQC wrapper
│   ├── qc_pdf.py              # QC report generation
│   ├── trim.py                # Cutadapt primer trimming
│   └── picrust2.py            # PICRUSt2 wrapper
├── dashboard/
│   ├── app.py                 # Dash app initialization
│   ├── layout.py              # Sidebar navigation + URL routing
│   └── pages/                 # 16 page modules (each registers its own callbacks)
│       ├── intro_page.py
│       ├── file_manager.py         # Upload FASTQ + attach metadata
│       ├── pipeline_status.py      # Launch + monitor pipeline (~1520 lines)
│       ├── datasets_page.py        # V-region extraction
│       ├── combine_page.py         # Merge BIOM files
│       ├── biom_browser_page.py    # Read-only BIOM inspection
│       ├── mothur_page.py          # BIOM ↔ MOTHUR format conversion
│       ├── subsampling_page.py     # Rarefaction + sample filtering
│       ├── rare_asv_page.py        # Low-prevalence ASV removal
│       ├── alpha_page.py           # Alpha diversity (Shannon, Simpson, Chao1, etc.)
│       ├── beta_page.py            # Beta diversity + ordination (PCoA, NMDS)
│       ├── taxonomy_page.py        # Stacked bar plots at all ranks
│       ├── diff_abundance_page.py  # 5-tool DA comparison with volcano plots
│       ├── pathways_page.py        # PICRUSt2 pathway DA comparison
│       ├── kegg_map_page.py        # Interactive KEGG pathway map viewer
│       └── picrust2_page.py        # Standalone PICRUSt2 runner
├── analysis/
│   ├── shared.py              # Helpers: BIOM↔DataFrame, metadata loading
│   ├── alpha.py               # Alpha diversity metrics
│   ├── beta.py                # Distance matrices, PCoA, NMDS, PERMANOVA
│   ├── taxonomy.py            # Taxonomy aggregation by rank
│   ├── diff_abundance.py      # Dispatcher for 5 DA tools (calls R scripts)
│   ├── pathways.py            # PICRUSt2 pathway comparison
│   ├── kegg_aggregation.py    # KO→EC→KEGG pathway aggregation
│   ├── kegg_map.py            # KEGG pathway map rendering
│   ├── pathway_plots.py       # Heatmap, error bar, PCA plots
│   └── r_runner.py            # R subprocess wrapper (streams output, parses JSON result)
├── data_manager/
│   ├── biom_ops.py            # Region detection/extraction, dataset combining
│   ├── mothur_convert.py      # BIOM ↔ MOTHUR shared/taxonomy conversion
│   ├── subsample.py           # Rarefaction to uniform depth
│   └── rare_asv.py            # Filter ASVs by prevalence/abundance
├── db/
│   ├── database.py            # Session management, init_db()
│   └── models.py              # 12 ORM tables (projects, uploads, fastq_files, datasets, samples, etc.)
└── utils/
    ├── file_handler.py        # Register local FASTQ files
    └── metadata_parser.py     # CSV/TSV metadata parsing

r_scripts/
├── run_dada2.R                # DADA2 pipeline (SE/PE/longread modes)
├── run_taxonomy.R             # SILVA taxonomy assignment
├── run_aldex2.R               # ALDEx2 differential abundance
├── run_ancombc.R              # ANCOM-BC2 differential abundance
├── run_deseq2.R               # DESeq2 differential abundance
├── run_linda.R                # LinDA differential abundance
├── run_maaslin2.R             # MaAsLin2 differential abundance
└── run_nmds.R                 # NMDS ordination

data/                          # All gitignored except .gitkeep files
├── uploads/                   # Raw FASTQ uploads
├── datasets/                  # Pipeline outputs (per dataset_id)
├── combined/                  # Merged datasets
├── picrust2_runs/             # PICRUSt2 outputs
├── references/                # SILVA 138.1 + E. coli ref
├── kegg_cache/                # Cached KEGG API responses (24h TTL)
└── exports/                   # User-downloaded files
```

## Pipeline Execution Flow

When a user launches a pipeline, these steps run sequentially in a background thread:

1. **FastQC** — Quality reports → `{dataset_dir}/qc/`
2. **Cutadapt** — Trim primers → `{dataset_dir}/trimmed/`
3. **Auto-truncation** — Analyze quality scores, find where 10bp sliding window mean < Q20. For PE: ensures `trunc_f + trunc_r >= insert_len + min_overlap`
4. **DADA2** — Denoise reads, remove chimeras → `asv_table.tsv`, `rep_seqs.fasta`, `track_reads.tsv`. Uses platform-specific error models for long reads.
5. **Taxonomy** — SILVA 138.1 assignment → `taxonomy.tsv`
6. **Phylogenetic tree** — MAFFT + FastTree → `tree.nwk`
7. **BIOM conversion** — Combine ASV table + taxonomy into HDF5 BIOM → `asv_table.biom`
8. **PICRUSt2** (optional) — Functional prediction in separate conda env

### Status Tracking

- Progress stored in `{dataset_dir}/status.json` with `current_step`, `progress_pct`, `steps_completed`, `pid`
- UI polls `GET /api/pipeline/status/{dataset_id}` via `dcc.Interval`
- Cancellation: sets `threading.Event` flag + kills subprocess via `os.killpg()`
- PID persists in `status.json` for recovery after server restart

## Data Flow Summary

```
Upload FASTQ → detect SE/PE + region + platform → attach metadata (CSV/TSV)
→ launch pipeline → DADA2 denoising → taxonomy → tree → BIOM
→ optional: filter ASVs, rarefy, combine datasets, extract V-regions
→ analysis: alpha/beta diversity, taxonomy plots, differential abundance, pathways
→ export: PNG/SVG/PDF plots, CSV tables
```

## Database (SQLite + SQLAlchemy)

Key tables in `app/db/models.py`:
- **projects** — Study grouping
- **uploads** — FASTQ upload batches (sequencing_type, variable_region, platform)
- **fastq_files** — Individual files (sample_name, read_direction, read_count)
- **upload_metadata** — Key-value per sample
- **datasets** — Pipeline outputs (status: pending/processing/complete/failed)
- **dataset_fastq_files** — M2M linking datasets to FASTQ files
- **samples** — Post-pipeline sample metadata (JSON)
- **analysis_results** — Cached analysis outputs
- **picrust2_runs** — Standalone PICRUSt2 jobs

Session pattern:
```python
from app.db.database import get_session
with get_session() as db:
    dataset = db.query(Dataset).filter(...).first()
```

## Key Code Patterns

### Dash Callbacks
```python
@dash_app.callback(
    Output("id-output", "children"),
    Input("id-button", "n_clicks"),
    [State("id-store", "data")],
)
def callback(n_clicks, store_data):
    if not n_clicks:
        return no_update
    return html.Div(...)
```

All page modules are eagerly imported in `main.py` so callbacks register before the server starts. Each page has a `get_layout()` function. URL routing is handled by a callback in `layout.py`.

### R Script Execution
```python
# In analysis/r_runner.py — spawns R in appropriate conda env
conda run -n {env_name} Rscript {script}.R --arg1 val1 --arg2 val2
```
Streams stdout/stderr to logger, parses JSON from last stdout line as result.

### DADA2 R Script (`r_scripts/run_dada2.R`)
- All output uses `log_msg()` (= `cat()` + `flush.console()`) to ensure output is flushed before potential crashes.
- Long-running steps (dereplication, denoising, merging, chimera removal) log elapsed time on completion.
- Auto-skips `filterAndTrim()` if `filtered/` directory already contains output files (also supported via explicit `--skip_filter` flag).
- The Python wrapper (`app/pipeline/dada2.py`) decodes signal names on crash (e.g., SIGSEGV, SIGKILL) and includes the last R output line in the error message.

### Pipeline Threading
- Global dicts: `_running_pipelines`, `_cancel_events`, `_active_procs`
- Each dataset_id gets its own cancel event
- Daemon threads auto-cleanup on shutdown

## Important Gotchas

1. **Truncation params (0, 0)** = auto-detect. The auto-detection in `quality.py` is critical for paired-end merge success.
2. **Long-read detection** triggers at avg read length > 900bp. Uses different DADA2 error models (PacBioErrfun vs loessErrfun).
3. **PICRUSt2 needs 11GB+ RAM**. Can fail independently without blocking the main pipeline.
4. **Primer detection threshold**: >30% of reads must contain a known primer for it to be detected.
5. **BIOM format is HDF5** (binary). Requires biom-format library. Rep seqs are embedded as observation metadata.
6. **KEGG API has rate limits** — responses cached 24h in `data/kegg_cache/`.
7. **Taxonomy can be None** — several places must handle missing taxonomy gracefully (past bug source).
8. **Multi-region combining** has two modes: by-sequence (same region) or by-taxonomy (cross-region, uses E. coli alignment positions).
9. **Callback import order matters** — all page modules must be imported before Dash starts.
10. **File uploads** use dash-uploader (chunked) rather than standard dcc.Upload (which blocks the browser for large files).
11. **Pipeline crash recovery**: If the R subprocess dies (segfault, OOM, signal), the daemon thread also dies silently. The UI detects this via PID liveness check and marks the dataset as "failed". On re-run, DADA2 automatically skips filtering if `filtered/` files already exist.
12. **Web app must be running** for the UI to work. Started via `bash run.sh` (uvicorn on port 7000+UID). If the app process dies, it must be manually restarted.

## Config Essentials (`app/config.py`)

- `DATA_DIR` = `{PROJECT_DIR}/data`
- `SILVA_TRAIN_SET` / `SILVA_SPECIES` = reference databases in `data/references/`
- `DADA2_DEFAULTS` = `{trim_left_f: 0, trim_left_r: 0, trunc_len_f: 0, trunc_len_r: 0, min_overlap: 12}`
- `LONGREAD_DADA2_DEFAULTS` = `{min_len: 1000, max_len: 1600, max_ee: 10, band_size: 32}`
- `conda_cmd(args, env_name)` = wrapper that builds `["conda", "run", "-n", env_name, "--no-capture-output"] + args`
- Threads default: `min(32, max(1, CPU_COUNT - 1))`

## Running the App

```bash
# One-command setup (installs all 4 conda envs + SILVA refs)
bash setup_ubuntu.sh

# Start the server
bash run.sh          # or: make run
# → Listens on http://0.0.0.0:{7000+UID}
```
