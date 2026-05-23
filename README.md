# 16S Pipeline -- 16S rRNA Microbiome Analysis Platform

A web-based tool for processing, managing, and visualizing 16S rRNA amplicon sequencing data. Built with Plotly Dash + FastAPI + SQLite.

**Supported platforms**: Docker (Windows/macOS/Linux), native Linux (Ubuntu)

**Supported input**: Illumina paired-end or single-end amplicon FASTQ files targeting specific 16S variable regions (V1-V2, V3-V4, V4, V4-V5, V5-V6). Full-length 16S long reads (PacBio HiFi, Nanopore) are also supported -- auto-detected at upload, processed with DADA2 using platform-appropriate error models.

**Three integrated tools:**

| Tool | Purpose |
|------|---------|
| **Pipeline Engine** | FASTQ.gz -> DADA2 -> Taxonomy -> Phylogeny -> BIOM -> PICRUSt2 |
| **Data Manager** | Browse, download, combine, subsample datasets across studies |
| **Analysis Dashboard** | Alpha/beta diversity, differential abundance, pathway analysis, KEGG maps |

**Analysis capabilities:**

- **Alpha diversity** -- Shannon, Simpson, observed OTUs with Kruskal-Wallis / Mann-Whitney tests
- **Beta diversity** -- Bray-Curtis / Jaccard distance, PCoA, NMDS, PERMANOVA (pairwise + global)
- **Taxonomy** -- Stacked bar plots at any taxonomic level
- **Differential abundance** -- 5 tools: ALDEx2, DESeq2, ANCOM-BC2, LinDA, MaAsLin2; all-pairwise mode; volcano plots
- **Pathway analysis** -- PICRUSt2 output analysis with multi-tool DA, KO-to-KEGG aggregation, errorbar/heatmap/PCA plots (ggpicrust2-inspired)
- **KEGG Map** -- Targeted pathway inspection with DA-colored KEGG maps
- **SRA Download** -- Fetch public datasets from NCBI SRA by accession
- **Analysis Report** -- Generate comprehensive PDF reports with all analysis results

---

## Table of Contents

- [Quick Start (Docker)](#quick-start-docker)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Running the App](#running-the-app)
- [Project Structure](#project-structure)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

---

## Quick Start (Docker)

The fastest way to run 16S Pipeline on **any operating system** (Windows, macOS, Linux). No conda, R, or system libraries needed -- everything is packaged in a single Docker image. See the [Tutorial](TUTORIAL.md) for a step-by-step guide with example data.

### Requirements

- [Docker Desktop](https://www.docker.com/products/docker-desktop/) (Windows/macOS) or Docker Engine (Linux)
- **RAM**: 8 GB minimum, 16 GB recommended (PICRUSt2 needs 11 GB+). Docker Desktop defaults to only ~4 GB — increase it in Docker Desktop → **Settings** → **Resources** → **Memory**
- **Disk**: ~15 GB for the Docker image
- **Note**: PICRUSt2-related features are not available on Apple Silicon (ARM64)

### Windows

1. Download and install [Docker Desktop for Windows](https://www.docker.com/products/docker-desktop/). It will enable WSL2 automatically if needed. Restart your PC when prompted.

2. Open Docker Desktop and wait until it shows **"Docker Desktop is running"** (green icon in the system tray).

3. Download this repository: click the green **Code** button on [GitHub](https://github.com/tatsu1207/16S-Pipeline) → **Download ZIP**, then extract it.

4. Open **PowerShell**, navigate to the extracted folder, and run:

```powershell
cd path\to\16S-Pipeline
docker compose up -d
```

5. Open **http://localhost:8016** in your browser.

> The first run builds the Docker image (~15 GB), which may take 15-30 minutes. From the second time, you can either run `docker compose up -d` again or open Docker Desktop and start the container from the **Containers** tab -- it will start instantly.

### macOS

1. Download and install [Docker Desktop for Mac](https://www.docker.com/products/docker-desktop/) (supports both Apple Silicon and Intel).

2. Download and extract this repository (or `git clone`).

3. Open **Terminal** and run:

```bash
cd path/to/16S-Pipeline
docker compose up -d
```

4. Open **http://localhost:8016** in your browser.

> The first run builds the Docker image, which may take 15-30 minutes. From the second time, you can either run `docker compose up -d` again or start the container from Docker Desktop.

> **Note**: PICRUSt2-related analysis (functional prediction, pathway analysis, KEGG maps) is not available on Apple Silicon (ARM64) due to lack of native support.

### Linux

1. Install [Docker Engine](https://docs.docker.com/engine/install/) if not already installed.

2. Clone and run:

```bash
git clone https://github.com/tatsu1207/16S-Pipeline.git
cd 16S-Pipeline
docker compose up -d
```

3. Open **http://localhost:8016** in your browser.

> The first run builds the Docker image, which may take 15-30 minutes. Subsequent starts with `docker compose up -d` are instant.

### Managing the container

```bash
docker compose up -d         # Start
docker compose down          # Stop
docker compose logs -f       # View logs

# Use a different port
PORT=9000 docker compose up -d
```

### Where is my data stored?

All data (uploads, pipeline outputs, database) is stored in a Docker volume called `pipeline-data`. Your data persists across container restarts and updates.

```bash
# Back up your data
docker compose down
docker run --rm -v pipeline-data:/data -v $(pwd):/backup alpine tar czf /backup/16s-backup.tar.gz /data

# View volume info
docker volume inspect 16s-pipeline_pipeline-data
```

> If you prefer a native installation without Docker (e.g., for development or HPC environments), see the sections below.

---

## Prerequisites

- **Operating System**: Ubuntu/Debian Linux. For Windows/macOS, use [Docker](#quick-start-docker).
- **Conda**: Miniforge recommended ([install guide](https://github.com/conda-forge/miniforge))
- **RAM**: 8 GB minimum, 16 GB recommended (16 GB required for PICRUSt2)
- **Disk Space**: ~10 GB for software + reference databases

---

## Installation

### Step 1: Clone the repository

```bash
git clone https://github.com/tatsu1207/16S-Pipeline.git
cd 16S-Pipeline
```

### Step 2: Run the setup script for your platform

The setup script uses a **5-environment architecture** to avoid dependency conflicts:

| Environment | Contents |
|-------------|----------|
| `microbiome_16S` | Python 3.11 + CLI tools (FastQC, Cutadapt, MAFFT, FastTree, vsearch, sra-tools). The web app runs here. |
| `dada2_16S` | R 4.3 + DADA2 (pre-built from bioconda, zero compilation) |
| `analysis_16S` | R 4.4 + ALDEx2, DESeq2, ANCOM-BC2 |
| `maaslin2_16S` | R 4.3 + MaAsLin2, LinDA, vegan (separate env due to R version conflicts) |
| `picrust2_16S` | PICRUSt2 functional prediction |

The script also:
- Installs all Python packages (FastAPI, Dash, scikit-bio, biom-format, etc.)
- Downloads SILVA 138.1 reference databases (optional, prompted)
- Generates `app/config.py` with auto-detected paths
- Skips any component that is already installed (safe to re-run)

```bash
chmod +x setup_ubuntu.sh
./setup_ubuntu.sh
```

> For Windows/macOS users: use the [Docker installation](#quick-start-docker) instead.

> Expected time: 15-30 minutes depending on internet speed and system.

### Step 3: Activate the environment

```bash
conda activate microbiome_16S
```

---

## Running the App

```bash
conda activate microbiome_16S
./run.sh
```

The app runs in the background. The port is auto-assigned based on your UID (7000 + UID). Open the URL shown in the terminal output.

To run manually in the foreground:

```bash
conda activate microbiome_16S
uvicorn app.main:app --reload --reload-exclude data --host 0.0.0.0 --port 8050
```

---

## Project Structure

```
16S-Pipeline/
├── app/
│   ├── main.py                  # FastAPI + Dash entry point
│   ├── config.py                # Auto-generated settings and paths
│   ├── api/                     # FastAPI REST endpoints
│   │   ├── pipeline.py          # Pipeline control API
│   │   └── upload.py            # File upload API
│   ├── pipeline/                # Pipeline Engine
│   │   ├── runner.py            # Pipeline orchestrator
│   │   ├── detect.py            # Auto-detect sequencing type + variable region
│   │   ├── quality.py           # Quality profiling + auto trunc_len detection
│   │   ├── qc.py                # FastQC quality control
│   │   ├── qc_pdf.py            # QC report PDF generation
│   │   ├── trim.py              # Cutadapt adapter trimming
│   │   ├── dada2.py             # DADA2 denoising (R wrapper)
│   │   ├── taxonomy.py          # Taxonomic assignment (R wrapper)
│   │   ├── phylogeny.py         # Phylogenetic tree building
│   │   ├── biom_convert.py      # BIOM format conversion
│   │   └── picrust2.py          # PICRUSt2 functional prediction
│   ├── data_manager/            # Data Management
│   │   ├── biom_ops.py          # BIOM region detection, extraction, combining
│   │   ├── mothur_convert.py    # Bidirectional BIOM/MOTHUR conversion
│   │   ├── rare_asv.py          # Rare ASV filtering
│   │   └── subsample.py         # Rarefaction subsampling
│   ├── analysis/                # Analysis Engine
│   │   ├── shared.py            # Shared BIOM/metadata helpers
│   │   ├── r_runner.py          # R subprocess wrapper
│   │   ├── alpha.py             # Alpha diversity (skbio + scipy)
│   │   ├── beta.py              # Beta diversity, PCoA, NMDS, PERMANOVA
│   │   ├── taxonomy.py          # Taxonomy aggregation
│   │   ├── diff_abundance.py    # Multi-tool DA dispatcher + volcano plots
│   │   ├── pathways.py          # PICRUSt2 pathway DA
│   │   ├── kegg_aggregation.py  # KO-to-KEGG aggregation + annotation
│   │   ├── kegg_map.py          # KEGG pathway map helpers
│   │   └── pathway_plots.py     # Errorbar, heatmap, PCA visualizations
│   ├── report/                  # Report Generation
│   │   ├── report_generator.py  # PDF report builder (matplotlib + fpdf2)
│   │   └── methods_text.py      # Auto-generate Materials & Methods text
│   ├── sra/                     # SRA Integration
│   │   └── downloader.py        # NCBI SRA download via prefetch + fasterq-dump
│   ├── dashboard/               # Plotly Dash UI
│   │   ├── app.py               # Dash app initialization
│   │   ├── layout.py            # Sidebar nav + page routing
│   │   └── pages/               # One file per page
│   ├── utils/                   # Utility modules
│   │   ├── file_handler.py      # Register local FASTQ files into DB
│   │   └── metadata_parser.py   # Metadata CSV/TSV parser + validator
│   └── db/                      # SQLAlchemy models + database
│       ├── database.py          # Session management
│       └── models.py            # ORM tables
├── r_scripts/                   # R analysis scripts
│   ├── run_dada2.R              # DADA2 pipeline
│   ├── run_taxonomy.R           # Taxonomy assignment
│   ├── run_nmds.R               # NMDS ordination (vegan)
│   ├── run_aldex2.R             # ALDEx2 DA
│   ├── run_deseq2.R             # DESeq2 DA
│   ├── run_ancombc.R            # ANCOM-BC2 DA
│   ├── run_linda.R              # LinDA DA
│   └── run_maaslin2.R           # MaAsLin2 DA
├── data/                        # Data storage (gitignored except placeholders)
│   ├── uploads/                 # User FASTQ uploads
│   ├── datasets/                # Processed pipeline outputs
│   ├── picrust2_runs/           # PICRUSt2 output directories
│   ├── kegg_cache/              # Cached KEGG API data (24h TTL)
│   ├── sra_cache/               # Cached SRA downloads
│   ├── references/              # SILVA databases + E. coli reference
│   ├── combined/                # Combined/merged datasets
│   └── exports/                 # User exports
├── Dockerfile                   # Docker image build
├── docker-compose.yml           # One-command Docker deployment
├── docker-entrypoint.sh         # Docker container startup script
├── setup_ubuntu.sh              # Setup script for Linux (Ubuntu/Debian/WSL2)
├── run.sh                       # Start the application (native install)
├── environment.yml              # Conda environment specification
├── requirements.txt             # Python dependencies (pip)
└── Makefile                     # Development commands
```

---

## Download Reference Databases

The setup script will offer to download SILVA automatically. If you skipped that step:

```bash
cd data/references

# SILVA 138.1 training set for DADA2 (~24 MB)
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz

# SILVA 138.1 species assignment (~77 MB)
wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz

cd ../..
```

The E. coli 16S reference (`ecoli_16S.fasta`) is included in the repository and used for primer orientation detection.

---

## Troubleshooting

### Docker: container exits immediately

```bash
# Check the logs for error messages
docker compose logs

# Ensure Docker Desktop is running (Windows/macOS)
# Ensure you have enough RAM allocated to Docker
```

On Windows, Docker Desktop defaults to using half your system RAM. To increase it: Docker Desktop > Settings > Resources > Memory.

### Docker: port 8016 already in use

```bash
# Use a different port (e.g., 9016)
PORT=9016 docker compose up -d
# Then open http://localhost:9016
```

### Docker: how to reset everything

```bash
docker compose down
docker volume rm 16s-pipeline_pipeline-data
docker compose up -d
```

This deletes all uploaded data, pipeline outputs, and the database.

### PICRUSt2 fails (exit code 137 / OOM)

PICRUSt2 requires ~11 GB RAM. Increase Docker's memory allocation: Docker Desktop > **Settings** > **Resources** > **Memory** (set to 16 GB).

> The main pipeline (FastQC, Cutadapt, DADA2, taxonomy, phylogeny) works fine with 8 GB. Only PICRUSt2 requires 16 GB. If PICRUSt2 fails, the rest of the pipeline still completes -- you can re-run PICRUSt2 later from the Pipeline Status page.

---

## Citation

If you use 16S Pipeline in your research, please cite:

> Unno T. (2026). 16S-Pipeline: A comprehensive web-based platform for end-to-end 16S rRNA amplicon sequencing analysis. *Journal of Microbiology*. Epub ahead of print. https://doi.org/10.71150/jm.2603014

```bibtex
@article{unno2026sixteenspipeline,
  title={16S-Pipeline: A comprehensive web-based platform for end-to-end 16S rRNA amplicon sequencing analysis},
  author={Unno, Tatsuya},
  journal={Journal of Microbiology},
  year={2026},
  doi={10.71150/jm.2603014},
  note={Epub ahead of print}
}
```

---

## License

MIT License -- see [LICENSE](LICENSE) for details.
