# 16S Analyzer -- 16S rRNA Microbiome Analysis Dashboard

A web-based tool for processing, managing, and visualizing 16S rRNA amplicon sequencing data. Built with Plotly Dash + FastAPI + SQLite.

**Supported input**: Illumina paired-end or single-end amplicon FASTQ files targeting specific 16S variable regions (V1-V2, V3-V4, V4, V4-V5, V5-V6). Full-length 16S long reads (PacBio HiFi, Nanopore) are also supported — auto-detected at upload, processed with DADA2 using platform-appropriate error models.

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

---

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Running the App](#running-the-app)
- [Project Structure](#project-structure)
- [Windows WSL Setup](#windows-wsl-setup)
- [Troubleshooting](#troubleshooting)
- [License](#license)

---

## Prerequisites

- **Operating System**: Ubuntu/Debian Linux (native or WSL2 on Windows)
- **Conda**: Miniconda or Miniforge installed ([install guide](https://docs.conda.io/en/latest/miniconda.html))
- **RAM**: 8 GB minimum, 16 GB recommended (16 GB required for PICRUSt2)
- **Disk Space**: ~10 GB for software + reference databases
- **System libraries** (install if missing):

```bash
sudo apt install -y build-essential curl git wget libxml2-dev \
    libcurl4-openssl-dev libssl-dev libfontconfig1-dev \
    libharfbuzz-dev libfribidi-dev libfreetype6-dev \
    libtiff5-dev libjpeg-dev libpng-dev
```

---

## Installation

### Step 1: Clone the repository

```bash
git clone https://github.com/tatsu1207/16S-Pipeline.git
cd 16S-Pipeline
```

### Step 2: Run the setup script

The setup script uses a **4-environment architecture** to avoid dependency conflicts:

| Environment | Contents |
|-------------|----------|
| `microbiome_16S` | Python 3.11 + CLI tools (FastQC, Cutadapt, MAFFT, FastTree, vsearch). The web app runs here. |
| `dada2_16S` | R 4.3 + DADA2 (pre-built from bioconda, zero compilation) |
| `analysis_16S` | R 4.3 + ALDEx2, DESeq2, ANCOM-BC2, MaAsLin2, LinDA, vegan |
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
│   ├── dashboard/               # Plotly Dash UI
│   │   ├── app.py               # Dash app initialization
│   │   ├── layout.py            # Sidebar nav + page routing
│   │   ├── components/
│   │   │   └── file_browser.py  # Server-side file/directory browser modal
│   │   └── pages/               # One file per page (16 pages)
│   │       ├── intro_page.py           # Landing page + quick-start guide
│   │       ├── file_manager.py         # Upload FASTQ, attach metadata
│   │       ├── pipeline_status.py      # Launch + monitor DADA2 pipeline
│   │       ├── datasets_page.py        # Inspect BIOM, extract sub-regions
│   │       ├── combine_page.py         # Merge BIOM files across studies
│   │       ├── biom_browser_page.py    # Read-only BIOM inspection
│   │       ├── subsampling_page.py     # Rarefy + filter samples
│   │       ├── rare_asv_page.py        # Remove low-prevalence ASVs
│   │       ├── mothur_page.py          # BIOM <-> MOTHUR conversion
│   │       ├── alpha_page.py           # Alpha diversity analysis
│   │       ├── beta_page.py            # Beta diversity + ordination
│   │       ├── taxonomy_page.py        # Taxonomy composition plots
│   │       ├── diff_abundance_page.py  # Differential abundance (5 tools)
│   │       ├── pathways_page.py        # PICRUSt2 pathway analysis
│   │       ├── picrust2_page.py        # Standalone PICRUSt2 runner
│   │       └── kegg_map_page.py        # KEGG pathway map viewer
│   ├── utils/                   # Utility modules
│   │   ├── file_handler.py      # Register local FASTQ files into DB
│   │   └── metadata_parser.py   # Metadata CSV/TSV parser + validator
│   └── db/                      # SQLAlchemy models + database
│       ├── database.py          # Session management
│       └── models.py            # 12 ORM tables
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
│   ├── references/              # SILVA databases + E. coli reference
│   ├── combined/                # Combined/merged datasets
│   └── exports/                 # User exports
├── setup_ubuntu.sh              # One-command installation script
├── setup_wsl2.sh                # WSL2-specific setup helper
├── run.sh                       # Start the application
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

## Windows WSL Setup

If you're on Windows, you need WSL2 with Ubuntu first. If you already have it, skip to [Installation](#installation).

### Step 1: Enable WSL

Open **PowerShell as Administrator** and run:

```powershell
wsl --install
```

This installs WSL2 with Ubuntu by default. If WSL is already installed but you need Ubuntu:

```powershell
wsl --install -d Ubuntu-24.04
```

Restart your computer when prompted.

### Step 2: Initial Ubuntu Setup

After restart, Ubuntu will open automatically (or search for "Ubuntu" in the Start menu). It will ask you to create a username and password.

### Step 3: Install system dependencies and Conda

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y build-essential curl git wget libxml2-dev \
    libcurl4-openssl-dev libssl-dev libfontconfig1-dev \
    libharfbuzz-dev libfribidi-dev libfreetype6-dev \
    libtiff5-dev libjpeg-dev libpng-dev

# Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
conda init bash
source ~/.bashrc
```

### Tips for WSL

- **Access WSL files from Windows**: Open File Explorer and go to `\\wsl$\Ubuntu\home\<your-username>`
- **Access Windows files from WSL**: Your C: drive is at `/mnt/c/`
- **VS Code integration**: Install the [Remote - WSL](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-wsl) extension

---

## Troubleshooting

### "conda: command not found"

```bash
eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
conda init bash
source ~/.bashrc
```

### R package installation fails (DADA2)

Make sure system libraries are installed (see [Prerequisites](#prerequisites)), then retry:

```bash
conda activate dada2_16S
Rscript -e 'BiocManager::install("dada2", force=TRUE)'
```

### "Permission denied" on setup_ubuntu.sh

```bash
chmod +x setup_ubuntu.sh
```

### Port already in use

```bash
# Find and kill the process using the port
lsof -ti:8050 | xargs kill -9

# Or use a different port
uvicorn app.main:app --reload --reload-exclude data --host 0.0.0.0 --port 8051
```

### PICRUSt2 installation fails

PICRUSt2 requires its own environment due to dependency conflicts:

```bash
conda create -n picrust2_16S -c bioconda -c conda-forge picrust2 -y
```

The setup script handles this automatically.

### WSL runs out of memory / PICRUSt2 OOM killed (exit code 137)

PICRUSt2's phylogenetic placement requires ~11 GB RAM. Create or edit `C:\Users\<YourUsername>\.wslconfig`:

```ini
[wsl2]
memory=16GB
swap=8GB
```

Then restart WSL from PowerShell: `wsl --shutdown`

> The pipeline (FastQC, Cutadapt, DADA2, taxonomy, phylogeny) works fine with 8 GB. Only PICRUSt2 requires 16 GB. If PICRUSt2 fails, the rest of the pipeline still completes -- you can re-run PICRUSt2 later from the Pipeline Status page.

---

## License

MIT License -- see [LICENSE](LICENSE) for details.
