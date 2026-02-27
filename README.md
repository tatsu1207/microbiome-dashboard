# 🧬 MicrobiomeDash — 16S rRNA Microbiome Analysis Dashboard

A web-based tool for processing, managing, and visualizing 16S rRNA amplicon sequencing data. Built with Plotly Dash + FastAPI + SQLite.

**Three integrated tools:**

| Tool | Purpose |
|------|---------|
| **Pipeline Engine** | FASTQ.gz → DADA2 → Taxonomy → PICRUSt2 → Microbiome Dataset |
| **Data Manager** | Browse, download, combine datasets across studies |
| **Analysis Dashboard** | Alpha/beta diversity, differential abundance, pathway analysis |

---

## Table of Contents

- [Prerequisites](#prerequisites)
- [Windows WSL Setup](#windows-wsl-setup)
- [Installation](#installation)
- [Download Reference Databases](#download-reference-databases)
- [Running the App](#running-the-app)
- [Project Structure](#project-structure)
- [Troubleshooting](#troubleshooting)
- [License](#license)

---

## Prerequisites

- **Operating System**: Windows 10 (version 2004+) or Windows 11 with WSL2
- **RAM**: 8 GB minimum, 16 GB recommended (16 GB required for PICRUSt2)
- **Disk Space**: ~10 GB for software + reference databases; more depending on your data
- **Internet**: Required for initial setup and reference database downloads

---

## Windows WSL Setup

If you already have WSL2 with Ubuntu installed, skip to [Installation](#installation).

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

After restart, Ubuntu will open automatically (or search for "Ubuntu" in the Start menu). It will ask you to create a username and password — these are for your Linux environment only.

### Step 3: Update Ubuntu

```bash
sudo apt update && sudo apt upgrade -y
```

### Step 4: Install system dependencies

```bash
sudo apt install -y build-essential curl git wget libxml2-dev \
    libcurl4-openssl-dev libssl-dev libfontconfig1-dev \
    libharfbuzz-dev libfribidi-dev libfreetype6-dev \
    libtiff5-dev libjpeg-dev libpng-dev
```

### Step 5: Install Miniconda

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
conda init bash
source ~/.bashrc
```

### Step 6: Verify

```bash
conda --version   # Should print conda version
python --version  # Should print Python version
R --version       # Will be installed by our setup script
```

You're ready for installation!

### Tips for WSL

- **Access WSL files from Windows**: Open File Explorer and go to `\\wsl$\Ubuntu\home\<your-username>`
- **Access Windows files from WSL**: Your C: drive is at `/mnt/c/`
- **Open WSL terminal**: Search "Ubuntu" in Start menu, or use Windows Terminal
- **VS Code integration**: Install the [Remote - WSL](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-wsl) extension for VS Code to edit files seamlessly

---

## Installation

### Step 1: Clone the repository

```bash
cd ~
git clone https://github.com/YOUR_USERNAME/microbiome-dashboard.git
cd microbiome-dashboard
```

### Step 2: Run the setup script

The setup script automatically:
- Detects and installs `mamba` (faster dependency resolver) if not present
- Creates the `microbiome` Conda environment with Python 3.11 + R 4.3
- Installs all Python, R, and bioinformatics dependencies
- Creates the project directory structure
- Initializes the SQLite database
- Downloads SILVA 138.1 reference databases (optional, prompted)

```bash
chmod +x setup.sh
./setup.sh
```

> ⏱ **Expected time**: 15–30 minutes depending on internet speed and system.
> The R package installation (especially DADA2) takes the longest.

### Step 3: Activate the environment

```bash
conda activate microbiome
```

> **Tip**: Add `conda activate microbiome` to your `~/.bashrc` if this is your primary project.

---

## Download Reference Databases

The setup script will offer to download these automatically. If you skipped that step:

```bash
mkdir -p data/references
cd data/references

# SILVA 138.1 training set for DADA2 (~24 MB)
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz

# SILVA 138.1 species assignment (~77 MB)
wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz

cd ../..
```

---

## Running the App

```bash
conda activate microbiome
./run.sh
```

Or manually:

```bash
conda activate microbiome
uvicorn app.main:app --reload --host 0.0.0.0 --port 8050
```

Open your **Windows browser** and go to: **http://localhost:8050**

---

## Project Structure

```
microbiome-dashboard/
├── app/                         # Application source code
│   ├── main.py                  # FastAPI + Dash entry point
│   ├── config.py                # Settings and paths
│   ├── api/                     # FastAPI REST endpoints
│   ├── pipeline/                # Tool 1: Pipeline Engine
│   ├── data_manager/            # Tool 2: Data Management Hub
│   ├── analysis/                # Tool 3: Analysis Engine
│   ├── dashboard/               # Plotly Dash UI (pages + components)
│   ├── db/                      # SQLAlchemy models + database
│   └── utils/                   # Shared utilities
├── r_scripts/                   # R scripts (DADA2, ANCOM-BC, ALDEx2)
├── data/                        # Data storage (uploads, datasets, references)
├── setup.sh                     # One-command installation script
├── run.sh                       # Start the application
├── environment.yml              # Conda environment specification
├── requirements.txt             # Python dependencies (pip)
├── Makefile                     # Common development commands
└── README.md                    # This file
```

---

## Development Commands

```bash
make run          # Start the app
make setup        # Run setup script
make db-reset     # Reset the database (WARNING: deletes all data)
make clean        # Remove temporary files and exports
make test         # Run tests (when available)
make help         # Show all available commands
```

---

## Troubleshooting

### "conda: command not found"

Your Conda installation is not in your PATH. Run:

```bash
eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
conda init bash
source ~/.bashrc
```

### R package installation fails (DADA2)

DADA2 requires system libraries. Make sure you installed them:

```bash
sudo apt install -y build-essential libxml2-dev libcurl4-openssl-dev \
    libssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev
```

Then retry:

```bash
Rscript -e 'BiocManager::install("dada2", force=TRUE)'
```

### "Permission denied" on setup.sh

```bash
chmod +x setup.sh
```

### Port 8050 already in use

```bash
# Find and kill the process using port 8050
lsof -ti:8050 | xargs kill -9

# Or use a different port
uvicorn app.main:app --reload --host 0.0.0.0 --port 8051
```

### PICRUSt2 installation fails

PICRUSt2 can be tricky with Conda. Try installing it separately:

```bash
conda create -n picrust2 -c bioconda -c conda-forge picrust2 -y
```

Then update `app/config.py` to point to the PICRUSt2 environment.

### WSL runs out of memory / PICRUSt2 OOM killed (exit code 137)

PICRUSt2's phylogenetic placement step loads a large reference tree that requires ~11 GB of RAM. With the default WSL2 memory limit, it will be killed by the OOM killer.

Create or edit `C:\Users\<YourUsername>\.wslconfig` in Windows:

```ini
[wsl2]
memory=16GB
swap=8GB
```

Then restart WSL from PowerShell:

```powershell
wsl --shutdown
```

Re-open your Ubuntu terminal — WSL will restart with the new limits. You can verify with `free -h` inside WSL.

> **Note**: The pipeline (FastQC, Cutadapt, DADA2, taxonomy, phylogeny) works fine with 8 GB. Only PICRUSt2 requires 16 GB. If PICRUSt2 fails, the rest of the pipeline still completes successfully — you can re-run PICRUSt2 later via the "Run PICRUSt2" button on the Pipeline Status page.

---

## License

MIT License — see [LICENSE](LICENSE) for details.
