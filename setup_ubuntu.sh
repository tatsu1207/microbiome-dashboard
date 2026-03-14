#!/usr/bin/env bash
# ============================================================================
# 16S Analyzer Setup Script
# One-command installation for the 16S Analyzer
#
# Architecture: 4 separate conda environments
#   microbiome — Python 3.11 + CLI bioinformatics tools (web app runs here)
#   dada2     — R 4.3 + bioconductor-dada2 (pre-built, zero compilation)
#   analysis  — R 4.3 + phyloseq/ANCOMBC/DESeq2/ALDEx2/Maaslin2/vegan + LinDA
#   picrust2  — PICRUSt2 (unchanged)
#
# Features:
#   - Skips any component that is already installed
#   - Auto-installs Miniforge3 (conda + mamba) if missing
#   - Uses bioconda pre-built R packages — no source compilation needed
#   - Re-runnable: safe to execute multiple times
# ============================================================================

set -eo pipefail  # Exit on any error; propagate failures through pipes

# --- Colors for output ---
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# --- Counters ---
INSTALLED_COUNT=0
SKIPPED_COUNT=0
FAILED_COUNT=0
FAILED_COMPONENTS=()

# --- Helper functions ---
info()    { echo -e "${BLUE}[INFO]${NC} $1"; }
success() { echo -e "${GREEN}[OK]${NC} $1"; }
warn()    { echo -e "${YELLOW}[WARN]${NC} $1"; }
skip()    { echo -e "${GREEN}[SKIP]${NC} $1 — already installed"; ((SKIPPED_COUNT++)) || true; }
error()   { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }
step()    { echo -e "\n${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"; \
            echo -e "${CYAN}  STEP $1: $2${NC}"; \
            echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"; }
installed() { ((INSTALLED_COUNT++)) || true; }
failed() { ((FAILED_COUNT++)) || true; FAILED_COMPONENTS+=("$1"); warn "Failed to install $1."; }

# Run a command while printing dots every 5 seconds to show it's still alive.
# Usage: run_with_dots "label" command arg1 arg2 ...
# The command's stdout/stderr are captured in RWD_OUTPUT. Exit code in RWD_EXIT.
RWD_OUTPUT=""
RWD_EXIT=0
run_with_dots() {
    local label="$1"; shift
    local logfile
    logfile=$(mktemp /tmp/rwd_XXXXXX.log)

    # Run command in background, capturing output
    "$@" > "${logfile}" 2>&1 &
    local cmd_pid=$!

    # Print dots while waiting
    echo -n "  ${label} "
    while kill -0 "${cmd_pid}" 2>/dev/null; do
        echo -n "."
        sleep 5
    done

    # Capture exit code without letting set -e abort the script
    RWD_EXIT=0
    wait "${cmd_pid}" || RWD_EXIT=$?
    RWD_OUTPUT=$(cat "${logfile}" 2>/dev/null || echo "")
    rm -f "${logfile}"
    echo ""  # newline after dots
}

# --- Cleanup trap for unexpected exits ---
cleanup_on_error() {
    local exit_code=$?
    if [[ ${exit_code} -ne 0 ]]; then
        echo ""
        echo -e "${RED}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
        echo -e "${RED}  Setup exited unexpectedly (exit code ${exit_code}).${NC}"
        echo -e "${RED}  The installation is incomplete. Re-run ./setup_ubuntu.sh${NC}"
        echo -e "${RED}  to resume — already-installed components will be skipped.${NC}"
        echo -e "${RED}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    fi
}
trap cleanup_on_error EXIT

# Remove stale (non-conda) env directories that block mamba create
clean_stale_env_dir() {
    local env_name="$1"
    local env_dir
    env_dir="$(conda info --base)/envs/${env_name}"
    if [[ -d "${env_dir}" ]] && [[ ! -f "${env_dir}/conda-meta/history" ]]; then
        warn "Removing stale non-conda directory at ${env_dir}"
        rm -rf "${env_dir}"
    fi
}

# Check if a command exists anywhere in PATH
has_cmd() { command -v "$1" &> /dev/null; }

# Check if a command exists specifically in a conda environment
has_env_cmd() {
    local env_prefix="$1" cmd="$2"
    [[ -n "${env_prefix}" ]] && [[ -x "${env_prefix}/bin/${cmd}" ]]
}

# Check if a Python package is importable in a given env
has_python_pkg() {
    conda run -n "$1" python -c "import $2" 2>/dev/null
}

# Check if an R package is installed in a given env
has_r_pkg() {
    conda run -n "$1" Rscript -e "if (!requireNamespace('$2', quietly=TRUE)) quit(status=1)" 2>/dev/null
}

# --- Configuration ---
ENV_NAME="microbiome_16S"
DADA2_ENV="dada2_16S"
ANALYSIS_ENV="analysis_16S"
PICRUST2_ENV="picrust2_16S"
PYTHON_VERSION="3.11"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${PROJECT_DIR}/data"
SILVA_DIR="${DATA_DIR}/references"

echo -e "${GREEN}"
echo "  ╔══════════════════════════════════════════════════════╗"
echo "  ║       🧬 16S Analyzer — Setup Script                  ║"
echo "  ║       16S rRNA Microbiome Analysis                   ║"
echo "  ╚══════════════════════════════════════════════════════╝"
echo -e "${NC}"

# ============================================================================
# STEP 0a: Install system dependencies
# ============================================================================
step "0a" "Checking system libraries"

SYS_DEPS=(libbz2-dev liblzma-dev libcurl4-openssl-dev zlib1g-dev libssl-dev libxml2-dev libpng-dev)
SYS_MISSING=()

for dep in "${SYS_DEPS[@]}"; do
    if dpkg -s "${dep}" &>/dev/null; then
        skip "System lib: ${dep}"
    else
        SYS_MISSING+=("${dep}")
    fi
done

if [ ${#SYS_MISSING[@]} -gt 0 ]; then
    info "Installing missing system libraries: ${SYS_MISSING[*]}"
    info "This requires sudo — you may be prompted for your password."
    if sudo apt-get install -y "${SYS_MISSING[@]}" 2>&1 | tail -5; then
        success "System libraries installed."
        installed
    else
        warn "Could not install system libraries."
        warn "Try manually: sudo apt-get install -y ${SYS_MISSING[*]}"
    fi
else
    info "All system libraries already present."
fi

# ============================================================================
# STEP 0b: Verify Conda is available
# ============================================================================
step "0b" "Checking Conda installation"

MINIFORGE_DIR="${HOME}/miniforge3"
MINIFORGE_INSTALLER="Miniforge3-Linux-$(uname -m).sh"

if ! has_cmd conda; then
    # Check if Miniforge3 exists but isn't in PATH
    if [[ -x "${MINIFORGE_DIR}/bin/conda" ]]; then
        info "Miniforge3 found at ${MINIFORGE_DIR} but not in PATH. Activating..."
        eval "$("${MINIFORGE_DIR}/bin/conda" shell.bash hook)"
        conda init bash 2>/dev/null
        success "Miniforge3 activated from ${MINIFORGE_DIR}"
    else
        info "Conda not found. Installing Miniforge3 (includes conda + mamba)..."

        if ! has_cmd wget && ! has_cmd curl; then
            error "Neither wget nor curl found. Install one first: sudo apt install wget"
        fi

        MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/${MINIFORGE_INSTALLER}"

        if [[ -f "/tmp/${MINIFORGE_INSTALLER}" ]]; then
            info "Using cached installer at /tmp/${MINIFORGE_INSTALLER}"
        elif has_cmd wget; then
            info "Downloading Miniforge3..."
            wget -q --show-progress -O "/tmp/${MINIFORGE_INSTALLER}" "${MINIFORGE_URL}"
        else
            info "Downloading Miniforge3..."
            curl -fsSL -o "/tmp/${MINIFORGE_INSTALLER}" "${MINIFORGE_URL}"
        fi

        bash "/tmp/${MINIFORGE_INSTALLER}" -b -p "${MINIFORGE_DIR}"
        eval "$("${MINIFORGE_DIR}/bin/conda" shell.bash hook)"
        conda init bash 2>/dev/null
        success "Miniforge3 installed to ${MINIFORGE_DIR}"
        installed
    fi
else
    success "Conda already installed."
fi

CONDA_VERSION=$(conda --version 2>&1)
info "Conda version: ${CONDA_VERSION}"

# Make sure conda commands work in this script
eval "$(conda shell.bash hook)"

# ============================================================================
# STEP 1: Verify Mamba (included with Miniforge3)
# ============================================================================
step "1" "Checking for Mamba"

if has_cmd mamba; then
    MAMBA_VERSION=$(mamba --version 2>&1 | head -1)
    success "Mamba found: ${MAMBA_VERSION}"
    SOLVER="mamba"
else
    warn "Mamba not found. Using conda (slower). Consider installing Miniforge3 for mamba support."
    SOLVER="conda"
fi

# ============================================================================
# STEP 2: Create 'microbiome' environment (Python + CLI tools)
# ============================================================================
step "2" "Checking Conda environment '${ENV_NAME}' (Python + CLI tools)"

if conda env list | grep -q "^${ENV_NAME} "; then
    skip "Conda environment '${ENV_NAME}'"
    info "To recreate from scratch, run: conda env remove -n ${ENV_NAME} && ./setup_ubuntu.sh"
else
    info "Creating '${ENV_NAME}' environment with Python ${PYTHON_VERSION}..."
    if ! ${SOLVER} create -n "${ENV_NAME}" -c conda-forge python=${PYTHON_VERSION} -y; then
        error "Failed to create conda environment '${ENV_NAME}'."
    fi
    success "Environment '${ENV_NAME}' created."
    installed
fi

# Activate for tool installs
if ! conda activate "${ENV_NAME}" 2>/dev/null; then
    error "Failed to activate conda environment '${ENV_NAME}'."
fi
success "Activated environment: ${ENV_NAME}"

ENV_PREFIX=$(conda info --envs | grep "^${ENV_NAME} " | awk '{print $NF}')

# Install bioinformatics CLI tools individually (avoids cross-channel conflicts)
info "Checking bioinformatics CLI tools..."

# Each entry: "conda_package:cmd1,cmd2"
BIOTOOLS_MAP=(
    "fastqc:fastqc"
    "cutadapt:cutadapt"
    "mafft:mafft"
    "fasttree:FastTree,fasttree"
    "bbmap:bbduk.sh"
    "sra-tools:prefetch,fasterq-dump"
)

for tool_pair in "${BIOTOOLS_MAP[@]}"; do
    pkg="${tool_pair%%:*}"
    cmds="${tool_pair##*:}"

    # Check if any of the command variants exist in the conda environment
    found=false
    IFS=',' read -ra CMD_ARRAY <<< "${cmds}"
    for cmd in "${CMD_ARRAY[@]}"; do
        if has_env_cmd "${ENV_PREFIX}" "${cmd}"; then
            found=true
            break
        fi
    done

    if ${found}; then
        skip "CLI tool: ${pkg}"
    else
        info "Installing ${pkg} ..."
        run_with_dots "conda: ${pkg}" \
            ${SOLVER} install -n "${ENV_NAME}" --override-channels -c conda-forge -c bioconda "${pkg}" -y
        echo "${RWD_OUTPUT:-}" | tail -5 || true

        # Re-check
        found=false
        for cmd in "${CMD_ARRAY[@]}"; do
            if has_env_cmd "${ENV_PREFIX}" "${cmd}"; then
                found=true
                break
            fi
        done

        if ${found}; then
            success "Installed: ${pkg}"
            installed
        else
            # Fallback: try pip (works for pure-Python tools like cutadapt)
            warn "Conda failed for ${pkg}, trying pip fallback..."
            pip install "${pkg}" 2>&1 | tail -5 || true
            found=false
            for cmd in "${CMD_ARRAY[@]}"; do
                if has_env_cmd "${ENV_PREFIX}" "${cmd}"; then
                    found=true
                    break
                fi
            done
            if ${found}; then
                success "Installed (via pip): ${pkg}"
                installed
            else
                failed "CLI tool: ${pkg}"
            fi
        fi
    fi
done

# Install Python packages via pip
info "Checking Python packages..."

PYTHON_PACKAGES=(
    "fastapi:fastapi"
    "uvicorn:uvicorn[standard]"
    "dash:dash"
    "dash_bootstrap_components:dash-bootstrap-components"
    "plotly:plotly"
    "sqlalchemy:sqlalchemy"
    "pandas:pandas"
    "numpy:numpy"
    "scipy:scipy"
    "skbio:scikit-bio"
    "multipart:python-multipart"
    "biom:biom-format"
    "matplotlib:matplotlib"
    "fpdf:fpdf2"
    "statsmodels:statsmodels"
    "dash_uploader:dash-uploader"
    "matplotlib_venn:matplotlib-venn"
)

PYTHON_MISSING=()

for pair in "${PYTHON_PACKAGES[@]}"; do
    module="${pair%%:*}"
    pip_name="${pair##*:}"

    if has_python_pkg "${ENV_NAME}" "${module}"; then
        skip "Python: ${pip_name}"
    else
        PYTHON_MISSING+=("${pip_name}")
    fi
done

if [ ${#PYTHON_MISSING[@]} -gt 0 ]; then
    info "Installing missing Python packages: ${PYTHON_MISSING[*]}"

    for pip_name in "${PYTHON_MISSING[@]}"; do
        for pair in "${PYTHON_PACKAGES[@]}"; do
            map_module="${pair%%:*}"
            map_pip="${pair##*:}"
            if [[ "${map_pip}" == "${pip_name}" ]]; then
                info "Installing ${pip_name} ..."
                pip install "${pip_name}" 2>&1 | tail -5 || true

                if has_python_pkg "${ENV_NAME}" "${map_module}"; then
                    success "Verified: ${pip_name}"
                    installed
                else
                    failed "Python: ${pip_name}"
                fi
                break
            fi
        done
    done
else
    info "All Python packages already present."
fi

# ============================================================================
# STEP 3: Create 'dada2' environment (R + DADA2, pre-built)
# ============================================================================
step "3" "Checking Conda environment '${DADA2_ENV}' (R + DADA2)"

if conda env list | grep -q "^${DADA2_ENV} "; then
    skip "Conda environment '${DADA2_ENV}'"
    info "To recreate: conda env remove -n ${DADA2_ENV} && ./setup_ubuntu.sh"
else
    clean_stale_env_dir "${DADA2_ENV}"
    info "Creating '${DADA2_ENV}' environment with pre-built bioconda packages..."
    info "(This uses pre-built binaries — no source compilation needed.)"
    run_with_dots "conda: ${DADA2_ENV}" \
        ${SOLVER} create -n "${DADA2_ENV}" --override-channels -c conda-forge -c bioconda \
            bioconductor-dada2 \
            r-optparse \
            r-jsonlite \
            -y
    echo "${RWD_OUTPUT:-}" | tail -5 || true
    if [[ ${RWD_EXIT} -eq 0 ]]; then
        success "Environment '${DADA2_ENV}' created (pre-built, zero compilation)."
        installed
    else
        warn "Failed to create '${DADA2_ENV}' environment:"
        echo "${RWD_OUTPUT:-}" | tail -20 || true
        failed "Environment: ${DADA2_ENV}"
    fi
fi

# ============================================================================
# STEP 4: Create 'analysis' environment (R + DA tools, pre-built)
# ============================================================================
step "4" "Checking Conda environment '${ANALYSIS_ENV}' (R + phyloseq/ANCOMBC/DESeq2/...)"

if conda env list | grep -q "^${ANALYSIS_ENV} "; then
    skip "Conda environment '${ANALYSIS_ENV}'"
    info "To recreate: conda env remove -n ${ANALYSIS_ENV} && ./setup_ubuntu.sh"
else
    clean_stale_env_dir "${ANALYSIS_ENV}"
    info "Creating '${ANALYSIS_ENV}' environment with pre-built bioconda packages..."
    info "(This uses pre-built binaries — no source compilation needed.)"
    run_with_dots "conda: ${ANALYSIS_ENV}" \
        ${SOLVER} create -n "${ANALYSIS_ENV}" --override-channels -c conda-forge -c bioconda \
            bioconductor-phyloseq \
            bioconductor-ancombc \
            bioconductor-deseq2 \
            bioconductor-aldex2 \
            bioconductor-maaslin2 \
            r-optparse \
            r-jsonlite \
            -y
    echo "${RWD_OUTPUT:-}" | tail -5 || true
    if [[ ${RWD_EXIT} -eq 0 ]]; then
        success "Environment '${ANALYSIS_ENV}' created (pre-built, zero compilation)."
        installed
    else
        warn "Failed to create '${ANALYSIS_ENV}' environment:"
        echo "${RWD_OUTPUT:-}" | tail -20 || true
        failed "Environment: ${ANALYSIS_ENV}"
    fi
fi

# --- Install vegan from CRAN (conda r-vegan has version conflicts with bioconda R) ---
if has_r_pkg "${ANALYSIS_ENV}" "vegan"; then
    skip "R package: vegan"
else
    info "Installing vegan from CRAN into '${ANALYSIS_ENV}' env..."
    run_with_dots "R: vegan" \
        conda run -n "${ANALYSIS_ENV}" --no-capture-output \
        Rscript -e "install.packages('vegan', repos='https://cloud.r-project.org', INSTALL_opts='--no-lock', Ncpus=4)"
    if has_r_pkg "${ANALYSIS_ENV}" "vegan"; then
        success "R package installed: vegan"
        installed
    else
        echo "${RWD_OUTPUT:-}" | tail -20 || true
        failed "R package: vegan"
    fi
fi

# --- Install LinDA from GitHub (only R package not on bioconda) ---
if has_r_pkg "${ANALYSIS_ENV}" "LinDA"; then
    skip "R package: LinDA"
else
    info "Installing LinDA from GitHub into '${ANALYSIS_ENV}' env..."
    run_with_dots "GitHub: LinDA" \
        conda run -n "${ANALYSIS_ENV}" --no-capture-output \
        Rscript -e "
        if (!requireNamespace('remotes', quietly=TRUE))
            install.packages('remotes', repos='https://cloud.r-project.org', INSTALL_opts='--no-lock')
        if (!requireNamespace('modeest', quietly=TRUE))
            install.packages('modeest', repos='https://cloud.r-project.org', INSTALL_opts='--no-lock')
        remotes::install_github('zhouhj1994/LinDA', upgrade='never', INSTALL_opts='--no-lock')
    "
    echo "${RWD_OUTPUT:-}" | tail -10 || true

    if has_r_pkg "${ANALYSIS_ENV}" "LinDA"; then
        success "LinDA installed from GitHub."
        installed
    else
        failed "R package: LinDA"
    fi
fi

# ============================================================================
# STEP 5: Install PICRUSt2
# ============================================================================
step "5" "Checking PICRUSt2 (separate environment)"

if conda env list | grep -q "^${PICRUST2_ENV} "; then
    skip "PICRUSt2 (${PICRUST2_ENV} environment)"
else
    info "Creating separate '${PICRUST2_ENV}' environment..."
    info "PICRUSt2 requires its own environment to avoid dependency conflicts."

    run_with_dots "conda: picrust2" \
        ${SOLVER} create -n "${PICRUST2_ENV}" --override-channels -c conda-forge -c bioconda picrust2 -y
    echo "${RWD_OUTPUT:-}" | tail -5 || true
    if [[ ${RWD_EXIT} -eq 0 ]]; then
        success "PICRUSt2 installed in separate '${PICRUST2_ENV}' environment."
        installed
    else
        warn "PICRUSt2 installation failed. You can install it manually later."
        warn "The pipeline will skip PICRUSt2 steps if it's not available."
    fi
fi

# ============================================================================
# STEP 6: Create project directory structure
# ============================================================================
step "6" "Checking project directory structure"

DIRS_CREATED=0

create_if_missing() {
    if [[ ! -d "$1" ]]; then
        mkdir -p "$1"
        ((DIRS_CREATED++)) || true
    fi
}

# App directories
create_if_missing "${PROJECT_DIR}/app/api"
create_if_missing "${PROJECT_DIR}/app/pipeline"
create_if_missing "${PROJECT_DIR}/app/data_manager"
create_if_missing "${PROJECT_DIR}/app/analysis"
create_if_missing "${PROJECT_DIR}/app/dashboard/pages"
create_if_missing "${PROJECT_DIR}/app/dashboard/components"
create_if_missing "${PROJECT_DIR}/app/db"
create_if_missing "${PROJECT_DIR}/app/utils"

# R scripts
create_if_missing "${PROJECT_DIR}/r_scripts"

# Data directories
create_if_missing "${DATA_DIR}/uploads"
create_if_missing "${DATA_DIR}/datasets"
create_if_missing "${DATA_DIR}/combined"
create_if_missing "${DATA_DIR}/exports"
create_if_missing "${DATA_DIR}/picrust2_runs"
create_if_missing "${DATA_DIR}/kegg_cache"
create_if_missing "${SILVA_DIR}"

# Create __init__.py files only where missing
INIT_CREATED=0
while IFS= read -r dir; do
    if [[ ! -f "${dir}/__init__.py" ]]; then
        touch "${dir}/__init__.py"
        ((INIT_CREATED++)) || true
    fi
done < <(find "${PROJECT_DIR}/app" -type d)

TOTAL_CREATED=$((DIRS_CREATED + INIT_CREATED))

if [[ ${TOTAL_CREATED} -gt 0 ]]; then
    success "Created ${DIRS_CREATED} directories and ${INIT_CREATED} __init__.py files."
    installed
else
    skip "Project directory structure"
fi

# ============================================================================
# STEP 7: Download SILVA 138.1 reference databases
# ============================================================================
step "7" "Checking SILVA 138.1 Reference Database"

SILVA_TRAIN="${SILVA_DIR}/silva_nr99_v138.1_train_set.fa.gz"
SILVA_SPECIES="${SILVA_DIR}/silva_species_assignment_v138.1.fa.gz"

# Minimum expected file sizes (bytes) to detect truncated downloads
SILVA_TRAIN_MIN_SIZE=20000000   # ~24 MB
SILVA_SPECIES_MIN_SIZE=70000000 # ~77 MB

# Helper: validate a downloaded file is not truncated
validate_download() {
    local file="$1" min_size="$2" label="$3"
    if [[ ! -f "${file}" ]]; then
        return 1
    fi
    local actual_size
    actual_size=$(stat --printf="%s" "${file}" 2>/dev/null || stat -f "%z" "${file}" 2>/dev/null || echo 0)
    if [[ ${actual_size} -lt ${min_size} ]]; then
        warn "${label} appears truncated (${actual_size} bytes, expected >${min_size})."
        warn "Removing corrupt file. Re-run setup to retry the download."
        rm -f "${file}"
        return 1
    fi
    return 0
}

# Helper: download a file with error handling (removes partial files on failure)
download_file() {
    local url="$1" dest="$2" label="$3"
    local tmp_dest="${dest}.part"
    local dl_ok=true
    if has_cmd wget; then
        wget -q --show-progress -O "${tmp_dest}" "${url}" || dl_ok=false
    else
        curl -fSL -o "${tmp_dest}" "${url}" || dl_ok=false
    fi

    if [[ "${dl_ok}" == "false" ]] || [[ ! -f "${tmp_dest}" ]]; then
        rm -f "${tmp_dest}"
        warn "${label} download failed. Re-run setup to retry."
        return 1
    fi

    mv "${tmp_dest}" "${dest}"
    return 0
}

SILVA_NEEDED=false

# Check existing files and validate they are not truncated
if [[ -f "${SILVA_TRAIN}" ]]; then
    if validate_download "${SILVA_TRAIN}" ${SILVA_TRAIN_MIN_SIZE} "SILVA training set"; then
        skip "SILVA training set"
    else
        SILVA_NEEDED=true
    fi
else
    SILVA_NEEDED=true
fi

if [[ -f "${SILVA_SPECIES}" ]]; then
    if validate_download "${SILVA_SPECIES}" ${SILVA_SPECIES_MIN_SIZE} "SILVA species assignment"; then
        skip "SILVA species assignment"
    else
        SILVA_NEEDED=true
    fi
else
    SILVA_NEEDED=true
fi

if [[ "${SILVA_NEEDED}" == "true" ]]; then
    echo ""
    info "SILVA 138.1 is required for taxonomic assignment (~100 MB total)."
    read -p "  Download missing files now? (Y/n): " DOWNLOAD_SILVA || DOWNLOAD_SILVA="Y"

    if [[ ! "${DOWNLOAD_SILVA}" =~ ^[Nn]$ ]]; then
        if ! has_cmd wget && ! has_cmd curl; then
            warn "Neither wget nor curl found. Skipping SILVA download."
            warn "Install wget (sudo apt install wget) and re-run setup."
        else
            SILVA_OK=true

            if [[ ! -f "${SILVA_TRAIN}" ]]; then
                info "Downloading SILVA training set (~24 MB)..."
                if download_file \
                    "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz" \
                    "${SILVA_TRAIN}" "SILVA training set"; then
                    if validate_download "${SILVA_TRAIN}" ${SILVA_TRAIN_MIN_SIZE} "SILVA training set"; then
                        installed
                    else
                        SILVA_OK=false
                        failed "SILVA training set"
                    fi
                else
                    SILVA_OK=false
                    failed "SILVA training set"
                fi
            fi

            if [[ ! -f "${SILVA_SPECIES}" ]]; then
                info "Downloading SILVA species assignment (~77 MB)..."
                if download_file \
                    "https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz" \
                    "${SILVA_SPECIES}" "SILVA species assignment"; then
                    if validate_download "${SILVA_SPECIES}" ${SILVA_SPECIES_MIN_SIZE} "SILVA species assignment"; then
                        installed
                    else
                        SILVA_OK=false
                        failed "SILVA species assignment"
                    fi
                else
                    SILVA_OK=false
                    failed "SILVA species assignment"
                fi
            fi

            if [[ "${SILVA_OK}" == "true" ]]; then
                success "SILVA 138.1 databases downloaded."
            fi
        fi
    else
        warn "Skipped SILVA download. You can download later — see README.md."
    fi
fi

# ============================================================================
# STEP 8: Create configuration file
# ============================================================================
step "8" "Checking configuration file"

CONFIG_FILE="${PROJECT_DIR}/app/config.py"

if [[ -f "${CONFIG_FILE}" ]]; then
    skip "Configuration file (app/config.py)"
    info "To regenerate, delete app/config.py and re-run setup."
else
    # Detect conda base path
    CONDA_BASE_PATH=$(conda info --base 2>/dev/null || echo "$HOME/miniforge3")

    cat > "${CONFIG_FILE}" << PYEOF
"""
16S Analyzer Configuration
Auto-generated by setup_ubuntu.sh — edit as needed.
"""
import os
from pathlib import Path

# --- Paths ---
PROJECT_DIR = Path("${PROJECT_DIR}")
DATA_DIR = PROJECT_DIR / "data"
UPLOAD_DIR = DATA_DIR / "uploads"
DATASET_DIR = DATA_DIR / "datasets"
COMBINED_DIR = DATA_DIR / "combined"
EXPORT_DIR = DATA_DIR / "exports"
PICRUST2_RUNS_DIR = DATA_DIR / "picrust2_runs"
REFERENCE_DIR = DATA_DIR / "references"
R_SCRIPTS_DIR = PROJECT_DIR / "r_scripts"

# --- Database ---
DATABASE_URL = f"sqlite:///{PROJECT_DIR / 'microbiome.db'}"

# --- SILVA 138.1 ---
SILVA_TRAIN_SET = REFERENCE_DIR / "silva_nr99_v138.1_train_set.fa.gz"
SILVA_SPECIES = REFERENCE_DIR / "silva_species_assignment_v138.1.fa.gz"

# --- Conda ---
CONDA_ENV_NAME = "${ENV_NAME}"
DADA2_ENV_NAME = "${DADA2_ENV}"
ANALYSIS_ENV_NAME = "${ANALYSIS_ENV}"
PICRUST2_ENV_NAME = "${PICRUST2_ENV}"
CONDA_BASE = Path(os.environ.get("CONDA_BASE", "${CONDA_BASE_PATH}"))

# Map R script filenames to the conda env that has their dependencies.
_R_SCRIPT_ENV_MAP = {
    "run_dada2.R": DADA2_ENV_NAME,
    "run_taxonomy.R": DADA2_ENV_NAME,
    "run_ancombc.R": ANALYSIS_ENV_NAME,
    "run_deseq2.R": ANALYSIS_ENV_NAME,
    "run_aldex2.R": ANALYSIS_ENV_NAME,
    "run_linda.R": ANALYSIS_ENV_NAME,
    "run_maaslin2.R": ANALYSIS_ENV_NAME,
    "run_nmds.R": ANALYSIS_ENV_NAME,
}


def r_script_env(script_name: str) -> str:
    """Return the conda env name for a given R script filename."""
    return _R_SCRIPT_ENV_MAP.get(script_name, CONDA_ENV_NAME)


def conda_cmd(args: list[str], env_name: str | None = None) -> list[str]:
    """Wrap a command to run inside a conda environment via conda run."""
    env = env_name or CONDA_ENV_NAME
    return ["conda", "run", "-n", env, "--no-capture-output"] + args

# --- Server ---
HOST = "0.0.0.0"
PORT = 7000 + os.getuid()  # e.g. UID 1000 -> port 8000, UID 1001 -> 8001
DEBUG = True

# --- Pipeline defaults ---
DADA2_DEFAULTS = {
    "trim_left_f": 0,
    "trim_left_r": 0,
    "trunc_len_f": 250,
    "trunc_len_r": 200,
    "min_overlap": 12,
    "threads": min(32, max(1, os.cpu_count() - 1)),  # Cap at 32, leave 1 core free
}
PYEOF

    success "Configuration file created at app/config.py"
    installed
fi

# ============================================================================
# STEP 9: Verify installation
# ============================================================================
step "9" "Verifying installation"

echo ""
info "Checking '${ENV_NAME}' environment (Python + CLI tools)..."

# Get env prefix for binary checks
MICROBIOME_PREFIX=$(conda info --envs | grep "^${ENV_NAME} " | awk '{print $NF}')

for tool_info in \
    "Python:python --version" \
    "FastQC:fastqc --version" \
    "Cutadapt:cutadapt --version" \
    "MAFFT:mafft --version 2>&1 | head -1"; do

    tool="${tool_info%%:*}"
    check_cmd="${tool_info#*:}"

    if eval "${check_cmd}" &> /dev/null; then
        VERSION=$(eval "${check_cmd}" 2>&1 | head -1)
        echo -e "  ${GREEN}✓${NC} ${tool}: ${VERSION}"
    else
        echo -e "  ${RED}✗${NC} ${tool}: not found"
    fi
done

# Python packages
info "Checking Python packages..."
python -c "
packages = {
    'fastapi': 'FastAPI',
    'uvicorn': 'Uvicorn',
    'dash': 'Dash',
    'dash_bootstrap_components': 'Dash Bootstrap',
    'dash_uploader': 'Dash Uploader',
    'plotly': 'Plotly',
    'sqlalchemy': 'SQLAlchemy',
    'pandas': 'Pandas',
    'numpy': 'NumPy',
    'scipy': 'SciPy',
    'skbio': 'scikit-bio',
    'biom': 'biom-format',
    'statsmodels': 'statsmodels',
    'matplotlib': 'Matplotlib',
    'fpdf': 'FPDF2',
    'multipart': 'python-multipart',
}
for module, name in packages.items():
    try:
        __import__(module)
        print(f'  \033[0;32m✓\033[0m {name}')
    except ImportError:
        print(f'  \033[0;31m✗\033[0m {name} — not installed')
"

# DADA2 env
echo ""
info "Checking '${DADA2_ENV}' environment (R + DADA2)..."
if conda env list | grep -q "^${DADA2_ENV} "; then
    for r_pkg in dada2 optparse jsonlite; do
        if has_r_pkg "${DADA2_ENV}" "${r_pkg}"; then
            echo -e "  ${GREEN}✓${NC} ${r_pkg}"
        else
            echo -e "  ${RED}✗${NC} ${r_pkg} — not installed"
        fi
    done
else
    echo -e "  ${RED}✗${NC} Environment '${DADA2_ENV}' not found"
fi

# Analysis env
echo ""
info "Checking '${ANALYSIS_ENV}' environment (R + DA tools)..."
if conda env list | grep -q "^${ANALYSIS_ENV} "; then
    for r_pkg in phyloseq ANCOMBC DESeq2 ALDEx2 Maaslin2 LinDA vegan optparse jsonlite; do
        if has_r_pkg "${ANALYSIS_ENV}" "${r_pkg}"; then
            echo -e "  ${GREEN}✓${NC} ${r_pkg}"
        else
            echo -e "  ${RED}✗${NC} ${r_pkg} — not installed"
        fi
    done
else
    echo -e "  ${RED}✗${NC} Environment '${ANALYSIS_ENV}' not found"
fi

# PICRUSt2
echo ""
if conda env list | grep -q "^${PICRUST2_ENV} "; then
    echo -e "  ${GREEN}✓${NC} PICRUSt2: available (${PICRUST2_ENV} env)"
else
    echo -e "  ${YELLOW}△${NC} PICRUSt2: not found (pathway analysis will be unavailable)"
fi

# SILVA
if [[ -f "${SILVA_TRAIN}" && -f "${SILVA_SPECIES}" ]]; then
    echo -e "  ${GREEN}✓${NC} SILVA 138.1: downloaded"
else
    echo -e "  ${YELLOW}△${NC} SILVA 138.1: not downloaded (see README)"
fi

# ============================================================================
# Summary
# ============================================================================
echo ""
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "  Summary: ${GREEN}${INSTALLED_COUNT} installed${NC} | ${BLUE}${SKIPPED_COUNT} skipped${NC} | ${RED}${FAILED_COUNT} failed${NC}"
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

if [[ ${FAILED_COUNT} -gt 0 ]]; then
    echo ""
    echo -e "${RED}  The following components failed to install:${NC}"
    for comp in "${FAILED_COMPONENTS[@]}"; do
        echo -e "    ${RED}✗${NC} ${comp}"
    done
    echo ""
    warn "Some features may not work until these are resolved."
    warn "You can re-run ./setup_ubuntu.sh to retry failed components."
    echo ""
    echo -e "${YELLOW}"
    echo "  ╔══════════════════════════════════════════════════════╗"
    echo "  ║       ⚠️  Setup completed with errors                ║"
    echo "  ╠══════════════════════════════════════════════════════╣"
    echo "  ║                                                      ║"
    echo "  ║  ${FAILED_COUNT} component(s) failed. See above for details.    ║"
    echo "  ║  Re-run ./setup_ubuntu.sh to retry.                  ║"
    echo "  ║                                                      ║"
    echo "  ║  To start the app anyway:                            ║"
    echo "  ║    conda activate ${ENV_NAME}                          ║"
    echo "  ║    ./run.sh                                          ║"
    echo "  ║                                                      ║"
    echo "  ╚══════════════════════════════════════════════════════╝"
    echo -e "${NC}"
elif [[ ${INSTALLED_COUNT} -eq 0 ]]; then
    echo ""
    echo -e "${GREEN}"
    echo "  ╔══════════════════════════════════════════════════════╗"
    echo "  ║       ✅ Everything already installed!                ║"
    echo "  ╠══════════════════════════════════════════════════════╣"
    echo "  ║                                                      ║"
    echo "  ║  To start the app:                                   ║"
    echo "  ║    conda activate ${ENV_NAME}                          ║"
    echo "  ║    ./run.sh                                          ║"
    echo "  ║                                                      ║"
    echo "  ║  Then open: http://localhost:\$(( 7000 + \$(id -u) ))      ║"
    echo "  ║                                                      ║"
    echo "  ╚══════════════════════════════════════════════════════╝"
    echo -e "${NC}"
else
    echo ""
    echo -e "${GREEN}"
    echo "  ╔══════════════════════════════════════════════════════╗"
    echo "  ║       🎉 Setup Complete!                             ║"
    echo "  ╠══════════════════════════════════════════════════════╣"
    echo "  ║                                                      ║"
    echo "  ║  To start the app:                                   ║"
    echo "  ║    conda activate ${ENV_NAME}                          ║"
    echo "  ║    ./run.sh                                          ║"
    echo "  ║                                                      ║"
    echo "  ║  Then open: http://localhost:\$(( 7000 + \$(id -u) ))      ║"
    echo "  ║                                                      ║"
    echo "  ╚══════════════════════════════════════════════════════╝"
    echo -e "${NC}"
fi
