#!/usr/bin/env bash
# ============================================================================
# MicrobiomeDash Setup Script
# One-command installation for the 16S rRNA Microbiome Analysis Dashboard
#
# Features:
#   - Skips any component that is already installed
#   - Auto-installs Miniforge3 (conda + mamba) if missing
#   - Handles PICRUSt2 dependency conflicts gracefully
#   - Re-runnable: safe to execute multiple times
# ============================================================================

set -e  # Exit on any error

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

# Check if a command exists anywhere in PATH
has_cmd() { command -v "$1" &> /dev/null; }

# Check if a command exists specifically in the conda environment
has_env_cmd() {
    [[ -n "${ENV_PREFIX}" ]] && [[ -x "${ENV_PREFIX}/bin/$1" ]]
}

# Check if a Python package is importable in the target env
has_python_pkg() {
    conda run -n "${ENV_NAME}" python -c "import $1" 2>/dev/null
}

# Check if an R package is installed in the target env
has_r_pkg() {
    conda run -n "${ENV_NAME}" Rscript -e "if (!requireNamespace('$1', quietly=TRUE)) quit(status=1)" 2>/dev/null
}

# --- Configuration ---
ENV_NAME="microbiome"
PYTHON_VERSION="3.11"
R_VERSION="4.3"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${PROJECT_DIR}/data"
SILVA_DIR="${DATA_DIR}/references"

echo -e "${GREEN}"
echo "  ╔══════════════════════════════════════════════════════╗"
echo "  ║       🧬 MicrobiomeDash — Setup Script              ║"
echo "  ║       16S rRNA Microbiome Analysis Dashboard         ║"
echo "  ╚══════════════════════════════════════════════════════╝"
echo -e "${NC}"

# ============================================================================
# STEP 0: Verify Conda is available
# ============================================================================
step "0" "Checking Conda installation"

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
# STEP 2: Create Conda environment
# ============================================================================
step "2" "Checking Conda environment '${ENV_NAME}'"

if conda env list | grep -q "^${ENV_NAME} "; then
    skip "Conda environment '${ENV_NAME}'"
    info "To recreate from scratch, run: conda env remove -n ${ENV_NAME} && ./setup_ubuntu.sh"
else
    info "Creating environment with Python ${PYTHON_VERSION} and R ${R_VERSION}..."
    ${SOLVER} create -n "${ENV_NAME}" python=${PYTHON_VERSION} r-base=${R_VERSION} -y
    success "Environment '${ENV_NAME}' created."
    installed
fi

# Activate environment
conda activate "${ENV_NAME}"
success "Activated environment: ${ENV_NAME}"

# Pin r-base to prevent accidental upgrades that break Bioconductor packages
ENV_PREFIX=$(conda info --envs | grep "^${ENV_NAME} " | awk '{print $NF}')
PINNED_FILE="${ENV_PREFIX}/conda-meta/pinned"
if ! grep -q "r-base" "${PINNED_FILE}" 2>/dev/null; then
    mkdir -p "${ENV_PREFIX}/conda-meta"
    echo "r-base ${R_VERSION}.*" >> "${PINNED_FILE}"
    success "Pinned r-base to ${R_VERSION}.* (prevents upgrade breaking Bioconductor)"
fi

# ============================================================================
# STEP 3: Install bioinformatics CLI tools via Conda
# ============================================================================
step "3" "Checking bioinformatics CLI tools"

BIOTOOLS_TO_INSTALL=()

# Each entry: "conda_package:cmd1,cmd2" (cmd variants to check)
BIOTOOLS_MAP=(
    "fastqc:fastqc"
    "cutadapt:cutadapt"
    "mafft:mafft"
)

for tool_pair in "${BIOTOOLS_MAP[@]}"; do
    pkg="${tool_pair%%:*}"
    cmds="${tool_pair##*:}"

    # Check if any of the command variants exist in the conda environment
    found=false
    IFS=',' read -ra CMD_ARRAY <<< "${cmds}"
    for cmd in "${CMD_ARRAY[@]}"; do
        if has_env_cmd "${cmd}"; then
            found=true
            break
        fi
    done

    if ${found}; then
        skip "${pkg}"
    else
        BIOTOOLS_TO_INSTALL+=("${pkg}")
    fi
done

if [ ${#BIOTOOLS_TO_INSTALL[@]} -gt 0 ]; then
    info "Installing missing tools: ${BIOTOOLS_TO_INSTALL[*]}"
    ${SOLVER} install -n "${ENV_NAME}" -c bioconda -c conda-forge \
        "${BIOTOOLS_TO_INSTALL[@]}" -y 2>&1 | tail -5
    success "Bioinformatics tools installed: ${BIOTOOLS_TO_INSTALL[*]}"
    installed
else
    info "All bioinformatics CLI tools already present."
fi

# ============================================================================
# STEP 4: Install PICRUSt2
# ============================================================================
step "4" "Checking PICRUSt2 (separate environment)"

PICRUST2_SEPARATE=true

if conda env list | grep -q "^picrust2 "; then
    skip "PICRUSt2 (picrust2 environment)"
else
    info "Creating separate 'picrust2' environment..."
    info "PICRUSt2 requires its own environment to avoid dependency conflicts."

    if ${SOLVER} create -n picrust2 -c bioconda -c conda-forge picrust2 -y 2>&1 | tail -5; then
        success "PICRUSt2 installed in separate 'picrust2' environment."
        installed
    else
        warn "PICRUSt2 installation failed. You can install it manually later."
        warn "The pipeline will skip PICRUSt2 steps if it's not available."
    fi
fi

# ============================================================================
# STEP 5: Install R packages
# ============================================================================
step "5" "Checking R packages"

# Make sure we're in the right env
conda activate "${ENV_NAME}"

# --- 5a-pre: Ensure BiocManager is installed and version-matched to R ---
# Conda may ship a BiocManager pinned to an older Bioconductor release that
# does not match the current R version, causing all BiocManager::install()
# calls to fail.  Detect the correct version and upgrade if needed.
info "Checking BiocManager version matches R..."
Rscript -e "
    if (!requireNamespace('BiocManager', quietly=TRUE))
        install.packages('BiocManager', repos='https://cloud.r-project.org')

    tryCatch({
        BiocManager::install(ask=FALSE, update=FALSE)
        message(sprintf('BiocManager %s OK for %s', BiocManager::version(), R.version.string))
    }, error = function(e) {
        message('BiocManager version mismatch detected, upgrading...')
        m <- regmatches(e\$message, regexpr(\"version = '[0-9.]+'\", e\$message))
        if (length(m) > 0) {
            ver <- gsub(\"version = '|'\", '', m)
            message(sprintf('Upgrading to Bioconductor %s', ver))
            BiocManager::install(version=ver, ask=FALSE)
        } else {
            message('Could not auto-detect version, trying fresh BiocManager...')
            install.packages('BiocManager', repos='https://cloud.r-project.org')
            BiocManager::install(ask=FALSE, update=FALSE)
        }
        message(sprintf('BiocManager now at %s', BiocManager::version()))
    })
" 2>&1 | tail -5

# --- 5a: Conda-installable R packages ---
# Install one-by-one so a conflict in one package doesn't block the rest.
# Each entry: "conda_package_name:r_namespace"
R_CONDA_MAP=(
    "r-optparse:optparse"
    "r-jsonlite:jsonlite"
    "bioconductor-dada2:dada2"
    "bioconductor-phyloseq:phyloseq"
)

for pair in "${R_CONDA_MAP[@]}"; do
    conda_pkg="${pair%%:*}"
    r_pkg="${pair##*:}"

    if has_r_pkg "${r_pkg}"; then
        skip "R package: ${r_pkg}"
    else
        info "Installing ${conda_pkg} ..."
        if ${SOLVER} install -n "${ENV_NAME}" -c bioconda -c conda-forge \
                "${conda_pkg}" -y 2>&1 | tail -5; then
            # Verify it actually installed
            if has_r_pkg "${r_pkg}"; then
                success "R package installed: ${r_pkg}"
                installed
            else
                warn "Conda reported success but ${r_pkg} is not loadable in R."
            fi
        else
            warn "Conda/Mamba could not install ${conda_pkg} (dependency conflict)."
            warn "Trying install via R install.packages / BiocManager..."
            # Fallback: install from R directly
            if [[ "${conda_pkg}" == bioconductor-* ]]; then
                Rscript -e "
                    if (!requireNamespace('BiocManager', quietly=TRUE))
                        install.packages('BiocManager', repos='https://cloud.r-project.org')
                    BiocManager::install('${r_pkg}', ask=FALSE, update=FALSE)
                " 2>&1 | tail -5
            else
                Rscript -e "install.packages('${r_pkg}', repos='https://cloud.r-project.org')" 2>&1 | tail -5
            fi

            if has_r_pkg "${r_pkg}"; then
                success "R package installed via fallback: ${r_pkg}"
                installed
            else
                warn "Failed to install ${r_pkg}. You may need to install it manually."
            fi
        fi
    fi
done

# --- 5b: Bioconductor-only R packages (ANCOM-BC, ALDEx2) ---
R_BIOC_NEEDED=()

if has_r_pkg "ANCOMBC"; then
    skip "R package: ANCOMBC"
else
    R_BIOC_NEEDED+=("ANCOMBC")
fi

if has_r_pkg "ALDEx2"; then
    skip "R package: ALDEx2"
else
    R_BIOC_NEEDED+=("ALDEx2")
fi

if has_r_pkg "DESeq2"; then
    skip "R package: DESeq2"
else
    R_BIOC_NEEDED+=("DESeq2")
fi

if has_r_pkg "Maaslin2"; then
    skip "R package: Maaslin2"
else
    R_BIOC_NEEDED+=("Maaslin2")
fi

if [ ${#R_BIOC_NEEDED[@]} -gt 0 ]; then
    info "Installing from Bioconductor: ${R_BIOC_NEEDED[*]}"
    info "This may take 10-15 minutes on first install..."

    # Install one-by-one to isolate failures
    for bioc_pkg in "${R_BIOC_NEEDED[@]}"; do
        info "Installing ${bioc_pkg} via BiocManager..."
        Rscript -e "
            if (!requireNamespace('BiocManager', quietly = TRUE))
                install.packages('BiocManager', repos='https://cloud.r-project.org')
            BiocManager::install('${bioc_pkg}', ask=FALSE, update=FALSE)
        " 2>&1 | tail -10

        if has_r_pkg "${bioc_pkg}"; then
            success "Bioconductor package installed: ${bioc_pkg}"
            installed
        else
            warn "Failed to install ${bioc_pkg}. You may need to install it manually."
        fi
    done
else
    info "All Bioconductor R packages already present."
fi

# --- 5c: LinDA (GitHub-only R package) ---
if has_r_pkg "LinDA"; then
    skip "R package: LinDA"
else
    info "Installing LinDA from GitHub..."
    Rscript -e "
        if (!requireNamespace('remotes', quietly=TRUE))
            install.packages('remotes', repos='https://cloud.r-project.org')
        remotes::install_github('zhouhj1994/LinDA', upgrade='never')
    " 2>&1 | tail -10

    if has_r_pkg "LinDA"; then
        success "LinDA installed from GitHub."
        installed
    else
        warn "Failed to install LinDA. You may need to install it manually."
    fi
fi

# --- 5d: vegan (for NMDS ordination) ---
if has_r_pkg "vegan"; then
    skip "R package: vegan"
else
    info "Installing vegan..."
    Rscript -e "install.packages('vegan', repos='https://cloud.r-project.org')" 2>&1 | tail -5

    if has_r_pkg "vegan"; then
        success "R package installed: vegan"
        installed
    else
        warn "Failed to install vegan."
    fi
fi

# ============================================================================
# STEP 6: Install Python packages
# ============================================================================
step "6" "Checking Python packages"

conda activate "${ENV_NAME}"

# Each entry: "import_module:pip_package_name"
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
)

PYTHON_MISSING=()

for pair in "${PYTHON_PACKAGES[@]}"; do
    module="${pair%%:*}"
    pip_name="${pair##*:}"
    
    if has_python_pkg "${module}"; then
        skip "Python: ${pip_name}"
    else
        PYTHON_MISSING+=("${pip_name}")
    fi
done

if [ ${#PYTHON_MISSING[@]} -gt 0 ]; then
    info "Installing missing Python packages: ${PYTHON_MISSING[*]}"
    pip install "${PYTHON_MISSING[@]}" 2>&1 | tail -10
    success "Python packages installed: ${PYTHON_MISSING[*]}"
    installed
else
    info "All Python packages already present."
fi

# ============================================================================
# STEP 7: Create project directory structure
# ============================================================================
step "7" "Checking project directory structure"

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
# STEP 8: Download SILVA 138.1 reference databases
# ============================================================================
step "8" "Checking SILVA 138.1 Reference Database"

SILVA_TRAIN="${SILVA_DIR}/silva_nr99_v138.1_train_set.fa.gz"
SILVA_SPECIES="${SILVA_DIR}/silva_species_assignment_v138.1.fa.gz"

SILVA_NEEDED=false

if [[ -f "${SILVA_TRAIN}" ]]; then
    skip "SILVA training set"
else
    SILVA_NEEDED=true
fi

if [[ -f "${SILVA_SPECIES}" ]]; then
    skip "SILVA species assignment"
else
    SILVA_NEEDED=true
fi

if [[ "${SILVA_NEEDED}" == "true" ]]; then
    echo ""
    info "SILVA 138.1 is required for taxonomic assignment (~100 MB total)."
    read -p "  Download missing files now? (Y/n): " DOWNLOAD_SILVA
    
    if [[ ! "${DOWNLOAD_SILVA}" =~ ^[Nn]$ ]]; then
        if [[ ! -f "${SILVA_TRAIN}" ]]; then
            info "Downloading SILVA training set (~24 MB)..."
            wget -q --show-progress -O "${SILVA_TRAIN}" \
                "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz"
            installed
        fi
        
        if [[ ! -f "${SILVA_SPECIES}" ]]; then
            info "Downloading SILVA species assignment (~77 MB)..."
            wget -q --show-progress -O "${SILVA_SPECIES}" \
                "https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz"
            installed
        fi
        
        success "SILVA 138.1 databases downloaded."
    else
        warn "Skipped SILVA download. You can download later — see README.md."
    fi
fi

# ============================================================================
# STEP 9: Create configuration file
# ============================================================================
step "9" "Checking configuration file"

CONFIG_FILE="${PROJECT_DIR}/app/config.py"

if [[ -f "${CONFIG_FILE}" ]]; then
    skip "Configuration file (app/config.py)"
    info "To regenerate, delete app/config.py and re-run setup."
else
    # Detect PICRUSt2 location
    if [[ "${PICRUST2_SEPARATE}" == "true" ]]; then
        PICRUST2_ENV="picrust2"
    else
        PICRUST2_ENV="${ENV_NAME}"
    fi

    # Detect conda base path
    CONDA_BASE_PATH=$(conda info --base 2>/dev/null || echo "$HOME/miniforge3")

    cat > "${CONFIG_FILE}" << PYEOF
"""
MicrobiomeDash Configuration
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
PICRUST2_ENV_NAME = "${PICRUST2_ENV}"
CONDA_BASE = Path(os.environ.get("CONDA_BASE", "${CONDA_BASE_PATH}"))


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
# STEP 10: Verify installation
# ============================================================================
step "10" "Verifying installation"

echo ""
info "Checking tools..."

# CLI tools verification
for tool_info in \
    "Python:python --version" \
    "R:R --version 2>&1 | head -1" \
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

# PICRUSt2
if has_cmd picrust2_pipeline.py; then
    echo -e "  ${GREEN}✓${NC} PICRUSt2: available (main env)"
elif conda env list | grep -q "^picrust2 "; then
    echo -e "  ${GREEN}✓${NC} PICRUSt2: available (separate env)"
else
    echo -e "  ${YELLOW}△${NC} PICRUSt2: not found (pathway analysis will be unavailable)"
fi

# SILVA
if [[ -f "${SILVA_TRAIN}" && -f "${SILVA_SPECIES}" ]]; then
    echo -e "  ${GREEN}✓${NC} SILVA 138.1: downloaded"
else
    echo -e "  ${YELLOW}△${NC} SILVA 138.1: not downloaded (see README)"
fi

# Python packages
echo ""
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

# R packages
echo ""
info "Checking R packages..."
Rscript -e '
    packages <- c("dada2", "ANCOMBC", "ALDEx2", "DESeq2", "Maaslin2", "LinDA", "vegan", "phyloseq", "optparse", "jsonlite")
    for (pkg in packages) {
        if (requireNamespace(pkg, quietly = TRUE)) {
            message("  \033[0;32m✓\033[0m ", pkg)
        } else {
            message("  \033[0;31m✗\033[0m ", pkg, " — not installed")
        }
    }
'

# ============================================================================
# Summary
# ============================================================================
echo ""
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "  Summary: ${GREEN}${INSTALLED_COUNT} newly installed${NC} | ${BLUE}${SKIPPED_COUNT} skipped (already present)${NC}"
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

if [[ ${INSTALLED_COUNT} -eq 0 ]]; then
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
