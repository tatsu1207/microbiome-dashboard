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

# Remove stale R library lock files that block subsequent installs
clean_r_locks() {
    if [[ -n "${ENV_PREFIX:-}" ]] && ls "${ENV_PREFIX}/lib/R/library"/00LOCK-* &>/dev/null; then
        rm -rf "${ENV_PREFIX}/lib/R/library"/00LOCK-*
    fi
}

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
# STEP 0a: Install system dependencies
# ============================================================================
step "0a" "Checking system libraries (needed for R package compilation)"

# These are required to compile dada2 (Rhtslib), curl/RCurl, openssl, xml2, and png from source.
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
        warn "Could not install system libraries. R packages that compile from source may fail."
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
# STEP 2: Create Conda environment
# ============================================================================
step "2" "Checking Conda environment '${ENV_NAME}'"

if conda env list | grep -q "^${ENV_NAME} "; then
    skip "Conda environment '${ENV_NAME}'"
    info "To recreate from scratch, run: conda env remove -n ${ENV_NAME} && ./setup_ubuntu.sh"
else
    info "Creating environment with Python ${PYTHON_VERSION} and R ${R_VERSION}..."
    if ! ${SOLVER} create -n "${ENV_NAME}" python=${PYTHON_VERSION} r-base=${R_VERSION} -y; then
        error "Failed to create conda environment '${ENV_NAME}'. Check disk space and network connectivity."
    fi
    success "Environment '${ENV_NAME}' created."
    installed
fi

# Activate environment
if ! conda activate "${ENV_NAME}" 2>/dev/null; then
    error "Failed to activate conda environment '${ENV_NAME}'. Try: conda env remove -n ${ENV_NAME} && ./setup_ubuntu.sh"
fi
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
# STEP 2b: Install C libraries inside conda env (needed by conda's GCC)
# ============================================================================
step "2b" "Checking C libraries in conda environment"

# Conda's R uses its own GCC which only searches conda lib paths, not system paths.
# These must be in the conda env for R packages to compile from source.
# Also includes R packages with C bindings (r-curl) that are hard to compile from source.
CONDA_CLIBS=(zlib libpng bzip2 xz openssl libcurl hdf5 glpk gmp r-curl r-ragg r-cairo r-microbenchmark)
CONDA_CLIBS_MISSING=()
for lib in "${CONDA_CLIBS[@]}"; do
    if ${SOLVER} list -n "${ENV_NAME}" "^${lib}$" 2>/dev/null | grep -q "${lib}"; then
        skip "Conda lib: ${lib}"
    else
        CONDA_CLIBS_MISSING+=("${lib}")
    fi
done

if [ ${#CONDA_CLIBS_MISSING[@]} -gt 0 ]; then
    info "Installing C libraries in conda env: ${CONDA_CLIBS_MISSING[*]}"
    for lib in "${CONDA_CLIBS_MISSING[@]}"; do
        ${SOLVER} install -n "${ENV_NAME}" -c conda-forge "${lib}" -y 2>&1 | tail -3 || \
            warn "Failed to install conda lib: ${lib}"
    done
    success "Conda C libraries installed."
    installed
else
    info "All C libraries already present in conda env."
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
    info "Installing missing tools individually: ${BIOTOOLS_TO_INSTALL[*]}"

    for pkg in "${BIOTOOLS_TO_INSTALL[@]}"; do
        info "Installing ${pkg} ..."
        run_with_dots "conda: ${pkg}" \
            ${SOLVER} install -n "${ENV_NAME}" -c bioconda -c conda-forge "${pkg}" -y
        echo "${RWD_OUTPUT:-}" | tail -5 || true
        if [[ ${RWD_EXIT} -ne 0 ]]; then
            warn "Conda/Mamba failed for ${pkg}:"
            echo "${RWD_OUTPUT:-}" | tail -20 || true
        fi

        # Verify the tool is actually available; fall back to pip if conda failed
        for tool_pair in "${BIOTOOLS_MAP[@]}"; do
            map_pkg="${tool_pair%%:*}"
            if [[ "${map_pkg}" == "${pkg}" ]]; then
                cmds="${tool_pair##*:}"
                IFS=',' read -ra CMD_ARRAY <<< "${cmds}"
                found=false
                for cmd in "${CMD_ARRAY[@]}"; do
                    if has_env_cmd "${cmd}"; then
                        found=true
                        break
                    fi
                done
                if ${found}; then
                    success "Verified: ${pkg}"
                    installed
                else
                    # Fallback: try pip install (works for pure-Python tools like cutadapt)
                    warn "Conda failed for ${pkg}, trying pip install as fallback..."
                    pip install "${pkg}" 2>&1 | tail -5 || true
                    # Re-check
                    for cmd in "${CMD_ARRAY[@]}"; do
                        if has_env_cmd "${cmd}"; then
                            found=true
                            break
                        fi
                    done
                    if ${found}; then
                        success "Verified (via pip): ${pkg}"
                        installed
                    else
                        failed "CLI tool: ${pkg}"
                    fi
                fi
                break
            fi
        done
    done
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

    run_with_dots "conda: picrust2" \
        ${SOLVER} create -n picrust2 -c bioconda -c conda-forge picrust2 -y
    echo "${RWD_OUTPUT:-}" | tail -5 || true
    if [[ ${RWD_EXIT} -eq 0 ]]; then
        success "PICRUSt2 installed in separate 'picrust2' environment."
        installed
    else
        warn "PICRUSt2 installation failed. You can install it manually later."
        warn "The pipeline will skip PICRUSt2 steps if it's not available."
    fi
fi

# ============================================================================
# STEP 5: Install R packages (BiocManager/CRAN primary, conda not used)
# ============================================================================
step "5" "Checking R packages"

# Ensure R can find packages in the conda env during byte-compilation.
# Without this, R's findpack() fails for packages like igraph during lazy loading.
R_LIB_DIR="${ENV_PREFIX}/lib/R/library"
export R_LIBS_USER="${R_LIB_DIR}"
export R_LIBS_SITE="${R_LIB_DIR}"

# Clean up stale R lock files from previously interrupted installs
if ls "${R_LIB_DIR}"/00LOCK-* &>/dev/null; then
    info "Removing stale R lock files from previous interrupted installs..."
    rm -rf "${R_LIB_DIR}"/00LOCK-*
    success "Lock files cleaned up."
fi

# --- 5a: Ensure BiocManager is installed and version-matched to R ---
info "Checking BiocManager version matches R..."
BIOCMGR_SCRIPT=$(mktemp /tmp/biocmgr_XXXXXX.R)
cat > "${BIOCMGR_SCRIPT}" << 'REOF'
if (!requireNamespace('BiocManager', quietly=TRUE))
    install.packages('BiocManager', repos='https://cloud.r-project.org', INSTALL_opts='--no-lock')

tryCatch({
    BiocManager::install(ask=FALSE, update=FALSE)
    message(sprintf('BiocManager %s OK for %s', BiocManager::version(), R.version.string))
}, error = function(e) {
    message('BiocManager version mismatch detected, upgrading...')
    m <- regmatches(e$message, regexpr("version = '[0-9.]+'", e$message))
    if (length(m) > 0) {
        ver <- gsub("version = '|'", '', m)
        message(sprintf('Upgrading to Bioconductor %s', ver))
        BiocManager::install(version=ver, ask=FALSE)
    } else {
        message('Could not auto-detect version, trying fresh BiocManager...')
        install.packages('BiocManager', repos='https://cloud.r-project.org', INSTALL_opts='--no-lock')
        BiocManager::install(ask=FALSE, update=FALSE)
    }
    message(sprintf('BiocManager now at %s', BiocManager::version()))
})
REOF
run_with_dots "BiocManager check" Rscript "${BIOCMGR_SCRIPT}"
echo "${RWD_OUTPUT:-}" | tail -5 || true
rm -f "${BIOCMGR_SCRIPT}"
if [[ ${RWD_EXIT} -ne 0 ]]; then
    warn "BiocManager version check had issues (may be OK)."
fi

# Update CRAN packages that conda ships as outdated versions.
# lifecycle >= 1.0.5 is required by treeio and many other modern R packages.
info "Updating core R dependencies from CRAN..."
run_with_dots "R: core deps" \
    Rscript -e "install.packages(c('lifecycle', 'rlang', 'cli', 'vctrs', 'pillar'), repos='https://cloud.r-project.org', INSTALL_opts='--no-lock', Ncpus=4)"
if [[ ${RWD_EXIT} -ne 0 ]]; then
    warn "Some core R deps may not have updated (continuing anyway)."
fi

# --- 5b-pre: treeio (Bioc 3.18's treeio 1.26.0 uses tidytree::random_ref
# which was removed from tidytree. Install latest from GitHub instead.) ---
if has_r_pkg "treeio"; then
    skip "R package: treeio"
else
    info "Installing treeio from GitHub (Bioc 3.18 version has a bug)..."
    clean_r_locks
    run_with_dots "R: treeio (GitHub)" \
        Rscript -e "
        if (!requireNamespace('remotes', quietly=TRUE))
            install.packages('remotes', repos='https://cloud.r-project.org', INSTALL_opts='--no-lock')
        remotes::install_github('YuLab-SMU/treeio', upgrade='never', INSTALL_opts='--no-lock')
    "
    if has_r_pkg "treeio"; then
        echo "${RWD_OUTPUT:-}" | tail -5 || true
        success "R package installed: treeio (GitHub)"
        installed
    else
        warn "Install output for treeio:"
        echo "${RWD_OUTPUT:-}" | tail -30 || true
        failed "R package: treeio"
    fi
fi

# --- 5b-pre2: ggrepel (needs older version for R 4.3; required by scater) ---
if has_r_pkg "ggrepel"; then
    skip "R package: ggrepel"
else
    info "Installing ggrepel 0.9.6 (0.9.7 requires R >= 4.5)..."
    clean_r_locks
    run_with_dots "R: ggrepel" \
        Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/ggrepel/ggrepel_0.9.6.tar.gz', repos=NULL, type='source', INSTALL_opts='--no-lock')"
    if has_r_pkg "ggrepel"; then
        success "R package installed: ggrepel 0.9.6"
        installed
    else
        warn "Install output for ggrepel:"
        echo "${RWD_OUTPUT:-}" | tail -20 || true
        failed "R package: ggrepel"
    fi
fi

# --- 5b: Install all R packages via BiocManager/CRAN ---
# BiocManager::install() handles both CRAN and Bioconductor packages.
# Packages are installed one-by-one to isolate failures.
# --no-lock avoids 00LOCK directory issues in conda R environments.
R_PACKAGES=(
    # CRAN packages (directly used by R scripts)
    optparse
    jsonlite
    vegan
    # Bioconductor packages (directly used by R scripts)
    dada2
    phyloseq
    # ANCOMBC dependency chain — installed individually before ANCOMBC so that
    # one failure doesn't cascade and block the rest.
    curl
    httr
    ggrastr
    TreeSummarizedExperiment
    scater
    mia
    ANCOMBC
    ALDEx2
    DESeq2
    Maaslin2
)

for r_pkg in "${R_PACKAGES[@]}"; do
    if has_r_pkg "${r_pkg}"; then
        skip "R package: ${r_pkg}"
    else
        info "Installing ${r_pkg} ..."
        clean_r_locks

        # Try normal install first. Fall back to --no-lock on lock errors.
        run_with_dots "R: ${r_pkg}" \
            Rscript -e "
            if (!requireNamespace('BiocManager', quietly=TRUE))
                install.packages('BiocManager', repos='https://cloud.r-project.org')
            BiocManager::install('${r_pkg}', ask=FALSE, update=FALSE)
        "
        if ! has_r_pkg "${r_pkg}"; then
            # Check if it was a lock failure and retry with --no-lock
            if echo "${RWD_OUTPUT:-}" | grep -q "failed to lock directory"; then
                warn "Lock error detected, retrying ${r_pkg} with --no-lock..."
                clean_r_locks
                run_with_dots "R: ${r_pkg} (no-lock)" \
                    Rscript -e "
                    if (!requireNamespace('BiocManager', quietly=TRUE))
                        install.packages('BiocManager', repos='https://cloud.r-project.org', INSTALL_opts='--no-lock')
                    BiocManager::install('${r_pkg}', ask=FALSE, update=FALSE, INSTALL_opts='--no-lock')
                "
            fi
        fi

        if has_r_pkg "${r_pkg}"; then
            echo "${RWD_OUTPUT:-}" | tail -5 || true
            success "R package installed: ${r_pkg}"
            installed
        else
            warn "Install output for ${r_pkg}:"
            echo "${RWD_OUTPUT:-}" | tail -30 || true
            failed "R package: ${r_pkg}"
        fi
    fi
done

# --- 5d: LinDA (GitHub-only R package) ---
if has_r_pkg "LinDA"; then
    skip "R package: LinDA"
else
    info "Installing LinDA from GitHub..."
    clean_r_locks
    run_with_dots "GitHub: LinDA" \
        Rscript -e "
        if (!requireNamespace('BiocManager', quietly=TRUE))
            install.packages('BiocManager', repos='https://cloud.r-project.org', INSTALL_opts='--no-lock')
        if (!requireNamespace('remotes', quietly=TRUE))
            install.packages('remotes', repos='https://cloud.r-project.org', INSTALL_opts='--no-lock')
        if (!requireNamespace('modeest', quietly=TRUE))
            install.packages('modeest', repos='https://cloud.r-project.org', INSTALL_opts='--no-lock')
        remotes::install_github('zhouhj1994/LinDA', upgrade='never', INSTALL_opts='--no-lock')
    "
    echo "${RWD_OUTPUT:-}" | tail -10 || true

    if has_r_pkg "LinDA"; then
        success "LinDA installed from GitHub."
        installed
    else
        failed "R package: LinDA"
    fi
fi

# ============================================================================
# STEP 6: Install Python packages
# ============================================================================
step "6" "Checking Python packages"

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
    info "Installing missing Python packages individually: ${PYTHON_MISSING[*]}"

    for pip_name in "${PYTHON_MISSING[@]}"; do
        # Find matching module from the map
        for pair in "${PYTHON_PACKAGES[@]}"; do
            map_module="${pair%%:*}"
            map_pip="${pair##*:}"
            if [[ "${map_pip}" == "${pip_name}" ]]; then
                info "Installing ${pip_name} ..."
                pip install "${pip_name}" 2>&1 | tail -5 || true

                if has_python_pkg "${map_module}"; then
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
