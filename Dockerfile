# ============================================================================
# 16S Pipeline — Dockerfile
# Builds all 5 conda environments + SILVA references into a single image.
# Image size: ~10-15 GB (5 conda envs + reference databases)
#
# Build:   docker build -t 16s-pipeline .
# Run:     docker compose up -d
# Open:    http://localhost:8016
# ============================================================================

FROM condaforge/mambaforge:latest

LABEL maintainer="16S Pipeline"
LABEL description="End-to-end 16S rRNA microbiome analysis platform"

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# ── System dependencies ─────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    zlib1g-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    wget \
    procps \
    lsof \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# ── Conda environment 1: microbiome_16S (Python + CLI tools) ───────────────
RUN mamba create -n microbiome_16S -c conda-forge python=3.11 -y && \
    mamba clean -afy

# Install bioinformatics CLI tools
RUN mamba install -n microbiome_16S --override-channels -c conda-forge -c bioconda \
    fastqc cutadapt mafft fasttree bbmap sra-tools -y && \
    mamba clean -afy

# Install Python packages via pip
RUN conda run -n microbiome_16S pip install --no-cache-dir \
    fastapi \
    "uvicorn[standard]" \
    dash \
    dash-bootstrap-components \
    plotly \
    sqlalchemy \
    pandas \
    numpy \
    scipy \
    scikit-bio \
    python-multipart \
    biom-format \
    matplotlib \
    fpdf2 \
    statsmodels \
    dash-uploader \
    matplotlib-venn

# ── Conda environment 2: dada2_16S (R + DADA2) ────────────────────────────
RUN mamba create -n dada2_16S --override-channels -c conda-forge -c bioconda \
    bioconductor-dada2 r-optparse r-jsonlite -y && \
    mamba clean -afy

# ── Conda environment 3: analysis_16S (R + DA tools) ──────────────────────
RUN mamba create -n analysis_16S --override-channels -c conda-forge -c bioconda \
    bioconductor-phyloseq bioconductor-ancombc bioconductor-deseq2 \
    bioconductor-aldex2 r-optparse r-jsonlite -y && \
    mamba clean -afy

# ── Conda environment 4: maaslin2_16S (MaAsLin2 + vegan + LinDA) ─────────
RUN mamba create -n maaslin2_16S --override-channels -c conda-forge -c bioconda \
    bioconductor-maaslin2 r-optparse r-jsonlite -y && \
    mamba clean -afy

# Install vegan from CRAN
RUN conda run -n maaslin2_16S Rscript -e \
    "install.packages('vegan', repos='https://cloud.r-project.org', INSTALL_opts='--no-lock', Ncpus=4)"

# Install LinDA from GitHub
RUN conda run -n maaslin2_16S Rscript -e " \
    install.packages('remotes', repos='https://cloud.r-project.org', INSTALL_opts='--no-lock'); \
    install.packages('modeest', repos='https://cloud.r-project.org', INSTALL_opts='--no-lock'); \
    remotes::install_github('zhouhj1994/LinDA', upgrade='never', INSTALL_opts='--no-lock')"

# ── Conda environment 5: picrust2_16S ─────────────────────────────────────
# PICRUSt2 has no linux-aarch64 package on bioconda; skip on arm64 builds
ARG TARGETARCH
RUN if [ "$TARGETARCH" = "amd64" ]; then \
        mamba create -n picrust2_16S --override-channels -c conda-forge -c bioconda \
            picrust2 -y && \
        mamba clean -afy; \
    else \
        echo "Skipping PICRUSt2 on $TARGETARCH (no bioconda package available)"; \
    fi

# ── Download SILVA 138.1 references to a staging location ─────────────────
# Stored in /opt/silva so they can be copied into the data volume at first run
RUN mkdir -p /opt/silva && \
    wget -q -O /opt/silva/silva_nr99_v138.1_train_set.fa.gz \
        "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz" && \
    wget -q -O /opt/silva/silva_species_assignment_v138.1.fa.gz \
        "https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz"

# ── Copy application code ─────────────────────────────────────────────────
COPY app/ /app/app/
COPY r_scripts/ /app/r_scripts/
COPY data/references/ecoli_16S.fasta /opt/silva/ecoli_16S.fasta
COPY docker-entrypoint.sh /app/docker-entrypoint.sh
RUN chmod +x /app/docker-entrypoint.sh

# ── Environment variables ─────────────────────────────────────────────────
ENV PORT=8016
ENV CONDA_BASE=/opt/conda
# Store database inside the data volume so it persists
ENV DATABASE_PATH=/app/data/microbiome.db

EXPOSE 8016

ENTRYPOINT ["/app/docker-entrypoint.sh"]
