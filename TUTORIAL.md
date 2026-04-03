# 16S-Pipeline: Step-by-Step Tutorial

This tutorial walks through a complete 16S rRNA amplicon analysis using 16S-Pipeline, from raw FASTQ files to publication-ready results. We use a published soil microbiome dataset (Veselovsky et al., 2025; BioProject [PRJNA1192699](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1192699)) comparing two soil types (SoilA and SoilB) across three sequencing platforms.

**What you will learn:**
- Download public datasets from NCBI SRA
- Run the automated DADA2 pipeline
- Perform alpha and beta diversity analysis
- Visualize taxonomic composition
- Run multi-method differential abundance testing
- Generate a PDF analysis report

**Time required:** ~30 minutes (excluding pipeline processing time)

> **Quick start option:** The repository includes a `test_samples/` folder containing subsampled versions of all datasets used in this tutorial (5,000 reads per sample, ~30 files total). You can skip the SRA download step and upload these files directly for a faster test run. A pre-made `metadata.tsv` is also included. The reduced read count means results will differ from full-depth analysis, but all pipeline features can be tested.

---

## Prerequisites

16S-Pipeline must be installed and running. See the [README](README.md) for installation instructions.

**Start the server:**

```bash
# Native installation
bash run.sh

# Docker
docker compose up -d
```

Open your browser and navigate to `http://localhost:8016` (Docker) or the URL shown in the terminal.

---

## Step 1: Download Data from NCBI SRA

We will download V4 Illumina paired-end data (6 samples: 3 SoilA + 3 SoilB).

1. Click **SRA Download** in the left sidebar.
2. Enter the following SRR accession numbers:

```
SRR31516611
SRR31516613
SRR31516614
SRR31516615
SRR31516622
SRR31516623
```

3. Enter a **Study name** (e.g., `Soil_V4`).
4. Click **Download**. The platform will:
   - Resolve each accession via the NCBI E-utilities API
   - Download FASTQ files using `prefetch` + `fasterq-dump`
   - Automatically detect paired-end sequencing, V4 variable region, and Illumina platform
   - Register files in the File Manager

> **Tip:** You can also enter a BioProject accession (e.g., `PRJNA1192699`) to download all runs at once. However, this project contains samples from multiple platforms, so downloading by individual SRR numbers gives you more control.

---

## Step 2: Upload Sample Metadata

Metadata associates each sample with experimental groups for downstream statistical analyses.

1. Click **File Manager** in the left sidebar.
2. In the **Upload Sample Metadata** section, click **Download Template** to get a pre-filled TSV template with your sample names.
3. Open the template in a spreadsheet editor and add a `group` column:

| sample_name   | group  |
|---------------|--------|
| SRR31516611   | soil-A |
| SRR31516613   | soil-A |
| SRR31516614   | soil-A |
| SRR31516615   | soil-B |
| SRR31516622   | soil-B |
| SRR31516623   | soil-B |

4. Save as TSV and upload it via the metadata upload area.
5. A preview table confirms the metadata was matched to your samples.

---

## Step 3: Run the DADA2 Pipeline

1. Click **Pipeline** in the left sidebar.
2. In the **Select Samples** table, check all 6 samples.
3. Enter a **Run name** (e.g., `V4_soil`).
4. Leave all parameters at their defaults:
   - **Truncation lengths**: Auto-detected (the platform analyzes quality scores using a 10 bp sliding window to find where mean quality drops below Q20)
   - **Threads**: Automatically set based on available CPUs
5. Click **Launch Pipeline**.

The pipeline executes the following steps automatically:

| Step | Tool | Output |
|------|------|--------|
| 1. Quality control | FastQC | Per-sample quality reports |
| 2. Primer trimming | Cutadapt | Trimmed FASTQ files |
| 3. Quality parameter selection | Auto-truncation | Optimal truncation lengths |
| 4. ASV inference | DADA2 | ASV count table + representative sequences |
| 5. Taxonomy assignment | SILVA v138.1 | Taxonomy classifications |
| 6. Phylogenetic tree | MAFFT + FastTree | Newick tree file |
| 7. BIOM generation | biom-format | HDF5 BIOM file |

Monitor progress in real-time through the status panel. Processing typically takes 5-15 minutes depending on your hardware.

> **What happens during auto-truncation?** For this V4 dataset, the algorithm typically selects forward truncation at ~231 bp and reverse at ~121 bp, yielding ~131 bp of overlap for successful read merging.

---

## Step 4: Data Management and Quality Control

Once the pipeline completes, 16S-Pipeline provides several tools to inspect, filter, and prepare your data before statistical analysis.

### 4a. BIOM Browser

1. Click **BIOM Browser** in the left sidebar.
2. Select your completed dataset from the dropdown.
3. The browser displays:
   - **Summary statistics**: Total samples, ASV count, total reads, sequencing depth range
   - **Per-sample read counts**: A table showing reads per sample
   - **Taxonomy overview**: Classification rates at each taxonomic rank

This is a quick way to verify your pipeline ran successfully.

### 4b. Outlier Detection

Before analysis, check for samples with atypical community compositions that may indicate technical issues.

1. Click **Beta Diversity** in the left sidebar.
2. Run a PCoA ordination — samples that cluster far from all others may be outliers.
3. The platform identifies samples with atypical community compositions based on beta diversity distances, helping you decide whether to exclude them from downstream analysis.

### 4c. Rare ASV Removal

Low-prevalence and low-abundance ASVs can introduce noise into statistical analyses.

1. Click **Rare ASV Removal** in the left sidebar.
2. Select your dataset.
3. Set filtering thresholds:
   - **Minimum prevalence**: Remove ASVs present in fewer than N% of samples (e.g., 5%)
   - **Minimum relative abundance**: Remove ASVs below a relative abundance threshold (e.g., 0.01%)
4. Click **Filter**. A new filtered dataset is created — the original is preserved.

> **When to use:** Recommended before differential abundance analysis to reduce false positives from rare, noisy ASVs.

### 4d. Subsampling (Rarefaction)

Uneven sequencing depth can bias diversity comparisons. Rarefaction normalizes all samples to the same depth.

1. Click **Subsampling** in the left sidebar.
2. Select your dataset. The interface shows per-sample read counts so you can choose an appropriate rarefaction depth.
3. Set the **rarefaction depth** (typically the minimum sample depth, or slightly below to retain more samples).
4. Optionally **exclude low-depth samples** that fall below your chosen threshold.
5. Click **Rarefy**. A new rarefied dataset is created.

> **When to use:** Recommended before alpha diversity analysis when sequencing depths vary substantially between samples. Not always necessary for beta diversity or differential abundance, which have their own normalization methods.

### 4e. Combine Datasets

If you have processed multiple sequencing runs or want to perform cross-study meta-analysis:

1. Click **Combine** in the left sidebar.
2. Select two or more completed datasets.
3. Choose the merge mode:
   - **By sequence** (same variable region): Merges ASV tables by matching identical sequences
   - **By taxonomy** (cross-region): Harmonizes different variable regions using E. coli 16S rRNA alignment positions to identify corresponding V-region boundaries
4. Click **Combine**. A new merged dataset is created.

> **Example use case:** Combine your V4 and V3-V4 datasets to compare results across sequencing strategies on the same samples.

### 4f. V-Region Extraction

Extract specific variable regions from processed datasets for targeted comparisons.

1. Click **Datasets** in the left sidebar.
2. Select a dataset and choose the target variable region to extract.
3. The platform identifies ASV sequences corresponding to the selected V-region using E. coli 16S reference alignment positions.

---

## Step 5: Alpha Diversity Analysis

Alpha diversity measures within-sample diversity.

1. Click **Alpha Diversity** in the left sidebar.
2. Select your dataset and choose `group` as the grouping variable.
3. The platform computes multiple metrics:
   - **Shannon entropy** — Accounts for both richness and evenness
   - **Simpson index** — Probability that two randomly chosen individuals belong to different species
   - **Observed ASVs** — Raw species richness
   - **Chao1 estimator** — Estimated total richness including unobserved species
   - **ACE estimator** — Abundance-based coverage estimator

4. Box plots are generated with Kruskal-Wallis test p-values and Dunn's post-hoc pairwise comparisons.
5. Click the camera icon in the plot toolbar to save figures as SVG. Download the diversity table as CSV using the download button.

> **Tip:** Plots are rendered with a dark background by default. To get white-background figures for publications, toggle the **Light/Dark** switch in the sidebar — all plots will update to white backgrounds, then use the camera icon to save.

> **Expected result:** SoilB shows significantly higher Shannon diversity than SoilA (Kruskal-Wallis p=0.0495).

---

## Step 6: Beta Diversity Analysis

Beta diversity measures between-sample differences in community composition.

1. Click **Beta Diversity** in the left sidebar.
2. Select your dataset and choose `group` as the grouping variable.
3. Choose a distance metric:
   - **Bray-Curtis** — Quantitative (considers abundance)
   - **Jaccard** — Qualitative (presence/absence only)

4. Choose an ordination method:
   - **PCoA** (Principal Coordinates Analysis) — Linear ordination
   - **NMDS** (Non-Metric Multidimensional Scaling) — Non-linear, preserves rank-order distances

5. The plot shows sample clustering with 95% confidence ellipses. PERMANOVA results (pseudo-F statistic and p-value) test whether groups differ significantly.

> **Expected result:** Clear separation of SoilA and SoilB along PC1 (explaining ~64% of variance), with PERMANOVA pseudo-F ≈ 6.9.

---

## Step 7: Taxonomic Composition

1. Click **Taxonomy** in the left sidebar.
2. Select your dataset.
3. Choose a taxonomic rank (Phylum, Class, Order, Family, Genus, or Species).
4. Adjust the number of top taxa to display (remaining are collapsed into "Others").

Interactive stacked bar plots show relative abundance per sample. At the phylum level, expect to see Actinobacteriota, Proteobacteria, and Firmicutes as dominant phyla in soil samples.

---

## Step 8: Differential Abundance Analysis

A key feature of 16S-Pipeline is the integration of five complementary differential abundance methods.

1. Click **Differential Abundance** in the left sidebar.
2. Select your dataset and choose `group` as the grouping variable.
3. Select one or more DA methods:

| Method | Approach | Characteristics |
|--------|----------|-----------------|
| **ALDEx2** | CLR transformation + Welch's t-test | Conservative, compositionally aware |
| **DESeq2** | Negative binomial GLM | Widely used, can be liberal with small n |
| **ANCOM-BC2** | Bias-corrected log-linear model | Very conservative, especially at small n |
| **LinDA** | Linear models on log-transformed data | Balanced sensitivity |
| **MaAsLin2** | Multivariable association | Handles complex designs |

4. Adjust **Filter Criteria** to focus on significant results:
   - **p-value**: Raw p-value threshold (e.g., < 0.05)
   - **q-value**: BH-adjusted p-value threshold (e.g., < 0.05) — recommended over raw p-value
   - **Abundance Diff (|log2FC|)**: Minimum fold-change magnitude (e.g., > 1.0 = 2-fold change)
   - **Effect Size (|effect|)**: ALDEx2-specific effect size filter (e.g., > 1.0)

5. Click **Run**. Each method produces:
   - A **volcano plot** showing effect size (x-axis) vs. significance (y-axis), with significant ASVs highlighted
   - A **results table** with log-fold changes, p-values, and adjusted p-values

6. The **Consensus** tab highlights ASVs detected as significant by multiple methods, providing a robust set of differentially abundant taxa.

> **Why use multiple methods?** Different DA methods can yield substantially different results on the same dataset. The consensus approach identifies the most robust findings — ASVs detected by 3+ methods are high-confidence candidates.

---

## Step 9: Functional Prediction with PICRUSt2

PICRUSt2 predicts functional potential (metabolic pathways, enzyme abundances) from 16S ASV data.

> **Note:** PICRUSt2 requires 11+ GB of RAM. This step is optional and does not affect other analyses.

1. Click **PICRUSt2** in the left sidebar.
2. Select your completed dataset from the dropdown.
3. Click **Run PICRUSt2**. The process takes 10-30 minutes depending on sample count and hardware.
4. Once complete, PICRUSt2 output files (KO metagenome, EC metagenome, MetaCyc pathways) are available for downstream analysis.

---

## Step 10: Pathway and KEGG Map Analysis

Once PICRUSt2 results are available, you can perform functional differential abundance analysis.

### Pathway Comparison

1. Click **Pathways** in the left sidebar.
2. Select the PICRUSt2 run and your metadata.
3. Choose `group` as the grouping variable and select DA methods (same five methods as Step 8).
4. Click **Run**. The platform applies the multi-method DA framework to predicted MetaCyc pathway abundances, generating:
   - **Error bar plots** showing pathway abundance differences between groups
   - **Heatmaps** of significantly different pathways
   - **PCA plots** of functional profiles

### KEGG Pathway Maps

1. Click **KEGG Map** in the left sidebar.
2. Select the PICRUSt2 run and upload your metadata.
3. Choose `group` as the grouping variable.
4. Select a group comparison (e.g., `soil-A vs soil-B`) for pairwise differential abundance coloring.
5. Choose a KEGG pathway (e.g., `map00190 — Oxidative phosphorylation`).
6. The platform:
   - Aggregates predicted KO and EC numbers to the selected KEGG pathway
   - Displays **KO coverage** (e.g., 106/224 KOs detected)
   - Shows a **pathway activity box plot** comparing the two selected groups
   - Renders the **KEGG pathway map** with detected enzymes highlighted — enzymes are colored by direction and significance of differential abundance (e.g., red = enriched in group A, blue = enriched in group B)
   - Provides a **per-KO abundance table** with Mann-Whitney U test statistics, p-values, and mean abundances per group

> **Tip:** KEGG API responses are cached for 24 hours to avoid rate limits. The first load of a pathway may take a few seconds. You can switch between different pairwise comparisons to see how enzyme enrichment patterns change.

---

## Step 11: Generate a PDF Report

1. Click **Analysis Report** in the left sidebar.
2. Select your dataset.
3. Choose `group` as the grouping variable.
4. Select which sections to include:
   - Dataset Summary
   - Materials & Methods (auto-generated text)
   - Alpha Diversity
   - Beta Diversity
   - Taxonomy Composition
   - Read Tracking
5. Click **Generate PDF Report**.

The report compiles all analyses into a single downloadable PDF, suitable for lab meetings, collaborator updates, or manuscript preparation.

---

## Optional: Repeat with Different Platforms

To compare results across sequencing platforms (as demonstrated in the manuscript), you can process additional datasets from the same study:

**V3-V4 Illumina paired-end:**
```
SRR31527673
SRR31527675
SRR31527676
SRR31527677
SRR31527684
SRR31527685
```

**PacBio HiFi full-length 16S:**
```
SRR31532541
SRR31532543
SRR31532544
SRR31532545
SRR31532552
SRR31532553
```

Each dataset is auto-detected and processed with platform-appropriate parameters. PacBio data uses long-read mode with PacBioErrfun error models and a 1,000-1,600 bp length filter. Full-length 16S sequences enable species-level taxonomic classification not achievable with short-read approaches.

---

## Additional Features

### Data Management
- **Subsampling**: Rarefy datasets to uniform sequencing depth (recommended before diversity analysis with uneven depth)
- **Rare ASV Filtering**: Remove low-prevalence ASVs to reduce noise
- **Dataset Combining**: Merge multiple datasets for cross-study meta-analysis
- **MOTHUR Conversion**: Export to MOTHUR-compatible shared/taxonomy format

### Functional Prediction
- **PICRUSt2**: Predict functional potential from 16S data (requires 11+ GB RAM)
- **Pathway Analysis**: Compare MetaCyc pathway abundances between groups using the same multi-method DA framework
- **KEGG Pathway Maps**: Interactive visualization of predicted KO/EC numbers mapped onto KEGG pathway diagrams

### SRA Submission Helper
- Generate NCBI SRA submission metadata spreadsheets from your registered files
- Auto-populates library strategy, source, layout, platform, and filenames

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Pipeline stuck at "processing" | Check if the R process is still running; the UI recovers via PID monitoring |
| DADA2 fails with SIGKILL | Insufficient RAM — try reducing the number of samples or threads |
| PICRUSt2 fails | Requires 11+ GB RAM; this step is optional and does not affect core results |
| "No files found" on upload | Ensure files end with `.fastq.gz` or `.fq.gz` |
| Differential abundance errors | Reduce thread count if you see "all connections in use" errors |

For additional help, visit the [GitHub repository](https://github.com/tatsu1207/16S-Pipeline).

---

## References

- Veselovsky V, et al. 2025. Comparative evaluation of sequencing platforms for 16S rRNA-based soil microbiome profiling. *Frontiers in Microbiology*. Volume 16.
- Callahan BJ, et al. 2016. DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*. 13: 581-583.
- Nearing JT, et al. 2022. Microbiome differential abundance methods produce different results across 38 datasets. *Nature Communications*. 13: 342.
