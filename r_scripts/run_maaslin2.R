#!/usr/bin/env Rscript
# ============================================================================
# MicrobiomeDash — MaAsLin2 differential abundance
#
# Standard interface: --counts, --metadata, --group_col, --ref_group,
#                     --test_group, --output
# Output: TSV with feature, log2fc, pvalue, qvalue
# ============================================================================

suppressPackageStartupMessages({
  library(Maaslin2)
  library(optparse)
  library(jsonlite)
})

option_list <- list(
  make_option("--counts",     type="character", help="Count matrix TSV"),
  make_option("--metadata",   type="character", help="Metadata TSV"),
  make_option("--group_col",  type="character", help="Grouping column name"),
  make_option("--ref_group",  type="character", help="Reference group"),
  make_option("--test_group", type="character", help="Test group"),
  make_option("--output",     type="character", help="Output TSV path"),
  make_option("--threads",    type="integer", default=1, help="Number of parallel cores")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Read data
counts <- read.csv(opt$counts, sep="\t", row.names=1, check.names=FALSE)
meta   <- read.csv(opt$metadata, sep="\t", check.names=FALSE)

# Align samples
common <- intersect(colnames(counts), meta$SampleID)
counts <- counts[, common, drop=FALSE]
meta   <- meta[match(common, meta$SampleID), ]
rownames(meta) <- meta$SampleID

# Set reference level
meta[[opt$group_col]] <- factor(meta[[opt$group_col]],
                                levels=c(opt$ref_group, opt$test_group))

# MaAsLin2 expects samples as rows, features as columns
counts_t <- as.data.frame(t(counts))

# Run MaAsLin2 (output to temp directory)
output_dir <- tempdir()
maaslin_dir <- file.path(output_dir, "maaslin2_output")

n_cores <- max(1L, opt$threads)
fit <- Maaslin2(
  input_data     = counts_t,
  input_metadata = meta,
  output         = maaslin_dir,
  fixed_effects  = c(opt$group_col),
  reference      = paste0(opt$group_col, ",", opt$ref_group),
  normalization  = "TSS",
  transform      = "LOG",
  analysis_method = "LM",
  min_abundance  = 0.0,
  min_prevalence = 0.1,
  plot_heatmap   = FALSE,
  plot_scatter   = FALSE,
  cores          = n_cores
)

# Read results
res <- fit$results

results_df <- data.frame(
  feature = res$feature,
  log2fc  = res$coef / log(2),  # MaAsLin2 coef is on log scale
  pvalue  = res$pval,
  qvalue  = res$qval,
  stringsAsFactors = FALSE
)

results_df <- results_df[order(results_df$qvalue), ]

write.table(results_df, file=opt$output, sep="\t", row.names=FALSE, quote=FALSE)

n_sig <- sum(results_df$qvalue < 0.05, na.rm=TRUE)
cat(toJSON(list(
  status = "success",
  n_features = nrow(results_df),
  n_significant = n_sig
), auto_unbox=TRUE))
cat("\n")
