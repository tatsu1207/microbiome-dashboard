#!/usr/bin/env Rscript
# ============================================================================
# MicrobiomeDash — ALDEx2 differential abundance
#
# Standard interface: --counts, --metadata, --group_col, --ref_group,
#                     --test_group, --output
# Output: TSV with feature, log2fc, pvalue, qvalue
# ============================================================================

suppressPackageStartupMessages({
  library(ALDEx2)
  library(optparse)
  library(jsonlite)
})

option_list <- list(
  make_option("--counts",     type="character", help="Count matrix TSV (feature_id + samples)"),
  make_option("--metadata",   type="character", help="Metadata TSV (SampleID + group col)"),
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

# Create conditions vector (ref=0, test=1)
conditions <- ifelse(meta[[opt$group_col]] == opt$ref_group, "ref", "test")

# Run ALDEx2
n_cores <- max(1L, opt$threads)
if (n_cores > 1) {
  suppressPackageStartupMessages(library(BiocParallel))
  bp <- MulticoreParam(workers=n_cores)
} else {
  bp <- NULL
}
aldex_out <- aldex(counts, conditions, mc.samples=128, test="t", effect=TRUE,
                   verbose=FALSE, BPPARAM=bp)

# Build output
results <- data.frame(
  feature = rownames(aldex_out),
  log2fc  = aldex_out$diff.btw,
  effect  = aldex_out$effect,
  pvalue  = aldex_out$wi.ep,
  qvalue  = aldex_out$wi.eBH,
  stringsAsFactors = FALSE
)

# Sort by qvalue
results <- results[order(results$qvalue), ]

write.table(results, file=opt$output, sep="\t", row.names=FALSE, quote=FALSE)

n_sig <- sum(results$qvalue < 0.05, na.rm=TRUE)
cat(toJSON(list(
  status = "success",
  n_features = nrow(results),
  n_significant = n_sig
), auto_unbox=TRUE))
cat("\n")
