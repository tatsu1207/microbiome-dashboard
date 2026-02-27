#!/usr/bin/env Rscript
# ============================================================================
# MicrobiomeDash — LinDA differential abundance
#
# Standard interface: --counts, --metadata, --group_col, --ref_group,
#                     --test_group, --output
# Output: TSV with feature, log2fc, pvalue, qvalue
# ============================================================================

suppressPackageStartupMessages({
  library(LinDA)
  library(optparse)
  library(jsonlite)
})

option_list <- list(
  make_option("--counts",     type="character", help="Count matrix TSV"),
  make_option("--metadata",   type="character", help="Metadata TSV"),
  make_option("--group_col",  type="character", help="Grouping column name"),
  make_option("--ref_group",  type="character", help="Reference group"),
  make_option("--test_group", type="character", help="Test group"),
  make_option("--output",     type="character", help="Output TSV path")
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

# Run LinDA
formula_str <- paste0("~ ", opt$group_col)
linda_out <- linda(otu.tab=as.data.frame(counts),
                   meta=meta,
                   formula=formula_str,
                   alpha=0.05)

# Extract results for the group comparison
# LinDA output structure: $output is a list keyed by variable name
var_name <- paste0(opt$group_col, opt$test_group)
res <- linda_out$output[[var_name]]

if (is.null(res)) {
  # Try first available result
  res <- linda_out$output[[1]]
}

results_df <- data.frame(
  feature = rownames(res),
  log2fc  = res$log2FoldChange,
  pvalue  = res$pvalue,
  qvalue  = res$padj,
  stringsAsFactors = FALSE
)

results_df <- results_df[!is.na(results_df$qvalue), ]
results_df <- results_df[order(results_df$qvalue), ]

write.table(results_df, file=opt$output, sep="\t", row.names=FALSE, quote=FALSE)

n_sig <- sum(results_df$qvalue < 0.05, na.rm=TRUE)
cat(toJSON(list(
  status = "success",
  n_features = nrow(results_df),
  n_significant = n_sig
), auto_unbox=TRUE))
cat("\n")
