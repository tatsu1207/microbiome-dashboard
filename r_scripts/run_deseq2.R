#!/usr/bin/env Rscript
# ============================================================================
# MicrobiomeDash — DESeq2 differential abundance
#
# Standard interface: --counts, --metadata, --group_col, --ref_group,
#                     --test_group, --output
# Output: TSV with feature, log2fc, pvalue, qvalue
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
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

# Ensure integer counts
counts <- round(counts)

# Set reference level
meta[[opt$group_col]] <- factor(meta[[opt$group_col]],
                                levels=c(opt$ref_group, opt$test_group))

# Create DESeq2 dataset
formula_str <- as.formula(paste0("~ ", opt$group_col))
dds <- DESeqDataSetFromMatrix(countData=as.matrix(counts),
                              colData=meta,
                              design=formula_str)

# Filter low-count features
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Run DESeq2 (poscounts handles sparse microbiome data where every gene has zeros)
n_cores <- max(1L, opt$threads)
if (n_cores > 1) {
  suppressPackageStartupMessages(library(BiocParallel))
  register(MulticoreParam(workers=n_cores))
  dds <- DESeq(dds, sfType="poscounts", quiet=TRUE, parallel=TRUE)
} else {
  dds <- DESeq(dds, sfType="poscounts", quiet=TRUE)
}
res <- results(dds, contrast=c(opt$group_col, opt$test_group, opt$ref_group))
res <- as.data.frame(res)

# Build output
results_df <- data.frame(
  feature = rownames(res),
  log2fc  = res$log2FoldChange,
  pvalue  = res$pvalue,
  qvalue  = res$padj,
  stringsAsFactors = FALSE
)

# Remove NA rows and sort
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
