#!/usr/bin/env Rscript
# ============================================================================
# MicrobiomeDash — NMDS ordination via vegan::metaMDS
#
# Input:  --dist_matrix  (TSV distance matrix with row/col labels)
# Output: --output       (TSV with NMDS1, NMDS2 columns)
# Prints JSON status to stdout.
# ============================================================================

suppressPackageStartupMessages({
  library(vegan)
  library(optparse)
  library(jsonlite)
})

option_list <- list(
  make_option("--dist_matrix", type="character", help="Path to distance matrix TSV"),
  make_option("--output",      type="character", help="Output coordinates TSV path"),
  make_option("--k",           type="integer",   default=2, help="Number of dimensions")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Read distance matrix
dm_df <- read.csv(opt$dist_matrix, sep="\t", row.names=1, check.names=FALSE)
dm <- as.dist(as.matrix(dm_df))

# Run NMDS
set.seed(42)
nmds <- metaMDS(dm, k=opt$k, trymax=100, trace=0)

# Write coordinates
coords <- as.data.frame(nmds$points)
colnames(coords) <- paste0("NMDS", seq_len(ncol(coords)))
write.table(coords, file=opt$output, sep="\t", quote=FALSE)

# Print JSON status
cat(toJSON(list(
  stress = nmds$stress,
  converged = nmds$converged
), auto_unbox=TRUE))
cat("\n")
