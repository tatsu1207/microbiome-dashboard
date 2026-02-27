#!/usr/bin/env Rscript
# ============================================================================
# MicrobiomeDash — ANCOM-BC2 differential abundance
#
# Standard interface: --counts, --metadata, --group_col, --ref_group,
#                     --test_group, --output
# Output: TSV with feature, log2fc, pvalue, qvalue
# ============================================================================

suppressPackageStartupMessages({
  library(ANCOMBC)
  library(phyloseq)
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

# Create phyloseq object
otu  <- otu_table(as.matrix(counts), taxa_are_rows=TRUE)
samp <- sample_data(meta)
ps   <- phyloseq(otu, samp)

# Run ANCOM-BC2
n_cores <- max(1L, opt$threads)
out <- ancombc2(data=ps, fix_formula=opt$group_col, p_adj_method="BH",
                group=opt$group_col, verbose=FALSE, n_cl=n_cores)

res <- out$res

# Extract the test group column (named like group_coltest_group)
lfc_col <- paste0("lfc_", opt$group_col, opt$test_group)
p_col   <- paste0("p_", opt$group_col, opt$test_group)
q_col   <- paste0("q_", opt$group_col, opt$test_group)

# Handle column name variations
all_cols <- colnames(res)
lfc_match <- grep("^lfc_", all_cols, value=TRUE)[1]
p_match   <- grep("^p_", all_cols, value=TRUE)[1]
q_match   <- grep("^q_", all_cols, value=TRUE)[1]

if (!is.na(lfc_match)) lfc_col <- lfc_match
if (!is.na(p_match))   p_col   <- p_match
if (!is.na(q_match))   q_col   <- q_match

results <- data.frame(
  feature = res$taxon,
  log2fc  = res[[lfc_col]] / log(2),  # ANCOM-BC uses natural log
  pvalue  = res[[p_col]],
  qvalue  = res[[q_col]],
  stringsAsFactors = FALSE
)

results <- results[order(results$qvalue), ]
write.table(results, file=opt$output, sep="\t", row.names=FALSE, quote=FALSE)

n_sig <- sum(results$qvalue < 0.05, na.rm=TRUE)
cat(toJSON(list(
  status = "success",
  n_features = nrow(results),
  n_significant = n_sig
), auto_unbox=TRUE))
cat("\n")
