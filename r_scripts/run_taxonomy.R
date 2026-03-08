#!/usr/bin/env Rscript
# ============================================================================
# MicrobiomeDash — SILVA 138.1 Taxonomy Assignment
#
# Uses DADA2's assignTaxonomy + addSpecies with SILVA reference databases.
# Outputs: taxonomy.tsv (ASV_ID, Kingdom, Phylum, ..., Species)
# ============================================================================

suppressPackageStartupMessages({
  library(dada2)
  library(optparse)
  library(jsonlite)
})

# ── CLI Arguments ─────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--rep_seqs",      type="character", help="Path to representative sequences FASTA"),
  make_option("--output",        type="character", help="Output taxonomy TSV path"),
  make_option("--silva_train",   type="character", help="SILVA training set path"),
  make_option("--silva_species", type="character", help="SILVA species assignment path"),
  make_option("--threads",       type="integer",   default=1, help="Number of threads"),
  make_option("--skip_species",  action="store_true", default=FALSE,
              help="Skip species-level assignment (faster)")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$rep_seqs) || is.null(opt$output) || is.null(opt$silva_train)) {
  stop("--rep_seqs, --output, and --silva_train are required")
}
if (!opt$skip_species && is.null(opt$silva_species)) {
  stop("--silva_species is required unless --skip_species is set")
}

# ── Main ──────────────────────────────────────────────────────────────────────

tryCatch({

  # Read representative sequences
  cat("Reading representative sequences...\n")
  lines <- readLines(opt$rep_seqs)
  header_idx <- grep("^>", lines)
  asv_ids <- sub("^>", "", lines[header_idx])
  seqs <- lines[header_idx + 1]
  names(seqs) <- asv_ids

  cat("Loaded", length(seqs), "ASV sequences\n")

  use_mt <- ifelse(opt$threads > 1, opt$threads, FALSE)

  # Assign taxonomy to genus level
  cat("Assigning taxonomy (this may take a while)...\n")
  taxa <- assignTaxonomy(seqs, opt$silva_train,
                         multithread=use_mt, tryRC=FALSE)

  # Add species-level assignment (single-threaded, slow)
  if (!opt$skip_species) {
    cat("Adding species assignments...\n")
    taxa <- addSpecies(taxa, opt$silva_species)
  } else {
    cat("Skipping species assignment (--skip_species)\n")
  }

  # Build output data frame
  tax_df <- data.frame(
    ASV_ID  = asv_ids,
    Kingdom = taxa[, "Kingdom"],
    Phylum  = taxa[, "Phylum"],
    Class   = taxa[, "Class"],
    Order   = taxa[, "Order"],
    Family  = taxa[, "Family"],
    Genus   = taxa[, "Genus"],
    Species = if ("Species" %in% colnames(taxa)) taxa[, "Species"] else NA,
    stringsAsFactors = FALSE
  )

  # Write output
  dir.create(dirname(opt$output), recursive=TRUE, showWarnings=FALSE)
  write.table(tax_df, opt$output, sep="\t", row.names=FALSE, quote=FALSE)
  cat("Wrote:", opt$output, "\n")

  # Success
  status <- toJSON(list(
    status    = "success",
    asv_count = length(asv_ids)
  ), auto_unbox=TRUE)
  cat("\n", status, "\n", sep="")
  quit(status=0)

}, error = function(e) {
  status <- toJSON(list(
    status  = "error",
    message = conditionMessage(e)
  ), auto_unbox=TRUE)
  cat("\n", status, "\n", sep="")
  quit(status=1)
})
