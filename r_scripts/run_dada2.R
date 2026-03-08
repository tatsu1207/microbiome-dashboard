#!/usr/bin/env Rscript
# ============================================================================
# MicrobiomeDash — DADA2 Denoising Pipeline
#
# Supports paired-end and single-end modes.
# Outputs: asv_table.tsv, rep_seqs.fasta, track_reads.tsv
# Prints JSON status to stdout on completion.
# ============================================================================

suppressPackageStartupMessages({
  library(dada2)
  library(optparse)
  library(jsonlite)
})

# Flushing log helper — ensures output reaches Python before a potential crash
log_msg <- function(...) {
  cat(...)
  flush.console()
}

# Decode signal number to name for crash diagnostics
signal_name <- function(sig) {
  names <- c("1"="SIGHUP","2"="SIGINT","4"="SIGILL","6"="SIGABRT",
             "7"="SIGBUS","8"="SIGFPE","9"="SIGKILL","11"="SIGSEGV",
             "13"="SIGPIPE","15"="SIGTERM")
  n <- as.character(sig)
  if (n %in% names(names)) names[[n]] else paste0("signal ", sig)
}

# ── CLI Arguments ─────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--input_dir",   type="character", help="Directory with FASTQ files"),
  make_option("--output_dir",  type="character", help="Output directory"),
  make_option("--mode",        type="character", default="paired", help="paired, single, or longread"),
  make_option("--trim_left_f", type="integer",   default=0,   help="Trim N bases from 5' of forward reads"),
  make_option("--trim_left_r", type="integer",   default=0,   help="Trim N bases from 5' of reverse reads"),
  make_option("--trunc_len_f", type="integer",   default=0,   help="Truncate forward reads (0=no truncation)"),
  make_option("--trunc_len_r", type="integer",   default=0,   help="Truncate reverse reads (0=no truncation)"),
  make_option("--min_overlap", type="integer",   default=12,  help="Min overlap for merging (PE only)"),
  make_option("--threads",     type="integer",   default=1,   help="Number of threads"),
  # Long-read options
  make_option("--fwd_primer",  type="character", default=NULL, help="Forward primer sequence (longread mode)"),
  make_option("--rev_primer",  type="character", default=NULL, help="Reverse primer sequence (longread mode)"),
  make_option("--min_len",     type="integer",   default=1000, help="Min read length after filtering (longread)"),
  make_option("--max_len",     type="integer",   default=1600, help="Max read length after filtering (longread)"),
  make_option("--max_ee",      type="integer",   default=10,   help="Max expected errors (longread)"),
  make_option("--band_size",   type="integer",   default=32,   help="DADA2 BAND_SIZE (longread)"),
  make_option("--platform",    type="character", default="pacbio", help="Platform: pacbio or nanopore (longread)"),
  make_option("--pool",        type="character", default="FALSE",  help="Pooling for dada(): FALSE, TRUE, or pseudo"),
  make_option("--skip_filter", action="store_true", default=FALSE, help="Skip filterAndTrim, use existing filtered/ files")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input_dir) || is.null(opt$output_dir)) {
  stop("--input_dir and --output_dir are required")
}

is_longread <- (opt$mode == "longread")
is_paired <- (opt$mode == "paired")

# Parse --pool: "FALSE" → FALSE, "TRUE" → TRUE, "pseudo" → "pseudo"
pool_opt <- if (opt$pool == "TRUE") TRUE else if (opt$pool == "pseudo") "pseudo" else FALSE

# ── Long-Read Pipeline (PacBio HiFi) ─────────────────────────────────────────

if (is_longread) {
tryCatch({

  dir.create(opt$output_dir, recursive=TRUE, showWarnings=FALSE)

  # Discover FASTQ files (all files in input directory, single-end only)
  fnFs <- sort(list.files(opt$input_dir, pattern="\\.(fastq|fq)(\\.gz)?$", full.names=TRUE))
  if (length(fnFs) == 0) stop("No FASTQ files found")
  sample.names <- sub("\\.(fastq|fq)(\\.gz)?$", "", basename(fnFs))
  sample.names <- sub("_R1[_.]?.*$|_1$", "", sample.names)
  log_msg("Found", length(fnFs), "long-read samples\n")

  use_mt <- ifelse(opt$threads > 1, opt$threads, FALSE)

  # ── Step 1: Remove Primers ──────────────────────────────────────────────

  if (!is.null(opt$fwd_primer) && !is.null(opt$rev_primer)) {
    log_msg("Removing primers (orient=TRUE)...\n")
    log_msg("  Forward:", opt$fwd_primer, "\n")
    log_msg("  Reverse:", opt$rev_primer, "\n")

    nop_dir <- file.path(opt$output_dir, "noprimers")
    dir.create(nop_dir, recursive=TRUE, showWarnings=FALSE)
    nopFs <- file.path(nop_dir, basename(fnFs))

    prim_out <- removePrimers(fnFs, nopFs,
                              primer.fwd=opt$fwd_primer,
                              primer.rev=dada2:::rc(opt$rev_primer),
                              orient=TRUE,
                              verbose=TRUE)

    # Track how many reads had primers removed
    prim_reads_in  <- prim_out[, "reads.in"]
    prim_reads_out <- prim_out[, "reads.out"]

    # Keep only samples with reads passing
    keep_prim <- prim_reads_out > 0
    if (sum(keep_prim) == 0) stop("No reads had primers removed for any sample")
    log_msg(sum(keep_prim), "of", length(keep_prim), "samples had primers removed\n")

    fnFs <- nopFs[keep_prim]
    sample.names <- sample.names[keep_prim]
    prim_reads_in  <- prim_reads_in[keep_prim]
    prim_reads_out <- prim_reads_out[keep_prim]
  } else {
    log_msg("No primers specified, skipping primer removal\n")
    prim_reads_in  <- NULL
    prim_reads_out <- NULL
  }

  # ── Step 2: Filter and Trim (length filtering, no truncation) ──────────

  filt_dir <- file.path(opt$output_dir, "filtered")
  dir.create(filt_dir, recursive=TRUE, showWarnings=FALSE)
  filtFs <- file.path(filt_dir, paste0(sample.names, "_filt.fastq.gz"))

  log_msg("Filtering (minLen=", opt$min_len, ", maxLen=", opt$max_len,
      ", maxEE=", opt$max_ee, ")...\n", sep="")
  filt_out <- filterAndTrim(
    fnFs, filtFs,
    minLen=opt$min_len, maxLen=opt$max_len,
    maxN=0, maxEE=opt$max_ee,
    rm.phix=FALSE, compress=TRUE,
    multithread=use_mt
  )

  keep <- filt_out[, "reads.out"] > 0
  if (sum(keep) == 0) stop("No reads passed filtering for any sample")
  log_msg(sum(keep), "of", length(keep), "samples passed filtering\n")

  sample.names <- sample.names[keep]
  filtFs <- filtFs[keep]
  if (!is.null(prim_reads_in)) {
    prim_reads_in  <- prim_reads_in[keep]
    prim_reads_out <- prim_reads_out[keep]
  }

  # ── Step 3: Learn Error Rates (PacBio error model) ─────────────────────

  if (opt$platform == "pacbio") {
    log_msg("Learning error rates (PacBioErrfun, BAND_SIZE=", opt$band_size, ")...\n", sep="")
    errF <- learnErrors(filtFs, errorEstimationFunction=PacBioErrfun,
                         BAND_SIZE=opt$band_size,
                         multithread=use_mt, randomize=TRUE)
  } else {
    log_msg("Learning error rates (loessErrfun, BAND_SIZE=", opt$band_size, ")...\n", sep="")
    errF <- learnErrors(filtFs, errorEstimationFunction=loessErrfun,
                         BAND_SIZE=opt$band_size,
                         multithread=use_mt, randomize=TRUE)
  }

  # ── Step 4: Dereplicate + Denoise ──────────────────────────────────────

  log_msg("Dereplicating (", length(filtFs), " samples)...\n", sep="")
  t0 <- proc.time()[["elapsed"]]
  derep <- derepFastq(filtFs)
  names(derep) <- sample.names
  log_msg("Dereplication done in ", round(proc.time()[["elapsed"]] - t0, 1), "s\n", sep="")

  log_msg("Denoising (BAND_SIZE=", opt$band_size, ", HOMOPOLYMER_GAP_PENALTY=-1, pool=", as.character(pool_opt), ")...\n", sep="")
  t0 <- proc.time()[["elapsed"]]
  dadaFs <- dada(derep, err=errF, multithread=use_mt,
                  BAND_SIZE=opt$band_size,
                  HOMOPOLYMER_GAP_PENALTY=-1,
                  pool=pool_opt)
  log_msg("Denoising done in ", round(proc.time()[["elapsed"]] - t0, 1), "s\n", sep="")

  # ── Step 5: Sequence Table + Chimera Removal ───────────────────────────

  seqtab <- makeSequenceTable(dadaFs)
  log_msg("Sequence table:", nrow(seqtab), "samples,", ncol(seqtab), "ASVs\n")

  log_msg("Removing chimeras...\n")
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                       multithread=use_mt)

  asv_count <- ncol(seqtab.nochim)
  sample_count <- nrow(seqtab.nochim)
  log_msg("After chimera removal:", sample_count, "samples,", asv_count, "ASVs\n")

  if (asv_count == 0) stop("No ASVs survived chimera removal")

  # ── Build Read Tracking Table ──────────────────────────────────────────

  if (length(sample.names) == 1) {
    dadaFs <- list(dadaFs)
  }

  get_n <- function(x) sum(getUniques(x))

  if (!is.null(prim_reads_in)) {
    track <- cbind(
      input=prim_reads_in,
      primers_removed=prim_reads_out,
      filtered=filt_out[keep, "reads.out"],
      denoised=sapply(dadaFs, get_n),
      nonchim=rowSums(seqtab.nochim)
    )
  } else {
    track <- cbind(
      input=filt_out[keep, "reads.in"],
      filtered=filt_out[keep, "reads.out"],
      denoised=sapply(dadaFs, get_n),
      nonchim=rowSums(seqtab.nochim)
    )
  }
  rownames(track) <- sample.names

  # ── Write Outputs ──────────────────────────────────────────────────────

  asv_ids <- paste0("ASV_", seq_len(asv_count))
  seqs <- colnames(seqtab.nochim)

  asv_table <- as.data.frame(t(seqtab.nochim))
  colnames(asv_table) <- sample.names
  asv_table <- cbind(ASV_ID = asv_ids, sequence = seqs, asv_table)

  asv_path <- file.path(opt$output_dir, "asv_table.tsv")
  write.table(asv_table, asv_path, sep="\t", row.names=FALSE, quote=FALSE)
  log_msg("Wrote:", asv_path, "\n")

  track_path <- file.path(opt$output_dir, "track_reads.tsv")
  track_df <- cbind(sample = rownames(track), as.data.frame(track))
  write.table(track_df, track_path, sep="\t", row.names=FALSE, quote=FALSE)
  log_msg("Wrote:", track_path, "\n")

  # ── Success ────────────────────────────────────────────────────────────

  status <- toJSON(list(
    status       = "success",
    asv_count    = asv_count,
    sample_count = sample_count
  ), auto_unbox=TRUE)
  log_msg("\n", status, "\n", sep="")
  quit(status=0)

}, error = function(e) {
  status <- toJSON(list(
    status  = "error",
    message = conditionMessage(e)
  ), auto_unbox=TRUE)
  log_msg("\n", status, "\n", sep="")
  quit(status=1)
})
}

# ── Short-Read Pipeline (Illumina PE/SE) ──────────────────────────────────────

tryCatch({

  dir.create(opt$output_dir, recursive=TRUE, showWarnings=FALSE)

  # Discover FASTQ files
  if (is_paired) {
    fnFs <- sort(list.files(opt$input_dir, pattern="_R1([_.].*)?\\.(fastq|fq)(\\.gz)?$|_1\\.(fastq|fq)(\\.gz)?$", full.names=TRUE))
    fnRs <- sort(list.files(opt$input_dir, pattern="_R2([_.].*)?\\.(fastq|fq)(\\.gz)?$|_2\\.(fastq|fq)(\\.gz)?$", full.names=TRUE))

    if (length(fnFs) == 0) stop("No forward (R1) FASTQ files found")
    if (length(fnFs) != length(fnRs)) stop("Mismatched number of forward and reverse files")

    # Extract sample names from forward reads
    sample.names <- sub("_R1([_.].*)?\\.(fastq|fq)(\\.gz)?$|_1\\.(fastq|fq)(\\.gz)?$", "", basename(fnFs))
    log_msg("Found", length(fnFs), "paired-end samples\n")
  } else {
    fnFs <- sort(list.files(opt$input_dir, pattern="\\.(fastq|fq)(\\.gz)?$", full.names=TRUE))
    if (length(fnFs) == 0) stop("No FASTQ files found")
    sample.names <- sub("\\.(fastq|fq)(\\.gz)?$", "", basename(fnFs))
    # Also strip _R1 suffixes if present in single-end mode
    sample.names <- sub("_R1[_.]?.*$|_1$", "", sample.names)
    log_msg("Found", length(fnFs), "single-end samples\n")
  }

  # ── Step 1: Filter and Trim ──────────────────────────────────────────────

  filt_dir <- file.path(opt$output_dir, "filtered")
  dir.create(filt_dir, recursive=TRUE, showWarnings=FALSE)

  # Use the --threads value for DADA2's multithread parameter.
  # When threads > 1 on Linux, DADA2 uses mclapply (forked parallelism).
  # (On Windows/WSL2 with fork issues, pass --threads 1 from the caller.)
  use_mt <- ifelse(opt$threads > 1, opt$threads, FALSE)

  # Auto-detect existing filtered files — skip filterAndTrim if they exist
  existing_filtFs <- sort(list.files(filt_dir, pattern="_F_filt\\.fastq\\.gz$", full.names=TRUE))
  skip_filter <- opt$skip_filter || length(existing_filtFs) > 0

  if (skip_filter) {
    filtFs <- existing_filtFs
    if (length(filtFs) == 0) stop("--skip_filter: no filtered forward files found in ", filt_dir)
    sample.names <- sub("_F_filt\\.fastq\\.gz$", "", basename(filtFs))
    if (is_paired) {
      filtRs <- sort(list.files(filt_dir, pattern="_R_filt\\.fastq\\.gz$", full.names=TRUE))
      if (length(filtRs) != length(filtFs)) stop("Mismatched F/R filtered files in ", filt_dir)
    }
    log_msg("Filtered files exist, skipping filterAndTrim: ", length(filtFs), " samples\n")

    # Build a synthetic filt_out matrix for the tracking table
    filt_out <- matrix(0L, nrow=length(filtFs), ncol=2)
    colnames(filt_out) <- c("reads.in", "reads.out")
    rownames(filt_out) <- basename(filtFs)
    keep <- rep(TRUE, length(filtFs))

  } else {
    filtFs <- file.path(filt_dir, paste0(sample.names, "_F_filt.fastq.gz"))

    if (is_paired) {
      filtRs <- file.path(filt_dir, paste0(sample.names, "_R_filt.fastq.gz"))
      log_msg("Filtering and trimming (paired-end)...\n")
      log_msg("  truncLen=c(", opt$trunc_len_f, ",", opt$trunc_len_r, ") trimLeft=c(",
          opt$trim_left_f, ",", opt$trim_left_r, ") truncQ=0 maxEE=c(5,5)\n", sep="")
      filt_out <- filterAndTrim(
        fnFs, filtFs, fnRs, filtRs,
        trimLeft  = c(opt$trim_left_f, opt$trim_left_r),
        truncLen  = c(opt$trunc_len_f, opt$trunc_len_r),
        maxN = 0, maxEE = c(5, 5), truncQ = 0,
        rm.phix = TRUE, compress = TRUE,
        multithread = use_mt
      )
      log_msg("  Filter results:\n")
      print(filt_out)
    } else {
      log_msg("Filtering and trimming (single-end)...\n")
      filt_out <- filterAndTrim(
        fnFs, filtFs,
        trimLeft  = opt$trim_left_f,
        truncLen  = opt$trunc_len_f,
        maxN = 0, maxEE = 5, truncQ = 0,
        rm.phix = TRUE, compress = TRUE,
        multithread = use_mt
      )
    }

    # Remove samples that had 0 reads pass filtering
    keep <- filt_out[, "reads.out"] > 0
    if (sum(keep) == 0) stop("No reads passed filtering for any sample")

    log_msg(sum(keep), "of", length(keep), "samples passed filtering\n")
    sample.names <- sample.names[keep]
    filtFs <- filtFs[keep]
    if (is_paired) filtRs <- filtRs[keep]
  }

  # ── Step 2: Learn Error Rates ────────────────────────────────────────────

  log_msg("Learning error rates (forward, randomize=TRUE)...\n")
  errF <- learnErrors(filtFs, multithread=use_mt, randomize=TRUE)

  if (is_paired) {
    log_msg("Learning error rates (reverse, randomize=TRUE)...\n")
    errR <- learnErrors(filtRs, multithread=use_mt, randomize=TRUE)
  }

  # ── Step 3: Dereplicate ────────────────────────────────────────────────

  log_msg("Dereplicating forward reads (", length(filtFs), " samples)...\n", sep="")
  t0 <- proc.time()[["elapsed"]]
  derepFs <- derepFastq(filtFs)
  names(derepFs) <- sample.names
  log_msg("Forward dereplication done in ", round(proc.time()[["elapsed"]] - t0, 1), "s\n", sep="")

  if (is_paired) {
    log_msg("Dereplicating reverse reads (", length(filtRs), " samples)...\n", sep="")
    t0 <- proc.time()[["elapsed"]]
    derepRs <- derepFastq(filtRs)
    names(derepRs) <- sample.names
    log_msg("Reverse dereplication done in ", round(proc.time()[["elapsed"]] - t0, 1), "s\n", sep="")
  }

  # ── Step 4: Denoise ──────────────────────────────────────────────────────

  log_msg("Denoising forward reads (pool=", as.character(pool_opt), ")...\n", sep="")
  t0 <- proc.time()[["elapsed"]]
  dadaFs <- dada(derepFs, err=errF, multithread=use_mt, pool=pool_opt)
  log_msg("Forward denoising done in ", round(proc.time()[["elapsed"]] - t0, 1), "s\n", sep="")

  if (is_paired) {
    log_msg("Denoising reverse reads (pool=", as.character(pool_opt), ")...\n", sep="")
    t0 <- proc.time()[["elapsed"]]
    dadaRs <- dada(derepRs, err=errR, multithread=use_mt, pool=pool_opt)
    log_msg("Reverse denoising done in ", round(proc.time()[["elapsed"]] - t0, 1), "s\n", sep="")
  }

  # ── Step 5: Merge (PE only) ──────────────────────────────────────────────

  if (is_paired) {
    log_msg("Merging paired reads...\n")
    t0 <- proc.time()[["elapsed"]]
    merged <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,
                         minOverlap=opt$min_overlap,
                         maxMismatch=1)
    seqtab <- makeSequenceTable(merged)
    log_msg("Merging done in ", round(proc.time()[["elapsed"]] - t0, 1), "s\n", sep="")
  } else {
    seqtab <- makeSequenceTable(dadaFs)
  }

  log_msg("Sequence table:", nrow(seqtab), "samples,", ncol(seqtab), "ASVs\n")

  # ── Step 6: Remove Chimeras ──────────────────────────────────────────────

  log_msg("Removing chimeras...\n")
  t0 <- proc.time()[["elapsed"]]
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                       multithread=use_mt)
  log_msg("Chimera removal done in ", round(proc.time()[["elapsed"]] - t0, 1), "s\n", sep="")

  asv_count <- ncol(seqtab.nochim)
  sample_count <- nrow(seqtab.nochim)
  log_msg("After chimera removal:", sample_count, "samples,", asv_count, "ASVs\n")

  if (asv_count == 0) stop("No ASVs survived chimera removal")

  # ── Build Read Tracking Table ────────────────────────────────────────────

  # When there's only 1 sample, dada()/mergePairs() return single objects
  # instead of lists.  Wrap them so sapply() works uniformly.
  if (length(sample.names) == 1) {
    dadaFs <- list(dadaFs)
    if (is_paired) {
      dadaRs <- list(dadaRs)
      merged <- list(merged)
    }
  }

  # Get denoised counts
  get_n <- function(x) sum(getUniques(x))

  if (is_paired) {
    track <- cbind(
      filt_out[keep, , drop=FALSE],
      sapply(dadaFs, get_n),
      sapply(dadaRs, get_n),
      sapply(merged, get_n),
      rowSums(seqtab.nochim)
    )
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  } else {
    track <- cbind(
      filt_out[keep, , drop=FALSE],
      sapply(dadaFs, get_n),
      rowSums(seqtab.nochim)
    )
    colnames(track) <- c("input", "filtered", "denoised", "nonchim")
  }

  rownames(track) <- sample.names

  # ── Write Outputs ────────────────────────────────────────────────────────

  # ASV table: rows = ASVs, columns = samples
  # First column is ASV_ID, second is the actual sequence, then one column per sample
  asv_ids <- paste0("ASV_", seq_len(asv_count))
  seqs <- colnames(seqtab.nochim)

  # Transpose: seqtab.nochim is samples x sequences → we want sequences x samples
  asv_table <- as.data.frame(t(seqtab.nochim))
  colnames(asv_table) <- sample.names
  asv_table <- cbind(ASV_ID = asv_ids, sequence = seqs, asv_table)

  asv_path <- file.path(opt$output_dir, "asv_table.tsv")
  write.table(asv_table, asv_path, sep="\t", row.names=FALSE, quote=FALSE)
  log_msg("Wrote:", asv_path, "\n")

  # Read tracking table
  track_path <- file.path(opt$output_dir, "track_reads.tsv")
  track_df <- cbind(sample = rownames(track), as.data.frame(track))
  write.table(track_df, track_path, sep="\t", row.names=FALSE, quote=FALSE)
  log_msg("Wrote:", track_path, "\n")

  # ── Success ──────────────────────────────────────────────────────────────

  status <- toJSON(list(
    status       = "success",
    asv_count    = asv_count,
    sample_count = sample_count
  ), auto_unbox=TRUE)
  log_msg("\n", status, "\n", sep="")
  quit(status=0)

}, error = function(e) {
  status <- toJSON(list(
    status  = "error",
    message = conditionMessage(e)
  ), auto_unbox=TRUE)
  log_msg("\n", status, "\n", sep="")
  quit(status=1)
})
