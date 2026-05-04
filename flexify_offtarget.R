# =============================================================================
# Flexify — Automated Off-Target Checking via BLAST
# =============================================================================
#
# These functions extend flexify_core.R with automated off-target checking of
# candidate probe sequences against a reference transcriptome.
#
# BIOLOGICAL RATIONALE:
#
#   FUSION PROBES:
#   Each 50 bp probe consists of two 25 bp halves. Exactly one half spans the
#   fusion breakpoint (the "junction half") — this is the specificity-determining
#   element unique to the fusion transcript. The other half lies entirely within
#   one of the two wild-type gene sequences and is expected to match that gene
#   perfectly by design.
#
#   Off-target checking for fusion probes is applied only to the junction half.
#   The junction half is identified from the fusion_point_displacement column:
#     - displacement > 0: junction falls in the RIGHT half → check right half
#     - displacement < 0: junction falls in the LEFT half  → check left half
#
#   A fusion probe fails if the junction half has a BLAST hit to any transcript
#   with fewer than min_mismatches effective mismatches.
#
#   NON-FUSION PROBES:
#   Non-fusion probes target wild-type transcripts — GFP, CRISPR constructs,
#   or any gene of interest. Both halves are designed to bind the intended
#   target, so both halves are queried against the database independently.
#
#   Because both halves bind wild-type sequence, the intended target transcript
#   WILL appear as a low-mismatch BLAST hit by design. Off-target concern arises
#   when a half also closely matches one or more ADDITIONAL transcripts. A probe
#   half is therefore flagged only if it has more than one unique database subject
#   with fewer than min_mismatches effective mismatches — the first close hit is
#   the expected on-target; any further ones indicate off-target binding risk.
#
#   Note: if the reference transcriptome contains multiple isoforms of the target
#   gene, these will each appear as close hits, potentially causing false positives.
#   Users should account for this when interpreting results.
#
# REQUIREMENTS:
#   - BLAST+ command-line tools installed and on PATH
#     Install via conda:  conda install -c bioconda blast
#     Or download from:  https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
#   - A pre-built BLAST nucleotide database of the reference transcriptome
#     Build with: makeblastdb -in transcriptome.fa -dbtype nucl -out transcriptome_db
#
# USAGE:
#   probe_df_filtered <- run_offtarget_check(
#     probe_df  = output_df,
#     blast_db  = "/path/to/transcriptome_db",
#     min_mismatches = 5,
#     n_threads = 4
#   )
#
# =============================================================================

library(tidyverse)
library(stringr)

# -----------------------------------------------------------------------------
# check_blast_available()
# Checks that BLAST+ is installed and accessible on the system PATH.
# -----------------------------------------------------------------------------

check_blast_available <- function() {
  result <- suppressWarnings(system("blastn -version", intern = TRUE, ignore.stderr = TRUE))
  if (length(result) == 0 || inherits(result, "try-error")) {
    stop(
      "BLAST+ does not appear to be installed or is not on your PATH.\n",
      "Install from: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/\n",
      "Then ensure 'blastn' is accessible from the command line."
    )
  }
  message("BLAST+ found: ", result[1])
  invisible(TRUE)
}


# -----------------------------------------------------------------------------
# extract_probe_halves()
# Given a probe string (with optional * and | markers), extracts the raw
# left and right 25 bp sequences for BLAST querying.
# -----------------------------------------------------------------------------

extract_probe_halves <- function(probe_str) {
  # Strip annotation characters (* and |) to get raw sequence
  raw <- gsub("[*|]", "", probe_str)
  raw <- toupper(trimws(raw))

  if (nchar(raw) < 50) {
    warning("Probe sequence shorter than expected 50 bp: ", probe_str)
    return(list(left = NA_character_, right = NA_character_))
  }

  left  <- substr(raw, 1, 25)
  right <- substr(raw, 26, 50)

  list(left = left, right = right)
}


# -----------------------------------------------------------------------------
# write_fasta()
# Writes a named character vector of sequences to a temporary FASTA file.
# Returns the file path.
# -----------------------------------------------------------------------------

write_fasta <- function(seqs, prefix = "probe") {
  tmp <- tempfile(pattern = paste0(prefix, "_"), fileext = ".fa")
  lines <- character(length(seqs) * 2)
  for (i in seq_along(seqs)) {
    lines[(i - 1) * 2 + 1] <- paste0(">", names(seqs)[i])
    lines[(i - 1) * 2 + 2] <- seqs[i]
  }
  writeLines(lines, tmp)
  tmp
}


# -----------------------------------------------------------------------------
# run_blastn()
# Runs blastn for a FASTA query file against a database.
# Returns a data frame of hits with: qseqid, sseqid, mismatch, length, qlen
# -----------------------------------------------------------------------------

run_blastn <- function(query_fa, blast_db, n_threads = 1, word_size = 7) {

  out_tmp <- tempfile(fileext = ".tsv")

  cmd <- sprintf(
    paste(
      "blastn",
      "-query %s",
      "-db %s",
      "-out %s",
      "-outfmt '6 qseqid sseqid mismatch length qlen'",
      "-task blastn-short",    # optimised for short sequences
      "-word_size %d",
      "-num_threads %d",
      "-max_target_seqs 500",  # retrieve enough hits to check off-targets
      "-evalue 1000"           # permissive e-value for short probes
    ),
    query_fa, blast_db, out_tmp, word_size, n_threads
  )

  exit_code <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

  if (exit_code != 0) {
    stop("BLAST failed with exit code: ", exit_code,
         "\nCommand was: ", cmd)
  }

  # Parse results
  if (file.info(out_tmp)$size == 0) {
    return(data.frame(
      qseqid = character(), sseqid = character(),
      mismatch = integer(), length = integer(), qlen = integer(),
      stringsAsFactors = FALSE
    ))
  }

  hits <- read.table(out_tmp, sep = "\t", header = FALSE,
                     stringsAsFactors = FALSE, quote = "",
                     col.names = c("qseqid", "sseqid", "mismatch", "length", "qlen"))
  unlink(out_tmp)
  hits
}


# -----------------------------------------------------------------------------
# is_offtarget_hit()
# Given a BLAST hits data frame for the junction half of a probe, returns TRUE
# if any hit has fewer than min_mismatches EFFECTIVE mismatches.
#
# Effective mismatches = mismatches within the aligned region + unaligned query
# bases. Treating unaligned bases as mismatches is essential for the junction
# half: BLAST will find partial alignments to the two fusion partner genes
# (because each gene contributes sequence to the junction half), but these
# partial alignments leave many query bases unaligned. The effective mismatch
# count correctly reflects that the junction half as a whole does not match
# any single reference transcript.
#
# A true off-target hit — where the junction sequence closely matches an
# unrelated transcript across most of its length — will have few unaligned
# bases and therefore a low effective mismatch count.
# -----------------------------------------------------------------------------

is_offtarget_hit <- function(hits, min_mismatches = 5) {
  if (nrow(hits) == 0) return(FALSE)
  # Effective mismatches: aligned mismatches + unaligned query bases
  effective_mm <- hits$mismatch + (hits$qlen - hits$length)
  any(effective_mm < min_mismatches)
}


# -----------------------------------------------------------------------------
# run_offtarget_check()
# Main function. Takes the ranked probe data frame output by
# create_probes_from_arriba() and adds an off-target check column,
# then optionally filters to passing probes only.
#
# Only the junction half of each probe is BLASTed (see module header for
# rationale). The junction half is determined from fusion_point_displacement.
#
# Args:
#   probe_df       — output of create_probes_from_arriba()
#   blast_db       — path to BLAST nucleotide database (no extension)
#   min_mismatches — minimum effective mismatches required to tolerate a hit (default 5)
#                    effective mismatches = aligned mismatches + unaligned query bases
#   n_threads      — number of BLAST threads (default 1)
#   filter_fails   — if TRUE, remove probes that fail off-target check (default TRUE)
#   batch_size     — number of sequences to BLAST per batch (default 50)
#
# Returns:
#   probe_df with added column:
#     offtarget_pass — logical: TRUE if the junction half passed off-target check
# -----------------------------------------------------------------------------

run_offtarget_check <- function(probe_df,
                                 blast_db,
                                 min_mismatches = 5,
                                 n_threads      = 1,
                                 filter_fails   = TRUE,
                                 batch_size     = 50) {

  check_blast_available()

  # Expand ~ so that system() receives an absolute path
  blast_db <- path.expand(blast_db)

  if (!file.exists(paste0(blast_db, ".nhr")) &&
      !file.exists(paste0(blast_db, ".nin")) &&
      !file.exists(paste0(blast_db, ".nsq")) &&
      !file.exists(paste0(blast_db, ".njs"))) {
    stop(
      "BLAST database not found at: ", blast_db, "\n",
      "Build with: makeblastdb -in transcriptome.fa -dbtype nucl -out ", blast_db
    )
  }

  n <- nrow(probe_df)
  message("Running off-target BLAST check on ", n, " probes (junction half only)...")

  # Determine which half is the junction half for each probe.
  # Positive displacement = junction in right half; negative = junction in left half.
  junction_side <- ifelse(probe_df$fusion_point_displacement > 0, "right", "left")

  # Extract the junction half sequence for each probe
  halves_list   <- lapply(probe_df$probe, extract_probe_halves)
  junction_seqs <- mapply(function(halves, side) {
    if (side == "left") halves$left else halves$right
  }, halves_list, junction_side, SIMPLIFY = TRUE)

  names(junction_seqs) <- paste0("probe_", seq_len(n), "_junction")

  # Remove any NA sequences (probes that failed half-extraction)
  valid_idx     <- !is.na(junction_seqs)
  junction_seqs <- junction_seqs[valid_idx]

  # Run BLAST in batches
  batches <- split(seq_along(junction_seqs),
                   ceiling(seq_along(junction_seqs) / batch_size))

  all_hits <- lapply(batches, function(idx) {
    batch_seqs <- junction_seqs[idx]
    query_fa   <- write_fasta(batch_seqs, prefix = "batch")
    hits       <- run_blastn(query_fa, blast_db, n_threads = n_threads)
    unlink(query_fa)
    hits
  })

  all_hits <- do.call(rbind, all_hits)

  # Evaluate each probe
  junction_fail <- logical(n)

  for (i in seq_len(n)) {
    jname <- paste0("probe_", i, "_junction")
    jhits <- if (nrow(all_hits) > 0) all_hits[all_hits$qseqid == jname, ] else all_hits
    junction_fail[i] <- is_offtarget_hit(jhits, min_mismatches)
  }

  probe_df$offtarget_pass <- !junction_fail

  n_pass <- sum(probe_df$offtarget_pass, na.rm = TRUE)
  n_fail <- n - n_pass
  message(sprintf("Off-target check complete: %d passed, %d failed.", n_pass, n_fail))

  if (filter_fails) {
    probe_df <- probe_df[probe_df$offtarget_pass, ]
    message(sprintf("%d probes retained after filtering.", nrow(probe_df)))
  }

  probe_df
}


# -----------------------------------------------------------------------------
# is_offtarget_hit_nonfusion()
# Variant of is_offtarget_hit() for non-fusion probe halves.
#
# Because the intended target transcript always produces a close BLAST hit by
# design, this function flags a probe half only when MORE THAN ONE unique
# subject sequence has fewer than min_mismatches effective mismatches. The
# single expected on-target hit passes; additional close hits indicate that
# the probe half may bind an unintended transcript.
#
# Args:
#   hits           — BLAST hit data frame for one 25 bp probe half
#   min_mismatches — effective mismatch threshold (default 5)
# Returns:
#   logical TRUE if the half has off-target risk (> 1 unique close subject)
# -----------------------------------------------------------------------------

is_offtarget_hit_nonfusion <- function(hits, min_mismatches = 5) {
  if (nrow(hits) == 0) return(FALSE)           # no hits at all — pass
  effective_mm      <- hits$mismatch + (hits$qlen - hits$length)
  close_hits        <- hits[effective_mm < min_mismatches, , drop = FALSE]
  n_unique_subjects <- length(unique(close_hits$sseqid))
  # > 1 unique subject = additional off-target hit beyond the expected on-target
  n_unique_subjects > 1
}


# -----------------------------------------------------------------------------
# run_offtarget_check_nonfusion()
# Off-target BLAST check for non-fusion probes.
#
# Both the LHS half (bases 1–25) and the RHS half (bases 26–50) of every probe
# are queried against the database independently. Each half is evaluated with
# is_offtarget_hit_nonfusion(): it is flagged only if more than one unique
# transcript matches it closely (the on-target hit is expected and does not
# count as a failure).
#
# A probe passes overall only if BOTH halves pass.
#
# Args:
#   probe_df       — output of create_nonfusion_probes()
#   blast_db       — path to BLAST nucleotide database (no extension)
#   min_mismatches — effective mismatch threshold (default 5)
#   n_threads      — number of BLAST threads (default 1)
#   filter_fails   — if TRUE, remove failing probes from output (default TRUE)
#   batch_size     — sequences per BLAST batch (default 50)
#
# Returns:
#   probe_df with added columns:
#     lhs_offtarget_pass — logical: TRUE if LHS half passed
#     rhs_offtarget_pass — logical: TRUE if RHS half passed
#     offtarget_pass     — logical: TRUE if both halves passed
# -----------------------------------------------------------------------------

run_offtarget_check_nonfusion <- function(probe_df,
                                           blast_db,
                                           min_mismatches = 5,
                                           n_threads      = 1,
                                           filter_fails   = TRUE,
                                           batch_size     = 50) {

  check_blast_available()

  blast_db <- path.expand(blast_db)

  if (!file.exists(paste0(blast_db, ".nhr")) &&
      !file.exists(paste0(blast_db, ".nin")) &&
      !file.exists(paste0(blast_db, ".nsq")) &&
      !file.exists(paste0(blast_db, ".njs"))) {
    stop(
      "BLAST database not found at: ", blast_db, "\n",
      "Build with: makeblastdb -in transcriptome.fa -dbtype nucl -out ", blast_db
    )
  }

  n <- nrow(probe_df)
  message("Running off-target BLAST check on ", n,
          " non-fusion probes (LHS and RHS halves independently)...")

  # Extract both halves for every probe
  halves_list <- lapply(probe_df$probe, extract_probe_halves)
  lhs_seqs    <- sapply(halves_list, `[[`, "left")
  rhs_seqs    <- sapply(halves_list, `[[`, "right")

  names(lhs_seqs) <- paste0("probe_", seq_len(n), "_lhs")
  names(rhs_seqs) <- paste0("probe_", seq_len(n), "_rhs")

  # Combine into a single query set; keep track of valid (non-NA) sequences
  all_seqs  <- c(lhs_seqs, rhs_seqs)
  valid_idx <- !is.na(all_seqs)
  all_seqs  <- all_seqs[valid_idx]

  # Run BLAST in batches over all halves
  batches  <- split(seq_along(all_seqs),
                    ceiling(seq_along(all_seqs) / batch_size))

  all_hits <- do.call(rbind, lapply(batches, function(idx) {
    batch_seqs <- all_seqs[idx]
    query_fa   <- write_fasta(batch_seqs, prefix = "batch_nf")
    hits       <- run_blastn(query_fa, blast_db, n_threads = n_threads)
    unlink(query_fa)
    hits
  }))

  lhs_fail <- logical(n)
  rhs_fail <- logical(n)

  for (i in seq_len(n)) {
    lname      <- paste0("probe_", i, "_lhs")
    rname      <- paste0("probe_", i, "_rhs")
    lhits      <- if (nrow(all_hits) > 0) all_hits[all_hits$qseqid == lname, ] else all_hits
    rhits      <- if (nrow(all_hits) > 0) all_hits[all_hits$qseqid == rname, ] else all_hits
    lhs_fail[i] <- is_offtarget_hit_nonfusion(lhits, min_mismatches)
    rhs_fail[i] <- is_offtarget_hit_nonfusion(rhits, min_mismatches)
  }

  probe_df$lhs_offtarget_pass <- !lhs_fail
  probe_df$rhs_offtarget_pass <- !rhs_fail
  probe_df$offtarget_pass     <- !lhs_fail & !rhs_fail

  n_pass <- sum(probe_df$offtarget_pass, na.rm = TRUE)
  n_fail <- n - n_pass
  message(sprintf("Off-target check complete: %d passed, %d failed.", n_pass, n_fail))

  if (filter_fails) {
    probe_df <- probe_df[probe_df$offtarget_pass, ]
    message(sprintf("%d probes retained after filtering.", nrow(probe_df)))
  }

  probe_df
}
