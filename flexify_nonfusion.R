# =============================================================================
# Flexify — Non-Fusion Probe Design
# =============================================================================
#
# This module designs candidate probe sequences for non-fusion (wild-type)
# transcripts and checks them for potential competition with the 10x Genomics
# standard Flex whole-transcriptome probe set.
#
# BIOLOGICAL RATIONALE:
#   For non-fusion targets, there is no junction to straddle. Instead, a
#   50 bp probe window is slid along the input transcript sequence and each
#   candidate is scored on the same quality metrics used for fusion probes:
#   GC content (each 25 bp half independently), ligation dinucleotide, and
#   homopolymer content.
#
#   The probe is still split 25|25 for LHS/RHS handle appending (identical
#   handle structure to fusion probes). There is no junction-position score.
#
# COMPETITION CHECK:
#   The 10x Genomics standard Flex probe set covers the whole transcriptome.
#   If a custom non-fusion probe sequence closely matches a standard probe
#   sequence, the two probes will compete for binding to the same mRNA
#   molecule, reducing the efficacy of the standard probe. The competition
#   check compares each half of the custom probe (25 bp LHS and 25 bp RHS
#   separately) against the corresponding halves of all standard probes.
#   A probe half is flagged if it matches any standard probe half with fewer
#   than max_mismatches mismatches (default: 2). A custom probe passes only
#   if BOTH halves pass.
#
# INPUT FORMAT for create_nonfusion_probes():
#   A data frame with columns:
#     - gene     : gene name
#     - sequence : mRNA transcript sequence (5'->3', at least 50 bp)
#
# DEPENDENCIES: tidyverse, stringr
# Requires flexify_core.R to be sourced first (uses check_GC, GC_rating,
# dinucleotide_check, check_for_homo_polymer, reverse_complement).
# =============================================================================

library(tidyverse)
library(stringr)


# =============================================================================
# TILING PROBE ENUMERATION
# =============================================================================

#' Design and rank all 50 bp candidate probes for a single non-fusion transcript.
#'
#' Slides a 50 bp window along the input mRNA sequence and scores every
#' position on GC content (each half independently), ligation dinucleotide,
#' and homopolymer content. The probe sequence is the reverse complement of
#' each mRNA window (antisense orientation, matching the standard Flex probes).
#'
#' There is no junction-position rating (unlike fusion probes).
#'
#' @param gene           character -- gene or target name
#' @param sequence       character -- mRNA sequence, 5'->3', >= 50 bp
#' @param PROBE_HALVES_FLAG logical -- mark probe half boundary with '|' (default TRUE)
#' @param MRNA_FLAG      logical -- include the mRNA target sequence column (default FALSE)
#' @return data frame of ranked probe candidates, or an empty data frame if
#'         the sequence is shorter than 50 bp
tile_sequence <- function(gene, sequence,
                          PROBE_HALVES_FLAG = TRUE,
                          MRNA_FLAG         = FALSE) {

  sequence <- toupper(trimws(sequence))
  n        <- nchar(sequence)

  if (n < 50) {
    warning("Sequence for '", gene, "' is shorter than 50 bp (", n, " bp). Skipping.")
    return(data.frame())
  }

  n_windows <- n - 49L

  # Generate probe sequences: RC of each 50 bp mRNA window
  probes    <- character(n_windows)
  positions <- integer(n_windows)
  for (i in seq_len(n_windows)) {
    probes[i]    <- reverse_complement(substr(sequence, i, i + 49L))
    positions[i] <- i
  }

  probe_df <- data.frame(
    probe            = probes,
    mRNA_position    = positions,
    stringsAsFactors = FALSE
  )

  # -- GC content of each half --
  first_half           <- lapply(probe_df$probe, str_sub, 1L, 25L)
  first_half_GC        <- unlist(lapply(first_half, check_GC))
  first_half_GC_rating <- unlist(lapply(first_half_GC, GC_rating))
  probe_df             <- cbind(probe_df, first_half_GC, first_half_GC_rating)

  second_half           <- lapply(probe_df$probe, str_sub, 26L, 50L)
  second_half_GC        <- unlist(lapply(second_half, check_GC))
  second_half_GC_rating <- unlist(lapply(second_half_GC, GC_rating))
  probe_df              <- cbind(probe_df, second_half_GC, second_half_GC_rating)

  # -- Ligation dinucleotide (positions 25-26) --
  di_results                  <- lapply(probe_df$probe, dinucleotide_check)
  probe_df$Dinucleotide        <- sapply(di_results, `[[`, 1)
  probe_df$Dinucleotide_Status <- sapply(di_results, `[[`, 2)
  probe_df$dinucleotide_rating <- as.numeric(sapply(di_results, `[[`, 3))

  # -- Homopolymer content --
  hp_results                 <- lapply(probe_df$probe, check_for_homo_polymer)
  probe_df$Homopolymer_Flag   <- sapply(hp_results, `[[`, 1)
  probe_df$homopolymer_rating <- as.numeric(sapply(hp_results, `[[`, 2))

  # -- Composite score: product of all rating columns --
  rating_cols  <- probe_df %>%
    select_if(grepl("rating", names(.))) %>%
    mutate_all(as.numeric)
  probe_df$Score <- apply(rating_cols, 1, prod)

  # -- Filter, rank, and tidy --
  probe_df <- probe_df %>%
    select_if(!grepl("rating", names(.))) %>%
    arrange(desc(Score)) %>%
    filter(Score > 0) %>%
    mutate(Ranking = seq_len(n()))

  # Prepend gene name
  probe_df <- cbind(data.frame(GENE = gene, stringsAsFactors = FALSE), probe_df)

  # Optional half-boundary marker
  if (PROBE_HALVES_FLAG) {
    probe_df$probe <- paste0(
      substr(probe_df$probe, 1L, 25L), "|",
      substr(probe_df$probe, 26L, 50L)
    )
  }

  # Optional mRNA target column
  if (MRNA_FLAG) {
    probe_df$mrna <- sapply(probe_df$probe, function(p) {
      reverse_complement(gsub("[|*]", "", p))
    })
  }

  probe_df
}


#' Design and rank probes for all entries in a non-fusion input data frame.
#'
#' @param input_df       data frame with columns 'gene' and 'sequence'
#' @param PROBE_HALVES_FLAG logical -- mark half boundary with '|' (default TRUE)
#' @param MRNA_FLAG      logical -- include mRNA target sequence (default FALSE)
#' @return combined data frame of all ranked probe candidates
create_nonfusion_probes <- function(input_df,
                                    PROBE_HALVES_FLAG = TRUE,
                                    MRNA_FLAG         = FALSE) {

  colnames(input_df) <- tolower(trimws(colnames(input_df)))

  required <- c("gene", "sequence")
  missing  <- setdiff(required, colnames(input_df))
  if (length(missing) > 0) {
    stop("Input CSV is missing required columns: ", paste(missing, collapse = ", "),
         "\nExpected columns: gene, sequence")
  }

  input_df$gene     <- trimws(input_df$gene)
  input_df$sequence <- toupper(trimws(input_df$sequence))

  output_list <- lapply(seq_len(nrow(input_df)), function(i) {
    tile_sequence(input_df$gene[i], input_df$sequence[i],
                  PROBE_HALVES_FLAG = PROBE_HALVES_FLAG,
                  MRNA_FLAG         = MRNA_FLAG)
  })

  output_list <- Filter(function(x) !is.null(x) && nrow(x) > 0, output_list)

  if (length(output_list) == 0) {
    warning("No probes could be designed. Check that all sequences are at least 50 bp.")
    return(data.frame())
  }

  do.call(rbind, output_list)
}


# =============================================================================
# FLEX PROBE COMPETITION CHECK
# =============================================================================

#' Load a 10x Genomics standard Flex probeset CSV.
#'
#' Skips the comment lines (starting with '#') at the top of the file and
#' reads the tabular data. Filters to included probes only if an 'included'
#' column is present.
#'
#' @param filepath character -- path to the 10x probeset CSV file
#' @return data frame with at minimum a 'probe_seq' column
load_flex_probeset <- function(filepath) {

  lines      <- readLines(filepath, warn = FALSE)
  data_start <- which(!startsWith(lines, "#"))[1]

  if (is.na(data_start)) {
    stop("No non-comment lines found in probeset file: ", filepath)
  }

  raw <- read.csv(
    text             = paste(lines[data_start:length(lines)], collapse = "\n"),
    stringsAsFactors = FALSE
  )

  if (!"probe_seq" %in% colnames(raw)) {
    stop("Probeset file must contain a 'probe_seq' column. Found: ",
         paste(colnames(raw), collapse = ", "))
  }

  # Filter to included probes only
  if ("included" %in% colnames(raw)) {
    n_before <- nrow(raw)
    raw      <- raw[which(raw$included == TRUE | raw$included == "TRUE"), ]
    message(sprintf("Probeset: %d of %d probes are marked 'included'.", nrow(raw), n_before))
  }

  raw$probe_seq <- toupper(trimws(raw$probe_seq))
  raw <- raw[nchar(raw$probe_seq) == 50, ]  # enforce 50 bp
  message("Probeset loaded: ", nrow(raw), " probes.")
  raw
}


#' Check custom non-fusion probes for competition with standard Flex probes.
#'
#' For each custom probe, the LHS half (bases 1-25) and RHS half (bases 26-50)
#' are compared independently against the corresponding halves of all standard
#' probes. Competition is reported if either half matches any standard probe's
#' corresponding half with fewer than max_mismatches mismatches (Hamming
#' distance over the aligned 25 bp).
#'
#' The comparison uses a vectorised character-matrix approach for efficiency:
#' all standard probe halves are pre-split into a matrix once, so each custom
#' probe requires only two matrix comparisons.
#'
#' @param probe_df     data frame -- output of create_nonfusion_probes()
#' @param probeset_df  data frame -- output of load_flex_probeset()
#' @param max_mismatches integer -- minimum mismatches required to NOT flag as
#'                       competing (default 2: probes with <= 2 mismatches are
#'                       flagged; i.e., >= 96% identical over 25 bp)
#' @return probe_df with four additional columns:
#'   - flex_lhs_min_mm  : minimum mismatches of LHS half vs all standard LHS halves
#'   - flex_rhs_min_mm  : minimum mismatches of RHS half vs all standard RHS halves
#'   - flex_competition : logical TRUE if either half has <= max_mismatches to any standard probe
#'   - flex_pass        : logical TRUE if probe passes (no competition detected)
#' Check fusion probes for competition with standard Flex probes (non-junction half only).
#'
#' For fusion probes, only the wild-type (non-junction) half is compared against
#' the standard probeset. The junction-spanning half binds a sequence that only
#' exists in the fusion transcript and will not be present in the standard probe
#' set, so it does not need checking.
#'
#' Which half is the non-junction half is determined by fusion_point_displacement:
#'   > 0  : junction falls in the right half  → check the LEFT half (LHS, bases 1-25)
#'   < 0  : junction falls in the left half   → check the RIGHT half (RHS, bases 26-50)
#'   == 0 : junction exactly at centre        → check BOTH halves
#'
#' @param probe_df     data frame -- output of create_probes_from_arriba() — must
#'                     contain a 'fusion_point_displacement' column
#' @param probeset_df  data frame -- output of load_flex_probeset()
#' @param max_mismatches integer -- minimum mismatches required to NOT flag as
#'                       competing (default 2: probes with <= 2 mismatches are
#'                       flagged; i.e., >= 96% identical over 25 bp)
#' @return probe_df with additional columns:
#'   - flex_nonjunction_mm  : minimum mismatches of the non-junction half vs all
#'                            standard probe halves (the relevant half only)
#'   - flex_competition     : logical TRUE if the non-junction half has <= max_mismatches
#'                            to any standard probe half
#'   - flex_pass            : logical TRUE if probe passes (no competition detected)
check_flex_competition_fusion <- function(probe_df, probeset_df, max_mismatches = 2L) {

  if (!"fusion_point_displacement" %in% colnames(probe_df)) {
    stop("probe_df must contain a 'fusion_point_displacement' column. ",
         "Ensure the input comes from create_probes_from_arriba().")
  }

  std_seqs <- unique(probeset_df$probe_seq)
  std_seqs <- std_seqs[nchar(std_seqs) == 50L]
  n_std    <- length(std_seqs)

  message(sprintf("Fusion competition check: %d custom probes vs %d standard probes...",
                  nrow(probe_df), n_std))

  # Pre-split standard probe halves into character matrices (once each)
  std_lhs_mat <- do.call(rbind, strsplit(substr(std_seqs, 1L, 25L), ""))   # n_std x 25
  std_rhs_mat <- do.call(rbind, strsplit(substr(std_seqs, 26L, 50L), ""))  # n_std x 25

  # Strip | and * markers from custom probe sequences
  clean_probes <- gsub("[|*]", "", probe_df$probe)

  # Vectorised Hamming distance: compare one 25bp query against all rows of a matrix
  min_hamming <- function(query_25, ref_mat) {
    if (nchar(query_25) != 25L) return(NA_integer_)
    q_mat <- matrix(strsplit(query_25, "")[[1]], nrow = nrow(ref_mat), ncol = 25L, byrow = TRUE)
    min(rowSums(ref_mat != q_mat))
  }

  nonjunction_mm <- integer(nrow(probe_df))

  for (i in seq_len(nrow(probe_df))) {
    p   <- clean_probes[i]
    fpd <- probe_df$fusion_point_displacement[i]

    if (fpd > 0) {
      # Junction is in the right half → left half is pure wild-type → check LHS
      nonjunction_mm[i] <- min_hamming(substr(p, 1L, 25L), std_lhs_mat)
    } else if (fpd < 0) {
      # Junction is in the left half → right half is pure wild-type → check RHS
      nonjunction_mm[i] <- min_hamming(substr(p, 26L, 50L), std_rhs_mat)
    } else {
      # Junction exactly at centre → check both halves, take minimum
      mm_lhs <- min_hamming(substr(p, 1L, 25L),  std_lhs_mat)
      mm_rhs <- min_hamming(substr(p, 26L, 50L), std_rhs_mat)
      nonjunction_mm[i] <- min(mm_lhs, mm_rhs)
    }
  }

  probe_df$flex_nonjunction_mm <- nonjunction_mm
  probe_df$flex_competition     <- (nonjunction_mm <= max_mismatches)
  probe_df$flex_pass            <- !probe_df$flex_competition

  n_pass <- sum(probe_df$flex_pass, na.rm = TRUE)
  n_fail <- nrow(probe_df) - n_pass
  message(sprintf("Fusion competition check complete: %d passed, %d flagged.", n_pass, n_fail))

  probe_df
}


check_flex_competition <- function(probe_df, probeset_df, max_mismatches = 2L) {

  std_seqs <- unique(probeset_df$probe_seq)
  std_seqs <- std_seqs[nchar(std_seqs) == 50L]
  n_std    <- length(std_seqs)

  message(sprintf("Competition check: %d custom probes vs %d standard probes...",
                  nrow(probe_df), n_std))

  # Pre-split standard probe halves into character matrices (once)
  std_lhs_mat <- do.call(rbind, strsplit(substr(std_seqs, 1L, 25L), ""))   # n_std x 25
  std_rhs_mat <- do.call(rbind, strsplit(substr(std_seqs, 26L, 50L), ""))  # n_std x 25

  # Strip | and * markers from custom probe sequences
  clean_probes <- gsub("[|*]", "", probe_df$probe)

  # Vectorised Hamming distance: compare one 25bp query against all rows of a matrix
  min_hamming <- function(query_25, ref_mat) {
    if (nchar(query_25) != 25L) return(NA_integer_)
    q_mat <- matrix(strsplit(query_25, "")[[1]], nrow = nrow(ref_mat), ncol = 25L, byrow = TRUE)
    min(rowSums(ref_mat != q_mat))
  }

  lhs_mm <- integer(nrow(probe_df))
  rhs_mm <- integer(nrow(probe_df))

  for (i in seq_len(nrow(probe_df))) {
    p         <- clean_probes[i]
    lhs_mm[i] <- min_hamming(substr(p, 1L, 25L),  std_lhs_mat)
    rhs_mm[i] <- min_hamming(substr(p, 26L, 50L), std_rhs_mat)
  }

  probe_df$flex_lhs_min_mm  <- lhs_mm
  probe_df$flex_rhs_min_mm  <- rhs_mm
  probe_df$flex_competition  <- (lhs_mm <= max_mismatches) | (rhs_mm <= max_mismatches)
  probe_df$flex_pass         <- !probe_df$flex_competition

  n_pass <- sum(probe_df$flex_pass, na.rm = TRUE)
  n_fail <- nrow(probe_df) - n_pass
  message(sprintf("Competition check complete: %d passed, %d flagged.", n_pass, n_fail))

  probe_df
}
