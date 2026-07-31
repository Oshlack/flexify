# =============================================================================
# Flexify — Handle & Barcode Appending
# =============================================================================
#
# This module appends the 10x Genomics handle sequences and Probe Barcodes
# to user-selected candidate probes, producing the full LHS and RHS probe
# sequences ready for submission to an oligonucleotide synthesis provider.
#
# PROBE STRUCTURE:
#
#   LHS probe (Left-Hand Side):
#     [21 bp constant handle] + [25 bp probe left half]
#     CCTTGGCACCCGAGAATTCCA + NNNNNNNNNNNNNNNNNNNNNNNNN
#
#   RHS probe (Right-Hand Side):
#     /5Phos/ + [25 bp probe right half] + [16 bp linker] + NN + [8 bp barcode] + [12 bp tail]
#     /5Phos/ + NNNNNNNNNNNNNNNNNNNNNNNNN + ACGCGGTTAGCACGTA + NN + XXXXXXXX + CGGTCCTAGCAA
#
# BARCODE ASSIGNMENT:
#   Each RHS probe is assigned one of 16 Probe Barcodes (BC001–BC016).
#   The barcode used must match the Probe Barcode of the corresponding whole
#   transcriptome probe in the hybridisation mix for that sample.
#   In a multiplexed experiment, each sample/pool gets a unique barcode.
#
# INPUT FORMAT for finalise_probes() / finalise_probes_v2():
#   A data frame with columns:
#     - GENE1    : first fusion partner gene name
#     - GENE2    : second fusion partner gene name
#     - probe    : 50 bp probe sequence (may contain | and * annotation markers)
#     - Barcode  : integer 1–16 indicating the barcode pool assignment (v1 only)
#
# INPUT FORMAT for finalise_nonfusion_probes() / finalise_nonfusion_probes_v2():
#   A data frame with columns:
#     - GENE     : target gene name
#     - probe    : 50 bp probe sequence (may contain | annotation marker)
#     - Barcode  : integer 1–16 indicating the barcode pool assignment (v1 only)
#
# DEPENDENCY: stringr
# =============================================================================

library(stringr)


# =============================================================================
# BARCODE LOOKUP TABLES
# =============================================================================

#' Probe Barcode sequences BC001–BC016 as specified by 10x Genomics.
#'
#' Each barcode is an 8-mer unique sequence appended to the RHS probe.
#' The barcode index (1–16) corresponds to the sample/pool position.
#' For multiplexed experiments, each sample must use the same barcode
#' as the corresponding whole transcriptome probe in the hybridisation mix.
PROBE_BARCODES <- c(
  BC001 = "ACTTTAGG",
  BC002 = "AACGGGAA",
  BC003 = "AGTAGGCT",
  BC004 = "ATGTTGAC",
  BC005 = "ACAGACCT",
  BC006 = "ATCCCAAC",
  BC007 = "AAGTAGAG",
  BC008 = "AGCTGTGA",
  BC009 = "ACAGTCTG",
  BC010 = "AGTGAGTG",
  BC011 = "AGAGGCAA",
  BC012 = "ACTACTCA",
  BC013 = "ATACGTCA",
  BC014 = "ATCATGTG",
  BC015 = "AACGCCGA",
  BC016 = "ATTCGGTT"
)

#' Pool names corresponding to each barcode index (1–16).
#' These match the naming convention used by 10x Genomics for sample pools.
POOL_NAMES <- c(
  "poolOne",     "poolTwo",     "poolThree",    "poolFour",
  "poolFive",    "poolSix",     "poolSeven",    "poolEight",
  "poolNine",    "poolTen",     "poolEleven",   "poolTwelve",
  "poolThirteen","poolFourteen","poolFifteen",  "poolSixteen"
)


# =============================================================================
# SEQUENCE CLEANING
# =============================================================================

#' Strip Flexify annotation markers from a probe sequence.
#'
#' Flexify optionally inserts '*' (marks the fusion breakpoint) and '|'
#' (marks the LHS/RHS half boundary) into probe sequences for display.
#' These characters must be removed before handle sequences are appended,
#' as they would otherwise be included in the final synthesis-ready sequence.
#'
#' @param probe_str character — probe sequence, possibly containing | and/or *
#' @return character — clean 50 bp uppercase sequence
clean_probe_sequence <- function(probe_str) {
  raw <- gsub("[*|]", "", probe_str)
  raw <- toupper(trimws(raw))

  if (nchar(raw) != 50) {
    warning(
      "Expected 50 bp after stripping markers but got ", nchar(raw),
      " bp. Check probe: ", probe_str
    )
  }
  return(raw)
}


#' Generate the full Right-Hand Side (RHS) probe sequence for Visium FFPE / CytAssist.
#'
#' Visium uses the same LHS handle as Chromium Flex but replaces the barcode
#' system with a 30-nucleotide poly-A tail on the RHS probe. No barcode is
#' embedded; spatial barcoding is handled by the Visium slide capture probes.
#'
#' RHS structure:
#'   /5Phos/ + [bases 26-50] + AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA  (30 × A)
#'
#' @param probe_str character — 50 bp probe sequence (markers will be stripped)
#' @return character — full Visium RHS probe sequence ready for ordering
add_rhs_handle_visium <- function(probe_str) {
  raw      <- clean_probe_sequence(probe_str)
  rhs_half <- substr(raw, 26, 50)
  paste0("/5Phos/", rhs_half, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
}


# =============================================================================
# HANDLE APPENDING FUNCTIONS
# =============================================================================

#' Generate the full Left-Hand Side (LHS) probe sequence.
#'
#' Prepends the 10x Genomics constant LHS handle to the first 25 bp of the
#' probe sequence. The handle sequence is identical for all probes.
#'
#' LHS structure:
#'   CCTTGGCACCCGAGAATTCCA  [21 bp constant handle]
#'   + [bases 1–25 of probe]
#'
#' @param probe_str character — 50 bp probe sequence (markers will be stripped)
#' @return character — full LHS probe sequence ready for ordering
add_lhs_handle <- function(probe_str) {
  raw      <- clean_probe_sequence(probe_str)
  lhs_half <- substr(raw, 1, 25)
  paste0("CCTTGGCACCCGAGAATTCCA", lhs_half)
}

#' Generate the full Right-Hand Side (RHS) probe sequence.
#'
#' Flanks the last 25 bp of the probe with the required 10x Genomics handle
#' sequences and the specified Probe Barcode.
#'
#' RHS structure:
#'   /5Phos/            [5-prime phosphorylation — required for ligation]
#'   + [bases 26–50 of probe]
#'   + ACGCGGTTAGCACGTA [16 bp linker / Constant Sequence]
#'   + NN               [2 bp spacer]
#'   + [8 bp Probe Barcode, unique per sample pool]
#'   + CGGTCCTAGCAA     [12 bp constant tail]
#'
#' @param probe_str   character — 50 bp probe sequence (markers will be stripped)
#' @param barcode_seq character — 8 bp barcode sequence (from PROBE_BARCODES)
#' @return character — full RHS probe sequence ready for ordering
add_rhs_handle <- function(probe_str, barcode_seq) {
  raw      <- clean_probe_sequence(probe_str)
  rhs_half <- substr(raw, 26, 50)
  paste0("/5Phos/", rhs_half, "ACGCGGTTAGCACGTA", "NN", barcode_seq, "CGGTCCTAGCAA")
}

#' Generate the full Right-Hand Side (RHS) probe sequence for GEM-X Flex v2.
#'
#' In GEM-X Flex v2, the barcode is NOT embedded in the custom probe sequence;
#' barcoding is handled separately by the kit reagents. The RHS probe therefore
#' ends with a short constant tail only.
#'
#' Two tail options are supported:
#'   - multiplex   (standard v2): /5Phos/ + [bases 26-50] + CCCATATAAGAAA  (13 bp)
#'   - singleplex  (4-sample kit): /5Phos/ + [bases 26-50] + CGGTCCTAGCAA  (12 bp)
#'
#' The LHS probe structure is identical to v1.
#'
#' @param probe_str character -- 50 bp probe sequence (markers will be stripped)
#' @param rhs_mode  character -- "multiplex" (default) or "singleplex"
#' @return character -- full v2 RHS probe sequence ready for ordering
add_rhs_handle_v2 <- function(probe_str, rhs_mode = "multiplex") {
  raw      <- clean_probe_sequence(probe_str)
  rhs_half <- substr(raw, 26, 50)
  tail     <- if (rhs_mode == "singleplex") "CGGTCCTAGCAA" else "CCCATATAAGAAA"
  paste0("/5Phos/", rhs_half, tail)
}


# =============================================================================
# MAIN FINALISE FUNCTION
# =============================================================================

#' Generate synthesis-ready LHS and RHS probe sequences for selected probes.
#'
#' Takes a data frame of user-selected candidate probes (one per fusion) and
#' returns a data frame with the complete LHS and RHS sequences, formatted
#' for submission to a commercial oligonucleotide synthesis provider.
#'
#' The Barcode column can be specified as:
#'   - Integer 1–16 (barcode index)
#'   - Character string "BC001"–"BC016" (barcode name)
#'
#' @param selected_df data frame with columns:
#'   - GENE1   : character — first fusion partner gene name
#'   - GENE2   : character — second fusion partner gene name
#'   - probe   : character — 50 bp probe sequence (| and * markers are stripped)
#'   - Barcode : integer 1–16 or character "BC001"–"BC016" — barcode assignment
#'
#' @return data frame with columns:
#'   - Fusion      : "GENE1::GENE2"
#'   - Barcode_ID  : e.g. "BC001"
#'   - Pool_Name   : e.g. "poolOne"
#'   - Barcode_Seq : 8 bp barcode sequence
#'   - LHS_Probe   : full LHS oligonucleotide sequence for ordering
#'   - RHS_Probe   : full RHS oligonucleotide sequence for ordering
finalise_probes <- function(selected_df) {

  # Resolve barcode index — accept either integer (1–16) or "BC001"–"BC016" strings
  barcode_raw <- selected_df$Barcode
  if (is.character(barcode_raw)) {
    # Convert "BC001" style to integer index
    barcode_idx <- match(barcode_raw, names(PROBE_BARCODES))
    if (any(is.na(barcode_idx))) {
      bad <- barcode_raw[is.na(barcode_idx)]
      stop("Unrecognised barcode name(s): ", paste(bad, collapse = ", "),
           "\nValid names are BC001–BC016.")
    }
  } else {
    barcode_idx <- as.integer(barcode_raw)
    if (any(is.na(barcode_idx)) || any(barcode_idx < 1) || any(barcode_idx > 16)) {
      stop("Barcode column must contain integers 1–16 or strings BC001–BC016.\n",
           "Check for missing or out-of-range values.")
    }
  }

  barcode_seqs <- PROBE_BARCODES[barcode_idx]
  barcode_ids  <- names(PROBE_BARCODES)[barcode_idx]
  pool_names   <- POOL_NAMES[barcode_idx]

  # Generate LHS and RHS sequences for each selected probe
  lhs_probes <- sapply(selected_df$probe, add_lhs_handle, USE.NAMES = FALSE)
  rhs_probes <- mapply(add_rhs_handle, selected_df$probe, barcode_seqs, USE.NAMES = FALSE)

  data.frame(
    Fusion      = paste(selected_df$GENE1, selected_df$GENE2, sep = "::"),
    Barcode_ID  = barcode_ids,
    Pool_Name   = pool_names,
    Barcode_Seq = as.character(barcode_seqs),
    LHS_Probe   = lhs_probes,
    RHS_Probe   = rhs_probes,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# GEM-X FLEX V2 FINALISE FUNCTION
# =============================================================================

#' Generate synthesis-ready LHS and RHS probe sequences for GEM-X Flex v2.
#'
#' In GEM-X Flex v2, barcoding is supplied by the kit reagents and is NOT
#' embedded in the custom probe sequence. This function therefore requires no
#' Barcode column in the input; only GENE1, GENE2, and probe are needed.
#'
#' @param selected_df data frame with columns:
#'   - GENE1 : character -- first fusion partner gene name
#'   - GENE2 : character -- second fusion partner gene name
#'   - probe : character -- 50 bp probe sequence (| and * markers are stripped)
#'
#' @param rhs_mode character -- RHS tail configuration:
#'   - "multiplex"  (default): appends CCCATATAAGAAA (13 bp, standard v2)
#'   - "singleplex"           : appends CGGTCCTAGCAA  (12 bp, 4-sample kit)
#'
#' @return data frame with columns:
#'   - Fusion    : "GENE1::GENE2"
#'   - RHS_Mode  : "multiplex" or "singleplex"
#'   - LHS_Probe : full LHS oligonucleotide sequence for ordering
#'   - RHS_Probe : full RHS oligonucleotide sequence for ordering
finalise_probes_v2 <- function(selected_df, rhs_mode = "multiplex") {

  if (!rhs_mode %in% c("multiplex", "singleplex")) {
    stop("rhs_mode must be 'multiplex' or 'singleplex'. Got: ", rhs_mode)
  }

  required <- c("GENE1", "GENE2", "probe")
  missing  <- setdiff(required, colnames(selected_df))
  if (length(missing) > 0) {
    stop("selected_df is missing required column(s): ", paste(missing, collapse = ", "))
  }

  lhs_probes <- sapply(selected_df$probe, add_lhs_handle,    USE.NAMES = FALSE)
  rhs_probes <- sapply(selected_df$probe, add_rhs_handle_v2,
                       rhs_mode = rhs_mode, USE.NAMES = FALSE)

  data.frame(
    Fusion    = paste(selected_df$GENE1, selected_df$GENE2, sep = "::"),
    RHS_Mode  = rhs_mode,
    LHS_Probe = lhs_probes,
    RHS_Probe = rhs_probes,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# NON-FUSION FINALISE FUNCTIONS
# =============================================================================

#' Generate synthesis-ready LHS and RHS sequences for selected non-fusion probes.
#'
#' Identical handle structure to the fusion probe finalise functions; the only
#' difference is that the gene label comes from a single GENE column rather
#' than GENE1/GENE2.
#'
#' @param selected_df data frame with columns:
#'   - GENE    : character — target gene name
#'   - probe   : character — 50 bp probe sequence (| markers are stripped)
#'   - Barcode : integer 1–16 or character "BC001"–"BC016" — barcode assignment
#'
#' @return data frame with columns:
#'   - Gene        : target gene name
#'   - Barcode_ID  : e.g. "BC001"
#'   - Pool_Name   : e.g. "poolOne"
#'   - Barcode_Seq : 8 bp barcode sequence
#'   - LHS_Probe   : full LHS oligonucleotide sequence for ordering
#'   - RHS_Probe   : full RHS oligonucleotide sequence for ordering
finalise_nonfusion_probes <- function(selected_df) {

  required <- c("GENE", "probe", "Barcode")
  missing  <- setdiff(required, colnames(selected_df))
  if (length(missing) > 0) {
    stop("selected_df is missing required column(s): ", paste(missing, collapse = ", "))
  }

  barcode_raw <- selected_df$Barcode
  if (is.character(barcode_raw)) {
    barcode_idx <- match(barcode_raw, names(PROBE_BARCODES))
    if (any(is.na(barcode_idx))) {
      bad <- barcode_raw[is.na(barcode_idx)]
      stop("Unrecognised barcode name(s): ", paste(bad, collapse = ", "),
           "\nValid names are BC001–BC016.")
    }
  } else {
    barcode_idx <- as.integer(barcode_raw)
    if (any(is.na(barcode_idx)) || any(barcode_idx < 1) || any(barcode_idx > 16)) {
      stop("Barcode column must contain integers 1–16 or strings BC001–BC016.")
    }
  }

  barcode_seqs <- PROBE_BARCODES[barcode_idx]
  barcode_ids  <- names(PROBE_BARCODES)[barcode_idx]
  pool_names   <- POOL_NAMES[barcode_idx]

  lhs_probes <- sapply(selected_df$probe, add_lhs_handle, USE.NAMES = FALSE)
  rhs_probes <- mapply(add_rhs_handle, selected_df$probe, barcode_seqs, USE.NAMES = FALSE)

  data.frame(
    Gene        = selected_df$GENE,
    Barcode_ID  = barcode_ids,
    Pool_Name   = pool_names,
    Barcode_Seq = as.character(barcode_seqs),
    LHS_Probe   = lhs_probes,
    RHS_Probe   = rhs_probes,
    stringsAsFactors = FALSE
  )
}


#' Generate synthesis-ready LHS and RHS sequences for non-fusion probes — GEM-X Flex v2.
#'
#' No barcode is embedded in the probe sequence; barcoding is handled by the kit.
#'
#' @param selected_df data frame with columns:
#'   - GENE  : character — target gene name
#'   - probe : character — 50 bp probe sequence (| markers are stripped)
#'
#' @param rhs_mode character — "multiplex" (default) or "singleplex"
#'
#' @return data frame with columns:
#'   - Gene      : target gene name
#'   - RHS_Mode  : "multiplex" or "singleplex"
#'   - LHS_Probe : full LHS oligonucleotide sequence for ordering
#'   - RHS_Probe : full RHS oligonucleotide sequence for ordering
finalise_nonfusion_probes_v2 <- function(selected_df, rhs_mode = "multiplex") {

  if (!rhs_mode %in% c("multiplex", "singleplex")) {
    stop("rhs_mode must be 'multiplex' or 'singleplex'. Got: ", rhs_mode)
  }

  required <- c("GENE", "probe")
  missing  <- setdiff(required, colnames(selected_df))
  if (length(missing) > 0) {
    stop("selected_df is missing required column(s): ", paste(missing, collapse = ", "))
  }

  lhs_probes <- sapply(selected_df$probe, add_lhs_handle, USE.NAMES = FALSE)
  rhs_probes <- sapply(selected_df$probe, add_rhs_handle_v2,
                       rhs_mode = rhs_mode, USE.NAMES = FALSE)

  data.frame(
    Gene      = selected_df$GENE,
    RHS_Mode  = rhs_mode,
    LHS_Probe = lhs_probes,
    RHS_Probe = rhs_probes,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# VISIUM FINALISE FUNCTIONS
# =============================================================================

#' Generate synthesis-ready LHS and RHS probe sequences for Visium FFPE / CytAssist.
#'
#' Visium probes share the same LHS handle as Chromium Flex but use a 30-nucleotide
#' poly-A tail on the RHS instead of the Flex barcode system. Spatial barcoding
#' is performed by the Visium slide capture probes, not by the custom probe sequence.
#'
#' @param selected_df data frame with columns:
#'   - GENE1 : character — first fusion partner gene name
#'   - GENE2 : character — second fusion partner gene name
#'   - probe : character — 50 bp probe sequence (| and * markers are stripped)
#'
#' @return data frame with columns:
#'   - Fusion    : "GENE1::GENE2"
#'   - LHS_Probe : full LHS oligonucleotide sequence for ordering
#'   - RHS_Probe : full RHS oligonucleotide sequence for ordering
finalise_probes_visium <- function(selected_df) {

  required <- c("GENE1", "GENE2", "probe")
  missing  <- setdiff(required, colnames(selected_df))
  if (length(missing) > 0) {
    stop("selected_df is missing required column(s): ", paste(missing, collapse = ", "))
  }

  lhs_probes <- sapply(selected_df$probe, add_lhs_handle,       USE.NAMES = FALSE)
  rhs_probes <- sapply(selected_df$probe, add_rhs_handle_visium, USE.NAMES = FALSE)

  data.frame(
    Fusion    = paste(selected_df$GENE1, selected_df$GENE2, sep = "::"),
    LHS_Probe = lhs_probes,
    RHS_Probe = rhs_probes,
    stringsAsFactors = FALSE
  )
}


#' Generate synthesis-ready LHS and RHS sequences for non-fusion probes — Visium FFPE / CytAssist.
#'
#' Same poly-A RHS tail as the fusion Visium function; gene label comes from
#' a single GENE column rather than GENE1/GENE2.
#'
#' @param selected_df data frame with columns:
#'   - GENE  : character — target gene name
#'   - probe : character — 50 bp probe sequence (| markers are stripped)
#'
#' @return data frame with columns:
#'   - Gene      : target gene name
#'   - LHS_Probe : full LHS oligonucleotide sequence for ordering
#'   - RHS_Probe : full RHS oligonucleotide sequence for ordering
finalise_nonfusion_probes_visium <- function(selected_df) {

  required <- c("GENE", "probe")
  missing  <- setdiff(required, colnames(selected_df))
  if (length(missing) > 0) {
    stop("selected_df is missing required column(s): ", paste(missing, collapse = ", "))
  }

  lhs_probes <- sapply(selected_df$probe, add_lhs_handle,       USE.NAMES = FALSE)
  rhs_probes <- sapply(selected_df$probe, add_rhs_handle_visium, USE.NAMES = FALSE)

  data.frame(
    Gene      = selected_df$GENE,
    LHS_Probe = lhs_probes,
    RHS_Probe = rhs_probes,
    stringsAsFactors = FALSE
  )
}
