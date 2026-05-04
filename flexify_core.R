# =============================================================================
# Flexify — Core Probe Design Functions
# =============================================================================
#
# This module contains all core functions for designing candidate probe
# sequences targeting fusion transcripts for the 10x Genomics Flex platform.
#
# WORKFLOW OVERVIEW:
#   1. Input: a data frame with columns gene1, gene2, gene1_transcript,
#             gene2_transcript (derived from Arriba fusion calls).
#   2. For each fusion, enumerate all 50 bp probe sequences spanning the
#      breakpoint across a range of junction positions.
#   3. Score each candidate probe on GC content, junction position,
#      ligation dinucleotide, and homopolymer content.
#   4. Return probes ranked by composite score; zero-scoring probes excluded.
#
# MAIN ENTRY POINT:
#   create_probes_from_arriba(input_df, ...)
#
# DEPENDENCIES: tidyverse, stringr
# =============================================================================

library(tidyverse)
library(stringr)


# =============================================================================
# SEQUENCE UTILITIES
# =============================================================================

#' Reverse complement of a DNA sequence.
#'
#' Reverses the input string and applies Watson-Crick base complementation
#' (A<->T, G<->C). Used to convert Arriba transcript sequences (which are
#' given in 5'->3' orientation relative to the fusion partner) into the
#' correct orientation for probe construction.
#'
#' @param input_string character — DNA sequence in uppercase (A/T/G/C)
#' @return character — reverse complement sequence
reverse_complement <- function(input_string) {
  reversed <- strsplit(input_string, NULL) %>%
    lapply(rev) %>%
    sapply(paste, collapse = "")
  updated <- chartr("ATGC", "TACG", reversed)
  return(updated)
}


# =============================================================================
# PROBE ENUMERATION
# =============================================================================

#' Enumerate all candidate 50 bp probe sequences spanning a fusion junction.
#'
#' For a given pair of transcript sequences flanking the breakpoint, generates
#' all possible 50 bp probes by sliding the junction position across the probe
#' from left to right, subject to a minimum contribution from each gene
#' (RESTRAINT_CONST bases minimum per half).
#'
#' The offset parameter encodes how far the fusion junction is displaced from
#' the probe centre:
#'   - Negative offset: junction falls in the left half (more bases from gene2)
#'   - Positive offset: junction falls in the right half (more bases from gene1)
#'   - Offsets within ±RESTRAINT_CONST are excluded (too few bases from one gene)
#'
#' @param left_string      character — left transcript sequence (gene2, RC'd)
#' @param right_string     character — right transcript sequence (gene1, RC'd)
#' @param RESTRAINT_CONST  integer — minimum bases required from each gene (default 5)
#' @return data frame with columns: probe (character), fusion_point_displacement (numeric)
produce_possible_probe_df <- function(left_string, right_string, RESTRAINT_CONST = 5) {
  list_with_probes_and_offset <- list()

  # Negative offsets: junction displaced into left half
  for (offset in (-25 + RESTRAINT_CONST):(-1 - RESTRAINT_CONST)) {
    probe <- str_c(
      str_sub(left_string, -25 - offset, -1),
      str_sub(right_string, 1, 25 - offset)
    )
    # Skip probes whose window overlaps a gap region (non-nucleotide characters
    # such as '.' or '_' from Arriba's gap markers indicate missing sequence).
    if (!grepl("[^ATGCNatgcn]", probe)) {
      list_with_probes_and_offset <- rbind(list_with_probes_and_offset, c(probe, offset))
    }
  }

  # Positive offsets: junction displaced into right half
  for (offset in (1 + RESTRAINT_CONST):(25 - RESTRAINT_CONST)) {
    probe <- str_c(
      str_sub(left_string, -25 - offset, -1),
      str_sub(right_string, 1, 25 - offset)
    )
    if (!grepl("[^ATGCNatgcn]", probe)) {
      list_with_probes_and_offset <- rbind(list_with_probes_and_offset, c(probe, offset))
    }
  }

  # Return an empty data frame if all candidates were excluded (e.g. gaps too
  # close to the breakpoint to allow any valid probe window).
  if (length(list_with_probes_and_offset) == 0) {
    return(data.frame(probe = character(), fusion_point_displacement = numeric(),
                      stringsAsFactors = FALSE))
  }

  probe_df <- as.data.frame(list_with_probes_and_offset)
  colnames(probe_df) <- c("probe", "fusion_point_displacement")
  probe_df <- transform(probe_df,
                        fusion_point_displacement = as.numeric(probe_df$fusion_point_displacement))
  return(probe_df)
}


# =============================================================================
# SCORING FUNCTIONS
# =============================================================================

#' Score the fusion junction position within a probe using a Gaussian function.
#'
#' Probes where the junction falls near the centre of either half score highest.
#' The Gaussian is centred at position 12.5 bp (centre of a 25 bp half), with
#' a standard deviation of 3 bp and a maximum score of 5.
#'
#' This reflects the expectation that probes with centrally positioned junctions
#' support more reliable hybridisation and ligation in the 10x Flex assay.
#'
#' @param x numeric — absolute fusion point displacement from probe centre
#' @return numeric — location score (0–5)
fusion_location_rating <- function(x) {
  x <- abs(x)
  rating <- 5 * exp(-(x - 12.5)^2 / (2 * (3^2)))
  return(rating)
}

#' Apply fusion location rating to all probes in a data frame.
#'
#' @param probe_df data frame with column fusion_point_displacement
#' @return probe_df with additional column location_rating
fusion_location_rating_df <- function(probe_df) {
  location_rating <- unlist(lapply(probe_df$fusion_point_displacement, fusion_location_rating))
  probe_df <- probe_df %>% cbind(location_rating)
  return(probe_df)
}

#' Compute GC content (%) for a 25 bp probe half.
#'
#' @param probe_half character — 25 bp nucleotide sequence
#' @return numeric — GC content as a percentage (0–100)
check_GC <- function(probe_half) {
  GC_content <- str_count(probe_half, pattern = c("C", "G"))
  GC_content <- sum(GC_content) / 25 * 100
  return(GC_content)
}

#' Convert GC content (%) to a tiered score.
#'
#' Based on 10x Genomics probe design guidelines:
#'   - Outside 44–72%: score 0 (probe excluded from output)
#'   - 44–49%: score 2
#'   - 50–54%: score 5 (optimal)
#'   - 55–59%: score 4
#'   - 60–64%: score 3
#'   - 65–69%: score 2
#'   - 70–72%: score 1
#'
#' @param x numeric — GC content percentage
#' @return integer — GC rating score (0, 1, 2, 3, 4, or 5)
GC_rating <- function(x) {
  if (x < 44 || x > 72) return(0)
  if (x >= 44 && x < 50) return(2)
  if (x >= 50 && x < 55) return(5)
  if (x >= 55 && x < 60) return(4)
  if (x >= 60 && x < 65) return(3)
  if (x >= 65 && x < 70) return(2)
  if (x >= 70 && x <= 72) return(1)
}

#' Check the ligation junction dinucleotide (positions 25–26 of the probe).
#'
#' In the 10x Flex design, the probe halves are ligated at the junction between
#' bases 25 and 26. Certain dinucleotides at this position are more common in
#' the standard transcriptome probe set and are preferred for reliable ligation.
#'
#' Preferred dinucleotides: AT, CA, CT, TA, TC, TG, TT
#'
#' @param probe character — full 50 bp probe sequence
#' @return list of: dinucleotide (character), status ("OK"/"Warning"), rating (numeric)
dinucleotide_check <- function(probe) {
  better_dinucleotides <- c("AT", "CA", "CT", "TA", "TC", "TG", "TT")
  probe_di <- str_sub(probe, 25, 26)
  dinucleotide_flag   <- "OK"
  dinucleotide_rating <- 3

  if (!(probe_di %in% better_dinucleotides)) {
    dinucleotide_flag   <- "Warning"
    dinucleotide_rating <- 1
  }
  return(list(probe_di, dinucleotide_flag, dinucleotide_rating))
}

#' Check for homopolymer runs within the probe sequence.
#'
#' Homopolymer runs of 4 or more identical consecutive bases can impair
#' hybridisation efficiency and are penalised. Runs of 4 receive a mild
#' penalty; runs of 5 or more receive a stronger penalty.
#'
#' @param probe character — full 50 bp probe sequence
#' @return list of: flag (character), rating (numeric)
check_for_homo_polymer <- function(probe) {
  homopolymer_flag   <- "OK"
  homopolymer_rating <- 2

  letter_runs        <- rle(unlist(str_split(probe, "")))
  longest_homopolymer <- max(letter_runs$lengths)

  if (longest_homopolymer == 4) {
    homopolymer_rating <- 1.5
    homopolymer_flag   <- "Warning (length 4)"
  }
  if (longest_homopolymer >= 5) {
    homopolymer_rating <- 1
    homopolymer_flag   <- "Warning (length ≥5)"
  }
  return(list(homopolymer_flag, homopolymer_rating))
}

#' Determine which probe half contains the fusion junction.
#'
#' A negative displacement means the junction falls in the left half;
#' positive means it falls in the right half.
#'
#' @param displacement numeric — fusion_point_displacement value
#' @return character — "left" or "right"
check_rhs <- function(displacement) {
  if (displacement < 0) return("left") else return("right")
}

#' Apply a mild penalty when the fusion junction falls in the left probe half.
#'
#' When PRIORITISE_RHS_FLAG is TRUE, probes where the junction falls in the
#' right half are preferred (they can be detected without ligation). Probes
#' with a left-half junction receive a 0.7x score multiplier.
#'
#' @param displacement numeric — fusion_point_displacement value
#' @return numeric — 0.7 if junction is in left half, 1.0 otherwise
check_lhs_rating <- function(displacement) {
  if (displacement < 0) return(0.7) else return(1)
}


# =============================================================================
# PROBE ANNOTATION HELPERS
# =============================================================================

#' Add gene name columns to the probe data frame.
#'
#' @param gene1    character — first fusion partner gene name
#' @param gene2    character — second fusion partner gene name
#' @param probe_df data frame of candidate probes
#' @return probe_df with GENE1 and GENE2 columns prepended
add_gene_names <- function(gene1, gene2, probe_df) {
  reps <- dim(probe_df)[1]
  probe_df$probe <- as.character(probe_df$probe)
  new_df <- data.frame(
    "GENE1" = rep(gene1, reps),
    "GENE2" = rep(gene2, reps),
    "probe" = unlist(probe_df$probe)
  )
  new_df <- left_join(new_df, probe_df)
  return(new_df)
}

#' Mark the fusion breakpoint within the probe sequence with an asterisk (*).
#'
#' Inserts a '*' at the exact position where the fusion junction occurs,
#' for visual identification. The position is determined by
#' fusion_point_displacement.
#'
#' @param probe_df single-row data frame with probe and fusion_point_displacement
#' @return character — probe sequence with '*' inserted at the junction
add_asterix <- function(probe_df) {
  # Use [[ ]] instead of $ because this function is called via apply(..., 1, ...)
  # which passes each row as a named character vector, not a data frame.
  disp      <- as.numeric(probe_df[["fusion_point_displacement"]])
  lhs <- paste0("^([A-Z]{", 25 + disp, "})([A-Z]+)$")
  rhs <- paste0("\\1", "*", "\\2")
  new_probe <- gsub(lhs, rhs, probe_df[["probe"]])
  return(new_probe)
}

#' Mark the boundary between the LHS and RHS probe halves with a pipe (|).
#'
#' Inserts a '|' at position 25|26 (or 26|27 when the junction has been
#' shifted by the asterisk insertion), to visually separate the two 25 bp
#' halves that will become the LHS and RHS probes.
#'
#' @param probe_df    single-row data frame with probe and fusion_point_displacement
#' @param ASTERIX_FLAG logical — whether an asterisk has already been inserted
#' @return character — probe sequence with '|' inserted at the half boundary
mark_halves <- function(probe_df, ASTERIX_FLAG) {
  # Use [[ ]] instead of $ because this function is called via apply(..., 1, ...)
  # which passes each row as a named character vector, not a data frame.
  disp  <- as.numeric(probe_df[["fusion_point_displacement"]])
  probe <- probe_df[["probe"]]
  if (ASTERIX_FLAG) {
    if (disp < 0) {
      return(str_c(str_sub(probe, 1, 26), "|", str_sub(probe, 27, -1)))
    }
  }
  return(str_c(str_sub(probe, 1, 25), "|", str_sub(probe, 26, -1)))
}

#' Compute the mRNA target sequence for a probe.
#'
#' Strips annotation markers (| and *) from the probe sequence and returns
#' the reverse complement, which corresponds to the mRNA sequence that the
#' probe pair will hybridise to.
#'
#' @param probe character — probe sequence (may contain | and * markers)
#' @return character — mRNA target sequence
mrna <- function(probe) {
  probe <- gsub("[[:punct:] ]+", "", probe)
  return(reverse_complement(probe))
}


# =============================================================================
# MAIN PROBE DESIGN FUNCTIONS
# =============================================================================

#' Design and rank all candidate probes for a single fusion.
#'
#' Enumerates all possible 50 bp probe sequences spanning the fusion junction,
#' scores each on GC content, junction position, ligation dinucleotide, and
#' homopolymer content, filters out zero-scoring probes, and returns the
#' remainder ranked by composite score (product of all rating components).
#'
#' @param gene1              character — first fusion partner gene name
#' @param gene2              character — second fusion partner gene name
#' @param fusion_half_1      character — Arriba transcript sequence for gene1
#' @param fusion_half_2      character — Arriba transcript sequence for gene2
#' @param RESTRAINT_CONST    integer — minimum bases per probe half from each gene (default 5)
#' @param PRIORITISE_RHS_FLAG logical — apply 0.7x penalty to left-half junction probes (default FALSE)
#' @param ASTERIX_FLAG       logical — mark fusion breakpoint with '*' (default TRUE)
#' @param PROBE_HALVES_FLAG  logical — mark probe half boundary with '|' (default TRUE)
#' @param MRNA_FLAG          logical — include mRNA target sequence column (default FALSE)
#' @return data frame of ranked candidate probes with quality metrics
create_probe <- function(gene1, gene2, fusion_half_1, fusion_half_2,
                         RESTRAINT_CONST    = 5,
                         PRIORITISE_RHS_FLAG = FALSE,
                         ASTERIX_FLAG       = TRUE,
                         PROBE_HALVES_FLAG  = TRUE,
                         MRNA_FLAG          = FALSE) {

  # Convert Arriba sequences to probe orientation (reverse complement)
  fusion_half_1_now_right <- reverse_complement(fusion_half_1)
  fusion_half_2_now_left  <- reverse_complement(fusion_half_2)

  # Enumerate all candidate probes across the allowed junction offset range.
  # If gaps are too close to the breakpoint, no valid candidates will be returned.
  probe_df <- produce_possible_probe_df(fusion_half_2_now_left, fusion_half_1_now_right, RESTRAINT_CONST)

  if (nrow(probe_df) == 0) {
    warning("No valid probe candidates for ", gene1, "::", gene2,
            " — gap regions may be too close to the breakpoint.")
    return(data.frame())
  }

  probe_df <- fusion_location_rating_df(probe_df)

  # Score GC content for each probe half independently
  first_half          <- lapply(probe_df$probe, str_sub, 1, 25)
  first_half_GC       <- unlist(lapply(first_half, check_GC))
  first_half_GC_rating <- unlist(lapply(first_half_GC, GC_rating))
  probe_df <- probe_df %>% cbind(first_half_GC) %>% cbind(first_half_GC_rating)

  second_half          <- lapply(probe_df$probe, str_sub, 26, 50)
  second_half_GC       <- unlist(lapply(second_half, check_GC))
  second_half_GC_rating <- unlist(lapply(second_half_GC, GC_rating))
  probe_df <- probe_df %>% cbind(second_half_GC) %>% cbind(second_half_GC_rating)

  # Score ligation junction dinucleotide (positions 25–26).
  # sapply with [[ indexing extracts each list element into a plain atomic vector,
  # avoiding the list-column problem that do.call(rbind, lapply(...)) would produce.
  di_results <- lapply(probe_df$probe, dinucleotide_check)
  probe_df$Dinucleotide        <- sapply(di_results, `[[`, 1)
  probe_df$Dinucleotide_Status <- sapply(di_results, `[[`, 2)
  probe_df$dinucleotide_rating <- as.numeric(sapply(di_results, `[[`, 3))

  # Score homopolymer content (same approach)
  hp_results <- lapply(probe_df$probe, check_for_homo_polymer)
  probe_df$Homopolymer_Flag   <- sapply(hp_results, `[[`, 1)
  probe_df$homopolymer_rating <- as.numeric(sapply(hp_results, `[[`, 2))

  # Record which half contains the fusion junction
  fusion_point_side <- unlist(lapply(probe_df$fusion_point_displacement, check_rhs))
  probe_df <- probe_df %>% cbind(fusion_point_side)

  # Optional: penalise probes where the junction falls in the left half
  if (PRIORITISE_RHS_FLAG) {
    lhs_rating <- unlist(lapply(probe_df$fusion_point_displacement, check_lhs_rating))
    probe_df <- probe_df %>% cbind(lhs_rating)
  }

  # Compute composite score as the product of all individual rating components.
  # Any zero-scoring property causes the composite to be zero, excluding that probe.
  rating_cols <- probe_df %>%
    select_if(grepl("rating", names(.))) %>%
    mutate_all(as.numeric)
  probe_df$Score <- apply(rating_cols, 1, prod)

  # Drop individual rating sub-columns; sort by score descending; filter zeros;
  # assign consecutive integer rankings starting from 1.
  probe_df_2 <- probe_df %>%
    select_if(!grepl("rating", names(.))) %>%
    arrange(desc(Score)) %>%
    filter(Score > 0) %>%
    mutate(Ranking = seq_len(n()))

  # Prepend gene name columns
  probe_df_2 <- add_gene_names(gene1, gene2, probe_df_2)

  # Optionally annotate the probe sequence with breakpoint and half-boundary markers.
  # mapply() operates directly on columns as vectors, avoiding the character-matrix
  # coercion that apply(..., 1, FUN) would introduce.
  if (ASTERIX_FLAG) {
    probe_df_2$probe <- mapply(function(probe, disp) {
      lhs_pattern <- paste0("^([A-Z]{", 25 + disp, "})([A-Z]+)$")
      gsub(lhs_pattern, "\\1*\\2", probe)
    }, probe_df_2$probe, probe_df_2$fusion_point_displacement)
  }
  if (PROBE_HALVES_FLAG) {
    probe_df_2$probe <- mapply(function(probe, disp) {
      if (ASTERIX_FLAG && disp < 0) {
        str_c(str_sub(probe, 1, 26), "|", str_sub(probe, 27, -1))
      } else {
        str_c(str_sub(probe, 1, 25), "|", str_sub(probe, 26, -1))
      }
    }, probe_df_2$probe, probe_df_2$fusion_point_displacement)
  }

  # Optionally include the mRNA target sequence.
  # sapply (not lapply) so the column is a plain character vector, not a list column.
  if (MRNA_FLAG) {
    probe_df_2$mrna <- sapply(probe_df_2$probe, mrna)
  }

  return(probe_df_2)
}

#' Wrapper to call create_probe() on a single row of the input data frame.
#'
#' Designed for use with apply() over rows of the Arriba input data frame.
#'
#' @param df a single-row data frame with columns gene1, gene2, gene1_transcript, gene2_transcript
#' @param ... additional arguments passed to create_probe()
#' @return data frame of ranked candidate probes for this fusion
data_frame_create_probe <- function(df, RESTRAINT_CONST = 5, PRIORITISE_RHS_FLAG = FALSE,
                                    ASTERIX_FLAG = TRUE, PROBE_HALVES_FLAG = TRUE, MRNA_FLAG = FALSE) {
  # Use [[ ]] instead of $ because apply(..., 1, ...) passes rows as named character vectors.
  create_probe(df[["gene1"]], df[["gene2"]], df[["gene1_transcript"]], df[["gene2_transcript"]],
               RESTRAINT_CONST, PRIORITISE_RHS_FLAG, ASTERIX_FLAG, PROBE_HALVES_FLAG, MRNA_FLAG)
}

#' Design and rank probes for all fusions in an Arriba input data frame.
#'
#' Applies create_probe() to each row of the input data frame and returns
#' a combined data frame of all candidate probes across all fusions, each
#' annotated with their gene names and ranked within their fusion.
#'
#' @param input_df data frame with columns: gene1, gene2, gene1_transcript, gene2_transcript
#' @param RESTRAINT_CONST    integer — minimum bases per probe half from each gene (default 5)
#' @param PRIORITISE_RHS_FLAG logical — penalise left-half junction probes (default FALSE)
#' @param ASTERIX_FLAG       logical — mark fusion breakpoint with '*' (default TRUE)
#' @param PROBE_HALVES_FLAG  logical — mark probe half boundary with '|' (default TRUE)
#' @param MRNA_FLAG          logical — include mRNA target sequence column (default FALSE)
#' @return data frame of all ranked candidate probes
create_probes_from_arriba <- function(input_df,
                                      RESTRAINT_CONST     = 5,
                                      PRIORITISE_RHS_FLAG = FALSE,
                                      ASTERIX_FLAG        = TRUE,
                                      PROBE_HALVES_FLAG   = TRUE,
                                      MRNA_FLAG           = FALSE) {
  # lapply over row indices keeps each row as a real single-row data frame,
  # so column access with $ works correctly inside create_probe().
  # (apply(..., 1, FUN) would coerce the data frame to a character matrix.)
  output_list <- lapply(seq_len(nrow(input_df)), function(i) {
    row <- input_df[i, , drop = FALSE]
    create_probe(row$gene1, row$gene2, row$gene1_transcript, row$gene2_transcript,
                 RESTRAINT_CONST, PRIORITISE_RHS_FLAG, ASTERIX_FLAG, PROBE_HALVES_FLAG, MRNA_FLAG)
  })

  # Remove NULL and empty data frames (e.g. fusions where gaps prevented probe design)
  # before rbind to avoid introducing NA rows into the combined output.
  output_list <- Filter(function(x) !is.null(x) && nrow(x) > 0, output_list)

  if (length(output_list) == 0) {
    warning("No probes could be designed for any fusion.")
    return(data.frame())
  }

  output_df <- do.call(rbind, output_list)
  return(output_df)
}


# =============================================================================
# INPUT PROCESSING
# =============================================================================

#' Validate and clean an Arriba fusion input data frame.
#'
#' Ensures that the required columns are present (gene1, gene2,
#' gene1_transcript, gene2_transcript), converts column names to lowercase,
#' trims whitespace, and converts sequences to uppercase.
#'
#' @param input_df data frame — raw CSV loaded from an Arriba-derived file
#' @return data frame — cleaned input ready for create_probes_from_arriba()
#' @export
process_arriba_transcript <- function(input_df) {
  # Normalise column names
  colnames(input_df) <- tolower(trimws(colnames(input_df)))

  required_cols <- c("gene1", "gene2", "gene1_transcript", "gene2_transcript")
  missing <- setdiff(required_cols, colnames(input_df))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "),
         "\nExpected columns: gene1, gene2, gene1_transcript, gene2_transcript")
  }

  # Clean sequence and name fields
  input_df$gene1_transcript <- toupper(trimws(input_df$gene1_transcript))
  input_df$gene2_transcript <- toupper(trimws(input_df$gene2_transcript))
  input_df$gene1            <- trimws(input_df$gene1)
  input_df$gene2            <- trimws(input_df$gene2)

  return(input_df)
}


# =============================================================================
# ARRIBA INPUT PARSING
# =============================================================================

#' Parse an Arriba fusion TSV file into the 4-column format required by
#' create_probes_from_arriba().
#'
#' Arriba outputs a tab-separated file where:
#'   - The first column is named "#gene1" (hash prefix from the comment-style header)
#'   - The second column is "gene2"
#'   - The "fusion_transcript" column contains the assembled fusion sequence with
#'     a "|" character marking the exact breakpoint between the two gene partners
#'   - "___" or "." in the fusion sequence marks low-coverage or intronic gaps
#'
#' This function extracts gene names and splits the fusion_transcript on "|" to
#' recover the two flanking sequences. Gap markers are preserved in the returned
#' sequences so that probe enumeration can identify and exclude any candidate
#' probe whose 50 bp window overlaps a gap position.
#'
#' Rows where the fusion_transcript is missing, ".", or lacks a "|" separator
#' are dropped with a warning.
#'
#' @param filepath character — path to the Arriba TSV file
#' @return data frame with columns gene1, gene2, gene1_transcript, gene2_transcript,
#'         ready for create_probes_from_arriba()
parse_arriba_tsv <- function(filepath) {
  # Read TSV; suppress warnings about incomplete final line which Arriba sometimes produces
  raw <- suppressWarnings(
    read.table(filepath, sep = "\t", header = TRUE,
               stringsAsFactors = FALSE, quote = "", comment.char = "")
  )

  # Arriba's first column is named "X.gene1" after R reads the "#gene1" header,
  # or sometimes "#gene1" depending on R version. Normalise to "gene1".
  colnames(raw) <- sub("^[#X\\.]+gene1$", "gene1", colnames(raw))
  colnames(raw) <- tolower(trimws(colnames(raw)))

  required <- c("gene1", "gene2", "fusion_transcript")
  missing  <- setdiff(required, colnames(raw))
  if (length(missing) > 0) {
    stop("Arriba TSV is missing expected columns: ", paste(missing, collapse = ", "),
         "\nEnsure this is a standard Arriba output file.")
  }

  n_input <- nrow(raw)

  # Drop rows where fusion_transcript is absent or uninformative
  valid <- raw$fusion_transcript != "." &
           !is.na(raw$fusion_transcript) &
           nchar(trimws(raw$fusion_transcript)) > 0 &
           grepl("|", raw$fusion_transcript, fixed = TRUE)

  n_dropped <- sum(!valid)
  if (n_dropped > 0) {
    warning(n_dropped, " row(s) dropped: fusion_transcript missing or lacks a '|' separator.")
  }
  raw <- raw[valid, ]

  if (nrow(raw) == 0) {
    stop("No valid fusions found in the Arriba TSV. ",
         "Check that the fusion_transcript column contains sequences separated by '|'.")
  }

  # Split fusion_transcript on "|" to recover the two flanking sequences.
  # Left of "|"  = gene1 sequence ending at the breakpoint (gene1_transcript)
  # Right of "|" = gene2 sequence starting at the breakpoint (gene2_transcript)
  split_seqs <- strsplit(raw$fusion_transcript, "|", fixed = TRUE)

  gene1_transcript <- sapply(split_seqs, function(x) trimws(x[1]))
  gene2_transcript <- sapply(split_seqs, function(x) trimws(x[2]))

  # Convert to uppercase and preserve gap markers (both "___" and "." from Arriba).
  # Gap characters are kept in the sequences so that probe enumeration can detect
  # and exclude any candidate whose 50 bp window overlaps a gap position.
  clean_gene1 <- toupper(trimws(gene1_transcript))
  clean_gene2 <- toupper(trimws(gene2_transcript))

  # Warn if either sequence has fewer than min_bp valid nucleotide bases adjacent
  # to the breakpoint. Gap characters do not count as usable sequence.
  # produce_possible_probe_df() needs at least 25 + RESTRAINT_CONST bp per side.
  min_bp <- 30L
  nucleotide_len1 <- nchar(gsub("[^ATGCNatgcn]", "", clean_gene1))
  nucleotide_len2 <- nchar(gsub("[^ATGCNatgcn]", "", clean_gene2))
  short1 <- nucleotide_len1 < min_bp
  short2 <- nucleotide_len2 < min_bp

  if (any(short1 | short2)) {
    bad_fusions <- paste(raw$gene1[short1 | short2], raw$gene2[short1 | short2], sep = "::")
    warning(
      length(bad_fusions), " fusion(s) have fewer than ", min_bp,
      " bp of valid nucleotide sequence adjacent to the breakpoint. ",
      "Probe design may produce no candidates for: ",
      paste(bad_fusions, collapse = ", ")
    )
  }

  data.frame(
    gene1             = trimws(raw$gene1),
    gene2             = trimws(raw$gene2),
    gene1_transcript  = clean_gene1,
    gene2_transcript  = clean_gene2,
    stringsAsFactors  = FALSE
  )
}
