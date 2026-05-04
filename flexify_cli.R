#!/usr/bin/env Rscript
# =============================================================================
# Flexify — Command-Line Interface
# =============================================================================
#
# A command-line wrapper for the Flexify probe design pipeline. Supports three
# modes: probe design from Arriba input, BLAST off-target filtering, and
# handle/barcode appending to generate final synthesis-ready sequences.
#
# USAGE:
#
#   Mode 1 — Design probes (default):
#     Rscript flexify_cli.R --input fusions.csv --output probes.csv
#     Rscript flexify_cli.R -i fusions.csv -o probes.csv --restraint 5 --mrna
#
#   Mode 2 — BLAST off-target check:
#     Rscript flexify_cli.R --mode blast \
#       --input probes.csv --blast-db /path/to/transcriptome_db \
#       --output probes_filtered.csv
#
#   Mode 3 — Append handles (v1, with barcode):
#     Rscript flexify_cli.R --mode finalise \
#       --input selected_probes.csv --output final_probes.csv
#
#     Input CSV columns required: GENE1, GENE2, probe, Barcode (integer 1-16 or "BC001"-"BC016")
#
#   Mode 3 — Append handles (v2, no barcode):
#     Rscript flexify_cli.R --mode finalise --assay-version v2 \
#       --input selected_probes.csv --output final_probes.csv
#
#     Optionally set --rhs-mode singleplex for the 4-sample kit tail (default: multiplex).
#     Input CSV columns required: GENE1, GENE2, probe (no Barcode column needed)
#
# DEPENDENCIES:
#   optparse, tidyverse, stringr
#   flexify_core.R, flexify_offtarget.R, flexify_handles.R must be in the
#   same directory as this script.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(stringr)
})

# Locate the directory containing this script so helper modules can be sourced
# regardless of where the script is called from.
script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)

source(file.path(script_dir, "flexify_core.R"))
source(file.path(script_dir, "flexify_offtarget.R"))
source(file.path(script_dir, "flexify_handles.R"))


# =============================================================================
# ARGUMENT DEFINITIONS
# =============================================================================

option_list <- list(

  # --- General ---------------------------------------------------------------
  make_option(c("-i", "--input"),
              type    = "character",
              default = NULL,
              help    = "Path to input CSV file [required]",
              metavar = "FILE"),

  make_option(c("-o", "--output"),
              type    = "character",
              default = NULL,
              help    = "Path for output CSV file [required]",
              metavar = "FILE"),

  make_option(c("-m", "--mode"),
              type    = "character",
              default = "design",
              help    = "Pipeline mode: 'design', 'blast', or 'finalise' [default: design]",
              metavar = "MODE"),

  make_option("--arriba",
              action  = "store_true",
              default = FALSE,
              help    = "Treat --input as an Arriba TSV file (auto-parses gene names and fusion_transcript column) [default: FALSE]"),

  # --- Design mode options ---------------------------------------------------
  make_option("--restraint",
              type    = "integer",
              default = 5L,
              help    = "Minimum bases each gene must contribute to a probe half [default: 5]",
              metavar = "INT"),

  make_option("--mrna",
              action  = "store_true",
              default = FALSE,
              help    = "Include mRNA target sequence column in output [default: FALSE]"),

  make_option("--prioritise-rhs",
              action  = "store_true",
              default = FALSE,
              help    = "Apply 0.7x score penalty to probes with junction in left half [default: FALSE]"),

  make_option("--no-asterix",
              action  = "store_true",
              default = FALSE,
              help    = "Do not mark the fusion breakpoint with '*' in output [default: FALSE]"),

  make_option("--no-halves",
              action  = "store_true",
              default = FALSE,
              help    = "Do not mark probe half boundary with '|' in output [default: FALSE]"),

  # --- BLAST mode options ----------------------------------------------------
  make_option("--blast-db",
              type    = "character",
              default = NULL,
              help    = "Path to BLAST nucleotide database (no extension) [required for blast mode]",
              metavar = "PATH"),

  make_option("--min-mismatches",
              type    = "integer",
              default = 5L,
              help    = "Minimum mismatches to off-target for a probe to pass [default: 5]",
              metavar = "INT"),

  make_option("--threads",
              type    = "integer",
              default = 1L,
              help    = "Number of BLAST threads [default: 1]",
              metavar = "INT"),

  make_option("--keep-fails",
              action  = "store_true",
              default = FALSE,
              help    = "Keep probes that fail off-target check (adds flag columns but does not filter) [default: FALSE]"),

  # --- Finalise mode options -------------------------------------------------
  make_option("--assay-version",
              type    = "character",
              default = "v1",
              help    = "Assay version for finalise mode: 'v1' (Chromium Flex, barcode embedded) or 'v2' (GEM-X Flex, no barcode) [default: v1]",
              metavar = "VERSION"),

  make_option("--rhs-mode",
              type    = "character",
              default = "multiplex",
              help    = "v2 RHS tail configuration: 'multiplex' (CCCATATAAGAAA) or 'singleplex' (CGGTCCTAGCAA) [default: multiplex]",
              metavar = "MODE")
)

# Parse arguments
opt_parser <- OptionParser(
  option_list  = option_list,
  description  = "\nFlexify: automated fusion probe design for 10x Genomics Flex.",
  epilogue     = paste(
    "Examples:",
    "  Rscript flexify_cli.R -i fusions.csv -o probes.csv",
    "  Rscript flexify_cli.R --mode blast -i probes.csv --blast-db /path/db -o filtered.csv",
    "  Rscript flexify_cli.R --mode finalise -i selected.csv -o final.csv",
    sep = "\n"
  )
)
opt <- parse_args(opt_parser)


# =============================================================================
# INPUT VALIDATION
# =============================================================================

if (is.null(opt$input)) {
  stop("--input is required. Run with --help for usage.")
}
if (!file.exists(opt$input)) {
  stop("Input file not found: ", opt$input)
}
if (is.null(opt$output)) {
  stop("--output is required. Run with --help for usage.")
}
if (!opt$mode %in% c("design", "blast", "finalise")) {
  stop("--mode must be one of: design, blast, finalise. Got: ", opt$mode)
}


# =============================================================================
# MODE: DESIGN
# Run probe design pipeline from an Arriba-derived fusion CSV.
# =============================================================================

if (opt$mode == "design") {
  message("Flexify | Mode: design")
  message("  Input:   ", opt$input)
  message("  Output:  ", opt$output)
  message("  Restraint constant: ", opt$restraint)

  # Parse input: Arriba TSV directly, or generic 4-column CSV
  if (isTRUE(opt$arriba)) {
    message("  Input format: Arriba TSV")
    cleaned_df <- tryCatch(
      parse_arriba_tsv(opt$input),
      error = function(e) stop("Arriba TSV parsing failed: ", conditionMessage(e))
    )
  } else {
    message("  Input format: generic CSV")
    raw_df <- read.csv(opt$input, stringsAsFactors = FALSE)
    cleaned_df <- tryCatch(
      process_arriba_transcript(raw_df),
      error = function(e) stop("Input validation failed: ", conditionMessage(e))
    )
  }

  n_fusions <- nrow(cleaned_df)
  message("  Fusions in input: ", n_fusions)

  # Run probe design
  output_df <- tryCatch(
    create_probes_from_arriba(
      cleaned_df,
      RESTRAINT_CONST     = opt$restraint,
      PRIORITISE_RHS_FLAG = isTRUE(opt[["prioritise-rhs"]]),
      ASTERIX_FLAG        = !isTRUE(opt[["no-asterix"]]),
      PROBE_HALVES_FLAG   = !isTRUE(opt[["no-halves"]]),
      MRNA_FLAG           = isTRUE(opt$mrna)
    ),
    error = function(e) stop("Probe design failed: ", conditionMessage(e))
  )

  write.csv(output_df, opt$output, row.names = FALSE)
  message("Done. ", nrow(output_df), " candidate probes written to: ", opt$output)
}


# =============================================================================
# MODE: BLAST
# Run BLAST off-target check on a ranked probe CSV.
# =============================================================================

if (opt$mode == "blast") {
  message("Flexify | Mode: blast off-target check")
  message("  Input:    ", opt$input)
  message("  BLAST db: ", opt[["blast-db"]])
  message("  Output:   ", opt$output)

  if (is.null(opt[["blast-db"]])) {
    stop("--blast-db is required in blast mode. Provide the path to your BLAST database.")
  }

  probe_df <- read.csv(opt$input, stringsAsFactors = FALSE)
  message("  Probes loaded: ", nrow(probe_df))

  filtered_df <- tryCatch(
    run_offtarget_check(
      probe_df       = probe_df,
      blast_db       = opt[["blast-db"]],
      min_mismatches = opt[["min-mismatches"]],
      n_threads      = opt$threads,
      filter_fails   = !isTRUE(opt[["keep-fails"]])
    ),
    error = function(e) stop("BLAST check failed: ", conditionMessage(e))
  )

  write.csv(filtered_df, opt$output, row.names = FALSE)
  message("Done. ", nrow(filtered_df), " probes written to: ", opt$output)
}


# =============================================================================
# MODE: FINALISE
# Append handle sequences and barcodes to a user-selected probe CSV.
# Input CSV must have columns: GENE1, GENE2, probe, Barcode (1-16 or BC001-BC016).
# If a 'Selected' column is present, only rows with Selected == TRUE are processed.
# =============================================================================

if (opt$mode == "finalise") {

  assay_ver <- opt[["assay-version"]]
  if (!assay_ver %in% c("v1", "v2")) {
    stop("--assay-version must be 'v1' or 'v2'. Got: ", assay_ver)
  }
  is_v2 <- assay_ver == "v2"

  rhs_mode <- opt[["rhs-mode"]]
  if (is_v2 && !rhs_mode %in% c("multiplex", "singleplex")) {
    stop("--rhs-mode must be 'multiplex' or 'singleplex'. Got: ", rhs_mode)
  }

  message("Flexify | Mode: finalise (append handles)")
  message("  Input:         ", opt$input)
  message("  Output:        ", opt$output)
  message("  Assay version: ", assay_ver)
  if (is_v2) message("  RHS mode:      ", rhs_mode)

  input_df <- read.csv(opt$input, stringsAsFactors = FALSE)

  # If a Selected column exists, filter to marked rows only
  if ("Selected" %in% colnames(input_df)) {
    n_before <- nrow(input_df)
    input_df <- input_df[which(input_df$Selected == TRUE | input_df$Selected == "TRUE"), ]
    message("  Filtered to Selected == TRUE: ", nrow(input_df), " of ", n_before, " rows")
  }

  if (nrow(input_df) == 0) {
    stop("No rows to process. If a 'Selected' column is present, ensure some rows are TRUE.")
  }

  # Column requirements differ by assay version
  required <- if (is_v2) c("GENE1", "GENE2", "probe") else c("GENE1", "GENE2", "probe", "Barcode")
  missing  <- setdiff(required, colnames(input_df))
  if (length(missing) > 0) {
    stop("Input CSV is missing required columns: ", paste(missing, collapse = ", "))
  }

  final_df <- tryCatch(
    if (is_v2) finalise_probes_v2(input_df, rhs_mode = rhs_mode)
    else        finalise_probes(input_df),
    error = function(e) stop("Finalise failed: ", conditionMessage(e))
  )

  write.csv(final_df, opt$output, row.names = FALSE)
  message("Done. ", nrow(final_df), " probe sequences written to: ", opt$output)
  message("  Each fusion produces one LHS and one RHS oligonucleotide.")
}
