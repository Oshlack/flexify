# =============================================================================
# Flexify — Testing Framework
# =============================================================================
#
# Run this file to validate every core function in Flexify using small
# synthetic inputs. No external dependencies (BLAST, Shiny, real data)
# are required for most tests. The BLAST test is skipped automatically
# if blastn is not on the PATH.
#
# USAGE:
#   Rscript test_flexify.R
#
# OUTPUT:
#   Each test prints  [PASS] or [FAIL] with a short description.
#   A summary line at the end reports the overall result.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})

# Locate the directory containing this script so that the Flexify modules
# can be sourced regardless of where the script is called from.
script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)

source(file.path(script_dir, "flexify_core.R"))
source(file.path(script_dir, "flexify_handles.R"))
source(file.path(script_dir, "flexify_offtarget.R"))
source(file.path(script_dir, "flexify_nonfusion.R"))


# =============================================================================
# TEST HARNESS
# =============================================================================

.pass_count <- 0L
.fail_count <- 0L
.test_log   <- character(0)

#' Run a single named test expression.
#'
#' @param name  character — short description of the test
#' @param expr  expression — should evaluate to TRUE for a pass
expect <- function(name, expr) {
  result <- tryCatch({
    val <- eval(expr, envir = parent.frame())
    if (isTRUE(val)) "PASS" else paste0("FAIL — expression returned: ", deparse(val))
  }, error = function(e) {
    paste0("FAIL — error: ", conditionMessage(e))
  }, warning = function(w) {
    # Warnings don't automatically fail a test; re-run without intercepting them
    val <- withCallingHandlers(
      eval(expr, envir = parent.frame()),
      warning = function(w) invokeRestart("muffleWarning")
    )
    if (isTRUE(val)) "PASS" else paste0("FAIL — expression returned: ", deparse(val))
  })

  status <- if (startsWith(result, "PASS")) "PASS" else "FAIL"
  if (status == "PASS") {
    .pass_count <<- .pass_count + 1L
    cat(sprintf("  [PASS] %s\n", name))
  } else {
    .fail_count <<- .fail_count + 1L
    cat(sprintf("  [FAIL] %s\n         %s\n", name, result))
  }
  .test_log <<- c(.test_log, sprintf("[%s] %s", status, name))
}

#' Print a section header.
section <- function(title) {
  cat(sprintf("\n── %s %s\n", title, strrep("─", max(0, 60 - nchar(title)))))
}

# =============================================================================
# SYNTHETIC TEST DATA
# =============================================================================

# A pair of 60 bp synthetic transcript sequences (gene1 = right of junction,
# gene2 = left of junction in Arriba convention).
# These are long enough for the probe enumerator and have ~50% GC.
GENE1 <- "BCR"
GENE2 <- "ABL1"
# 60 bp: designed so GC ~ 50% and no long homopolymers
SEQ1 <- "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"  # 60 bp
SEQ2 <- "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"  # 61 bp – trimmed in use

# 80 bp non-fusion sequence for tile_sequence / create_nonfusion_probes tests.
# ~50% GC, no homopolymer runs >= 4.
NONFUSION_GENE <- "GFP"
NONFUSION_SEQ  <- paste0(rep("ATGCATGC", 10), collapse = "")  # 80 bp

# Synthetic standard probeset data frame for competition check tests.
# Mimics the output of load_flex_probeset().
# Probe 1: LHS = 25xA, RHS = 25xT
# Probe 2: LHS = 25xG, RHS = 25xC
synth_probeset_df <- data.frame(
  probe_seq = c(
    paste0(strrep("A", 25), strrep("T", 25)),
    paste0(strrep("G", 25), strrep("C", 25))
  ),
  included = TRUE,
  stringsAsFactors = FALSE
)


# =============================================================================
# SECTION 1: SEQUENCE UTILITIES
# =============================================================================
section("1. reverse_complement()")

expect("RC of 'ATGC' is 'GCAT'",
  quote(reverse_complement("ATGC") == "GCAT"))

expect("RC of 'AAAA' is 'TTTT'",
  quote(reverse_complement("AAAA") == "TTTT"))

expect("RC of 'GCTAGC' is 'GCTAGC' (palindrome)",
  quote(reverse_complement("GCTAGC") == "GCTAGC"))

expect("RC is its own inverse (double complement = identity)",
  quote(reverse_complement(reverse_complement(SEQ1)) == SEQ1))


# =============================================================================
# SECTION 2: PROBE ENUMERATION
# =============================================================================
section("2. produce_possible_probe_df()")

left_str  <- reverse_complement(SEQ2)
right_str <- reverse_complement(SEQ1)

probe_df <- tryCatch(
  produce_possible_probe_df(left_str, right_str, RESTRAINT_CONST = 5),
  error = function(e) NULL
)

expect("Returns a data frame",
  quote(is.data.frame(probe_df)))

expect("Has 'probe' and 'fusion_point_displacement' columns",
  quote(all(c("probe", "fusion_point_displacement") %in% colnames(probe_df))))

expect("All probes are exactly 50 bp",
  quote(all(nchar(probe_df$probe) == 50)))

expect("No probes have |displacement| < RESTRAINT_CONST (5)",
  quote(all(abs(probe_df$fusion_point_displacement) >= 5)))

expect("Displacements range from -(25-5) to +(25-5), excluding 0",
  quote({
    d <- probe_df$fusion_point_displacement
    min(d) >= -20 && max(d) <= 20 && !any(d == 0)
  }))


# =============================================================================
# SECTION 3: SCORING FUNCTIONS
# =============================================================================
section("3a. GC scoring — check_GC() and GC_rating()")

expect("check_GC on 'ATATATATATAT ATATAT ATATAT' (0% GC) = 0",
  quote(check_GC("AAAAAAAAAAAAAAAAAAAAAAAAAA") == 0))  # 26A = 0% GC (but check_GC uses 25bp)

expect("check_GC on all-GC 25-mer = 100",
  quote(check_GC(strrep("G", 25)) == 100))

expect("GC_rating(50) = 5 (optimal)",
  quote(GC_rating(50) == 5))

expect("GC_rating(43) = 0 (below minimum)",
  quote(GC_rating(43) == 0))

expect("GC_rating(73) = 0 (above maximum)",
  quote(GC_rating(73) == 0))

expect("GC_rating(44) = 2",
  quote(GC_rating(44) == 2))

expect("GC_rating(70) = 1",
  quote(GC_rating(70) == 1))


section("3b. Junction position scoring — fusion_location_rating()")

expect("Score at displacement 12.5 = 5 (optimum)",
  quote(abs(fusion_location_rating(12.5) - 5) < 0.001))

expect("Score at displacement 0 is lower than at 12.5",
  quote(fusion_location_rating(0) < fusion_location_rating(12.5)))

expect("Score is symmetric (same for +12 and -12)",
  quote(abs(fusion_location_rating(12) - fusion_location_rating(-12)) < 0.001))

expect("Score is non-negative",
  quote(fusion_location_rating(25) >= 0))


section("3c. Dinucleotide check — dinucleotide_check()")

# Build a probe with 'AT' at positions 25-26 (preferred)
probe_AT <- paste0(strrep("A", 24), "AT", strrep("G", 24))  # 50 bp, positions 25-26 = "AT"
di_AT    <- dinucleotide_check(probe_AT)

expect("dinucleotide_check returns a list of length 3",
  quote(length(di_AT) == 3))

expect("Preferred dinucleotide 'AT' gets status 'OK'",
  quote(di_AT[[2]] == "OK"))

expect("Preferred dinucleotide 'AT' gets rating 3",
  quote(as.numeric(di_AT[[3]]) == 3))

probe_GG <- paste0(strrep("A", 24), "GG", strrep("C", 24))  # positions 25-26 = "GG" (not preferred)
di_GG    <- dinucleotide_check(probe_GG)

expect("Non-preferred dinucleotide 'GG' gets status 'Warning'",
  quote(di_GG[[2]] == "Warning"))

expect("Non-preferred dinucleotide 'GG' gets rating 1",
  quote(as.numeric(di_GG[[3]]) == 1))


section("3d. Homopolymer check — check_for_homo_polymer()")

probe_ok     <- paste0(strrep("ATGC", 12), "AT")   # no homopolymer runs ≥4
probe_run4   <- paste0("AAAACCGGTTCCAATTGGCCAATTG", "CCAATTGGCCAATTGGCCAATTTG")  # build one with run of 4
probe_run5   <- paste0(strrep("A", 5), strrep("GCGC", 11), "GCGCG")             # 5-A at start

hp_ok    <- check_for_homo_polymer(probe_ok)
hp_run5  <- check_for_homo_polymer(probe_run5)

expect("check_for_homo_polymer returns a list of length 2",
  quote(length(hp_ok) == 2))

expect("Clean probe (no run ≥4) gets rating 2",
  quote(as.numeric(hp_ok[[2]]) == 2))

expect("Probe with run of 5 gets rating 1",
  quote(as.numeric(hp_run5[[2]]) == 1))

expect("Probe with run of 5 gets a Warning flag",
  quote(grepl("Warning", hp_run5[[1]])))


# =============================================================================
# SECTION 4: SINGLE-FUSION PROBE DESIGN
# =============================================================================
section("4. create_probe() — single fusion")

result <- tryCatch(
  create_probe(GENE1, GENE2, SEQ1, SEQ2,
               RESTRAINT_CONST = 5, ASTERIX_FLAG = TRUE, PROBE_HALVES_FLAG = TRUE),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)

expect("create_probe() returns a data frame",
  quote(is.data.frame(result)))

expect("Output has GENE1 and GENE2 columns",
  quote(all(c("GENE1", "GENE2") %in% colnames(result))))

expect("GENE1 values are correct",
  quote(all(result$GENE1 == GENE1)))

expect("GENE2 values are correct",
  quote(all(result$GENE2 == GENE2)))

expect("Ranking column is present and starts at 1",
  quote("Ranking" %in% colnames(result) && result$Ranking[1] == 1))

expect("Rankings are consecutive integers starting at 1",
  quote(all(result$Ranking == seq_len(nrow(result)))))

expect("Score column is present and all > 0",
  quote("Score" %in% colnames(result) && all(result$Score > 0)))

expect("Scores are in descending order",
  quote(all(diff(result$Score) <= 0)))

expect("Asterisk (*) is present in at least one probe",
  quote(any(grepl("\\*", result$probe))))

expect("Pipe (|) is present in at least one probe",
  quote(any(grepl("\\|", result$probe))))

# Test with ASTERIX and HALVES off — clean 50-char sequences
result_clean <- tryCatch(
  create_probe(GENE1, GENE2, SEQ1, SEQ2,
               ASTERIX_FLAG = FALSE, PROBE_HALVES_FLAG = FALSE),
  error = function(e) NULL
)

expect("With markers off, probes are exactly 50 characters",
  quote(!is.null(result_clean) && all(nchar(result_clean$probe) == 50)))

# Test mRNA flag
result_mrna <- tryCatch(
  create_probe(GENE1, GENE2, SEQ1, SEQ2, MRNA_FLAG = TRUE,
               ASTERIX_FLAG = FALSE, PROBE_HALVES_FLAG = FALSE),
  error = function(e) NULL
)

expect("With MRNA_FLAG = TRUE, 'mrna' column is present",
  quote(!is.null(result_mrna) && "mrna" %in% colnames(result_mrna)))

expect("mRNA sequence is RC of probe (for first row)",
  quote({
    p    <- result_mrna$probe[1]
    m    <- as.character(result_mrna$mrna[1])
    reverse_complement(p) == m
  }))


# =============================================================================
# SECTION 5: MULTI-FUSION PROBE DESIGN
# =============================================================================
section("5. create_probes_from_arriba() — multiple fusions")

# Build a 3-row input data frame (two real-ish fusions + one duplicate)
input_multi <- data.frame(
  gene1             = c("BCR",  "EML4",  "BCR"),
  gene2             = c("ABL1", "ALK",   "ABL1"),
  gene1_transcript  = c(SEQ1,   SEQ1,    SEQ1),
  gene2_transcript  = c(SEQ2,   SEQ2,    SEQ2),
  stringsAsFactors  = FALSE
)

multi_result <- tryCatch(
  create_probes_from_arriba(input_multi, RESTRAINT_CONST = 5),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)

expect("create_probes_from_arriba() returns a data frame",
  quote(is.data.frame(multi_result)))

expect("Output contains rows from all 3 input fusions",
  quote({
    genes <- paste(multi_result$GENE1, multi_result$GENE2, sep = "::")
    all(c("BCR::ABL1", "EML4::ALK") %in% genes)
  }))

expect("Ranking resets to 1 for each fusion",
  quote({
    min_rank_per_fusion <- multi_result %>%
      group_by(GENE1, GENE2) %>%
      summarise(min_rank = min(Ranking), .groups = "drop")
    all(min_rank_per_fusion$min_rank == 1)
  }))


# =============================================================================
# SECTION 6: ARRIBA TSV PARSING
# =============================================================================
section("6. parse_arriba_tsv()")

# Write a minimal synthetic Arriba-style TSV to a temp file
arriba_tsv_content <- paste(
  "#gene1\tgene2\tfusion_transcript\tconfidence",
  paste0("BCR\tABL1\t", SEQ1, "|", SEQ2, "\thigh"),
  paste0("EML4\tALK\t",  SEQ1, "|", SEQ2, "\thigh"),
  paste0("SKIP\tME\t.\thigh"),                           # should be dropped (no |)
  sep = "\n"
)

tmp_tsv <- tempfile(fileext = ".tsv")
writeLines(arriba_tsv_content, tmp_tsv)

parsed <- tryCatch(
  suppressWarnings(parse_arriba_tsv(tmp_tsv)),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)

expect("parse_arriba_tsv() returns a data frame",
  quote(is.data.frame(parsed)))

expect("Correct number of valid rows returned (2 of 3)",
  quote(nrow(parsed) == 2))

expect("Output has expected columns",
  quote(all(c("gene1", "gene2", "gene1_transcript", "gene2_transcript") %in% colnames(parsed))))

expect("gene1 column contains correct values",
  quote(all(parsed$gene1 %in% c("BCR", "EML4"))))

expect("Sequences are uppercase",
  quote(all(parsed$gene1_transcript == toupper(parsed$gene1_transcript))))

# Test with a TSV that has '___' gap markers.
# Gaps are now preserved in the returned sequences so that probe enumeration
# can exclude any candidate whose window overlaps a gap position.
arriba_gap_content <- paste(
  "#gene1\tgene2\tfusion_transcript",
  paste0("PTEN\tTP53\t", "AAAA___", SEQ1, "|", SEQ2, "___TTTT"),
  sep = "\n"
)
tmp_gap <- tempfile(fileext = ".tsv")
writeLines(arriba_gap_content, tmp_gap)

parsed_gap <- tryCatch(
  suppressWarnings(parse_arriba_tsv(tmp_gap)),
  error = function(e) NULL
)

expect("'___' gaps are preserved in gene1_transcript (not trimmed)",
  quote(!is.null(parsed_gap) && grepl("___", parsed_gap$gene1_transcript[1], fixed = TRUE)))

expect("'___' gaps are preserved in gene2_transcript (not trimmed)",
  quote(!is.null(parsed_gap) && grepl("___", parsed_gap$gene2_transcript[1], fixed = TRUE)))

# Probes spanning the gap should be excluded during enumeration.
# The gap is far from the breakpoint here (AAAA___ prefix / ___TTTT suffix),
# so valid probes should still be produced from the clean sequence near the junction.
probes_gap <- tryCatch(
  suppressWarnings(
    create_probes_from_arriba(parsed_gap, RESTRAINT_CONST = 5,
                              ASTERIX_FLAG = FALSE, PROBE_HALVES_FLAG = FALSE)
  ),
  error = function(e) NULL
)

expect("Probes are still designed when gaps are far from the breakpoint",
  quote(!is.null(probes_gap) && nrow(probes_gap) > 0))

expect("No designed probe contains a gap character",
  quote(!is.null(probes_gap) && !any(grepl("[^ATGCNatgcn*|]", probes_gap$probe))))

# Test that missing fusion_transcript column raises an error
bad_tsv_content <- "col1\tcol2\nA\tB\n"
tmp_bad <- tempfile(fileext = ".tsv")
writeLines(bad_tsv_content, tmp_bad)

parsed_bad <- tryCatch(
  parse_arriba_tsv(tmp_bad),
  error = function(e) e
)
expect("Missing fusion_transcript column raises an error",
  quote(inherits(parsed_bad, "error")))


# =============================================================================
# SECTION 7: GENERIC CSV INPUT VALIDATION
# =============================================================================
section("7. process_arriba_transcript()")

valid_df <- data.frame(
  Gene1             = c("BCR"),
  GENE2             = c("ABL1"),
  gene1_transcript  = c(SEQ1),
  gene2_transcript  = c(SEQ2),
  stringsAsFactors  = FALSE
)

processed <- tryCatch(
  process_arriba_transcript(valid_df),
  error = function(e) NULL
)

expect("process_arriba_transcript() returns a data frame",
  quote(is.data.frame(processed)))

expect("Column names are lowercased",
  quote(all(colnames(processed) == tolower(colnames(processed)))))

expect("Sequences are uppercase",
  quote(processed$gene1_transcript[1] == toupper(SEQ1)))

# Missing column should raise an error
bad_df <- data.frame(gene1 = "BCR", gene2 = "ABL1", stringsAsFactors = FALSE)
err_result <- tryCatch(process_arriba_transcript(bad_df), error = function(e) e)

expect("Missing columns raise an error",
  quote(inherits(err_result, "error")))

expect("Error message mentions missing column names",
  quote(grepl("gene1_transcript", conditionMessage(err_result))))


# =============================================================================
# SECTION 8: HANDLE & BARCODE APPENDING
# =============================================================================
section("8. finalise_probes()")

LHS_HANDLE <- "CCTTGGCACCCGAGAATTCCA"   # 21 bp
RHS_TAIL   <- "CGGTCCTAGCAA"             # 12 bp
RHS_LINKER <- "ACGCGGTTAGCACGTA"         # 16 bp
BC001_SEQ  <- "ACTTTAGG"                 # 8 bp

# A clean 50 bp probe (no markers)
test_probe <- paste0(strrep("A", 25), strrep("T", 25))

sel_df <- data.frame(
  GENE1   = "BCR",
  GENE2   = "ABL1",
  probe   = test_probe,
  Barcode = 1L,
  stringsAsFactors = FALSE
)

final <- tryCatch(finalise_probes(sel_df), error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })

expect("finalise_probes() returns a data frame",
  quote(is.data.frame(final)))

expect("Output has expected columns",
  quote(all(c("Fusion", "Barcode_ID", "Pool_Name", "Barcode_Seq", "LHS_Probe", "RHS_Probe") %in% colnames(final))))

expect("Fusion name formatted as GENE1::GENE2",
  quote(final$Fusion[1] == "BCR::ABL1"))

expect("Barcode_ID is 'BC001' for barcode 1",
  quote(final$Barcode_ID[1] == "BC001"))

expect("Pool_Name is 'poolOne' for barcode 1",
  quote(final$Pool_Name[1] == "poolOne"))

expect("LHS_Probe starts with the 21 bp constant handle",
  quote(startsWith(final$LHS_Probe[1], LHS_HANDLE)))

expect("LHS_Probe has correct total length (21 + 25 = 46 bp)",
  quote(nchar(final$LHS_Probe[1]) == 46))

expect("RHS_Probe starts with /5Phos/",
  quote(startsWith(final$RHS_Probe[1], "/5Phos/")))

expect("RHS_Probe contains BC001 barcode sequence",
  quote(grepl(BC001_SEQ, final$RHS_Probe[1])))

expect("RHS_Probe ends with the 12 bp constant tail",
  quote(endsWith(final$RHS_Probe[1], RHS_TAIL)))

expect("RHS_Probe has correct total length (/5Phos/ + 25 + 16 + 2 + 8 + 12 = 70 characters)",
  quote(nchar(final$RHS_Probe[1]) == nchar("/5Phos/") + 25 + 16 + 2 + 8 + 12))

# Test with '|' and '*' markers in probe — should be stripped
marked_probe <- paste0(strrep("A", 12), "*", strrep("A", 13), "|", strrep("T", 25))
sel_marked <- data.frame(GENE1="X", GENE2="Y", probe=marked_probe, Barcode=2L,
                         stringsAsFactors=FALSE)
final_marked <- tryCatch(finalise_probes(sel_marked), error = function(e) NULL)

expect("Annotation markers (*  |) are stripped before handle appending",
  quote(!is.null(final_marked) && !grepl("[*|]", final_marked$LHS_Probe[1])))

# Test barcode string format "BC003"
sel_str_bc <- data.frame(GENE1="A", GENE2="B", probe=test_probe, Barcode="BC003",
                         stringsAsFactors=FALSE)
final_str <- tryCatch(finalise_probes(sel_str_bc), error = function(e) NULL)

expect("String barcode 'BC003' resolves correctly",
  quote(!is.null(final_str) && final_str$Barcode_ID[1] == "BC003"))

# Out-of-range barcode should raise an error
sel_bad_bc <- data.frame(GENE1="A", GENE2="B", probe=test_probe, Barcode=99L,
                         stringsAsFactors=FALSE)
err_bc <- tryCatch(finalise_probes(sel_bad_bc), error = function(e) e)
expect("Out-of-range barcode (99) raises an error",
  quote(inherits(err_bc, "error")))


# =============================================================================
# SECTION 8B: GEM-X FLEX V2 HANDLE APPENDING
# =============================================================================
section("8b. finalise_probes_v2()")

V2_MULTIPLEX_TAIL  <- "CCCATATAAGAAA"   # 13 bp — standard v2 multiplex
V2_SINGLEPLEX_TAIL <- "CGGTCCTAGCAA"    # 12 bp — 4-sample singleplex kit

# Reuse test_probe, LHS_HANDLE, and marked_probe constants from section 8

# --- add_rhs_handle_v2(): multiplex ---
rhs_v2_multi <- tryCatch(add_rhs_handle_v2(test_probe, rhs_mode = "multiplex"),
                          error = function(e) NULL)

expect("add_rhs_handle_v2() multiplex starts with /5Phos/",
  quote(!is.null(rhs_v2_multi) && startsWith(rhs_v2_multi, "/5Phos/")))

expect("add_rhs_handle_v2() multiplex ends with CCCATATAAGAAA",
  quote(!is.null(rhs_v2_multi) && endsWith(rhs_v2_multi, V2_MULTIPLEX_TAIL)))

expect("add_rhs_handle_v2() multiplex has correct total length (/5Phos/ + 25 + 13 = 45 chars)",
  quote(!is.null(rhs_v2_multi) && nchar(rhs_v2_multi) == nchar("/5Phos/") + 25 + 13))

# --- add_rhs_handle_v2(): singleplex ---
rhs_v2_single <- tryCatch(add_rhs_handle_v2(test_probe, rhs_mode = "singleplex"),
                           error = function(e) NULL)

expect("add_rhs_handle_v2() singleplex ends with CGGTCCTAGCAA",
  quote(!is.null(rhs_v2_single) && endsWith(rhs_v2_single, V2_SINGLEPLEX_TAIL)))

expect("add_rhs_handle_v2() singleplex has correct total length (/5Phos/ + 25 + 12 = 44 chars)",
  quote(!is.null(rhs_v2_single) && nchar(rhs_v2_single) == nchar("/5Phos/") + 25 + 12))

# --- finalise_probes_v2(): multiplex ---
sel_v2 <- data.frame(GENE1 = "BCR", GENE2 = "ABL1", probe = test_probe,
                     stringsAsFactors = FALSE)
final_v2 <- tryCatch(finalise_probes_v2(sel_v2, rhs_mode = "multiplex"),
                     error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })

expect("finalise_probes_v2() returns a data frame",
  quote(is.data.frame(final_v2)))

expect("finalise_probes_v2() has expected columns (no Barcode columns)",
  quote(!is.null(final_v2) &&
        all(c("Fusion", "RHS_Mode", "LHS_Probe", "RHS_Probe") %in% colnames(final_v2)) &&
        !any(c("Barcode_ID", "Barcode_Seq", "Pool_Name") %in% colnames(final_v2))))

expect("finalise_probes_v2() LHS_Probe starts with the 21 bp constant handle",
  quote(!is.null(final_v2) && startsWith(final_v2$LHS_Probe[1], LHS_HANDLE)))

expect("finalise_probes_v2() LHS_Probe has correct length (46 bp)",
  quote(!is.null(final_v2) && nchar(final_v2$LHS_Probe[1]) == 46))

expect("finalise_probes_v2() multiplex RHS_Probe ends with CCCATATAAGAAA",
  quote(!is.null(final_v2) && endsWith(final_v2$RHS_Probe[1], V2_MULTIPLEX_TAIL)))

expect("finalise_probes_v2() RHS_Mode column records 'multiplex'",
  quote(!is.null(final_v2) && final_v2$RHS_Mode[1] == "multiplex"))

# --- finalise_probes_v2(): singleplex ---
final_v2_single <- tryCatch(finalise_probes_v2(sel_v2, rhs_mode = "singleplex"),
                             error = function(e) NULL)

expect("finalise_probes_v2() singleplex RHS_Probe ends with CGGTCCTAGCAA",
  quote(!is.null(final_v2_single) && endsWith(final_v2_single$RHS_Probe[1], V2_SINGLEPLEX_TAIL)))

# --- Marker stripping ---
sel_v2_marked <- data.frame(GENE1 = "X", GENE2 = "Y", probe = marked_probe,
                             stringsAsFactors = FALSE)
final_v2_marked <- tryCatch(finalise_probes_v2(sel_v2_marked), error = function(e) NULL)

expect("finalise_probes_v2() strips * and | markers before appending handles",
  quote(!is.null(final_v2_marked) && !grepl("[*|]", final_v2_marked$LHS_Probe[1])))

# --- Error cases ---
err_v2_cols <- tryCatch(
  finalise_probes_v2(data.frame(GENE1 = "A", GENE2 = "B", stringsAsFactors = FALSE)),
  error = function(e) e)
expect("finalise_probes_v2() raises an error when 'probe' column is missing",
  quote(inherits(err_v2_cols, "error")))

err_v2_mode <- tryCatch(finalise_probes_v2(sel_v2, rhs_mode = "badmode"),
                         error = function(e) e)
expect("finalise_probes_v2() raises an error for an invalid rhs_mode",
  quote(inherits(err_v2_mode, "error")))


# =============================================================================
# SECTION 8C: NON-FUSION HANDLE APPENDING (v1)
# =============================================================================
section("8c. finalise_nonfusion_probes()")

nf_sel_v1 <- data.frame(
  GENE    = "GFP",
  probe   = test_probe,
  Barcode = 1L,
  stringsAsFactors = FALSE
)

final_nf_v1 <- tryCatch(
  finalise_nonfusion_probes(nf_sel_v1),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)

expect("finalise_nonfusion_probes() returns a data frame",
  quote(is.data.frame(final_nf_v1)))

expect("Output has expected columns (Gene, Barcode_ID, Pool_Name, Barcode_Seq, LHS_Probe, RHS_Probe)",
  quote(!is.null(final_nf_v1) &&
        all(c("Gene", "Barcode_ID", "Pool_Name", "Barcode_Seq",
              "LHS_Probe", "RHS_Probe") %in% colnames(final_nf_v1))))

expect("Gene column contains the target gene name",
  quote(!is.null(final_nf_v1) && final_nf_v1$Gene[1] == "GFP"))

expect("Barcode_ID is 'BC001' for barcode 1",
  quote(!is.null(final_nf_v1) && final_nf_v1$Barcode_ID[1] == "BC001"))

expect("LHS_Probe starts with the 21 bp constant handle",
  quote(!is.null(final_nf_v1) && startsWith(final_nf_v1$LHS_Probe[1], LHS_HANDLE)))

expect("LHS_Probe has correct total length (21 + 25 = 46 bp)",
  quote(!is.null(final_nf_v1) && nchar(final_nf_v1$LHS_Probe[1]) == 46))

expect("RHS_Probe starts with /5Phos/",
  quote(!is.null(final_nf_v1) && startsWith(final_nf_v1$RHS_Probe[1], "/5Phos/")))

expect("RHS_Probe contains BC001 barcode sequence",
  quote(!is.null(final_nf_v1) && grepl(BC001_SEQ, final_nf_v1$RHS_Probe[1])))

expect("RHS_Probe ends with the 12 bp constant tail",
  quote(!is.null(final_nf_v1) && endsWith(final_nf_v1$RHS_Probe[1], RHS_TAIL)))

# | markers in probe should be stripped
nf_marked_probe <- paste0(strrep("A", 25), "|", strrep("T", 25))
nf_sel_marked <- data.frame(GENE = "GFP", probe = nf_marked_probe, Barcode = 1L,
                              stringsAsFactors = FALSE)
final_nf_marked <- tryCatch(finalise_nonfusion_probes(nf_sel_marked), error = function(e) NULL)
expect("| marker is stripped before handle appending",
  quote(!is.null(final_nf_marked) && !grepl("|", final_nf_marked$LHS_Probe[1], fixed = TRUE)))

# Missing GENE column raises an error
err_nf_gene <- tryCatch(
  finalise_nonfusion_probes(data.frame(probe = test_probe, Barcode = 1L,
                                        stringsAsFactors = FALSE)),
  error = function(e) e)
expect("Missing GENE column raises an error",
  quote(inherits(err_nf_gene, "error")))


# =============================================================================
# SECTION 8D: NON-FUSION HANDLE APPENDING (v2)
# =============================================================================
section("8d. finalise_nonfusion_probes_v2()")

nf_sel_v2 <- data.frame(
  GENE  = "GFP",
  probe = test_probe,
  stringsAsFactors = FALSE
)

final_nf_v2 <- tryCatch(
  finalise_nonfusion_probes_v2(nf_sel_v2, rhs_mode = "multiplex"),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)

expect("finalise_nonfusion_probes_v2() returns a data frame",
  quote(is.data.frame(final_nf_v2)))

expect("Output has expected columns (Gene, RHS_Mode, LHS_Probe, RHS_Probe) — no Barcode columns",
  quote(!is.null(final_nf_v2) &&
        all(c("Gene", "RHS_Mode", "LHS_Probe", "RHS_Probe") %in% colnames(final_nf_v2)) &&
        !any(c("Barcode_ID", "Barcode_Seq", "Pool_Name") %in% colnames(final_nf_v2))))

expect("Gene column is correct",
  quote(!is.null(final_nf_v2) && final_nf_v2$Gene[1] == "GFP"))

expect("LHS_Probe has correct length (46 bp)",
  quote(!is.null(final_nf_v2) && nchar(final_nf_v2$LHS_Probe[1]) == 46))

expect("Multiplex RHS_Probe ends with CCCATATAAGAAA",
  quote(!is.null(final_nf_v2) && endsWith(final_nf_v2$RHS_Probe[1], V2_MULTIPLEX_TAIL)))

expect("RHS_Mode column records 'multiplex'",
  quote(!is.null(final_nf_v2) && final_nf_v2$RHS_Mode[1] == "multiplex"))

final_nf_v2_single <- tryCatch(
  finalise_nonfusion_probes_v2(nf_sel_v2, rhs_mode = "singleplex"),
  error = function(e) NULL)
expect("Singleplex RHS_Probe ends with CGGTCCTAGCAA",
  quote(!is.null(final_nf_v2_single) && endsWith(final_nf_v2_single$RHS_Probe[1], V2_SINGLEPLEX_TAIL)))

err_nf_v2_mode <- tryCatch(finalise_nonfusion_probes_v2(nf_sel_v2, rhs_mode = "badmode"),
                             error = function(e) e)
expect("Invalid rhs_mode raises an error",
  quote(inherits(err_nf_v2_mode, "error")))

err_nf_v2_col <- tryCatch(
  finalise_nonfusion_probes_v2(data.frame(probe = test_probe, stringsAsFactors = FALSE)),
  error = function(e) e)
expect("Missing GENE column raises an error",
  quote(inherits(err_nf_v2_col, "error")))


# =============================================================================
# SECTION 9: END-TO-END PIPELINE
# =============================================================================
section("9. End-to-end: parse_arriba_tsv → create_probes_from_arriba → finalise_probes")

e2e_result <- tryCatch({
  # Step 1: parse
  parsed_e2e <- suppressWarnings(parse_arriba_tsv(tmp_tsv))  # reuse temp file from section 6

  # Step 2: design
  probes_e2e  <- create_probes_from_arriba(parsed_e2e, RESTRAINT_CONST = 5,
                                            ASTERIX_FLAG = TRUE, PROBE_HALVES_FLAG = TRUE)

  # Step 3: pick top probe per fusion; assign barcode 1
  selected_e2e <- probes_e2e %>%
    group_by(GENE1, GENE2) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(Barcode = 1L)

  # Step 4: finalise
  final_e2e <- finalise_probes(selected_e2e)
  final_e2e
}, error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })

expect("End-to-end pipeline returns a data frame",
  quote(is.data.frame(e2e_result)))

expect("End-to-end pipeline produces one row per fusion",
  quote(!is.null(e2e_result) && nrow(e2e_result) == 2))

expect("All LHS probes are the correct length (46 bp)",
  quote(!is.null(e2e_result) && all(nchar(e2e_result$LHS_Probe) == 46)))


# =============================================================================
# SECTION 10: BLAST OFF-TARGET CHECK (skipped if blastn not on PATH)
# =============================================================================
section("10. BLAST off-target check — check_blast_available()")

blast_available <- tryCatch({
  check_blast_available()
  TRUE
}, error = function(e) FALSE)

if (blast_available) {
  cat("  [INFO] blastn found on PATH — running BLAST tests.\n")
  # (BLAST tests require a real database; skipped here)
  cat("  [SKIP] Full BLAST test requires a pre-built database (--blast-db).\n")
} else {
  cat("  [SKIP] blastn not found on PATH — BLAST tests skipped.\n")
  cat("         Install BLAST+ and add it to PATH to enable these tests.\n")
}


# =============================================================================
# SECTION 11: NON-FUSION PROBE DESIGN
# =============================================================================
section("11. tile_sequence() and create_nonfusion_probes()")

nf_tiled <- tryCatch(
  tile_sequence(NONFUSION_GENE, NONFUSION_SEQ, PROBE_HALVES_FLAG = TRUE, MRNA_FLAG = FALSE),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)

expect("tile_sequence() returns a data frame",
  quote(is.data.frame(nf_tiled)))

expect("Output has GENE, probe, Ranking, Score columns",
  quote(!is.null(nf_tiled) &&
        all(c("GENE", "probe", "Ranking", "Score") %in% colnames(nf_tiled))))

expect("GENE column contains the correct gene name",
  quote(!is.null(nf_tiled) && all(nf_tiled$GENE == NONFUSION_GENE)))

expect("Ranking starts at 1 and is consecutive",
  quote(!is.null(nf_tiled) &&
        nf_tiled$Ranking[1] == 1 &&
        all(nf_tiled$Ranking == seq_len(nrow(nf_tiled)))))

expect("Scores are in descending order",
  quote(!is.null(nf_tiled) && all(diff(nf_tiled$Score) <= 0)))

expect("All scores are > 0 (filter removed zero-score probes)",
  quote(!is.null(nf_tiled) && all(nf_tiled$Score > 0)))

expect("Probes with PROBE_HALVES_FLAG contain '|' at position 26",
  quote(!is.null(nf_tiled) && all(substr(nf_tiled$probe, 26, 26) == "|")))

expect("Probe sequences (stripped of |) are exactly 50 bp",
  quote(!is.null(nf_tiled) &&
        all(nchar(gsub("[|]", "", nf_tiled$probe)) == 50)))

# Short sequence should warn and return empty data frame
nf_short <- tryCatch(
  suppressWarnings(tile_sequence("SHORT", "ATGCATGC", PROBE_HALVES_FLAG = FALSE)),
  error = function(e) NULL)
expect("tile_sequence() returns empty data frame for sequence < 50 bp",
  quote(!is.null(nf_short) && nrow(nf_short) == 0))

# MRNA_FLAG adds mrna column
nf_mrna <- tryCatch(
  tile_sequence(NONFUSION_GENE, NONFUSION_SEQ, PROBE_HALVES_FLAG = FALSE, MRNA_FLAG = TRUE),
  error = function(e) NULL)
expect("MRNA_FLAG = TRUE adds 'mrna' column",
  quote(!is.null(nf_mrna) && "mrna" %in% colnames(nf_mrna)))
expect("mrna column contains the RC of each probe (first row check)",
  quote(!is.null(nf_mrna) &&
        nf_mrna$mrna[1] == reverse_complement(nf_mrna$probe[1])))

# create_nonfusion_probes() — multi-gene input
nf_input_df <- data.frame(
  gene     = c("GFP", "CRISPR_TARGET"),
  sequence = c(NONFUSION_SEQ, NONFUSION_SEQ),
  stringsAsFactors = FALSE
)

nf_multi <- tryCatch(
  create_nonfusion_probes(nf_input_df, PROBE_HALVES_FLAG = FALSE),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)

expect("create_nonfusion_probes() returns a data frame",
  quote(is.data.frame(nf_multi)))

expect("Output contains rows for both input genes",
  quote(!is.null(nf_multi) &&
        all(c("GFP", "CRISPR_TARGET") %in% nf_multi$GENE)))

expect("Ranking resets to 1 for each gene",
  quote(!is.null(nf_multi) && {
    min_ranks <- tapply(nf_multi$Ranking, nf_multi$GENE, min)
    all(min_ranks == 1)
  }))

# Missing required columns should raise an error
err_nf_cols <- tryCatch(
  create_nonfusion_probes(data.frame(gene = "GFP", stringsAsFactors = FALSE)),
  error = function(e) e)
expect("Missing 'sequence' column raises an error",
  quote(inherits(err_nf_cols, "error")))


# =============================================================================
# SECTION 12: FLEX COMPETITION CHECK — NON-FUSION PROBES
# =============================================================================
section("12. check_flex_competition() — non-fusion probes")

# Probe that matches the synthetic standard probe (LHS = 25xA, RHS = 25xT)
nf_comp_match <- data.frame(
  GENE          = "MATCH",
  probe         = paste0(strrep("A", 25), strrep("T", 25)),  # identical to synth_probeset_df probe 1
  mRNA_position = 1L, Score = 5, Ranking = 1L,
  first_half_GC = 0, second_half_GC = 100,
  Dinucleotide = "AT", Dinucleotide_Status = "OK", Homopolymer_Flag = "OK",
  stringsAsFactors = FALSE
)

# Probe that does NOT match any standard probe (many mismatches on both halves)
nf_comp_pass <- data.frame(
  GENE          = "PASS",
  probe         = paste0(strrep("C", 25), strrep("G", 25)),  # C/G probe — differs from A/T and G/C std probes
  mRNA_position = 1L, Score = 5, Ranking = 1L,
  first_half_GC = 100, second_half_GC = 100,
  Dinucleotide = "CG", Dinucleotide_Status = "Warning", Homopolymer_Flag = "OK",
  stringsAsFactors = FALSE
)

nf_comp_df <- rbind(nf_comp_match, nf_comp_pass)
nf_comp_result <- tryCatch(
  check_flex_competition(nf_comp_df, synth_probeset_df, max_mismatches = 2L),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)

expect("check_flex_competition() returns a data frame",
  quote(is.data.frame(nf_comp_result)))

expect("Adds flex_lhs_min_mm, flex_rhs_min_mm, flex_competition, flex_pass columns",
  quote(!is.null(nf_comp_result) &&
        all(c("flex_lhs_min_mm", "flex_rhs_min_mm", "flex_competition", "flex_pass") %in%
            colnames(nf_comp_result))))

expect("Probe identical to a standard probe is flagged (flex_competition = TRUE)",
  quote(!is.null(nf_comp_result) &&
        isTRUE(nf_comp_result$flex_competition[nf_comp_result$GENE == "MATCH"])))

expect("Probe identical to standard probe does not pass (flex_pass = FALSE)",
  quote(!is.null(nf_comp_result) &&
        !isTRUE(nf_comp_result$flex_pass[nf_comp_result$GENE == "MATCH"])))

expect("Probe with many mismatches to all standard probes passes (flex_pass = TRUE)",
  quote(!is.null(nf_comp_result) &&
        isTRUE(nf_comp_result$flex_pass[nf_comp_result$GENE == "PASS"])))

expect("flex_pass and flex_competition are complementary",
  quote(!is.null(nf_comp_result) &&
        all(nf_comp_result$flex_pass == !nf_comp_result$flex_competition)))

expect("LHS min mismatches are non-negative integers",
  quote(!is.null(nf_comp_result) && all(nf_comp_result$flex_lhs_min_mm >= 0)))


# =============================================================================
# SECTION 13: FLEX COMPETITION CHECK — FUSION PROBES
# =============================================================================
section("13. check_flex_competition_fusion() — fusion probes")

# Fusion probes: displacement > 0 means junction falls in RIGHT half,
# so the LEFT half is the wild-type (non-junction) half that gets checked.

# Probe A: LHS = 25xA (matches synth standard probe 1 LHS) — should be FLAGGED
# Probe B: LHS = 25xC (no match to any standard probe) — should PASS
comp_fusion_df <- data.frame(
  GENE1 = c("BCR", "EML4"),
  GENE2 = c("ABL1", "ALK"),
  probe = c(
    paste0(strrep("A", 25), strrep("T", 25)),  # LHS matches std probe 1 → flagged
    paste0(strrep("C", 25), strrep("G", 25))   # LHS = 25xC → no match → pass
  ),
  fusion_point_displacement = c(5L, 5L),  # junction in right half → check left half
  Score = c(5, 5), Ranking = c(1L, 1L),
  first_half_GC = c(0, 100), second_half_GC = c(100, 100),
  Dinucleotide = c("AT", "CG"), Dinucleotide_Status = c("OK", "Warning"),
  Homopolymer_Flag = c("OK", "OK"),
  stringsAsFactors = FALSE
)

comp_fusion_result <- tryCatch(
  check_flex_competition_fusion(comp_fusion_df, synth_probeset_df, max_mismatches = 2L),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)

expect("check_flex_competition_fusion() returns a data frame",
  quote(is.data.frame(comp_fusion_result)))

expect("Adds flex_nonjunction_mm, flex_competition, flex_pass columns",
  quote(!is.null(comp_fusion_result) &&
        all(c("flex_nonjunction_mm", "flex_competition", "flex_pass") %in%
            colnames(comp_fusion_result))))

expect("Probe whose non-junction half matches a standard probe is flagged",
  quote(!is.null(comp_fusion_result) &&
        isTRUE(comp_fusion_result$flex_competition[comp_fusion_result$GENE1 == "BCR"])))

expect("Probe whose non-junction half differs from all standard probes passes",
  quote(!is.null(comp_fusion_result) &&
        isTRUE(comp_fusion_result$flex_pass[comp_fusion_result$GENE1 == "EML4"])))

expect("flex_pass and flex_competition are complementary",
  quote(!is.null(comp_fusion_result) &&
        all(comp_fusion_result$flex_pass == !comp_fusion_result$flex_competition)))

# displacement < 0 → junction in left half → check RIGHT half
comp_fusion_rhs_df <- data.frame(
  GENE1 = "BCR", GENE2 = "ABL1",
  probe = paste0(strrep("A", 25), strrep("T", 25)),  # RHS = 25xT, matches std probe 1 RHS
  fusion_point_displacement = -5L,  # junction in left → check right
  Score = 5, Ranking = 1L,
  first_half_GC = 0, second_half_GC = 100,
  Dinucleotide = "AT", Dinucleotide_Status = "OK", Homopolymer_Flag = "OK",
  stringsAsFactors = FALSE
)
comp_rhs_result <- tryCatch(
  check_flex_competition_fusion(comp_fusion_rhs_df, synth_probeset_df, max_mismatches = 2L),
  error = function(e) NULL)
expect("Negative displacement: right (wild-type) half is correctly checked",
  quote(!is.null(comp_rhs_result) && isTRUE(comp_rhs_result$flex_competition[1])))

# Missing fusion_point_displacement raises an error
err_no_fpd <- tryCatch(
  check_flex_competition_fusion(
    data.frame(GENE1 = "A", GENE2 = "B", probe = test_probe, stringsAsFactors = FALSE),
    synth_probeset_df
  ),
  error = function(e) e)
expect("Missing fusion_point_displacement column raises an error",
  quote(inherits(err_no_fpd, "error")))


# =============================================================================
# SECTION 14: is_offtarget_hit_nonfusion() HELPER
# =============================================================================
section("14. is_offtarget_hit_nonfusion()")

no_hits <- data.frame(qseqid=character(), sseqid=character(),
                      mismatch=integer(), length=integer(), qlen=integer(),
                      stringsAsFactors=FALSE)

one_subject_hits <- data.frame(
  qseqid   = c("probe_1_lhs", "probe_1_lhs"),
  sseqid   = c("transcript_1", "transcript_1"),  # same subject, two alignments
  mismatch = c(0L, 1L),
  length   = c(25L, 24L),
  qlen     = c(25L, 25L),
  stringsAsFactors = FALSE
)

two_subject_hits <- data.frame(
  qseqid   = c("probe_1_lhs", "probe_1_lhs"),
  sseqid   = c("transcript_1", "transcript_2"),  # two DIFFERENT subjects
  mismatch = c(0L, 1L),
  length   = c(25L, 25L),
  qlen     = c(25L, 25L),
  stringsAsFactors = FALSE
)

high_mm_hits <- data.frame(
  qseqid   = c("probe_1_lhs", "probe_1_lhs"),
  sseqid   = c("transcript_1", "transcript_2"),
  mismatch = c(10L, 12L),  # effective mm >> min_mismatches → not close hits
  length   = c(25L, 25L),
  qlen     = c(25L, 25L),
  stringsAsFactors = FALSE
)

expect("No hits at all → returns FALSE (pass)",
  quote(!isTRUE(is_offtarget_hit_nonfusion(no_hits, min_mismatches = 5))))

expect("Single unique subject with close hits → returns FALSE (on-target only)",
  quote(!isTRUE(is_offtarget_hit_nonfusion(one_subject_hits, min_mismatches = 5))))

expect("Two unique subjects with close hits → returns TRUE (off-target found)",
  quote(isTRUE(is_offtarget_hit_nonfusion(two_subject_hits, min_mismatches = 5))))

expect("Hits with high effective mismatches do not trigger flag",
  quote(!isTRUE(is_offtarget_hit_nonfusion(high_mm_hits, min_mismatches = 5))))


# =============================================================================
# SUMMARY
# =============================================================================
total <- .pass_count + .fail_count
cat(sprintf(
  "\n══════════════════════════════════════════════════════════════\n  Results: %d / %d tests passed",
  .pass_count, total
))

if (.fail_count == 0) {
  cat("  ✓ All tests passed!\n")
} else {
  cat(sprintf("  ✗ %d test(s) FAILED — see details above.\n", .fail_count))
}
cat("══════════════════════════════════════════════════════════════\n\n")

# Return non-zero exit code if any tests failed (useful for CI).
# Guard against interactive sessions (e.g. RStudio) where quit() kills the process.
if (.fail_count > 0 && !interactive()) quit(status = 1L)
