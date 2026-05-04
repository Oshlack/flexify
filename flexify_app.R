# =============================================================================
# Flexify — Shiny App
# =============================================================================
#
# A three-tab Shiny application for designing, checking, and finalising
# custom probes for the 10x Genomics Flex platform.
#
# TAB 1 — Design Probes
#   Mode selector: Fusion probes or Non-fusion (wild-type) probes.
#   Fusion: upload Arriba TSV or generic 4-column CSV, set design parameters,
#     generate ranked candidate 50 bp probe sequences spanning each junction.
#   Non-fusion: upload a gene/sequence CSV, generate tiled probes ranked by
#     GC content, ligation dinucleotide, and homopolymer content.
#
# TAB 2 — Off-Target & Competition Check (optional)
#   BLAST off-target check (fusion probes only): queries the junction-spanning
#     half of each fusion probe against a reference transcriptome database.
#   Flex competition check (both probe types): compares the non-junction half
#     of fusion probes, or both halves of non-fusion probes, against the
#     standard 10x Flex whole-transcriptome probe set.
#
# TAB 3 — Select & Finalise
#   Select one probe per fusion/gene using radio buttons (or upload a CSV).
#   Assign Probe Barcodes (v1) or omit them (v2). Generates full LHS and RHS
#   oligonucleotide sequences ready for synthesis ordering.
#   Works for both fusion and non-fusion probe workflows.
#
# USAGE:
#   Rscript -e "shiny::runApp('flexify_app.R')"
#   or open flexify_app.R in RStudio and click Run App.
#
# DEPENDENCIES:
#   shiny, tidyverse, DT, stringr
#   flexify_core.R, flexify_offtarget.R, flexify_handles.R,
#   flexify_nonfusion.R must be in the same directory as this file.
# =============================================================================

library(shiny)
library(tidyverse)
library(DT)
library(stringr)

source("flexify_core.R")
source("flexify_offtarget.R")
source("flexify_handles.R")
source("flexify_nonfusion.R")

# Paths to bundled 10x Genomics standard probeset files.
# Shiny sets the working directory to the app folder, so relative paths work.
BUNDLED_PROBESETS <- list(
  "v1" = "data/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv",
  "v2" = "data/Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv"
)


# =============================================================================
# SHARED STYLES
# =============================================================================

app_css <- "
  @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600&family=JetBrains+Mono:wght@400;500&display=swap');

  body {
    font-family: 'Inter', system-ui, -apple-system, sans-serif;
    font-size: 14px;
    line-height: 1.55;
    background-color: #f0eaff;
    color: #1a0a2e;
  }

  .navbar { background-color: #1a0a2e !important; border-bottom: 2px solid #e040fb; }
  .navbar-brand { color: #e040fb !important; font-weight: 600; letter-spacing: 0.03em; }
  .navbar-nav > li > a { color: #d4b8ff !important; }
  .navbar-nav > li.active > a { background-color: #e040fb !important; color: #fff !important; }
  .navbar-nav > li > a:hover { background-color: #9c27b0 !important; color: #fff !important; }

  .well {
    background-color: #ffffff;
    border: 1px solid #d4b8ff;
    border-radius: 8px;
    box-shadow: 0 1px 6px rgba(224, 64, 251, 0.07);
  }

  .btn-primary {
    background-color: #e040fb; border-color: #c61ddf; color: #fff; font-weight: 500;
  }
  .btn-primary:hover { background-color: #c61ddf; border-color: #a800c0; }
  .btn-success {
    background-color: #00bcd4; border-color: #0097a7; color: #fff; font-weight: 500;
  }
  .btn-success:hover { background-color: #0097a7; border-color: #00788a; }

  h3 {
    color: #6a0080;
    border-bottom: 2px solid #e040fb;
    padding-bottom: 6px;
    font-weight: 600;
    letter-spacing: -0.01em;
  }
  h4 { color: #4a0080; font-weight: 600; letter-spacing: -0.01em; }

  .top-probe-box {
    background: #1a0a2e;
    border-left: 4px solid #e040fb;
    padding: 12px 16px;
    margin-bottom: 10px;
    border-radius: 6px;
    font-family: 'JetBrains Mono', 'Courier New', monospace;
    font-size: 12.5px;
    color: #f0e6ff;
    line-height: 1.6;
  }

  .fusion-label {
    font-weight: 600; color: #6a0080; font-size: 14px;
    margin-bottom: 4px; letter-spacing: -0.01em;
  }

  .info-box {
    background: #f3e5ff; border-left: 4px solid #ce93d8;
    padding: 10px 14px; border-radius: 6px;
    font-size: 13px; margin-bottom: 14px; color: #3a006a;
  }

  .error-box {
    background: #fce4ec; border-left: 4px solid #f06292;
    padding: 10px 14px; border-radius: 6px;
    font-size: 13px; margin-top: 14px; color: #880e4f;
  }

  .success-box {
    background: #e0f7fa; border-left: 4px solid #00bcd4;
    padding: 10px 14px; border-radius: 6px;
    font-size: 13px; margin-bottom: 14px; color: #004d5e;
  }

  .fusion-panel {
    background: #ffffff; border: 1px solid #d4b8ff;
    border-radius: 8px; padding: 14px 18px; margin-bottom: 16px;
    box-shadow: 0 1px 6px rgba(224, 64, 251, 0.06);
  }
  .fusion-panel h4 { margin-top: 0; }

  .step-banner {
    background: #1a0a2e; color: #00e5ff;
    padding: 9px 16px; border-radius: 6px; margin-bottom: 16px;
    font-size: 13px; font-weight: 500;
    border-left: 3px solid #e040fb; letter-spacing: 0.01em;
  }
"


# =============================================================================
# USER INTERFACE
# =============================================================================

ui <- navbarPage(
  title  = "Flexify - Probe Designer",
  header = tags$head(
    tags$link(rel  = "stylesheet",
              href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600&family=JetBrains+Mono:wght@400;500&display=swap"),
    tags$style(HTML(app_css))
  ),

  # ---------------------------------------------------------------------------
  # TAB 1: DESIGN PROBES
  # ---------------------------------------------------------------------------
  tabPanel("Design Probes",
    sidebarLayout(
      sidebarPanel(
        width = 3,

        h3("Probe Type"),
        radioButtons("probe_design_mode", NULL,
                     choices = list(
                       "Fusion probes"              = "fusion",
                       "Non-fusion (wild-type) probes" = "nonfusion"
                     ),
                     selected = "fusion"),

        hr(),

        # ── Fusion inputs ─────────────────────────────────────────────────────
        conditionalPanel(
          condition = "input.probe_design_mode == 'fusion'",
          h3("Input"),
          radioButtons("input_mode", "Input format:",
                       choices = list(
                         "Arriba TSV (direct output)"      = "arriba",
                         "Generic CSV (any fusion caller)" = "generic"
                       ),
                       selected = "arriba"),
          conditionalPanel(
            condition = "input.input_mode == 'arriba'",
            div(class = "info-box",
                "Upload an Arriba fusion TSV file. Gene names and sequences are ",
                "extracted from the ", tags$code("fusion_transcript"), " column."
            ),
            fileInput("arriba_file", "Upload Arriba TSV",
                      accept = c(".tsv", ".txt"),
                      buttonLabel = "Browse...", placeholder = "No file selected")
          ),
          conditionalPanel(
            condition = "input.input_mode == 'generic'",
            div(class = "info-box",
                "Upload a CSV with columns: ",
                tags$code("gene1"), ", ", tags$code("gene2"), ", ",
                tags$code("gene1_transcript"),
                " (sequence ending at breakpoint), ",
                tags$code("gene2_transcript"),
                " (sequence starting at breakpoint)."
            ),
            fileInput("csv_file", "Upload fusion CSV",
                      accept = ".csv",
                      buttonLabel = "Browse...", placeholder = "No file selected")
          ),
          hr(),
          h3("Parameters"),
          sliderInput("restraint_const", "Min. bases per probe half",
                      min = 1, max = 12, value = 5, step = 1),
          tags$small(tags$i(
            "Minimum bases each fusion partner must contribute to a candidate probe."
          )),
          br(), br(),
          checkboxInput("asterix_flag",        "Mark fusion point with *",           value = TRUE),
          checkboxInput("probe_halves_flag",   "Mark probe halves with |",           value = TRUE),
          checkboxInput("mrna_flag",           "Include mRNA target sequence",       value = FALSE),
          checkboxInput("prioritise_rhs_flag", "Penalise left-half junction probes", value = FALSE)
        ),

        # ── Non-fusion inputs ─────────────────────────────────────────────────
        conditionalPanel(
          condition = "input.probe_design_mode == 'nonfusion'",
          h3("Input"),
          div(class = "info-box",
              "Upload a CSV with two columns: ",
              tags$code("gene"), " (target name) and ",
              tags$code("sequence"), " (mRNA sequence, 5'→3', at least 50 bp)."
          ),
          fileInput("nonfusion_csv", "Upload gene/sequence CSV",
                    accept = ".csv",
                    buttonLabel = "Browse...", placeholder = "No file selected"),
          hr(),
          h3("Parameters"),
          checkboxInput("nonfusion_halves_flag", "Mark probe halves with |",    value = TRUE),
          checkboxInput("nonfusion_mrna_flag",   "Include mRNA target sequence", value = FALSE)
        ),

        hr(),
        actionButton("run_btn", "Design Probes",
                     class = "btn-primary btn-block",
                     style = "width:100%; font-size:15px; padding:10px;"),
        br(), br(),
        conditionalPanel(
          condition = "output.tab1_ready",
          downloadButton("download_all_btn", "Download All Probes (CSV)",
                         class = "btn-success btn-block",
                         style = "width:100%; font-size:13px; padding:8px;"),
          br(), br(),
          downloadButton("download_selection_template_btn", "Download Selection Template",
                         class = "btn-success btn-block",
                         style = "width:100%; font-size:13px; padding:8px;"),
          tags$small(tags$i(
            "Adds 'Selected' and 'Barcode' columns. Fill these in and re-upload in Tab 3."
          ))
        )
      ),

      mainPanel(
        width = 9,
        conditionalPanel(
          condition = "!output.tab1_ready",
          div(style = "text-align:center; margin-top:80px; color:#888;",
              tags$h4("Select a probe type, upload your input, and click 'Design Probes'."),
              tags$p("Flexify will generate and rank candidate 50 bp probe sequences.")
          )
        ),
        conditionalPanel(
          condition = "output.tab1_ready",
          uiOutput("tab1_warnings_ui"),
          h3("Top-Ranked Probe per Target"),
          uiOutput("top_probes_ui"),
          br(),
          h3("All Candidate Probes"),
          div(style = "font-size:13px;", DTOutput("results_table"))
        ),
        uiOutput("tab1_error_ui")
      )
    )
  ),

  # ---------------------------------------------------------------------------
  # TAB 2: OFF-TARGET & COMPETITION CHECK
  # ---------------------------------------------------------------------------
  tabPanel("Off-Target & Competition Check",
    sidebarLayout(
      sidebarPanel(
        width = 3,

        # ── BLAST section (fusion probes only) ───────────────────────────────
        h3("BLAST Off-Target Check"),
        div(class = "info-box",
            tags$b("Fusion probes:"), " queries the junction-spanning half only.",
            br(),
            tags$b("Non-fusion probes:"), " queries each 25 bp half independently; ",
            "a half is flagged only if it closely matches more than one unique transcript ",
            "(the on-target hit is expected and does not count as a failure).",
            br(), br(),
            "Requires BLAST+ on PATH and a pre-built database:",
            br(),
            tags$code("makeblastdb -in transcriptome.fa -dbtype nucl -out transcriptome_db")
        ),
        radioButtons("blast_source", "Probe input:",
                     choices = list(
                       "Use results from Design tab" = "app",
                       "Upload probe CSV"            = "upload"
                     ),
                     selected = "app"),
        conditionalPanel(
          condition = "input.blast_source == 'upload'",
          fileInput("blast_csv_file", "Upload probe CSV",
                    accept = ".csv",
                    buttonLabel = "Browse...", placeholder = "No file selected")
        ),
        textInput("blast_db", "BLAST database path",
                  placeholder = "/path/to/transcriptome_db"),
        tags$small(tags$i("Path to BLAST database, without file extension.")),
        br(), br(),
        numericInput("min_mismatches", "Min. mismatches for off-target",
                     value = 5, min = 1, max = 25),
        numericInput("blast_threads", "BLAST threads",
                     value = 1, min = 1, max = 32),
        checkboxInput("filter_fails", "Remove probes that fail off-target check", value = TRUE),
        actionButton("blast_run_btn", "Run BLAST Check",
                     class = "btn-primary btn-block",
                     style = "width:100%; font-size:15px; padding:10px;"),
        br(), br(),
        conditionalPanel(
          condition = "output.tab2_ready",
          downloadButton("blast_download_btn", "Download BLAST-Filtered Probes (CSV)",
                         class = "btn-success btn-block",
                         style = "width:100%; font-size:13px; padding:8px;")
        ),

        hr(),

        # ── Flex Competition section (both probe types) ──────────────────────
        h3("Flex Competition Check"),
        div(class = "info-box",
            tags$b("Fusion probes:"), " checks the non-junction (wild-type) half against ",
            "the standard probeset.", br(),
            tags$b("Non-fusion probes:"), " checks both halves independently."
        ),
        radioButtons("comp_probeset_source", "Probeset source:",
                     choices = list(
                       "Use bundled 10x probeset" = "bundled",
                       "Upload custom probeset"   = "upload"
                     ),
                     selected = "bundled"),
        conditionalPanel(
          condition = "input.comp_probeset_source == 'bundled'",
          selectInput("comp_bundled_probeset", "Select probeset:",
                      choices = list(
                        "Chromium Flex v1 (GRCh38-2024-A)" = "v1",
                        "GEM-X Flex v2 (GRCh38-2024-A)"   = "v2"
                      ),
                      selected = "v1")
        ),
        conditionalPanel(
          condition = "input.comp_probeset_source == 'upload'",
          fileInput("comp_probeset_file", "Upload probeset CSV",
                    accept = ".csv",
                    buttonLabel = "Browse...", placeholder = "No file selected")
        ),
        numericInput("comp_max_mm", "Max mismatches to flag competition",
                     value = 2L, min = 0L, max = 10L),
        tags$small(tags$i(
          "Probes with ≤ this many mismatches to any standard probe half are flagged."
        )),
        br(), br(),
        actionButton("comp_run_btn", "Run Flex Competition Check",
                     class = "btn-primary btn-block",
                     style = "width:100%; font-size:15px; padding:10px;"),
        br(), br(),
        conditionalPanel(
          condition = "output.tab_competition_ready",
          downloadButton("comp_download_btn", "Download Competition Results (CSV)",
                         class = "btn-success btn-block",
                         style = "width:100%; font-size:13px; padding:8px;")
        )
      ),

      mainPanel(
        width = 9,
        conditionalPanel(
          condition = "!output.tab2_ready && !output.tab_competition_ready",
          div(style = "text-align:center; margin-top:80px; color:#888;",
              tags$h4("Run BLAST off-target check or Flex competition check."),
              tags$p("Design tab results are used automatically, or upload a probe CSV.")
          )
        ),

        # BLAST results
        conditionalPanel(
          condition = "output.tab2_ready",
          uiOutput("blast_summary_ui"),
          h3("BLAST-Filtered Probes"),
          div(style = "font-size:13px;", DTOutput("blast_results_table"))
        ),
        uiOutput("tab2_error_ui"),

        # Competition results
        conditionalPanel(
          condition = "output.tab_competition_ready",
          br(),
          h3("Flex Competition Check Results"),
          uiOutput("competition_summary_ui"),
          div(style = "font-size:13px;", DTOutput("competition_results_table"))
        ),
        uiOutput("tab_competition_error_ui")
      )
    )
  ),

  # ---------------------------------------------------------------------------
  # TAB 3: SELECT & FINALISE
  # ---------------------------------------------------------------------------
  tabPanel("Select & Finalise",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h3("Input"),
        radioButtons("finalise_source", "Probe source:",
                     choices = list(
                       "Select from app results" = "app",
                       "Upload selection CSV"    = "upload"
                     ),
                     selected = "app"),

        conditionalPanel(
          condition = "input.finalise_source == 'upload'",
          conditionalPanel(
            condition = "input.assay_version == 'v1' || !input.assay_version",
            div(class = "info-box",
                "Upload a CSV with columns: ",
                tags$code("GENE1"), " + ", tags$code("GENE2"),
                " (fusion) or ", tags$code("GENE"), " (non-fusion), ",
                tags$code("probe"), ", ", tags$code("Barcode"), " (integer 1–16).",
                br(), br(),
                "Use the 'Selection Template' from Tab 1 as a starting point."
            )
          ),
          conditionalPanel(
            condition = "input.assay_version == 'v2'",
            div(class = "info-box",
                "Upload a CSV with columns: ",
                tags$code("GENE1"), " + ", tags$code("GENE2"),
                " (fusion) or ", tags$code("GENE"), " (non-fusion), ",
                tags$code("probe"), ".",
                br(), br(),
                tags$b("GEM-X Flex v2:"), " no Barcode column required."
            )
          ),
          fileInput("finalise_csv_file", "Upload selection CSV",
                    accept = ".csv",
                    buttonLabel = "Browse...", placeholder = "No file selected")
        ),

        conditionalPanel(
          condition = "input.finalise_source == 'app'",
          uiOutput("app_results_status_ui")
        ),

        hr(),
        h3("Assay Version"),
        radioButtons("assay_version", NULL,
                     choices = list(
                       "Chromium Flex (v1)" = "v1",
                       "GEM-X Flex (v2)"   = "v2"
                     ),
                     selected = "v1"),

        conditionalPanel(
          condition = "input.assay_version == 'v2'",
          radioButtons("v2_rhs_mode", "RHS configuration:",
                       choices = list(
                         "Multiplex (CCCATATAAGAAA)"  = "multiplex",
                         "Singleplex (CGGTCCTAGCAA)"  = "singleplex"
                       ),
                       selected = "multiplex"),
          div(class = "info-box",
              tags$b("GEM-X Flex v2:"), " no barcode embedded in probe. ",
              "Barcode is supplied separately by kit reagents.",
              br(), br(),
              tags$b("Multiplex:"), " standard v2 tail (CCCATATAAGAAA).", br(),
              tags$b("Singleplex:"), " 4-sample kit tail (CGGTCCTAGCAA)."
          )
        ),

        hr(),
        actionButton("generate_btn", "Generate Final Probes",
                     class = "btn-primary btn-block",
                     style = "width:100%; font-size:15px; padding:10px;"),
        br(), br(),
        conditionalPanel(
          condition = "output.tab3_ready",
          downloadButton("final_download_btn", "Download Final Probes (CSV)",
                         class = "btn-success btn-block",
                         style = "width:100%; font-size:13px; padding:8px;")
        )
      ),

      mainPanel(
        width = 9,
        conditionalPanel(
          condition = "input.finalise_source == 'app'",
          uiOutput("selection_ui")
        ),
        conditionalPanel(
          condition = "input.finalise_source == 'upload'",
          uiOutput("upload_preview_ui")
        ),
        conditionalPanel(
          condition = "output.tab3_ready",
          br(),
          h3("Final Probe Sequences"),
          div(class = "info-box",
              "These sequences are ready for submission to an oligonucleotide synthesis provider. ",
              "Each target produces one LHS and one RHS probe oligonucleotide."
          ),
          div(style = "font-size:13px;", DTOutput("final_table"))
        ),
        uiOutput("tab3_error_ui")
      )
    )
  )
)


# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {

  # ---------------------------------------------------------------------------
  # Shared reactive state
  # ---------------------------------------------------------------------------

  # Combined design results (fusion or non-fusion depending on mode)
  design_results <- reactiveVal(NULL)

  # "fusion" or "nonfusion" — set whenever design is run
  probe_type     <- reactiveVal(NULL)

  # BLAST-filtered results (fusion only)
  blast_results  <- reactiveVal(NULL)

  # Flex competition check results (fusion or non-fusion)
  competition_results <- reactiveVal(NULL)

  # Synthesis-ready final output
  final_results  <- reactiveVal(NULL)

  # Per-tab error and warning accumulators
  tab1_error            <- reactiveVal(NULL)
  tab2_error            <- reactiveVal(NULL)
  tab3_error            <- reactiveVal(NULL)
  tab_competition_error <- reactiveVal(NULL)
  tab1_warnings         <- reactiveVal(NULL)

  # Three-tier priority: competition results > BLAST results > design results
  working_probes <- reactive({
    if (!is.null(competition_results())) competition_results()
    else if (!is.null(blast_results()))  blast_results()
    else                                 design_results()
  })


  # ---------------------------------------------------------------------------
  # TAB 1: Design Probes
  # ---------------------------------------------------------------------------

  observeEvent(input$run_btn, {
    design_results(NULL)
    probe_type(NULL)
    blast_results(NULL)
    competition_results(NULL)
    final_results(NULL)
    tab1_error(NULL)
    tab1_warnings(NULL)

    tryCatch({
      warnings_collected <- character(0)
      mode <- input$probe_design_mode

      if (mode == "fusion") {

        if (input$input_mode == "arriba") {
          req(input$arriba_file)
          cleaned_df <- withCallingHandlers(
            parse_arriba_tsv(input$arriba_file$datapath),
            warning = function(w) {
              warnings_collected <<- c(warnings_collected, conditionMessage(w))
              invokeRestart("muffleWarning")
            }
          )
        } else {
          req(input$csv_file)
          raw_df     <- read.csv(input$csv_file$datapath, stringsAsFactors = FALSE)
          cleaned_df <- process_arriba_transcript(raw_df)
        }

        output_df <- withCallingHandlers(
          create_probes_from_arriba(
            cleaned_df,
            RESTRAINT_CONST     = input$restraint_const,
            PRIORITISE_RHS_FLAG = input$prioritise_rhs_flag,
            ASTERIX_FLAG        = input$asterix_flag,
            PROBE_HALVES_FLAG   = input$probe_halves_flag,
            MRNA_FLAG           = input$mrna_flag
          ),
          warning = function(w) {
            warnings_collected <<- c(warnings_collected, conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        )

      } else {  # nonfusion

        req(input$nonfusion_csv)
        raw_df <- read.csv(input$nonfusion_csv$datapath, stringsAsFactors = FALSE)

        output_df <- withCallingHandlers(
          create_nonfusion_probes(
            raw_df,
            PROBE_HALVES_FLAG = input$nonfusion_halves_flag,
            MRNA_FLAG         = input$nonfusion_mrna_flag
          ),
          warning = function(w) {
            warnings_collected <<- c(warnings_collected, conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        )

        if (nrow(output_df) == 0) {
          stop("No probes could be designed. Check that sequences are at least 50 bp.")
        }
      }

      if (length(warnings_collected) > 0) tab1_warnings(warnings_collected)
      probe_type(mode)
      design_results(output_df)

    }, error = function(e) {
      tab1_error(conditionMessage(e))
    })
  })

  output$tab1_ready <- reactive({ !is.null(design_results()) })
  outputOptions(output, "tab1_ready", suspendWhenHidden = FALSE)

  # Top-ranked probe summary boxes — works for both fusion and non-fusion
  output$top_probes_ui <- renderUI({
    req(design_results())
    df    <- design_results()
    ptype <- probe_type()

    if (ptype == "fusion") {
      groups     <- unique(paste(df$GENE1, df$GENE2, sep = "::"))
      get_sub_df <- function(g) {
        parts <- strsplit(g, "::")[[1]]
        df[which(df$GENE1 == parts[1] & df$GENE2 == parts[2]), ]
      }
    } else {
      groups     <- unique(df$GENE)
      get_sub_df <- function(g) df[which(df$GENE == g), ]
    }

    boxes <- lapply(groups, function(g) {
      sub_df <- get_sub_df(g)
      top    <- sub_df[which.min(sub_df$Ranking), ]

      div(
        div(class = "fusion-label",
            paste0(g, "  —  ", nrow(sub_df), " candidate probe(s)")),
        div(class = "top-probe-box",
            tags$b("Top probe: "),    top$probe, br(),
            tags$b("Score: "),        round(as.numeric(top$Score), 3), "  |  ",
            tags$b("LGC: "),          round(as.numeric(top$first_half_GC), 1), "%  |  ",
            tags$b("RGC: "),          round(as.numeric(top$second_half_GC), 1), "%  |  ",
            tags$b("Dinucleotide: "), top$Dinucleotide, " (", top$Dinucleotide_Status, ")  |  ",
            tags$b("Homopolymer: "),  top$Homopolymer_Flag
        )
      )
    })
    do.call(tagList, boxes)
  })

  output$results_table <- renderDT({
    req(design_results())
    display_df <- apply(design_results(), 2, as.character) %>% as.data.frame()
    datatable(
      display_df,
      options = list(
        pageLength = 15, scrollX = TRUE, dom = "lfrtip",
        columnDefs = list(list(targets = "_all", className = "dt-left"))
      ),
      rownames = FALSE, filter = "top"
    )
  })

  output$download_all_btn <- downloadHandler(
    filename = function() paste0("flexify_all_probes_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content  = function(file) write.csv(design_results(), file, row.names = FALSE)
  )

  output$download_selection_template_btn <- downloadHandler(
    filename = function() paste0("flexify_selection_template_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content  = function(file) {
      df          <- design_results()
      df$Selected <- FALSE
      df$Barcode  <- NA_integer_
      write.csv(df, file, row.names = FALSE)
    }
  )

  output$tab1_error_ui <- renderUI({
    req(tab1_error())
    div(class = "error-box", tags$b("Error: "), tab1_error())
  })

  output$tab1_warnings_ui <- renderUI({
    req(tab1_warnings())
    msgs <- tab1_warnings()
    div(class = "info-box",
        tags$b(paste0("⚠  ", length(msgs), " warning(s) from input parsing:")),
        tags$ul(lapply(msgs, tags$li))
    )
  })


  # ---------------------------------------------------------------------------
  # TAB 2: BLAST Off-Target Check (fusion probes only)
  # ---------------------------------------------------------------------------

  observeEvent(input$blast_run_btn, {
    blast_results(NULL)
    competition_results(NULL)
    final_results(NULL)
    tab2_error(NULL)

    tryCatch({

      if (trimws(input$blast_db) == "") stop("Please enter the path to your BLAST database.")

      # Resolve probe input and probe type
      if (input$blast_source == "app") {
        probe_df <- design_results()
        if (is.null(probe_df)) {
          stop("No design results available. Run Tab 1 first, or switch to 'Upload probe CSV'.")
        }
        ptype <- probe_type()
      } else {
        req(input$blast_csv_file)
        probe_df <- read.csv(input$blast_csv_file$datapath, stringsAsFactors = FALSE)
        # Infer type from columns: non-fusion has GENE but not fusion_point_displacement
        ptype <- if ("fusion_point_displacement" %in% colnames(probe_df)) "fusion" else "nonfusion"
      }

      # Route to the appropriate BLAST function
      filtered_df <- if (!is.null(ptype) && ptype == "nonfusion") {
        run_offtarget_check_nonfusion(
          probe_df       = probe_df,
          blast_db       = trimws(input$blast_db),
          min_mismatches = input$min_mismatches,
          n_threads      = input$blast_threads,
          filter_fails   = input$filter_fails
        )
      } else {
        run_offtarget_check(
          probe_df       = probe_df,
          blast_db       = trimws(input$blast_db),
          min_mismatches = input$min_mismatches,
          n_threads      = input$blast_threads,
          filter_fails   = input$filter_fails
        )
      }

      blast_results(filtered_df)

    }, error = function(e) {
      tab2_error(conditionMessage(e))
    })
  })

  output$tab2_ready <- reactive({ !is.null(blast_results()) })
  outputOptions(output, "tab2_ready", suspendWhenHidden = FALSE)

  output$blast_summary_ui <- renderUI({
    req(blast_results())
    div(class = "success-box",
        tags$b("BLAST check complete. "),
        nrow(blast_results()), " probe(s) retained after off-target filtering."
    )
  })

  output$blast_results_table <- renderDT({
    req(blast_results())
    display_df <- apply(blast_results(), 2, as.character) %>% as.data.frame()
    datatable(
      display_df,
      options = list(
        pageLength = 15, scrollX = TRUE, dom = "lfrtip",
        columnDefs = list(list(targets = "_all", className = "dt-left"))
      ),
      rownames = FALSE, filter = "top"
    )
  })

  output$blast_download_btn <- downloadHandler(
    filename = function() paste0("flexify_blast_filtered_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content  = function(file) write.csv(blast_results(), file, row.names = FALSE)
  )

  output$tab2_error_ui <- renderUI({
    req(tab2_error())
    div(class = "error-box", tags$b("Error: "), tab2_error())
  })


  # ---------------------------------------------------------------------------
  # TAB 2: Flex Competition Check (both probe types)
  # ---------------------------------------------------------------------------

  observeEvent(input$comp_run_btn, {
    competition_results(NULL)
    tab_competition_error(NULL)

    tryCatch({
      # Use BLAST results if available; otherwise design results
      probe_df <- if (!is.null(blast_results())) blast_results()
                  else design_results()

      if (is.null(probe_df)) {
        stop("No probe results available. Run Tab 1 (Design Probes) first.")
      }

      ptype <- probe_type()

      # Load probeset
      probeset_df <- if (input$comp_probeset_source == "bundled") {
        ps_path <- BUNDLED_PROBESETS[[input$comp_bundled_probeset]]
        if (!file.exists(ps_path)) stop("Bundled probeset file not found: ", ps_path)
        load_flex_probeset(ps_path)
      } else {
        req(input$comp_probeset_file)
        load_flex_probeset(input$comp_probeset_file$datapath)
      }

      # Route to correct competition check function based on probe type
      result_df <- if (!is.null(ptype) && ptype == "nonfusion") {
        # Non-fusion: check both halves independently
        check_flex_competition(
          probe_df, probeset_df,
          max_mismatches = as.integer(input$comp_max_mm)
        )
      } else {
        # Fusion: check non-junction half only
        if (!"fusion_point_displacement" %in% colnames(probe_df)) {
          stop("Probe data does not contain a 'fusion_point_displacement' column. ",
               "Ensure probes were designed using Flexify's fusion design workflow.")
        }
        check_flex_competition_fusion(
          probe_df, probeset_df,
          max_mismatches = as.integer(input$comp_max_mm)
        )
      }

      competition_results(result_df)

    }, error = function(e) {
      tab_competition_error(conditionMessage(e))
      showNotification(conditionMessage(e), type = "error", duration = 10)
    })
  })

  output$tab_competition_ready <- reactive({ !is.null(competition_results()) })
  outputOptions(output, "tab_competition_ready", suspendWhenHidden = FALSE)

  output$competition_summary_ui <- renderUI({
    req(competition_results())
    df     <- competition_results()
    n_pass <- sum(df$flex_pass, na.rm = TRUE)
    n_fail <- nrow(df) - n_pass
    div(class = "success-box",
        tags$b("Competition check complete. "),
        n_pass, " probe(s) passed, ", n_fail,
        " flagged for potential competition with standard Flex probes."
    )
  })

  output$competition_results_table <- renderDT({
    req(competition_results())
    display_df <- apply(competition_results(), 2, as.character) %>% as.data.frame()
    datatable(
      display_df,
      options = list(
        pageLength = 15, scrollX = TRUE, dom = "lfrtip",
        columnDefs = list(list(targets = "_all", className = "dt-left"))
      ),
      rownames = FALSE, filter = "top"
    )
  })

  output$comp_download_btn <- downloadHandler(
    filename = function() paste0("flexify_competition_checked_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content  = function(file) write.csv(competition_results(), file, row.names = FALSE)
  )

  output$tab_competition_error_ui <- renderUI({
    req(tab_competition_error())
    div(class = "error-box", tags$b("Error: "), tab_competition_error())
  })


  # ---------------------------------------------------------------------------
  # TAB 3: Select & Finalise
  # ---------------------------------------------------------------------------

  # Status banner — three-tier description of what's loaded
  output$app_results_status_ui <- renderUI({
    if (!is.null(competition_results())) {
      df     <- competition_results()
      n_pass <- sum(df$flex_pass, na.rm = TRUE)
      div(class = "success-box",
          "Using competition-checked results (", n_pass, " of ", nrow(df),
          " probes passed).")
    } else if (!is.null(blast_results())) {
      div(class = "success-box",
          "Using BLAST-filtered results (", nrow(blast_results()), " probes). ",
          "No competition check applied.")
    } else if (!is.null(design_results())) {
      div(class = "info-box",
          "Using Design tab results (", nrow(design_results()), " probes). ",
          "No BLAST or competition check applied.")
    } else {
      div(class = "info-box",
          "No results available yet. Run Tab 1 (Design Probes) first.")
    }
  })

  # Dynamic selection UI — handles both fusion and non-fusion probe data frames
  output$selection_ui <- renderUI({
    df    <- working_probes()
    ptype <- probe_type()

    if (is.null(df)) {
      return(div(style = "text-align:center; margin-top:60px; color:#888;",
                 tags$h4("No probe results available."),
                 tags$p("Run Tab 1 to design probes, then return here to select and finalise.")))
    }

    is_v1      <- is.null(input$assay_version) || input$assay_version == "v1"
    has_blast  <- "offtarget_pass" %in% colnames(df)
    has_comp   <- "flex_pass"      %in% colnames(df)
    is_nonfusion <- !is.null(ptype) && ptype == "nonfusion"

    # Build group labels and sub-df accessor
    if (is_nonfusion) {
      groups     <- unique(df$GENE)
      get_sub    <- function(g) df[which(df$GENE == g), , drop = FALSE]
    } else {
      groups     <- unique(paste(df$GENE1, df$GENE2, sep = "::"))
      get_sub    <- function(g) {
        parts <- strsplit(g, "::")[[1]]
        df[which(df$GENE1 == parts[1] & df$GENE2 == parts[2]), , drop = FALSE]
      }
    }

    panels <- lapply(groups, function(g) {
      sub_df  <- get_sub(g)
      sub_df  <- sub_df[order(as.numeric(sub_df$Ranking)), ]
      safe_id <- gsub("[^A-Za-z0-9]", "_", g)

      choices <- setNames(
        as.character(seq_len(nrow(sub_df))),
        sapply(seq_len(nrow(sub_df)), function(i) {
          row       <- sub_df[i, ]
          blast_str <- if (has_blast) {
            if (isTRUE(as.logical(row$offtarget_pass))) " | BLAST: pass" else " | BLAST: fail"
          } else ""
          comp_str <- if (has_comp) {
            if (isTRUE(as.logical(row$flex_pass))) " | Flex: pass" else " | Flex: FLAGGED"
          } else ""
          paste0(
            "Rank ", row$Ranking,
            " | Score: ", round(as.numeric(row$Score), 2),
            " | ", row$probe,
            " | LGC: ", round(as.numeric(row$first_half_GC), 1), "%",
            " RGC: ", round(as.numeric(row$second_half_GC), 1), "%",
            " | ", row$Dinucleotide, " (", row$Dinucleotide_Status, ")",
            " | ", row$Homopolymer_Flag,
            blast_str, comp_str
          )
        })
      )

      barcode_col <- if (is_v1) {
        column(3,
               selectInput(
                 inputId  = paste0("barcode_sel_", safe_id),
                 label    = "Assign barcode:",
                 choices  = setNames(as.character(1:16), names(PROBE_BARCODES)),
                 selected = "1"
               )
        )
      } else NULL

      div(class = "fusion-panel",
          h4(g),
          fluidRow(
            column(if (is_v1) 9L else 12L,
                   radioButtons(
                     inputId  = paste0("probe_sel_", safe_id),
                     label    = "Select probe:",
                     choices  = choices,
                     selected = "1"
                   )
            ),
            barcode_col
          )
      )
    })

    is_v1_banner <- is.null(input$assay_version) || input$assay_version == "v1"
    banner_text  <- if (is_v1_banner) {
      "Select one probe per target and assign a Probe Barcode, then click 'Generate Final Probes'."
    } else {
      "Select one probe per target, then click 'Generate Final Probes'. No barcode needed for GEM-X Flex v2."
    }

    tagList(
      div(class = "step-banner", banner_text),
      do.call(tagList, panels)
    )
  })

  # Upload preview
  output$upload_preview_ui <- renderUI({
    req(input$finalise_csv_file)
    tryCatch({
      df <- read.csv(input$finalise_csv_file$datapath, stringsAsFactors = FALSE)
      selected_df <- if ("Selected" %in% colnames(df)) {
        df[which(df$Selected == TRUE | df$Selected == "TRUE"), ]
      } else df

      # Detect type from column names
      is_nonfusion_upload <- "GENE" %in% colnames(selected_df) &&
                             !("GENE1" %in% colnames(selected_df))
      group_col    <- if (is_nonfusion_upload) "GENE" else c("GENE1", "GENE2")
      n_targets    <- if (is_nonfusion_upload) {
        length(unique(selected_df$GENE))
      } else {
        length(unique(paste(selected_df$GENE1, selected_df$GENE2)))
      }

      tagList(
        div(class = "success-box",
            tags$b("File loaded: "), nrow(selected_df),
            " probe(s) selected across ", n_targets, " target(s)."
        ),
        h4("Preview (first 10 rows)"),
        div(style = "font-size:12px; overflow-x:auto;",
            tableOutput("upload_preview_table"))
      )
    }, error = function(e) {
      div(class = "error-box", tags$b("Error reading file: "), conditionMessage(e))
    })
  })

  output$upload_preview_table <- renderTable({
    req(input$finalise_csv_file)
    df <- read.csv(input$finalise_csv_file$datapath, stringsAsFactors = FALSE)
    if ("Selected" %in% colnames(df)) df <- df[which(df$Selected == TRUE | df$Selected == "TRUE"), ]
    head(df, 10)
  })

  # Generate Final Probes
  observeEvent(input$generate_btn, {
    final_results(NULL)
    tab3_error(NULL)

    tryCatch({
      assay_ver    <- if (is.null(input$assay_version)) "v1" else input$assay_version
      is_v2        <- assay_ver == "v2"
      rhs_mode     <- if (is.null(input$v2_rhs_mode)) "multiplex" else input$v2_rhs_mode
      ptype        <- probe_type()

      if (input$finalise_source == "app") {
        df <- working_probes()
        if (is.null(df)) stop("No probe results available in the app. Run Tab 1 first.")

        is_nonfusion <- !is.null(ptype) && ptype == "nonfusion"

        if (is_nonfusion) {
          groups  <- unique(df$GENE)
          get_sub <- function(g) df[which(df$GENE == g), , drop = FALSE]
        } else {
          groups  <- unique(paste(df$GENE1, df$GENE2, sep = "::"))
          get_sub <- function(g) {
            parts <- strsplit(g, "::")[[1]]
            df[which(df$GENE1 == parts[1] & df$GENE2 == parts[2]), , drop = FALSE]
          }
        }

        selected_rows <- lapply(groups, function(g) {
          sub_df  <- get_sub(g)
          sub_df  <- sub_df[order(as.numeric(sub_df$Ranking)), ]
          safe_id <- gsub("[^A-Za-z0-9]", "_", g)

          sel_idx <- as.integer(input[[paste0("probe_sel_", safe_id)]])
          if (length(sel_idx) == 0 || is.na(sel_idx) || sel_idx > nrow(sub_df)) sel_idx <- 1L

          row <- sub_df[sel_idx, , drop = FALSE]

          if (!is_v2) {
            barcode_idx <- as.integer(input[[paste0("barcode_sel_", safe_id)]])
            if (length(barcode_idx) == 0 || is.na(barcode_idx)) barcode_idx <- 1L
            row$Barcode <- barcode_idx
          }
          row
        })

        selected_df <- do.call(rbind, selected_rows)

      } else {
        # Upload path
        req(input$finalise_csv_file)
        selected_df <- read.csv(input$finalise_csv_file$datapath, stringsAsFactors = FALSE)
        if ("Selected" %in% colnames(selected_df)) {
          selected_df <- selected_df[which(selected_df$Selected == TRUE |
                                           selected_df$Selected == "TRUE"), ]
        }
        if (nrow(selected_df) == 0) {
          stop("No rows with Selected == TRUE found in the uploaded file.")
        }

        # Infer probe type from column names
        is_nonfusion <- "GENE" %in% colnames(selected_df) &&
                        !("GENE1" %in% colnames(selected_df))
      }

      # Call appropriate finalise function based on probe type and assay version
      final_df <- if (is_nonfusion && is_v2) {
        finalise_nonfusion_probes_v2(selected_df, rhs_mode = rhs_mode)
      } else if (is_nonfusion) {
        finalise_nonfusion_probes(selected_df)
      } else if (is_v2) {
        finalise_probes_v2(selected_df, rhs_mode = rhs_mode)
      } else {
        finalise_probes(selected_df)
      }

      final_results(final_df)
      showNotification(
        paste0("✓ Final probes generated for ", nrow(final_df),
               " target(s). Scroll down to see the table."),
        type = "message", duration = 6
      )

    }, error = function(e) {
      tab3_error(conditionMessage(e))
      showNotification(conditionMessage(e), type = "error", duration = 10)
    })
  })

  output$tab3_ready <- reactive({ !is.null(final_results()) })
  outputOptions(output, "tab3_ready", suspendWhenHidden = FALSE)

  output$final_table <- renderDT({
    req(final_results())
    datatable(
      final_results(),
      options = list(
        pageLength = 20, scrollX = TRUE, dom = "lfrtip",
        columnDefs = list(list(targets = "_all", className = "dt-left"))
      ),
      rownames = FALSE
    )
  })

  output$final_download_btn <- downloadHandler(
    filename = function() paste0("flexify_final_probes_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content  = function(file) write.csv(final_results(), file, row.names = FALSE)
  )

  output$tab3_error_ui <- renderUI({
    req(tab3_error())
    div(class = "error-box", tags$b("Error: "), tab3_error())
  })

}


# =============================================================================
# LAUNCH
# =============================================================================

shinyApp(ui = ui, server = server)
