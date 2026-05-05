# Flexify

An R tool for automated probe design for the 10x Genomics Flex platform. Flexify supports two probe design workflows:

- **Fusion probes** — 50 bp probes spanning a fusion gene junction, scored and ranked across all possible junction offsets.
- **Non-fusion probes** — 50 bp antisense probes tiled along a wild-type transcript sequence (e.g. GFP, CRISPR reporter, or any exogenous/custom gene).

Both workflows support **Chromium Flex (v1)** and **GEM-X Flex (v2)** assay formats, automated BLAST-based off-target screening, Flex probeset competition checking, and generate full synthesis-ready LHS and RHS probe sequences including the required 10x Genomics handle sequences.

---

## Quick Start (Shiny App)

This section walks you through launching the Flexify app with no command-line experience required.

### Step 1 — Install R

Download and install R from [https://cran.r-project.org](https://cran.r-project.org). Choose the version for your operating system (Windows, macOS, or Linux) and follow the installer prompts.

### Step 2 — Install RStudio

RStudio is a user-friendly interface for running R. Download the free Desktop version from [https://posit.co/download/rstudio-desktop](https://posit.co/download/rstudio-desktop) and install it.

### Step 3 — Install required R packages

Open RStudio. In the **Console** panel at the bottom, paste the following and press Enter:

```r
install.packages(c("shiny", "tidyverse", "DT", "stringr", "optparse"))
```

This only needs to be done once. It may take a few minutes to complete.

### Step 4 — Download Flexify

Click the green **Code** button at the top of this page and select **Download ZIP**. Unzip the downloaded folder somewhere convenient on your computer.

### Step 5 — Open and run the app

1. In RStudio, go to **File → Open File** and navigate to the unzipped Flexify folder.
2. Open `flexify_app.R`.
3. Click the **Run App** button that appears in the top-right corner of the editor panel.

The Flexify app will open in a new window (or in your browser). You can now upload your input file and start designing probes.

> **Note:** The off-target BLAST check (Tab 2) requires additional software (BLAST+) and a pre-built reference transcriptome database. This step is optional — see the [Installation](#installation) section below if you need it. Most users can skip it and proceed directly from Tab 1 to Tab 3.

---

## Installation

### R dependencies

```r
install.packages(c("shiny", "tidyverse", "DT", "stringr", "optparse"))
```

### BLAST+ (required for off-target checking only)

The easiest way to install BLAST+ is via conda or mamba, which handles all platforms automatically:

```bash
# With conda:
conda install -c bioconda blast

# With mamba (faster):
mamba install -c bioconda blast
```

Alternatively, platform-specific binaries are available from NCBI:
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

Verify installation:
```bash
blastn -version
```

Build a nucleotide database from your reference transcriptome FASTA:
```bash
makeblastdb -in transcriptome.fa -dbtype nucl -out transcriptome_db
```

---

## File Structure

```
Flexify/
├── flexify_app.R          # Shiny app (main entry point for interactive use)
├── flexify_cli.R          # Command-line interface
├── flexify_core.R         # Core fusion probe design functions
├── flexify_nonfusion.R    # Non-fusion probe design and competition check functions
├── flexify_offtarget.R    # BLAST off-target checking functions (fusion and non-fusion)
├── flexify_handles.R      # Handle/barcode appending functions (fusion and non-fusion)
└── README.md
```

All `.R` files must be in the same directory.

---

## Input Formats

### Fusion probes

Flexify supports two fusion input formats in both the Shiny app and the CLI.

#### Option 1: Arriba TSV (direct output)

Upload the TSV file produced directly by [Arriba](https://github.com/suhrig/arriba). Flexify parses the `fusion_transcript` column automatically, extracting gene names and splitting the sequence at the `|` breakpoint marker. Rows where the fusion transcript is absent or lacks a `|` separator are skipped with a warning.

#### Option 2: Generic CSV (any fusion caller)

For output from any other fusion caller (e.g. STAR-Fusion, FusionCatcher), prepare a CSV with the following four columns:

| Column | Description |
|---|---|
| `gene1` | Name of the first fusion partner gene |
| `gene2` | Name of the second fusion partner gene |
| `gene1_transcript` | mRNA sequence of gene1 ending at the breakpoint (5'→3', at least 25 bp) |
| `gene2_transcript` | mRNA sequence of gene2 starting at the breakpoint (5'→3', at least 25 bp) |

Column names are case-insensitive. Sequences should be in DNA alphabet (A/T/G/C). Both sequences must contribute at least 25 nucleotides to allow probe enumeration across the full junction offset range.

**Example:**

```
gene1,gene2,gene1_transcript,gene2_transcript
BCR,ABL1,ATGCGT...CCAGTA,TTAGCC...GAATTC
EML4,ALK,CCGTAA...TTAGCA,AACGGT...CCTGAA
```

### Non-fusion probes

Prepare a CSV with the following two columns:

| Column | Description |
|---|---|
| `gene` | Gene or construct name (used to label outputs) |
| `sequence` | Full mRNA/transcript sequence in DNA alphabet (A/T/G/C); at least 50 bp |

Column names are case-insensitive. Probes are generated as 50 bp reverse-complement windows tiled across the sequence. Windows with GC content outside 44–72% on either half are excluded.

**Example:**

```
gene,sequence
GFP,ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAAC...
TagBFP,ATGGTGAGCGTGATCAAGCCCGACATGCCCATCGTGGAGGGCGGCATGGACGGCTACGTGCTGGAGCCCTTC...
```

---

## Shiny App

Launch the interactive app from RStudio or the command line:

```r
# From RStudio: open flexify_app.R and click 'Run App'

# From the command line:
Rscript -e "shiny::runApp('flexify_app.R')"
```

The app has three tabs:

### Tab 1 — Design Probes

Select the probe design mode at the top of the sidebar:

- **Fusion probes**: Upload an Arriba TSV or a generic 4-column CSV and configure design parameters (restraint, asterix/halves/mRNA markers, penalise left-half junction).
- **Non-fusion probes**: Upload a 2-column gene/sequence CSV and configure output options (halves marker, mRNA sequence).

Click **Design Probes** to run. Results are shown as a summary box (top-ranked probe per fusion or per gene) and a full sortable/filterable table. Two downloads are available:

- **Download All Probes**: complete ranked output as CSV.
- **Download Selection Template**: same output with empty `Selected` (FALSE) and `Barcode` (NA) columns — fill these in externally and re-upload in Tab 3.

#### Fusion design parameters

- **Min. bases per probe half**: the minimum number of bases each fusion partner must contribute to a candidate probe (default: 5).
- **Mark fusion point with \***: inserts an asterisk at the exact fusion breakpoint position.
- **Mark probe halves with |**: inserts a pipe character at the LHS/RHS boundary (position 25|26).
- **Include mRNA target sequence**: appends the reverse complement of the probe.
- **Penalise left-half junction probes**: applies a 0.7× score multiplier to probes where the junction falls in the left half.

#### Non-fusion design parameters

- **Mark probe halves with |**: inserts a pipe character at the LHS/RHS boundary.
- **Include mRNA target sequence**: appends the mRNA window the probe targets.

### Tab 2 — Off-Target & Competition Check

Both checks are optional and can be run in either order after Tab 1.

#### BLAST off-target check

Screens probe sequences against a reference transcriptome database using BLAST.

**Fusion probes**: only the junction-spanning half (the half crossing the fusion breakpoint) is queried. The non-junction half binds a wild-type sequence by design and is not screened. A probe is flagged as off-target if its junction half has a hit with fewer than `min_mismatches` effective mismatches to any transcript.

**Non-fusion probes**: both the LHS (bases 1–25) and RHS (bases 26–50) halves are queried independently. Because each half will always have a perfect match to its intended target gene, a half is flagged as off-target only if close hits (fewer than `min_mismatches` effective mismatches) are found to **more than one unique transcript**. This flags probes where either half may bind unintended targets while accepting the expected on-target hit.

Effective mismatches = `aligned_mismatches + (query_length − alignment_length)`. Provide the BLAST database path (the path prefix used with `makeblastdb`, without file extension).

#### Flex probeset competition check

Checks whether any probe half closely matches a sequence already present in the standard 10x Genomics Flex whole-transcriptome probeset, which could compete for the same target and reduce signal.

**Fusion probes**: only the non-junction (wild-type) half is checked, since the junction half spans a novel sequence not present in the standard probeset.

**Non-fusion probes**: both halves are checked independently, since both bind wild-type sequence.

A probe half is flagged as competing if its Hamming distance to any standard probeset sequence is ≤ `max_mismatches` (default: 2). Provide the bundled probeset CSV or select from the pre-loaded v1/v2 probesets.

Results from both checks are carried forward to Tab 3.

### Tab 3 — Select & Finalise

**Select your assay version** (v1 or v2) in the sidebar before generating final probes.

**Fusion probes**: each fusion (GENE1::GENE2) is shown as a panel with radio buttons listing all ranked candidate probes and a barcode dropdown (v1 only). Select one probe per fusion, assign barcodes (v1), then click **Generate Final Probes**.

**Non-fusion probes**: each gene is shown as a panel with radio buttons. Select one probe per gene, assign barcodes (v1 only), then click **Generate Final Probes**.

**CSV upload path**: upload a CSV with columns `GENE1`/`GENE2` (fusion) or `GENE` (non-fusion), `probe`, and `Barcode` (v1 only, integer 1–16 or string BC001–BC016). If a `Selected` column is present, only rows with `Selected == TRUE` are processed.

The output table contains the full **LHS** and **RHS** probe sequences ready for submission to an oligonucleotide synthesis provider.

#### Chromium Flex (v1)

Barcode-embedded format. Each RHS probe encodes one of 16 Probe Barcodes (BC001–BC016). The barcode must match the corresponding whole transcriptome probe in the hybridisation mix for that sample.

#### GEM-X Flex (v2)

Barcoding is handled by kit reagents and is **not** embedded in the custom probe sequence. Select the RHS configuration:

- **Multiplex (CCCATATAAGAAA)**: standard v2 tail for multiplexed experiments.
- **Singleplex (CGGTCCTAGCAA)**: tail for the 4-sample singleplex kit.

---

## Command-Line Interface

### Mode 1: Design Probes (default)

```bash
# Fusion probes from generic CSV:
Rscript flexify_cli.R --input fusions.csv --output probes.csv

# Fusion probes from Arriba TSV:
Rscript flexify_cli.R --arriba --input fusions.tsv --output probes.csv

# With optional flags:
Rscript flexify_cli.R \
  --input fusions.csv \
  --output probes.csv \
  --restraint 5 \
  --mrna \
  --prioritise-rhs
```

**Flags:**

| Flag | Default | Description |
|---|---|---|
| `--input / -i` | required | Input fusion CSV or Arriba TSV (with `--arriba`) |
| `--output / -o` | required | Output ranked probe CSV |
| `--arriba` | FALSE | Parse input as an Arriba TSV file |
| `--restraint` | 5 | Minimum bases per probe half from each gene |
| `--mrna` | FALSE | Include mRNA target sequence column |
| `--prioritise-rhs` | FALSE | Penalise left-half junction probes |
| `--no-asterix` | FALSE | Omit fusion breakpoint marker (*) |
| `--no-halves` | FALSE | Omit probe half boundary marker (|) |

### Mode 2: BLAST Off-Target Check

```bash
Rscript flexify_cli.R \
  --mode blast \
  --input probes.csv \
  --blast-db /path/to/transcriptome_db \
  --output probes_filtered.csv

# Keep failed probes in output:
Rscript flexify_cli.R --mode blast -i probes.csv \
  --blast-db /path/db --output probes_flagged.csv --keep-fails
```

**Additional flags:**

| Flag | Default | Description |
|---|---|---|
| `--blast-db` | required | Path to BLAST database (no extension) |
| `--min-mismatches` | 5 | Minimum effective mismatches for a hit to be considered off-target |
| `--threads` | 1 | Number of BLAST threads |
| `--keep-fails` | FALSE | Retain failed probes with flag columns rather than removing |

### Mode 3: Finalise Probes

```bash
# Chromium Flex v1 (default) — barcode embedded in probe:
Rscript flexify_cli.R \
  --mode finalise \
  --input selected_probes.csv \
  --output final_probes.csv

# GEM-X Flex v2 multiplex — no barcode in probe:
Rscript flexify_cli.R \
  --mode finalise \
  --assay-version v2 \
  --input selected_probes.csv \
  --output final_probes.csv

# GEM-X Flex v2 singleplex (4-sample kit):
Rscript flexify_cli.R \
  --mode finalise \
  --assay-version v2 --rhs-mode singleplex \
  --input selected_probes.csv \
  --output final_probes.csv
```

**v1:** Input CSV must contain `GENE1`, `GENE2`, `probe`, and `Barcode` (integer 1–16 or string BC001–BC016).

**v2:** Input CSV requires only `GENE1`, `GENE2`, and `probe` — no `Barcode` column.

If a `Selected` column is present in either case, only rows marked `TRUE` are processed.

**Additional flags:**

| Flag | Default | Description |
|---|---|---|
| `--assay-version` | v1 | Assay version: `v1` (Chromium Flex) or `v2` (GEM-X Flex) |
| `--rhs-mode` | multiplex | v2 RHS tail: `multiplex` (CCCATATAAGAAA) or `singleplex` (CGGTCCTAGCAA) |

---

## Scoring System

### Fusion probes

Each candidate probe is scored on four criteria. The composite score is the product of all four components; any zero score excludes the probe from the output.

| Criterion | Method | Score range |
|---|---|---|
| GC content (LHS half) | Tiered, based on 10x guidelines | 0–5 (0 if outside 44–72%) |
| GC content (RHS half) | Tiered, based on 10x guidelines | 0–5 (0 if outside 44–72%) |
| Junction position | Gaussian, optimum at 12.5 bp from centre of each half | 0–5 |
| Ligation dinucleotide | Preferred set: AT, CA, CT, TA, TC, TG, TT | 1 or 3 |
| Homopolymer runs | Penalised for runs of 4+ identical bases | 1–2 |

### Non-fusion probes

Non-fusion probes are scored on GC content and homopolymer content only (there is no junction position score, since all probes tile uniformly across the transcript). Probes with GC content outside 44–72% on either half receive a score of zero and are excluded.

---

## Probe Structure

The LHS handle sequence is identical for both assay versions and both probe types. The RHS structure depends on assay version.

### LHS probe (v1 and v2)
```
CCTTGGCACCCGAGAATTCCA  [21 bp constant handle]
+ [bases 1–25 of the 50 bp probe]
```
Total length: 46 bp.

### RHS probe — Chromium Flex (v1)
```
/5Phos/                [5-prime phosphorylation]
+ [bases 26–50 of the 50 bp probe]
+ ACGCGGTTAGCACGTA     [16 bp linker / Constant Sequence]
+ NN                   [2 bp spacer]
+ [8 bp Probe Barcode] [BC001–BC016, unique per sample pool]
+ CGGTCCTAGCAA         [12 bp constant tail]
```
Total length: 70 characters (excluding `/5Phos/`).

### RHS probe — GEM-X Flex v2 (multiplex)
```
/5Phos/                [5-prime phosphorylation]
+ [bases 26–50 of the 50 bp probe]
+ CCCATATAAGAAA        [13 bp constant tail]
```
Total length: 38 characters (excluding `/5Phos/`). No barcode in probe.

### RHS probe — GEM-X Flex v2 (singleplex, 4-sample kit)
```
/5Phos/                [5-prime phosphorylation]
+ [bases 26–50 of the 50 bp probe]
+ CGGTCCTAGCAA         [12 bp constant tail]
```
Total length: 37 characters (excluding `/5Phos/`). No barcode in probe.

### Probe Barcodes (v1 only)

| Barcode | Sequence | Pool |
|---|---|---|
| BC001 | ACTTTAGG | poolOne |
| BC002 | AACGGGAA | poolTwo |
| BC003 | AGTAGGCT | poolThree |
| BC004 | ATGTTGAC | poolFour |
| BC005 | ACAGACCT | poolFive |
| BC006 | ATCCCAAC | poolSix |
| BC007 | AAGTAGAG | poolSeven |
| BC008 | AGCTGTGA | poolEight |
| BC009 | ACAGTCTG | poolNine |
| BC010 | AGTGAGTG | poolTen |
| BC011 | AGAGGCAA | poolEleven |
| BC012 | ACTACTCA | poolTwelve |
| BC013 | ATACGTCA | poolThirteen |
| BC014 | ATCATGTG | poolFourteen |
| BC015 | AACGCCGA | poolFifteen |
| BC016 | ATTCGGTT | poolSixteen |

Each v1 RHS probe must use the same barcode as the corresponding whole transcriptome probe in the hybridisation mix for that sample.

---

## Off-Target Assessment

Off-target specificity is assessed by BLAST alignment against the reference transcriptome. The strategy differs by probe type.

**Fusion probes**: only the junction-spanning half is screened. The non-junction half binds a wild-type sequence by design and is expected to have a perfect match to one of the fusion partner genes. A probe is flagged as off-target if the junction half has a hit with fewer than `min_mismatches` effective mismatches to any transcript.

**Non-fusion probes**: both halves (LHS: bases 1–25; RHS: bases 26–50) are screened independently. Because each half will always have a 0-mismatch hit to its intended target gene, a half is only flagged if close hits are found to more than one unique transcript. A probe passes if both halves pass. Note: multiple isoforms of the same gene may generate multiple hits that are each acceptable; the check is intentionally conservative and targets truly off-target binding to unrelated transcripts.

Effective mismatches = `aligned_mismatches + (query_length − alignment_length)`.

This can be run:
- **In the app**: using Tab 2 with a locally installed BLAST+ and pre-built database.
- **Via the CLI**: using `--mode blast`.
- **Manually**: by exporting the probe CSV from Tab 1, running BLAST externally, and filtering before re-importing in Tab 3.

---

## Citation

If you use Flexify in your research, please cite:

> [Citation to be added upon publication]

---

## Contact

For issues or questions, please open an issue on GitHub:
https://github.com/victoriasc/FLEX-fusion-probe-design
