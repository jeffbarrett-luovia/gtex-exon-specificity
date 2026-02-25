# analysis.R
# GTEx Tissue-Specific Exon Analysis — entry point and example workflow.
#
# Usage: source("analysis.R") from this directory, or run interactively.
# Requires an internet connection to reach the GTEx Portal API.

# ── 1. Dependencies ──────────────────────────────────────────────────────────

if (!requireNamespace("gtexr", quietly = TRUE)) {
  install.packages("gtexr")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}

library(gtexr)
library(dplyr)
library(tidyr)

# ── 2. Source project functions ───────────────────────────────────────────────

source("R/query.R")
source("R/specificity.R")

# ── 3. Example: tissue-specific genes ────────────────────────────────────────
# TNNI3 — cardiac troponin I, heart-specific
# MBP   — myelin basic protein, brain-specific
# ACTB  — beta-actin, housekeeping (expect low Tau)
# GAPDH — glycolysis enzyme, housekeeping (expect low Tau)

genes_of_interest <- c("TNNI3", "MBP", "ACTB", "GAPDH")

message("=== Fetching GTEx exon data for: ", paste(genes_of_interest, collapse = ", "))

long_data <- fetch_tissue_exon_data(
  gene_symbols = genes_of_interest,
  dataset_id   = "gtex_v10",
  tissues      = NULL          # NULL = all GTEx v10 tissues
)

message("=== Computing tissue specificity scores...")

exon_scores <- score_exon_specificity(long_data)

# ── 4. Inspect results ────────────────────────────────────────────────────────

message("\n--- Top 10 most tissue-specific exons (Tau > 0.8) ---")
specific_exons <- exon_scores |>
  filter(tau > 0.8) |>
  arrange(desc(tau))

print(specific_exons, n = 10)

message("\n--- Per-gene Tau summary ---")
gene_summary <- exon_scores |>
  group_by(gene_symbol) |>
  summarise(
    n_exons      = n(),
    mean_tau     = round(mean(tau, na.rm = TRUE), 3),
    max_tau      = round(max(tau, na.rm = TRUE), 3),
    median_tau   = round(median(tau, na.rm = TRUE), 3),
    .groups = "drop"
  ) |>
  arrange(desc(mean_tau))

print(gene_summary)

# ── 5. Verification checks ────────────────────────────────────────────────────

message("\n--- Verification ---")

# 5a. TNNI3 and MBP should have high Tau
high_tau_genes <- gene_summary |>
  filter(gene_symbol %in% c("TNNI3", "MBP")) |>
  pull(mean_tau)

if (all(high_tau_genes > 0.8)) {
  message("PASS: TNNI3 and MBP have mean Tau > 0.8 (tissue-specific as expected)")
} else {
  message("WARN: TNNI3/MBP Tau lower than expected — check API response")
}

# 5b. ACTB and GAPDH should have low Tau
low_tau_genes <- gene_summary |>
  filter(gene_symbol %in% c("ACTB", "GAPDH")) |>
  pull(mean_tau)

if (all(low_tau_genes < 0.3)) {
  message("PASS: ACTB and GAPDH have mean Tau < 0.3 (ubiquitous as expected)")
} else {
  message("WARN: ACTB/GAPDH Tau higher than expected — review results")
}

# 5c. Coordinate join: all exons should have chromosome populated
missing_coords <- sum(is.na(exon_scores$chromosome))
if (missing_coords == 0) {
  message("PASS: All exons have chromosome coordinates populated")
} else {
  message("WARN: ", missing_coords, " exon(s) missing coordinates — inspect join")
}

# ── 6. Save results ───────────────────────────────────────────────────────────

output_file <- "exon_specificity_scores.tsv"
readr_available <- requireNamespace("readr", quietly = TRUE)

if (readr_available) {
  readr::write_tsv(exon_scores, output_file)
} else {
  write.table(exon_scores, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
}

message("\nResults written to: ", output_file)
message("Rows: ", nrow(exon_scores), "  Columns: ", ncol(exon_scores))
