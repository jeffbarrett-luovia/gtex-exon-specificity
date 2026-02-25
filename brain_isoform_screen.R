# brain_isoform_screen.R
#
# Identifies genes that are brain-dominant and have at least one exon with
# differential expression between Brain and Muscle (either direction).
#
# Rationale: a muscle-enriched exon flags a peripheral isoform; an antibody
# targeting an exon absent from that isoform would be selective for
# brain-derived protein in blood (cf. MAPT exon 8).

suppressPackageStartupMessages({
  library(gtexr)
  library(dplyr)
  library(tidyr)
})

source("R/query.R")
source("R/specificity.R")

# ── Parameters ────────────────────────────────────────────────────────────────

# Set to a positive integer to run on only the first N genes (for testing).
# Set to NULL (or Inf) for the full run.
TRIAL_N <- NULL

# Minimum reads-per-base in an exon's max tissue to treat it as expressed.
# Used for cross-gene comparisons (expr_fraction is not suitable here as its
# denominator scales with the most-expressed gene in the dataset).
RPB_FLOOR <- 0.5

# Minimum fraction of expressed exons with Brain as max organ
# for a gene to pass the brain-dominance filter
BRAIN_DOMINANCE_THRESHOLD <- 0.5

# log2(brain_rpb / muscle_rpb) threshold for calling an exon
# brain- or muscle-enriched. log2(2) = 1 means 2-fold minimum.
ENRICH_LOG2FC <- 1.0

# Small pseudocount added to RPB before log2FC to handle zeros
PSEUDOCOUNT <- 0.001

# Seconds to sleep between gene fetches — keeps load on GTEx API reasonable
SLEEP_BETWEEN_GENES <- 0.5

# ── Gene list: all well-characterised protein-coding genes ───────────────────
# Use HPA "Evidence at protein level" — ~18,500 genes with confirmed protein
# expression, covering the well-characterised portion of the proteome.
# The brain-dominance filter downstream selects the relevant subset.

message("Fetching gene list from Human Protein Atlas...")

hpa_all <- read.delim(
  paste0("https://www.proteinatlas.org/api/search_download.php",
         "?search=&columns=g,pe&compress=no&format=tsv"),
  stringsAsFactors = FALSE
)

hpa_genes <- hpa_all |>
  dplyr::filter(Evidence == "Evidence at protein level") |>
  dplyr::pull(Gene)

message("Genes with protein evidence: ", length(hpa_genes))

genes <- hpa_genes
if (!is.null(TRIAL_N) && is.finite(TRIAL_N)) {
  genes <- head(genes, TRIAL_N)
  message("TRIAL MODE: using first ", length(genes), " genes")
}

# ── Fetch data (with per-gene checkpointing) ─────────────────────────────────
# Each gene is saved to cache/ as it completes. Re-running the script skips
# already-fetched genes, so a failure can be resumed without starting over.

CACHE_DIR <- "cache/brain_screen"
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

message("=== Fetching GTEx v10 exon data for ", length(genes), " genes ===")
message("Cache: ", CACHE_DIR)

long_data_list <- vector("list", length(genes))

for (i in seq_along(genes)) {
  gene <- genes[[i]]
  cache_file <- file.path(CACHE_DIR, paste0(gene, ".rds"))

  if (file.exists(cache_file)) {
    long_data_list[[i]] <- readRDS(cache_file)
    next
  }

  message("[", i, "/", length(genes), "] Fetching ", gene, "...")
  result <- tryCatch(
    fetch_tissue_exon_data(gene),
    error = function(e) {
      message("  WARN: failed for ", gene, " — ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(result) && nrow(result) > 0) {
    saveRDS(result, cache_file)
    long_data_list[[i]] <- result
  }

  Sys.sleep(SLEEP_BETWEEN_GENES)
}

long_data <- dplyr::bind_rows(long_data_list)
message("Total rows fetched: ", nrow(long_data))

# ── Score ─────────────────────────────────────────────────────────────────────

message("\n=== Scoring exon specificity (organ level) ===")
scores <- score_exon_specificity(long_data, level = "organ")

# ── Organ expression pivot ────────────────────────────────────────────────────

organ_expr <- pivot_organ_expression(long_data)

# Join Brain and Muscle RPB onto scores
# (use any_of in case a gene has zero expression in one organ)
exon_detail <- scores |>
  dplyr::left_join(
    organ_expr |> dplyr::select(gencode_id, exon_id,
                                brain_rpb  = dplyr::any_of("Brain"),
                                muscle_rpb = dplyr::any_of("Muscle")),
    by = c("gencode_id", "exon_id")
  ) |>
  dplyr::mutate(
    brain_rpb  = tidyr::replace_na(brain_rpb, 0),
    muscle_rpb = tidyr::replace_na(muscle_rpb, 0),
    brain_muscle_log2fc = log2((brain_rpb + PSEUDOCOUNT) /
                               (muscle_rpb + PSEUDOCOUNT))
  )

# ── Gene-level filter: brain-dominant ─────────────────────────────────────────

gene_summary <- exon_detail |>
  dplyr::group_by(gene_symbol) |>
  dplyr::summarise(
    n_expressed       = sum(max_median > RPB_FLOOR, na.rm = TRUE),
    n_brain_max       = sum(max_tissue == "Brain" &
                            max_median > RPB_FLOOR, na.rm = TRUE),
    pct_brain         = n_brain_max / pmax(n_expressed, 1),
    n_brain_enriched  = sum(max_median > RPB_FLOOR &
                            brain_muscle_log2fc >  ENRICH_LOG2FC, na.rm = TRUE),
    n_muscle_enriched = sum(max_median > RPB_FLOOR &
                            brain_muscle_log2fc < -ENRICH_LOG2FC, na.rm = TRUE),
    .groups = "drop"
  )

brain_dominant_genes <- gene_summary |>
  dplyr::filter(pct_brain >= BRAIN_DOMINANCE_THRESHOLD)

message("\nGenes passing brain-dominance filter (>= ",
        round(BRAIN_DOMINANCE_THRESHOLD * 100), "% of expressed exons peak in Brain):")
print(brain_dominant_genes)

# ── Exon classification ───────────────────────────────────────────────────────

exon_classified <- exon_detail |>
  dplyr::filter(gene_symbol %in% brain_dominant_genes$gene_symbol) |>
  dplyr::mutate(
    exon_class = dplyr::case_when(
      max_median < RPB_FLOOR                            ~ "low_expression",
      brain_muscle_log2fc >  ENRICH_LOG2FC              ~ "brain_enriched",
      brain_muscle_log2fc < -ENRICH_LOG2FC              ~ "muscle_enriched",
      TRUE                                              ~ "shared"
    )
  ) |>
  dplyr::left_join(
    gene_summary |> dplyr::select(gene_symbol, pct_brain,
                                  n_brain_enriched, n_muscle_enriched),
    by = "gene_symbol"
  )

# ── Results ───────────────────────────────────────────────────────────────────

message("\n=== Exon classification for brain-dominant genes ===")

exon_classified |>
  dplyr::select(gene_symbol, exon_number, chromosome, start, end,
                brain_rpb, muscle_rpb, brain_muscle_log2fc,
                exon_class, expr_fraction, tau) |>
  print(n = Inf, width = 120)

message("\n=== Genes with differential brain/muscle exons ===")

gene_contrast_summary <- exon_classified |>
  dplyr::filter(exon_class %in% c("brain_enriched", "muscle_enriched")) |>
  dplyr::count(gene_symbol, exon_class) |>
  tidyr::pivot_wider(names_from = exon_class, values_from = n, values_fill = 0) |>
  dplyr::left_join(
    gene_summary |> dplyr::select(gene_symbol, n_expressed, pct_brain),
    by = "gene_symbol"
  ) |>
  dplyr::arrange(dplyr::desc(pct_brain))

print(gene_contrast_summary, n = Inf)

# ── Save outputs ──────────────────────────────────────────────────────────────

write.table(exon_classified,       "brain_screen_exons.tsv",  sep = "\t", row.names = FALSE, quote = FALSE)
write.table(gene_contrast_summary, "brain_screen_genes.tsv",  sep = "\t", row.names = FALSE, quote = FALSE)
message("\nSaved: brain_screen_exons.tsv, brain_screen_genes.tsv")
