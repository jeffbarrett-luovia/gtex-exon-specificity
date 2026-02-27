# R/specificity.R
# Tissue specificity scoring for exon-level GTEx expression data.

#' Compute Tau tissue specificity score
#'
#' Tau ranges from 0 (ubiquitously expressed) to 1 (tissue-specific).
#' Formula: tau = sum(1 - x_hat) / (n - 1), where x_hat = x / max(x).
#'
#' Reference: Yanai et al. (2005) Bioinformatics 21(5):650-659.
#'
#' @param x Numeric vector of expression values across tissues (non-negative).
#' @return Scalar Tau value, or NA if all values are zero.
compute_tau <- function(x) {
  if (all(x == 0, na.rm = TRUE)) return(NA_real_)

  x_max <- max(x, na.rm = TRUE)
  x_hat <- x / x_max
  n <- length(x)

  sum(1 - x_hat, na.rm = TRUE) / (n - 1)
}

#' Score tissue specificity for every exon in a long-format expression tibble
#'
#' @param long_df Long-format tibble as returned by fetch_tissue_exon_data().
#'   Required columns: gencode_id, gene_symbol, exon_id, tissue, median,
#'   and optionally chromosome, start, end, strand, exon_number, organ.
#' @param level "organ" (default) collapses sub-tissues to their broad organ
#'   group (max RPB per organ) before computing Tau, correcting for unequal
#'   tissue representation (e.g. 13 brain regions vs 1 heart). "tissue"
#'   computes Tau on all individual tissues as-is.
#' @return One-row-per-exon summary tibble with Tau and related metrics.
score_exon_specificity <- function(long_df, level = c("organ", "tissue")) {
  level <- match.arg(level)

  # Coordinate columns that may be present
  coord_cols <- c("chromosome", "start", "end", "strand", "exon_number")
  has_coords  <- all(coord_cols %in% names(long_df))

  # Collapse to organ level if requested
  scoring_df <- if (level == "organ" && "organ" %in% names(long_df)) {
    long_df |>
      dplyr::group_by(gencode_id, gene_symbol, exon_id, organ) |>
      dplyr::summarise(median = max(median, na.rm = TRUE), .groups = "drop") |>
      dplyr::rename(tissue = organ)
  } else {
    long_df
  }

  # Core per-exon summary
  scored <- scoring_df |>
    dplyr::group_by(gencode_id, gene_symbol, exon_id) |>
    dplyr::summarise(
      tau                = compute_tau(median),
      max_tissue         = tissue[which.max(median)],
      max_median         = max(median, na.rm = TRUE),
      mean_all           = mean(median, na.rm = TRUE),
      mean_other         = {
        max_idx <- which.max(median)
        others  <- median[-max_idx]
        if (length(others) == 0) NA_real_ else mean(others, na.rm = TRUE)
      },
      n_tissues          = dplyr::n(),
      n_tissues_expressed = sum(median > 0, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      fold_change = dplyr::if_else(mean_other == 0, NA_real_, max_median / mean_other)
    )

  # expr_fraction: exon's peak expression relative to the gene's top exon
  scored <- scored |>
    dplyr::group_by(gene_symbol) |>
    dplyr::mutate(
      expr_fraction = max_median / max(max_median, na.rm = TRUE)
    ) |>
    dplyr::ungroup()

  if (!has_coords) {
    return(scored)
  }

  # Pull one coordinate row per exon (they're constant across tissues)
  coords <- long_df |>
    dplyr::select(gencode_id, exon_id, dplyr::all_of(coord_cols)) |>
    dplyr::distinct()

  dplyr::left_join(scored, coords, by = c("gencode_id", "exon_id")) |>
    dplyr::select(
      gene_symbol, gencode_id, exon_id, exon_number,
      chromosome, start, end, strand,
      n_tissues, n_tissues_expressed,
      max_tissue, max_median, expr_fraction, mean_all, mean_other, fold_change,
      tau
    ) |>
    dplyr::arrange(gene_symbol, exon_number)
}

#' Pivot long expression data to one column per organ
#'
#' Returns max RPB per organ per exon, suitable for pairwise organ contrasts.
#' Organ names are normalised (spaces → underscores) for clean column names.
#'
#' @param long_df Long-format tibble as returned by fetch_tissue_exon_data().
#'   Must contain an 'organ' column.
#' @return Wide tibble: one row per (gencode_id, gene_symbol, exon_id),
#'   one column per organ containing max RPB across its sub-tissues.
pivot_organ_expression <- function(long_df) {
  if (!"organ" %in% names(long_df)) {
    stop("long_df must contain an 'organ' column; use fetch_tissue_exon_data().")
  }

  long_df |>
    dplyr::mutate(organ = gsub(" ", "_", organ)) |>
    dplyr::group_by(gencode_id, gene_symbol, exon_id, organ) |>
    dplyr::summarise(max_rpb = max(median, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from  = organ,
      values_from = max_rpb,
      values_fill = 0
    )
}
