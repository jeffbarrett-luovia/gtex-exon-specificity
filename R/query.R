# R/query.R
# GTEx Portal API v2 query functions.
# Wraps gtexr package calls with consistent column naming and serial iteration.

# Map GTEx dataset IDs to their GENCODE versions
.dataset_gencode_version <- c(
  gtex_v7  = "v19",
  gtex_v8  = "v26",
  gtex_v10 = "v39"
)

#' Resolve gene symbols to versioned GENCODE IDs
#'
#' @param gene_symbols Character vector of HGNC gene symbols.
#' @param gencode_version GENCODE version string (e.g. "v39"). Defaults to "v39".
#' @return Tibble with columns: gene_symbol, gencode_id.
#'   Symbols that failed to resolve are dropped with a warning.
lookup_gencode_ids <- function(gene_symbols, gencode_version = "v39") {
  result <- gtexr::get_genes(geneIds = gene_symbols, gencodeVersion = gencode_version)

  resolved <- dplyr::select(result, gene_symbol = geneSymbol, gencode_id = gencodeId)

  missing <- setdiff(tolower(gene_symbols), tolower(resolved$gene_symbol))
  if (length(missing) > 0) {
    warning("Could not resolve gene symbol(s): ", paste(missing, collapse = ", "))
  }

  resolved
}

#' Fetch median exon expression across tissues for a set of GENCODE IDs
#'
#' Queries GTEx serially (one gene at a time) as required by the API.
#'
#' @param gencode_ids Character vector of versioned GENCODE IDs.
#' @param dataset_id GTEx dataset identifier. Default: "gtex_v8".
#' @param tissues Character vector of tissue site detail IDs to restrict to,
#'   or NULL (default) to fetch all 54 tissues.
#' @return Long-format tibble: one row per (exon x tissue).
#'   Columns: gencode_id, gene_symbol, exon_id, tissue, median, unit, dataset_id.
fetch_exon_expression <- function(gencode_ids,
                                  dataset_id = "gtex_v10",
                                  tissues = NULL) {
  results <- vector("list", length(gencode_ids))

  for (i in seq_along(gencode_ids)) {
    gid <- gencode_ids[[i]]
    message("  Fetching expression for ", gid, " (", i, "/", length(gencode_ids), ")")

    raw <- gtexr::get_median_exon_expression(
      gencodeIds          = gid,
      datasetId           = dataset_id,
      tissueSiteDetailIds = tissues,
      itemsPerPage        = 10000L,
      .return_raw         = TRUE
    )

    results[[i]] <- dplyr::bind_rows(raw$data) |>
      dplyr::transmute(
        gencode_id  = gencodeId,
        gene_symbol = geneSymbol,
        exon_id     = exonId,
        tissue      = tissueSiteDetailId,
        median      = median,
        unit        = unit,
        dataset_id  = datasetId
      )
  }

  dplyr::bind_rows(results)
}

#' Fetch collapsed gene model exon coordinates for a set of GENCODE IDs
#'
#' @param gencode_ids Character vector of versioned GENCODE IDs.
#' @return Tibble with columns: gencode_id, exon_id, exon_number,
#'   chromosome, start, end, strand.
fetch_exon_coordinates <- function(gencode_ids, dataset_id = "gtex_v10") {
  raw <- gtexr::get_collapsed_gene_model_exon(
    gencodeId    = gencode_ids,
    datasetId    = dataset_id,
    itemsPerPage = 10000L
  )

  dplyr::transmute(
    raw,
    gencode_id  = gencodeId,
    gene_symbol = geneSymbol,
    exon_id     = exonId,
    exon_number = as.integer(exonNumber),
    chromosome  = chromosome,
    start       = start,
    end         = end,
    strand      = strand
  )
}

#' Convenience wrapper: fetch expression + coordinates for a list of gene symbols
#'
#' Orchestrates lookup_gencode_ids(), fetch_exon_expression(), and
#' fetch_exon_coordinates(), then left-joins coordinates onto expression.
#'
# Tissues excluded by default: cell lines and testis
.default_exclude_tissues <- c(
  "Cells_EBV-transformed_lymphocytes",
  "Cells_Cultured_fibroblasts",
  "Testis"
)

#' @param gene_symbols Character vector of HGNC gene symbols.
#' @param dataset_id GTEx dataset identifier. Default: "gtex_v10".
#' @param tissues Character vector of tissue site detail IDs, or NULL for all.
#' @param exclude_tissues Tissues to drop before scoring. Defaults to cell lines
#'   and testis. Set to character(0) to keep all tissues.
#' @return Combined long-format tibble ready for score_exon_specificity().
fetch_tissue_exon_data <- function(gene_symbols,
                                   dataset_id = "gtex_v10",
                                   tissues = NULL,
                                   exclude_tissues = .default_exclude_tissues) {
  gencode_version <- .dataset_gencode_version[[dataset_id]]
  if (is.null(gencode_version)) {
    stop("Unknown dataset_id '", dataset_id, "'. Known: ",
         paste(names(.dataset_gencode_version), collapse = ", "))
  }

  message("Step 1/3: Resolving gene symbols (GENCODE ", gencode_version, ")...")
  ids <- lookup_gencode_ids(gene_symbols, gencode_version = gencode_version)

  message("Step 2/3: Fetching exon expression data...")
  expr <- fetch_exon_expression(
    gencode_ids = ids$gencode_id,
    dataset_id  = dataset_id,
    tissues     = tissues
  )

  message("Step 3/3: Fetching exon coordinates...")
  coords <- fetch_exon_coordinates(gencode_ids = ids$gencode_id, dataset_id = dataset_id)

  # Drop gene_symbol from coords to avoid duplication; it comes from expr
  coords_slim <- dplyr::select(coords, -gene_symbol)

  message("Fetching tissue-to-organ mapping...")
  tissue_map <- gtexr::get_tissue_site_detail() |>
    dplyr::select(tissue = tissueSiteDetailId, organ = tissueSite)

  message("Joining coordinates onto expression data...")
  result <- dplyr::left_join(expr, coords_slim, by = c("gencode_id", "exon_id")) |>
    dplyr::mutate(
      median = median / (end - start),
      unit   = "reads_per_base"
    ) |>
    dplyr::left_join(tissue_map, by = "tissue")

  if (length(exclude_tissues) > 0) {
    result <- dplyr::filter(result, !tissue %in% exclude_tissues)
  }

  result
}
