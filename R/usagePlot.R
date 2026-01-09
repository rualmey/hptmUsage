#' Generate Usage Plots
#'
#' @description
#' This function generates line plots showing the usage of histone variants,
#' precursors, or PTMs across samples. It visualizes the relevant usage
#' normalization factors, corrected usages, and model estimates.
#'
#' @param object A `QFeatures` object containing the histone data.
#' @param features A named `list()` of feature names (hPTM or variant group) to
#'   plot, where names should correspond to the assay names in `object`.
#' @param include_precursors Logical. Whether to generate lineplots for all
#'   precursors corresponding to hPTMs or variants to plot.
#' @param design_formula A formula used for the linear model. Defaults to `~group`.
#'   Used to reconstruct model estimates.
#' @param show_legend Logical. Whether to show the legend. Defaults to `TRUE`.
#' @returns A list of `ggplot` objects.
#' @export
#' @examples
#' # TODO
usagePlot <- function(
  object,
  features,
  include_precursors = TRUE,
  design_formula = ~group,
  show_legend = TRUE
) {
  stopifnot("'object' must be a QFeatures object." = is(object, "QFeatures"))

  # retrieve relevant precursors
  precursor_features <- purrr::imap(features, \(x, idx) .get_precursor_features(object, idx, x))
  # reduce to only non-duplicate assays and features
  precursor_features <- tibble::tibble(
    assay = purrr::map_chr(precursor_features, "assay"),
    feats = purrr::map(precursor_features, "feats")
  ) |>
    dplyr::group_by(.data[["assay"]]) |>
    dplyr::summarise(
      feats = list(unique(unlist(feats))),
      .groups = "drop"
    ) |>
    tibble::deframe()
  features <- c(features, precursor_features)

  plot_list <- purrr::imap(features, function(feats, assay) {
    # retrieve normalization factors
    normalization_factors <- .fetch_normalization_factors(object, assay, feats)
    level <- .detect_feature_level(object, assay)

    plots <- vector("list", length(feats))
    names(plots) <- feats
    for (i in seq_along(feats)) {
      target <- feats[i]
      normalization <- normalization_factors$scaling_factors[i, ]

      # TODO HERE
    }

    # Filter out NULLs if any skips happened
    # TODO needed?
    plots <- plots[!purrr::map_lgl(plots, is.null)]
  })

  return(plot_list)
}

# ------------------------------------------------------------------------------
# Internal Helper Functions
# ------------------------------------------------------------------------------

.get_precursor_features <- function(object, assay_name, features) {
  link <- QFeatures::assayLink(object, assay_name)
  precursor_assay <- link@from
  link_col <- link@fcol

  rd <- SummarizedExperiment::rowData(object[[precursor_assay]])
  feat_mask <- rd[[link_col]] %in% features
  precursor_features <- SummarizedExperiment::rowData(object[[precursor_assay]][feat_mask])$precursor |>
    unique()

  # PTM assays have a deconvoluted assay "in-between" precursor and hPTM
  if (stringr::str_ends(precursor_assay, "Deconv")) {
    precursor_assay <- QFeatures::assayLink(object, precursor_assay)@from
  }

  return(list(
    assay = precursor_assay,
    feats = precursor_features
  ))
}

.fetch_normalization_factors <- function(object, assay, feats) {
  metadata <- metadata(object[[assay]])
  if (length(metadata) == 0) {
    # PTM or variant assay, so need to look back
    metadata <- metadata(object[[QFeatures::assayLink(object, assay)@from]])
  }

  factor_name <- rowData(object[[assay]][feats])[[metadata$usage_level]]
  if (metadata$usage_level == "histone") {
    factor_name <- as.integer(factor_name)
  }

  scaling_factors <- metadata$scaling_factors[factor_name, ]

  return(list(
    scaling_factors = scaling_factors,
    usage_level = metadata$usage_level
  ))
}

.detect_feature_level <- function(object, assay) {
  fcol <- QFeatures::assayLink(object, assay)@fcol
  if (fcol == "histone_group") {
    return("variant")
  } else if (fcol == "hptm") {
    return("ptm")
  } else {
    return("precursor")
  }
}
