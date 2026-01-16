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
  assays,
  model = "factor",
  significant_only = TRUE
) {
  stopifnot("'object' must be a QFeatures object." = is(object, "QFeatures"))
  stopifnot("'assays' must only contain names 'ptm' and/or 'variant'" = all(names(assays) %in% c("ptm", "variant")))
  purrr::imap(assays, \(x, idx) list(QFeatures:::.normIndex(object, x)) |> purrr::set_names(idx))

  # retrieve feature data = summarized and estimated PTM or variant usage
  features <- purrr::imap(assays, \(x, idx) .get_features(x, idx, object, model, significant_only))
  summarized_usage <- purrr::imap(features, \(x, idx) .get_summarized_usage(x, idx, object))
  estimated_usage <- purrr::imap(features, \(x, idx) .get_estimated_usage(x, idx, object, model))

  # retrieve precursor data = abundance and usage of corresponding precursors
  precursor_features <- purrr::imap(features, \(x, idx) .get_precursor_features(x, idx, object))
  # TODO
  return(precursor_features)

  # retrieve normalization factors
  norm_factors <- purrr::imap(features, \(x, idx) .get_norm_factors(x, idx, object))
}

# ------------------------------------------------------------------------------
# Internal Helper Functions
# ------------------------------------------------------------------------------

.get_features <- function(assay_names, assay_level, object, model, significant_only) {
  purrr::map(assay_names, function(x) {
    rd <- rowData(object[[x]])
    # rowData columns with data.frames containing adjPval
    contrast_df <- paste0(names(model), "/") |>
      paste(collapse = "|") |>
      stringr::str_starts(names(rd), pattern = _)
    contrast_df <- purrr::map(names(rd)[contrast_df], \(x) tibble::as_tibble(rd[[x]], rownames = "name")) |>
      dplyr::bind_rows()
    if (significant_only) {
      contrast_df <- contrast_df |>
        dplyr::filter(.data[["adjPval"]] <= 0.05)
    }
    contrast_df <- contrast_df |>
      dplyr::arrange(.data[["adjPval"]]) |>
      dplyr::pull(.data[["name"]]) |>
      unique()
  }) |>
    # store assay name for later retrieval
    purrr::set_names(assay_names)
}

.get_summarized_usage <- function(features, assay_level, object) {
  purrr::imap(features, function(x, idx) {
    SummarizedExperiment::assay(object[[idx]][x]) |>
      tibble::as_tibble(rownames = "feature") |>
      tidyr::pivot_longer(!feature, names_to = "sample", values_to = "value") |>
      dplyr::mutate(level = if (assay_level == "ptm") "PTM\nUsage" else "Variant\nUsage")
  })
}

.get_estimated_usage <- function(features, assay_level, object, model) {
  purrr::imap(features, function(x, idx) {
    purrr::imap(model, function(model_formula, model_name) {
      design_mat <- model.matrix(as.formula(model_formula), data = droplevels(colData(object)))
      purrr::map(rowData(object[[idx]][x])[[model_name]], function(model) {
        coef <- msqrob2::getCoef(model)
        # fix ridge regression as this prepends "ridge" to the factor name
        names(coef) <- stringr::str_remove(names(coef), "^ridge")
        fixed_coef <- intersect(names(coef), colnames(design_mat))
        # only calculate fixed effects as random effects are of no interest
        design_mat[, fixed_coef, drop = FALSE] %*%
          coef[fixed_coef] |>
          as.numeric() |>
          purrr::set_names(rownames(colData(object)))
      }) |>
        tibble::as_tibble() |>
        dplyr::mutate(sample = rownames(colData(object))) |>
        tidyr::pivot_longer(!sample, names_to = "feature", values_to = "value") |>
        dplyr::relocate(feature) |>
        dplyr::mutate(
          level = if (assay_level == "ptm") {
            paste0("PTM Estimated\nModel '", model_name, "'")
          } else {
            paste0("Variant Estimated\nModel '", model_name, "'")
          }
        )
    }) |>
      dplyr::bind_rows()
  })
}

.link_data <- function(object, assay, assay_level = NULL) {
  if (is.null(assay_level) || assay_level == "variant") {
    link <- QFeatures::assayLink(object, assay)
    precursor_assay <- link@from
    link_fcol <- link@fcol
  } else if (assay_level == "ptm") {
    links <- QFeatures::assayLinks(object, assay)
    precursor_assay <- names(links)[[3]]
    link_fcol <- links[[1]]@fcol
  } else {
    # incorrect names?
    stop()
  }

  list(precursor_assay = precursor_assay, fcol = link_fcol)
}

.get_precursor_features <- function(assay_features, assay_level, object) {
  purrr::imap(assay_features, function(x, idx) {
    link_data <- .link_data(object, idx)
    rd <- rowData(object[[link_data$precursor_assay]])
    purrr::map(x, \(y) rd[rd[[link_data$fcol]] == y, ]$precursor) |>
      purrr::set_names(x)
  })
}

.get_norm_factors <- function(assay_features, assay_level, object) {
  purrr::imap(assay_features, function(x, idx) {
    # retrieve metadata
    precursor_assay <- .link_data(object, idx, assay_level)$precursor_assay
    metadata <- metadata(object[[precursor_assay]])

    # match factors to features
    factor_names <- rowData(object[[idx]][x])[[metadata$usage_level]]
    if (metadata$usage_level == "histone") {
      factor_names <- as.integer(factor_names)
    }
    scaling_factors <- metadata$scaling_factors[factor_names, ]

    list(
      scaling_factors = scaling_factors,
      usage_level = metadata$usage_level
    )
  })
}


# OLD TODO

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
