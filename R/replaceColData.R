#' Safely Replace `colData` of a QFeatures Object
#'
#' TODO Description here. Included error checking...
#'
#' @param object The QFeatures object in which to replace the colData.
#' @param colData TODO The path (`character(1)`) pointing to a colData `.csv` file. This file should contain at least the
#'  columns of quantCols, original_name, group, include, and outlier.
#' @param coltypes Column types used by `readr::read_csv()`, also see `vignette("readr")` for more details.
#' @returns A `QFeatures` object with colData replaced.
#' @export
#' @name replaceColData
#' @include hptmUsage-package.R
#' @examples
#' # TODO
#' NULL
setGeneric("replaceColData", function(object, colData, ...) standardGeneric("replaceColData"))

#' @rdname replaceColData
setMethod(
  "replaceColData",
  c("QFeatures_OR_SummarizedExperiment", "character"),
  function(
    object,
    colData,
    custom_coltypes = NULL
  ) {
    # check arguments
    stopifnot(file.exists(colData))

    # read colData csv and check required columns
    coltypes <- list(
      quantCols = readr::col_character(),
      original_name = readr::col_character(),
      group = readr::col_factor(),
      include = readr::col_logical(),
      outlier = readr::col_logical(),
      .default = readr::col_factor()
    ) |>
      c(custom_coltypes) |>
      do.call(what = readr::cols, args = _)
    col_data <- readr::read_csv(colData, col_types = coltypes, progress = FALSE)
    stopifnot(all(c("quantCols", "original_name", "group", "include", "outlier") %in% colnames(col_data)))

    # replace colData and return
    replace_col_data(object, col_data)
  }
)

#' @rdname replaceColData
setMethod(
  "replaceColData",
  c("QFeatures_OR_SummarizedExperiment", "tbl_df"),
  function(
    object,
    colData,
    custom_coltypes = NULL
  ) {
    # check arguments TODO
    # stopifnot(file.exists(colData))

    # # read colData csv and check required columns
    # coltypes <- list(
    #   quantCols = readr::col_character(),
    #   original_name = readr::col_character(),
    #   group = readr::col_factor(),
    #   include = readr::col_logical(),
    #   outlier = readr::col_logical(),
    #   .default = readr::col_factor()
    # ) |>
    #   c(custom_coltypes) |>
    #   do.call(what = readr::cols, args = _)
    # col_data <- readr::read_csv(colData, col_types = coltypes, progress = FALSE)
    # stopifnot(all(c("quantCols", "original_name", "group", "include", "outlier") %in% colnames(col_data)))

    # replace colData and return
    replace_col_data(object, colData)
  }
)

replace_col_data <- function(object, col_data) {
  # make sure the same quantCols are present
  existing_quantcols <- rownames(colData(object))
  quantcols_not_meta <- setdiff(existing_quantcols, col_data$quantCols)
  meta_not_quantcols <- setdiff(col_data$quantCols, existing_quantcols)
  if (length(quantcols_not_meta) >= 1 || length(meta_not_quantcols) >= 1) {
    error_str <- "Samples do not match between the `object` and `colData`."
    if (length(quantcols_not_meta) >= 1) {
      error_str <- paste(
        error_str,
        paste("  Samples not in `quantCols` column of the colData:", paste(quantcols_not_meta, collapse = ", ")),
        sep = "\n"
      )
    }
    if (length(meta_not_quantcols) >= 1) {
      error_str <- paste(
        error_str,
        paste("  Samples not in present in the `object`:", paste(meta_not_quantcols, collapse = ", ")),
        sep = "\n"
      )
    }
    stop(error_str)
  }

  # replace colData in QFeatures object, preserve ordering from QFeatures object
  col_data <- col_data |>
    dplyr::arrange(match(.data[["quantCols"]], existing_quantcols)) |>
    as("DataFrame")
  rownames(col_data) <- col_data$quantCols
  colData(object) <- col_data
  object
}
