#' Safely Replace `colData`
#'
#' Replaces `colData` of the given `QFeatures` or `SummarizedExperiment`,
#' including checking for common errors.
#'
#' @param object The `QFeatures` or `SummarizedExperiment` object of which the
#'   `colData` should be replaced.
#' @param colData New colData to be added to `object`. Can be either a path
#'   (`character(1)`) pointing to a colData `.csv` file, for example generated
#'   using [readProgenesis()], or a `tbl_df` where rows are samples and columns
#'   are properties of the samples. The samples need to match those in the
#'   `object` and the columns should contain at least the "quantCols",
#'   "original_name", "group", "include", and "outlier".
#' @param ... Additional arguments passed to specific methods.
#' @returns A `QFeatures` or `SummarizedExperiment` (same as supplied) with the
#'   new colData.
#' @export
#' @name replaceColData
#' @include hptmUsage-package.R
setGeneric("replaceColData", function(object, colData, ...) standardGeneric("replaceColData"))

#' @rdname replaceColData
#' @param custom_coltypes A named `list(...)` of types for non-standard columns,
#'   defaulting to [readr::col_factor()]. See [readr::cols()] for more
#'   information. There is no need to specify the types for "quantCols",
#'   "original_name", "group", "include", and "outlier".
#' @examples
#' replaceColData(
#'   ncbtoy,
#'   hptmUsageData("ncbtoy_coldata.csv"),
#'   list(
#'     ms_run = readr::col_integer(),
#'     date_collected = readr::col_date(format = "%y%m%d"),
#'     prep_batch = readr::col_factor()
#'   )
#' )
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
    stopifnot(all(
      c("quantCols", "original_name", "group", "include", "outlier") %in%
        colnames(readr::read_csv(colData, show_col_types = FALSE, progress = FALSE))
    ))

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

    # replace colData and return
    .replace_col_data(object, col_data)
  }
)

#' @rdname replaceColData
#' @examples
#' readr::read_csv(hptmUsageData("ncbtoy_coldata.csv")) |>
#'   dplyr::mutate(group = as.factor(group)) |>
#'   replaceColData(
#'     ncbtoy,
#'     colData = _
#'   )
setMethod(
  "replaceColData",
  c("QFeatures_OR_SummarizedExperiment", "tbl_df"),
  function(
    object,
    colData
  ) {
    # check arguments
    stopifnot(all(c("quantCols", "original_name", "group", "include", "outlier") %in% colnames(colData)))

    # check column types
    expected_coltypes <- c(
      quantCols = "character",
      original_name = "character",
      group = "factor",
      include = "logical",
      outlier = "logical"
    )
    # fmt: skip
    coltypes <- colData |>
      dplyr::summarize_all(class) |>
      _[names(expected_coltypes)] == expected_coltypes
    coltypes <- colnames(coltypes)[!coltypes]
    if (length(coltypes) >= 1) {
      warning(
        "Column types do not match expectation: ",
        paste(coltypes, expected_coltypes[coltypes], sep = "~", collapse = ", ")
      )
    }

    .replace_col_data(object, colData)
  }
)

.replace_col_data <- function(object, col_data) {
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
