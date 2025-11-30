#' Deconvolute (hPTM) Features
#'
#' @description
#' This is primarily an internal function required before usage aggregation.
#'
#' @details
#' One feature can contain multiple hPTMs. Such features will be replicated and
#' assigned one hPTM each. This is needed to be able to use [calculateUsage()],
#' which relies on [QFeatures::aggregateFeatures()].
#'
#' Different hPTMs can have the same underlying information. These "degenerate"
#' PTMs cannot be discerned from one another, so they are added into "hPTM
#' groups".
#'
#' @param object The `QFeatures` or `SummarizedExperiment` object to be
#'    deconvoluted.
#' @param ... Additional arguments passed to specific methods.
#' @param deconv A `character(1)` pointing to the rowData column defining the
#'   features to deconvolute; for example, hPTMs after alignment mapping.
#' @param sep The `character(1)` used to separate different hPTMs/... in the
#'   `deconv` column.
#' @param group Optional, a `character(1)` pointing to the rowData column
#'   defining the level at which deconvolution features (e.g., hPTMs) will be
#'   grouped; for example, the histone family or histone variant group.
#' @returns
#' A `QFeatures` with deconvoluted assay(s) added or a deconvoluted
#' `SummarizedExperiment` (same as supplied).
#' @export
#' @name deconvolute
setGeneric(
  "deconvolute",
  function(
    object,
    ...,
    deconv = "mods_ref",
    sep = ";",
    group = "histone_group"
  ) {
    standardGeneric("deconvolute")
  },
  signature = "object"
)

#' @rdname deconvolute
#' @param i The index (`integer()`) or name (`character()`) of the assay(s) to
#'   be processed.
#' @param name Name(s) of the new assay(s) to add to the QFeatures object. Must
#'   have the same length as i.
#' @param identifier A `character(1)` naming a rowData variable uniquely
#'   defining the features.
#' @examples
#' new_qf <- hptm_benchmark |>
#'   deconvolute(i = "precursor")
#' # By default, a hPTM will be defined as the unique combination of histone
#' # group and PTM, e.g., H33K27Me3
#' rowData(new_qf[["precursorDeconv"]])[938, ]
#' # but if two hPTMs are only quantified from the same features, then these
#' # will form one "hPTM group".
#' rowData(new_qf[["precursorDeconv"]])[910, ]
#' # and if the unaligned amino acid index is to be used instead
#' new_qf_unaligned <- hptm_benchmark |>
#'   deconvolute(i = "precursor", deconv = "mods_var")
setMethod(
  "deconvolute",
  "QFeatures",
  function(
    object,
    i,
    name = "precursorDeconv",
    identifier = "feature_number",
    deconv = "mods_ref",
    sep = ";",
    group = "histone_group"
  ) {
    # Check arguments
    i <- QFeatures:::.normIndex(object, i)
    stopifnot(length(i) == length(name))
    # Create and add new assays with the transformed data
    for (j in seq_along(i)) {
      from <- i[[j]]
      to <- name[[j]]
      object <- QFeatures::addAssay(
        object,
        deconvolute(object[[from]], deconv = deconv, group = group, sep = sep),
        to
      )
      object <- QFeatures::addAssayLink(
        object,
        from = from,
        to = to,
        varFrom = identifier,
        varTo = identifier
      )
    }

    object
  }
)

#' @rdname deconvolute
#' @examples
#' new_se <- hptm_benchmark[[2]] |>
#'   deconvolute()
#' # By default, a hPTM will be defined as the unique combination of histone
#' # group and PTM, e.g., H33K27Me3
#' rowData(new_se)[938, ]
#' # but if two hPTMs are only quantified from the same features, then these
#' # will form one "hPTM group".
#' rowData(new_se)[910, ]
#' # and if the unaligned amino acid index is to be used instead
#' new_se_unaligned <- hptm_benchmark[[2]] |>
#'   deconvolute(deconv = "mods_var")
setMethod(
  "deconvolute",
  "SummarizedExperiment",
  function(
    object,
    deconv = "mods_ref",
    sep = ";",
    group = "histone_group"
  ) {
    # filter out features that carry no hPTMs
    object <- object[!is.na(rowData(object)[[deconv]]), ]
    deconv_map <- .deconv_grouper(object, deconv, sep, group)
    new_se <- object[deconv_map$rows_same_ptm, ]
    rowData(new_se)$hptm <- deconv_map$ptms
    # rownames need to be unique so just reset
    rownames(new_se) <- seq_len(nrow(new_se))

    new_se
  }
)

.deconv_grouper <- function(object, deconv, sep, group) {
  degenerate_groups <- rowData(object) |>
    tibble::as_tibble(rownames = ".rowname") |>
    dplyr::select(dplyr::all_of(c(".rowname", deconv, group))) |>
    tidyr::separate_longer_delim(dplyr::all_of(deconv), sep) |>
    dplyr::summarise(rows_same_ptm = list(sort(.data[[".rowname"]])), .by = dplyr::all_of(c(deconv, group)))
  if (!is.null(group)) {
    processed_groups <- degenerate_groups |>
      dplyr::summarise(
        ptms = stringr::str_flatten(!!rlang::sym(deconv), collapse = sep),
        # by grouping earlier, this will only contain 1 unique histone family
        grouper = dplyr::first(!!rlang::sym(group)),
        .by = dplyr::all_of("rows_same_ptm")
      ) |>
      dplyr::mutate(
        ptms = stringr::str_c(.data[["grouper"]], .data[["ptms"]], sep = "#"),
        .keep = "unused"
      )
  } else {
    processed_groups <- degenerate_groups |>
      dplyr::summarise(
        ptms = stringr::str_flatten(!!rlang::sym(deconv), collapse = sep),
        .by = dplyr::all_of("rows_same_ptm")
      )
  }

  processed_groups |> tidyr::unnest(dplyr::all_of("rows_same_ptm"))
}
