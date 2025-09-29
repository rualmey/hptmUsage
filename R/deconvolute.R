#' Deconvolute (hPTM) Features WIP
#'
#' @description
#' This is primarily an internal function required before usage aggregation and
#' normalization.
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
#' @param object The `QFeatures` or `SummarizedExperiment` object in which the
#'   contaminants should be tagged.
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
    group = NULL
  ) {
    standardGeneric("deconvolute")
  },
  signature = "object"
)

#' @rdname deconvolute
#' @param i The index (`integer()`) or name (`character()`) of the assay(s) to
#'   be processed.
#' @param fcol A `character(1)` naming a rowData variable defining the unique
#'   features to deconvolute.
#' @param name Name(s) of the new assay(s) to add to the QFeatures object. Must
#'   have the same length as i.
#' @param filter Either an instance of class `AnnotationFilter` or a formula.
#'   See [QFeatures::filterFeatures()]
#' @examples
#' # TODO
setMethod(
  "deconvolute",
  "QFeatures",
  function(
    object,
    i,
    fcol,
    name = "precursorDeconv",
    filter = NULL,
    deconv = "mods_ref",
    sep = ";",
    group = NULL
  ) {
    # Check arguments
    i <- QFeatures:::.normIndex(object, i)
    stopifnot(length(i) == length(name))
    # Create and add new assays with the transformed data
    for (j in seq_along(i)) {
      from <- i[[j]]
      to <- name[[j]]
      # drop if NA in deconv column, but only for a new assay
      .object <- substitute(object)
      .deconv <- as.name(deconv)
      .filter <- substitute(~ !is.na(deconv), list(deconv = .deconv))
      .from <- substitute(from)
      no_na_se <- eval.parent(
        substitute(
          filterFeatures(object, filter, from, keep = TRUE),
          list(object = .object, filter = .filter, from = .from)
        )
      )[[from]]
      from <- paste0(from, "Filtered")
      object <- addAssay(object, no_na_se, from)
      object <- addAssayLink(object, from = .from, to = from, varFrom = fcol, varTo = fcol)
      # prepend top-level to deconvolution value
      if (!is.null(group)) {
        # reshape deconv column in the form of "group + deconv"
        rowData(object[[from]])$deconvoluted <- mapply(
          function(x, y) {
            deconv_split <- strsplit(y, split = sep, fixed = TRUE)
            top_paste <- sapply(deconv_split, \(i) paste(x, i, sep = "#"))
            paste(top_paste, collapse = sep)
          },
          rowData(object[[from]])[[group]],
          rowData(object[[from]])[[deconv]]
        )
      }
      if (!is.null(filter)) {
        object <- filterFeatures(object, filter, from, keep = TRUE)
      }
      object <- addAssay(
        object,
        deconvolute(object[[from]], deconv = "deconvoluted", group = group, sep = sep),
        to
      )
      object <- addAssayLink(object, from = from, to = to, varFrom = fcol, varTo = fcol)
    }
    if (validObject(object)) {
      # necessary or already checked during calls?
      return(object)
    } else {
      return(validObject(object))
    }
  }
)

# setMethod(
#   "deconvolute",
#   "MultiAssayExperiment",
#   function(
#     object,
#     deconv = "mods_ref",
#     sep = ";",
#     group = NULL
#   ) {
#     # get all unique PTM groups (collapse 'degenerate' PTMs)
#     deconv_groups <- .deconv_grouper(object, deconv = deconv, group = group, sep = sep)
#     deconv_groups <- strsplit(deconv_groups, ",", fixed = TRUE)
#     deconv_se <- list()
#     for (i in seq_along(deconv_groups)) {
#       se <- object[deconv_groups[[i]], ]
#       rowData(se)[[deconv]] <- names(deconv_groups)[[i]]
#       deconv_se[[names(deconv_groups)[[i]]]] <- se
#     }
#     se_bind <- do.call(rbind, deconv_se)
#     rownames(se_bind) <- seq_len(dim(se_bind)[[1]])
#     se_bind
#   }
# )

.deconv_grouper2 <- function(object, deconv, group, sep) {
  # get all unique groups and collapse degenerated groups that come from the same data
  all_groups <- strsplit(rowData(object)[[deconv]], split = sep, fixed = TRUE)
  unique_groups <- unique(unlist(all_groups))
  deconv_fcol_map <- sapply(unique_groups, function(i) {
    row_matches <- sapply(all_groups, \(j) i %in% j)
    rownames(object[row_matches, ])
  })
  deconv_fcol_map <- split(names(deconv_fcol_map), sapply(deconv_fcol_map, paste, collapse = ","))
  if (!is.null(group)) {
    deconv_groups <- sapply(deconv_fcol_map, stringr::str_split, pattern = "#") |>
      lapply(function(x) {
        top_level_str <- sapply(x, \(x) x[1])
        deconv_str <- paste(sapply(x, \(x) x[2]), collapse = ";")
        paste(top_level_str[[1]], deconv_str, sep = "#")
      })
  } else {
    deconv_groups <- lapply(deconv_fcol_map, paste, collapse = ";")
  }
  setNames(names(deconv_groups), deconv_groups)
}

#' @rdname deconvolute
#' @examples
#' # TODO
setMethod(
  "deconvolute",
  "SummarizedExperiment",
  function(
    object,
    deconv = "mods_ref",
    sep = ";",
    group = NULL
  ) {
    deconv_map <- .deconv_grouper(object, deconv, group, sep)
    reconstruction_plan <- deconv_map |>
      tidyr::unnest(.rows)
    new_se <- object[reconstruction_plan$.rows, ]
    new_rd <- SummarizedExperiment::rowData(new_se)
    new_rd[[deconv]] <- reconstruction_plan$.new_deconv
    SummarizedExperiment::rowData(new_se) <- S4Vectors::DataFrame(new_rd)
    rownames(new_se) <- seq_len(nrow(new_se))
    new_se
  }
)

.deconv_grouper <- function(object, deconv, group, sep) {
  degenerate_groups <- SummarizedExperiment::rowData(object) |>
    tibble::as_tibble(rownames = ".rowname") |>
    dplyr::select(dplyr::all_of(c(".rowname", deconv))) |>
    tidyr::separate_rows(!!rlang::sym(deconv), sep = sep) |>
    dplyr::group_by(!!rlang::sym(deconv)) |>
    dplyr::summarise(.rows = list(sort(.rowname)), .groups = "drop") |>
    dplyr::mutate(.origin_id = purrr::map_chr(.rows, paste, collapse = ",")) |>
    dplyr::group_by(.origin_id) |>
    dplyr::summarise(
      .ptms = list(!!rlang::sym(deconv)),
      .rows = dplyr::first(.rows),
      .groups = "drop"
    )
  if (!is.null(group)) {
    processed_groups <- degenerate_groups |>
      dplyr::mutate(
        .new_deconv = purrr::map_chr(
          .ptms,
          ~ {
            parts <- stringr::str_split(.x, "#", n = 2)
            group_val <- purrr::map_chr(parts, 1) |> purrr::first()
            ptm_vals <- purrr::map_chr(parts, 2) |> paste(collapse = ";")
            paste(group_val, ptm_vals, sep = "#")
          }
        )
      )
  } else {
    processed_groups <- degenerate_groups |>
      dplyr::mutate(.new_deconv = purrr::map_chr(.ptms, paste, collapse = ";"))
  }
  processed_groups |>
    dplyr::select(.new_deconv, .rows)
}
