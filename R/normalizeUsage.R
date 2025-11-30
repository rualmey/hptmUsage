#' Usage Normalization of hPTMs
#'
#' @description
#' This is primarily an internal function required before usage aggregation.
#'
#' @details
#' Histone precursors are normalized by subtracting the parent protein abundance
#' from the (globally normalized) feature abundances. See also Demeulemeester et
#' al. (1).
#'
#' @param object The `QFeatures` or `SummarizedExperiment` object to be
#'    usage normalized.
#' @param ... Additional arguments passed to specific methods.
#' @param usage_level A `character(1)` defining the level at which usage will be
#'    calculated, i.e., what is considered the "parent protein". Must be one of:
#'   * `"histone"` = the total histone abundance (default)
#'   * `"histone_family"` = the total family abundance
#'   * `"histone_group"` = the histone group abundance
#' @returns
#' A `QFeatures` with usage normalized assay(s) added or a deconvoluted
#' `SummarizedExperiment` (same as supplied).
#' @export
#' @references
#' (1) Demeulemeester N, Gébelin M, Gomes LC, Lingor P, Carapito C, Martens L,
#' Clement L. msqrob2PTM: Differential Abundance and Differential Usage Analysis
#' of MS-Based Proteomics Data at the Posttranslational Modification and
#' Peptidoform Level. Molecular & Cellular Proteomics. 2024 Feb 1;23(2):100708.
#' <https://pubs.acs.org/doi/full/10.1021/acs.jproteome.2c00145>.
#' @name normalizeUsage
setGeneric(
  "normalizeUsage",
  function(
    object,
    ...,
    usage_level = c("histone", "histone_family", "histone_group")
  ) {
    standardGeneric("normalizeUsage")
  },
  signature = "object"
)

#' @rdname normalizeUsage
#' @param i The index (`integer()`) or name (`character()`) of the assay(s) to
#'   be processed.
#' @param name Name(s) of the new assay(s) to add to the QFeatures object. Must
#'   have the same length as i.
#' @param ... Additional arguments passed to [QFeatures::aggregateFeatures()].
#' @examples
#' \dontrun{
#' # By default, the entire chromatosome abundance (i.e., aggregated value of
#' # all histone features) is used as normalization factor
#' new_qf <- hptm_benchmark |>
#'   normalizeUsage(i = "precursorHistone")
#' # It is also possible to to correct for the histone family abundance
#' new_qf_hf <- hptm_benchmark |>
#'   normalizeUsage(i = "precursorHistone", usage_level = "histone_family")
#' # or to do a "variant-corrected" analysis
#' new_qf_vc <- hptm_benchmark |>
#'   normalizeUsage(i = "precursorHistone", usage_level = "histone_group")
#'
#' # The default aggregation function (`MsCoreUtils::robustSummary()`) can be
#' # slow if feature groups are large, for example when `usage_level = "histone`
#' # In this case, it can be faster to use `MsCoreUtils::medianPolish()` instead
#' new_qf_mp <- hptm_benchmark |>
#'   normalizeUsage(
#'     i = "precursorHistone",
#'     fun = MsCoreUtils::medianPolish,
#'     na.rm = TRUE
#'   )
#'}
setMethod(
  "normalizeUsage",
  "QFeatures",
  function(
    object,
    i,
    name = "precursorNorm",
    ...,
    usage_level = c("histone", "histone_family", "histone_group")
  ) {
    # Check arguments
    i <- QFeatures:::.normIndex(object, i)
    stopifnot(length(i) == length(name))
    usage_level <- match.arg(usage_level)

    # Create and add new assays with the transformed data
    for (j in seq_along(i)) {
      from <- i[[j]]
      to <- name[[j]]
      stopifnot(
        "'usage_level' should not contain any `NA`, filter to only histones" = !any(is.na(rowData(object[[from]])[[
          usage_level
        ]]))
      )
      object <- QFeatures::addAssay(
        object,
        normalizeUsage(object[[from]], usage_level = usage_level, ...),
        to
      )
      object <- QFeatures::addAssayLinkOneToOne(object, from = from, to = to)
    }

    object
  }
)

#' @rdname normalizeUsage
#' @param ... Additional arguments passed to [QFeatures::aggregateFeatures()].
#' @examples
#' \dontrun{
#' # By default, the entire chromatosome abundance (i.e., aggregated value of
#' # all histone features) is used as normalization factor
#' new_se <- hptm_benchmark[[5]] |>
#'   normalizeUsage()
#' # It is also possible to to correct for the histone family abundance
#' new_se_hf <- hptm_benchmark[[5]] |>
#'   normalizeUsage(usage_level = "histone_family")
#' # or to do a "variant-corrected" analysis
#' new_se_vc <- hptm_benchmark[[5]] |>
#'   normalizeUsage(usage_level = "histone_group")
#'
#' # The default aggregation function (`MsCoreUtils::robustSummary()`) can be
#' # slow if feature groups are large, for example when `usage_level = "histone`
#' # In this case, it can be faster to use `MsCoreUtils::medianPolish()` instead
#' new_se_mp <- hptm_benchmark[[5]] |>
#'   normalizeUsage(fun = MsCoreUtils::medianPolish, na.rm = TRUE)
#'}
setMethod(
  "normalizeUsage",
  "SummarizedExperiment",
  function(
    object,
    ...,
    usage_level = c("histone", "histone_family", "histone_group")
  ) {
    usage_level <- match.arg(usage_level)
    stopifnot(
      "'usage_level' should not contain any `NA`, filter to only histones" = !any(is.na(rowData(object)[[usage_level]]))
    )

    agg <- QFeatures::aggregateFeatures(
      object,
      fcol = usage_level,
      ...
    )

    # extract protein abundances as numeric scaling factors
    scaling_factors <- SummarizedExperiment::assay(agg, "assay", withDimnames = FALSE)
    # and match to original assay dimensions
    usage_level_identifier <- rowData(object)[[usage_level]]
    # workaround for boolean "histone" column
    if (usage_level == "histone") {
      usage_level_identifier <- as.character(usage_level_identifier)
    }
    scaling_factors <- scaling_factors[usage_level_identifier, ]

    SummarizedExperiment::assay(object) <- SummarizedExperiment::assay(object) - scaling_factors

    object
  }
)
