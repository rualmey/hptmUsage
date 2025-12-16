#' Calculate hPTM/Variant Usage
#'
#' @description
#' hPTM or histone variant usage is calculated by first normalizing histone
#' peptidoforms (i.e., subtracting the "parent protein" abundance) and then
#' summarizing to the hPTM/variant level.
#'
#' @param object The `QFeatures` or `SummarizedExperiment` object from which
#'   usage needs to be calculated.
#' @param ... Additional arguments passed to specific methods.
#' @param target A `character(1)` defining whether PTM or histone variant usage
#'   will be calculated.
#' @param usage_level A `character(1)` defining the level at which usage will be
#'    calculated, i.e., what is considered the "parent protein". Must be one of:
#'   * `"histone"` = the total histone abundance (default)
#'   * `"histone_family"` = the total family abundance
#'   * `"histone_group"` = the histone group abundance, this should not be used if `target = "variant"`
#' @param deconv A `character(1)` pointing to the rowData column defining the
#'   features to deconvolute; for example, hPTMs after alignment mapping.
#' @param sep The `character(1)` used to separate different hPTMs/... in the
#'   `deconv` column.
#' @param group A `character(1)` pointing to the rowData column defining the
#'   level at which deconvolution features (e.g., hPTMs) will be grouped; for
#'   example, the histone family or histone variant group.
#' @returns
#' A `QFeatures` with usage assay(s) added or a usage `SummarizedExperiment`
#' (same as supplied).
#' @export
#' @references
#' (1) Demeulemeester N, Gébelin M, Gomes LC, Lingor P, Carapito C, Martens L,
#' Clement L. msqrob2PTM: Differential Abundance and Differential Usage Analysis
#' of MS-Based Proteomics Data at the Posttranslational Modification and
#' Peptidoform Level. Molecular & Cellular Proteomics. 2024 Feb 1;23(2):100708.
#' <https://pubs.acs.org/doi/full/10.1021/acs.jproteome.2c00145>.
#' @name calculateUsage
setGeneric(
  "calculateUsage",
  function(
    object,
    ...,
    target = c("ptm", "variant"),
    usage_level = c("histone", "histone_family", "histone_group"),
    deconv = "mods_ref",
    sep = ";",
    group = "histone_group"
  ) {
    standardGeneric("calculateUsage")
  },
  signature = "object"
)

#' @rdname calculateUsage
#' @param i The index (`integer()`) or name (`character()`) of the assay(s) to
#'   be processed.
#' @param name Name(s) of the new assay(s) to add to the QFeatures object. Must
#'   have the same length as i.
#' @param identifier A `character(1)` naming a rowData variable uniquely
#'   defining the features.
#' @param ... Additional arguments passed to [QFeatures::aggregateFeatures()].
#' @examples
#' \dontrun{
#' # By default, the entire chromatosome abundance (i.e., aggregated value of
#' # all histone features) is used as normalization factor
#' new_qf <- hptm_benchmark |>
#'   calculateUsage(i = "precursorHistone")
#' #' # A hPTM will be defined as the unique combination of histone
#' # group and PTM, e.g., H33K27Me3
#' rowData(new_qf[["ptm"]])[258, ]
#' # but if two hPTMs are only quantified from the same features, then these
#' # will form one "hPTM group".
#' rowData(new_qf[["ptm"]])[257, ]
#'
#' # It is also possible to correct for the histone family abundance
#' new_qf <- new_qf |>
#'   calculateUsage(i = "precursorHistone", name = "ptm_family", usage_level = "histone_family") |>
#' # or to do a "variant-corrected" analysis
#'   calculateUsage(i = "precursorHistone", name = "ptm_group", usage_level = "histone_group")
#'
#' # The default aggregation function (`MsCoreUtils::robustSummary()`) can be
#' # slow if feature groups are large, for example when `usage_level = "histone`
#' # In this case, it can be faster to use `MsCoreUtils::medianPolish()` instead
#' new_qf_mp <- hptm_benchmark |>
#'   calculateUsage(
#'     i = "precursorHistone",
#'     fun = MsCoreUtils::medianPolish,
#'     na.rm = TRUE
#'   )
#'
#' # Finally, if the unaligned amino acid index is to be used instead
#' new_qf_unaligned <- hptm_benchmark |>
#'   calculateUsage(i = "precursorHistone", deconv = "mods_var")
#'
#' # Histone variant usage is also easily calculated
#' new_qf_variant <- hptm_benchmark |>
#'   calculateUsage(i = "precursorHistone", target = "variant")
#'}
setMethod(
  "calculateUsage",
  "QFeatures",
  function(
    object,
    i,
    name = NULL,
    identifier = "feature_number",
    ...,
    target = c("ptm", "variant"),
    usage_level = c("histone", "histone_family", "histone_group"),
    deconv = "mods_ref",
    sep = ";",
    group = "histone_group"
  ) {
    # Check arguments
    i <- QFeatures:::.normIndex(object, i)
    target <- match.arg(target)
    if (is.null(name)) {
      name <- target
    }
    stopifnot(length(i) == length(name))
    usage_level <- match.arg(usage_level)

    object <- object |>
      normalizeUsage(
        i = i,
        name = paste0(i, "Norm", stringr::str_to_sentence(usage_level)),
        usage_level = usage_level,
        ...
      )
    if (target == "ptm") {
      object <- object |>
        deconvolute(
          i = paste0(i, "Norm", stringr::str_to_sentence(usage_level)),
          name = paste0(i, "Norm", stringr::str_to_sentence(usage_level), "Deconv"),
          identifier = identifier,
          deconv = deconv,
          sep = sep,
          group = group
        ) |>
        QFeatures::aggregateFeatures(
          i = paste0(i, "Norm", stringr::str_to_sentence(usage_level), "Deconv"),
          fcol = "hptm",
          name = name,
          ...
        )
    } else {
      object <- object |>
        QFeatures::aggregateFeatures(
          i = paste0(i, "Norm", stringr::str_to_sentence(usage_level)),
          fcol = "histone_group",
          name = name,
          ...
        )
    }
  }
)

#' @rdname calculateUsage
#' @param ... Additional arguments passed to [QFeatures::aggregateFeatures()].
#' @examples
#' \dontrun{
#' # By default, the entire chromatosome abundance (i.e., aggregated value of
#' # all histone features) is used as normalization factor
#' new_se <- hptm_benchmark[[5]] |>
#'   calculateUsage()
#' # A hPTM will be defined as the unique combination of histone
#' # group and PTM, e.g., H33K27Unmod
#' rowData(new_se)[258, ]
#' # but if two hPTMs are only quantified from the same features, then these
#' # will form one "hPTM group".
#' rowData(new_se)[257, ]
#'
#' # It is also possible to correct for the histone family abundance
#' new_se_hf <- hptm_benchmark[[5]] |>
#'   calculateUsage(usage_level = "histone_family")
#' # or to do a "variant-corrected" analysis
#' new_se_vc <- hptm_benchmark[[5]] |>
#'   calculateUsage(usage_level = "histone_group")
#'
#' # The default aggregation function (`MsCoreUtils::robustSummary()`) can be
#' # slow if feature groups are large, for example when `usage_level = "histone`
#' # In this case, it can be faster to use `MsCoreUtils::medianPolish()` instead
#' new_se_mp <- hptm_benchmark[[5]] |>
#'   calculateUsage(fun = MsCoreUtils::medianPolish, na.rm = TRUE)
#'
#' # Finally, if the unaligned amino acid index is to be used instead
#' new_se_unaligned <- hptm_benchmark[[5]] |>
#'   calculateUsage(deconv = "mods_var")
#'
#' # Histone variant usage is also easily calculated
#' new_se_variant <- hptm_benchmark[[5]] |>
#'   calculateUsage(target = "variant")
#'}
setMethod(
  "calculateUsage",
  "SummarizedExperiment",
  function(
    object,
    ...,
    target = c("ptm", "variant"),
    usage_level = c("histone", "histone_family", "histone_group"),
    deconv = "mods_ref",
    sep = ";",
    group = "histone_group"
  ) {
    target <- match.arg(target)
    usage_level <- match.arg(usage_level)

    object <- object |>
      normalizeUsage(usage_level = usage_level, ...)
    if (target == "ptm") {
      object <- object |>
        deconvolute(deconv = deconv, sep = sep, group = group) |>
        # consider MsCoreUtils::medianPolish() if fcol = histone?; this is faster and should perform similarly here
        QFeatures::aggregateFeatures(fcol = "hptm", ...)
    } else {
      object <- object |>
        QFeatures::aggregateFeatures(fcol = "histone_group", ...)
    }

    return(object)
  }
)
