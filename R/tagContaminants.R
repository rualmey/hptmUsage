#' Tag Contaminant Proteins
#'
#' @description
#' Use custom contaminant libraries from Frankenfield et al. (1), where histone
#' entries (Cont_A1A4R1, Cont_P62803, Cont_Q64475) were manually removed, to tag
#' contaminant proteins for later filtering. See
#' <https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics>
#' for more information.
#'
#' @param object The `QFeatures` or `SummarizedExperiment` object in which the
#'   contaminants should be tagged.
#' @param ... Additional arguments passed to specific methods.
#' @param library The name (`character(1)`) of the contaminant library to use.
#'   Must be one of:
#'   * `"universal"` (default)
#'   * `"cell_culture"`
#'   * `"mouse_tissue"`
#'   * `"rat_tissue"`
#'   * `"neuron_culture"`
#'   * `"stem_cell_culture"`
#' @param sequence The name (`character(1)`) of the rowData column containing
#'   (stripped) peptide sequences.
#' @returns
#' A `QFeatures` or `SummarizedExperiment` (same as supplied) with the rowData
#' containing one new `logical()` column "contaminant".
#' @export
#' @references
#' (1) Frankenfield AM, Ni J, Ahmed M, Hao L. Protein contaminants matter:
#' building universal protein contaminant libraries for DDA and DIA proteomics.
#' Journal of proteome research. 2022 Jul 6;21(9):2104-13.
#' <https://pubs.acs.org/doi/full/10.1021/acs.jproteome.2c00145>.
#' @name tagContaminants
setGeneric(
  "tagContaminants",
  function(
    object,
    ...,
    library = c("universal", "cell_culture", "mouse_tissue", "rat_tissue", "neuron_culture", "stem_cell_culture"),
    sequence = "sequence"
  ) {
    standardGeneric("tagContaminants")
  },
  signature = "object"
)

#' @rdname tagContaminants
#' @param i The index (`integer()`) or name (`character()`) of the assay(s) to be processed.
#' @examples
#' # By default, contaminants will be tagged from the "universal" library
#' qf <- ncbtoy |>
#'   tagContaminants(i = "precursorRaw")
#' rowData(qf)
#' # but we can easily specify other libraries
#' qf_alt <- ncbtoy |>
#'   tagContaminants(i = "precursorRaw", library = "stem_cell_culture")
#' rowData(qf_alt)
setMethod(
  "tagContaminants",
  "QFeatures",
  function(
    object,
    i,
    library = c("universal", "cell_culture", "mouse_tissue", "rat_tissue", "neuron_culture", "stem_cell_culture"),
    sequence = "sequence"
  ) {
    i <- QFeatures:::.normIndex(object, i)
    for (j in i) {
      object <- QFeatures::replaceAssay(
        object,
        tagContaminants(object[[j]], library = library, sequence = sequence),
        j
      )
    }
    object
  }
)

#' @rdname tagContaminants
#' @examples
#' # By default, contaminants will be tagged from the "universal" library
#' se <- ncbtoy[[1]] |>
#'   tagContaminants()
#' rowData(se)
#' # but we can easily specify other libraries
#' se_alt <- ncbtoy[[1]] |>
#'   tagContaminants(library = "stem_cell_culture")
#' rowData(se_alt)
setMethod(
  "tagContaminants",
  "SummarizedExperiment",
  function(
    object,
    library = c("universal", "cell_culture", "mouse_tissue", "rat_tissue", "neuron_culture", "stem_cell_culture"),
    sequence = "sequence"
  ) {
    # Check arguments
    library <- match.arg(library)
    stopifnot(
      "`library` must be a single valid library name" = is.character(library) &&
        length(library) == 1 &&
        library %in% names(.contaminants),
      "`sequence` must be a single valid column name in `rowData(object)`" = is.character(sequence) &&
        length(sequence) == 1 &&
        sequence %in% colnames(rowData(object))
    )
    # match contaminant sequences and tag them
    contam_idx <- Biostrings::AAStringSet(rowData(object)[[sequence]]) |>
      Biostrings::vwhichPDict(.contaminants[[library]]) |>
      unlist() |>
      unique()

    rowData(object)$contaminant <- seq_len(nrow(object)) %in% contam_idx

    object
  }
)
