#' Match Sequences to Histone Variants
#'
#' @description
#' Adds columns to the `rowData` of the given `QFeatures` or `SummarizedExperiment`,
#' including checking for common errors.
#'
#' @param object The `QFeatures` or `SummarizedExperiment` object in which the
#'   sequences should be matched to histone variants.
#' @param matching_subject An `AAStringSetList` of histone families and
#'   corresponding sequences. For example, the "unaligned" element resulting
#'   from [alignHistones()].
#' @param sequence_col The column name (`character(1)`) containing the sequences.
#' @param ... Additional arguments passed to specific methods.
#' @returns A `QFeatures` or `SummarizedExperiment` (same as supplied) with
#'   modified rowData of histone variant matches and the location.
#' @export
#' @name matchHistones
setGeneric(
  "matchHistones",
  function(object, matching_subject, ..., sequence_col = "sequence") {
    standardGeneric("matchHistones")
  },
  signature = c("object", "matching_subject")
)

#' @rdname matchHistones
#' @param i The index (`integer(1)`) or name (`character(1)`) of the assay to be processed.
#' @examples
#' # TODO
#' NULL
setMethod(
  "matchHistones",
  c("QFeatures", "AAStringSetList"),
  function(
    object,
    i,
    matching_subject,
    sequence_col = "sequence"
  ) {
    i <- QFeatures:::.normIndex(object, i)
    purrr::map(
      i,
      \(x) `rowData<-`(object[[x]], matchHistones(object[[x]], matching_subject, sequence_col = sequence_col)),
      .progress = TRUE
    )
  }
)

#' @rdname matchHistones
#' @examples
#' # TODO
#' NULL
# TODO rework this part!
setMethod(
  "matchHistones",
  c("SummarizedExperiment", "AAStringSetList"),
  function(
    object,
    matching_subject,
    sequence_col = "sequence"
  ) {
    which_pattern_matches <- rowData(object)[[sequence_col]] |>
      Biostrings::AAStringSet() |>
      lapply(matching_subject, Biostrings::vwhichPDict, pdict = _) |>
      lapply(\(x) unlist(x) |> unique())
    # check for family ambiguity and drop these rows if applicable
    if (length(matching_subject) > 1) {
      family_counts <- unlist(which_pattern_matches) |> table()
      family_ambiguous_row_idx <- names(family_counts)[family_counts > 1] |> as.integer()
      if (length(family_ambiguous_row_idx) > 0) {
        print(paste(
          "Dropping sequences with ambiguous family assignment:",
          paste(rowData(object)[family_ambiguous_row_idx, sequence_col], collapse = ";")
        ))
        object <- object[-family_ambiguous_row_idx, ]
      }
    }

    # get matching histone family from MIndex objects
    histone_family <- character(length = nrow(object))
    for (match_family in names(which_pattern_matches)) {
      histone_family[which_pattern_matches[[match_family]]] <- match_family
    }
    histone_family[histone_family == ""] <- NA
    rowData(object)[["histone_family"]] <- histone_family
    # flag core histones vs. H1
    core_histones <- rep.int(TRUE, length(histone_family))
    core_histones[histone_family == "H1" | is.na(histone_family)] <- FALSE
    rowData(object)[["core_histone"]] <- core_histones
    # flag histones with logical
    histone <- rep.int(FALSE, length(histone_family))
    histone[!is.na(histone_family)] <- TRUE
    rowData(object)[["histone"]] <- histone

    # match peptides to all supplied histone family AAStringSet objects
    # FUTURE TODO use vmatchPDict() instead for faster matching, but not yet implemented in Biostrings (03/2025)
    pattern_matches <- lapply(
      matching_subject,
      \(x) sapply(rowData(object)[[sequence_col]], Biostrings::vmatchPattern, x)
    )
    # check for multiple matches in one variant
    position_ambiguous_row_idx <- sapply(
      pattern_matches,
      \(x) as.logical(lapply(x, \(x) any(S4Vectors::elementNROWS(x) > 1)))
    )
    # combine information from multiple families
    position_ambiguous_row_idx <- apply(position_ambiguous_row_idx, 1, any)
    if (any(unlist(position_ambiguous_row_idx))) {
      print(paste(
        "Dropping sequences with multiple possible positions per variant:",
        paste(rowData(object)[position_ambiguous_row_idx, sequence_col], collapse = ";")
      ))
      object <- object[!position_ambiguous_row_idx, ]
      pattern_matches <- lapply(pattern_matches, \(x) x[!position_ambiguous_row_idx])
    }

    # create histone groups from MIndex objects
    histone_group <- lapply(pattern_matches, .mindex_list_to_protein_group)
    # merge matches from different families into one vector if applicable
    if (length(histone_group) > 1) {
      histone_group <- do.call(\(...) paste0(...), histone_group) |>
        sapply(\(x) if (x == "") NA else x)
    } else {
      histone_group <- sapply(histone_group[[1]], \(x) if (x == "") NA else x)
    }
    rowData(object)[["histone_group"]] <- histone_group

    # create index positions from MIndex objects
    start_idx <- lapply(pattern_matches, .mindex_list_to_idx, start = TRUE)
    # merge start indexes from different families into one vector if applicable
    if (length(start_idx) > 1) {
      start_idx <- do.call(\(...) mapply(c, ..., MoreArgs = list(use.names = FALSE)), start_idx) |>
        lapply(\(x) if (length(x) == 0) NA else x)
    } else {
      start_idx <- lapply(start_idx[[1]], \(x) if (length(x) == 0) NA else x)
    }
    rowData(object)[["start_index"]] <- start_idx
    # create end positions from MIndex objects
    end_idx <- lapply(pattern_matches, .mindex_list_to_idx, start = FALSE)
    # merge end indexes from different families into one vector if applicable
    if (length(end_idx) > 1) {
      end_idx <- do.call(\(...) mapply(c, ..., MoreArgs = list(use.names = FALSE)), end_idx) |>
        lapply(\(x) if (length(x) == 0) NA else x)
    } else {
      end_idx <- lapply(end_idx[[1]], \(x) if (length(x) == 0) NA else x)
    }
    rowData(object)[["end_index"]] <- end_idx

    object
  }
)

# Convert a list of MIndex object into a character vector of all proteins matching the
# sequence
#
# @param pattern_matches A list of MIndex objects containing the match between
# @importMethodsFrom Biostrings names
.mindex_list_to_protein_group <- function(pattern_matches) {
  mapply(
    \(x, y) paste(x[y], collapse = "/"),
    lapply(pattern_matches, names),
    lapply(pattern_matches, \(x) S4Vectors::elementNROWS(x) |> as.logical())
  )
}

# Convert a list of MIndex object into a list of integer vectors of the start/end
# indexes of the sequence in the matching proteins
#
# @param pattern_matches A list of MIndex objects containing the match between
# @param start If `TRUE`, returns the start indexes, otherwise returns the end indexes
.mindex_list_to_idx <- function(pattern_matches, start = TRUE) {
  mapply(
    \(x, y) as.integer(x[y]) - 1L, # - 1 because N-term methionine should get index 0
    lapply(pattern_matches, if (start) Biostrings::startIndex else Biostrings::endIndex),
    lapply(pattern_matches, \(x) S4Vectors::elementNROWS(x) |> as.logical())
  )
}
