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
#' @param progress Show a progress bar? Defaults to `TRUE`.
#' @examples
#' ncbtoy |>
#'   matchHistones(aligned_histones, 1)
setMethod(
  "matchHistones",
  c("QFeatures", "AAStringSetList"),
  function(
    object,
    matching_subject,
    i,
    progress = TRUE,
    sequence_col = "sequence"
  ) {
    i <- QFeatures:::.normIndex(object, i)
    if (progress) {
      pb <- progress::progress_bar$new(
        total = length(i),
        show_after = 0,
        format = "(:spin) [:bar] assay :current of :total"
      )
      pb$message("Matching sequences to histones...")
      # to ensure the bar only is removed once the loop is completed
      pb$tick(-1)
    }
    for (j in i) {
      if (progress) {
        pb$tick()
      }
      QFeatures::replaceAssay(
        object,
        matchHistones(
          object[[j]],
          matching_subject,
          progress = if (progress) pb else FALSE,
          sequence_col = sequence_col
        ),
        j
      )
    }
    # to complete the progress and remove the bar
    if (progress) {
      pb$tick()
    }
    object
  }
)

setMethod(
  "matchHistones",
  c("SummarizedExperiment", "AAStringSetList"),
  function(
    object,
    matching_subject,
    progress = TRUE,
    sequence_col = "sequence"
  ) {
    rd <- rowData(object)

    # workaround for fast matching as Biostrings::vmatchPDict() is not yet implemented (07/2025)
    # FUTURE TODO use vmatchPDict() instead
    sequences <- rd[[sequence_col]]
    # first use vwhichPDict to quickly get histone family -> feature match
    match_df <- purrr::imap(
      matching_subject,
      function(family_sequences, family) {
        pattern_matches <- Biostrings::vwhichPDict(pdict = Biostrings::AAStringSet(sequences), family_sequences) |>
          unlist() |>
          unique()
        # then do a nested vmatchPattern instead of one top-level vmatchPDict to get histone variant -> feature match
        # this is quite slow, so limit the features to match and only match families that match the feature
        # preserve feature number as name so that list_rbind can retrieve this
        df <- purrr::set_names(sequences[pattern_matches], pattern_matches) |>
          # here is the nested vmatchPattern
          purrr::map(\(pattern) Biostrings::vmatchPattern(pattern, family_sequences) |> as.data.frame()) |>
          purrr::list_rbind(names_to = "idx") |>
          mutate(idx = as.integer(.data[["idx"]]))
        if (nrow(df) == 0) {
          return(NULL)
        }
        df |>
          tibble::as_tibble() |>
          mutate(
            histone = TRUE,
            histone_family = family,
            core_histone = if (family == "H1") FALSE else TRUE,
            histone_variant = purrr::map_chr(.data[["group"]], \(x) names(family_sequences)[[x]]),
            start_index = .data[["start"]] - 1L, # do not count initiator methionine
            end_index = .data[["end"]] - 1L, # do not count initiator methionine
            .keep = "unused"
          ) |>
          dplyr::select(-c(.data[["group_name"]], .data[["width"]]))
      }
    ) |>
      dplyr::bind_rows()

    # "error" checking
    if (nrow(match_df) == 0) {
      stop("No histone sequences could be matched, please check your arguments")
    }
    fam_ambiguous_rows <- match_df |>
      dplyr::distinct(.data[["idx"]], .data[["histone_family"]]) |>
      dplyr::count(.data[["idx"]]) |>
      dplyr::filter(.data[["n"]] > 1) |>
      dplyr::pull(.data[["idx"]])
    pos_ambiguous_rows <- match_df |>
      dplyr::group_by(.data[["idx"]], .data[["histone_variant"]]) |>
      dplyr::filter(dplyr::n() > 1) |>
      dplyr::pull(.data[["idx"]]) |>
      unique()
    rows_to_drop <- c()
    if (length(fam_ambiguous_rows) > 0) {
      message(
        "Dropping sequences with ambiguous family assignment: ",
        paste(rd[fam_ambiguous_rows, sequence_col], collapse = ", ")
      )
      rows_to_drop <- c(rows_to_drop, fam_ambiguous_rows)
    }
    if (length(pos_ambiguous_rows) > 0) {
      message(
        "Dropping sequences with multiple possible positions per variant: ",
        paste(rd[pos_ambiguous_rows, sequence_col], collapse = ";")
      )
      rows_to_drop <- c(rows_to_drop, pos_ambiguous_rows)
    }
    rows_to_drop <- unique(rows_to_drop)

    # combine variants of the same family
    match_summary <- match_df |>
      dplyr::group_by(.data[["idx"]]) |>
      dplyr::summarise(
        histone = dplyr::first(.data[["histone"]]),
        histone_family = dplyr::first(.data[["histone_family"]]),
        core_histone = dplyr::first(.data[["core_histone"]]),
        histone_group = paste(.data[["histone_variant"]], collapse = "/"),
        start_index = list(.data[["start_index"]]),
        end_index = list(.data[["end_index"]]),
        .groups = "drop"
      ) |>
      as("DataFrame")

    # combine with original rowData, make sure to preserve order
    rd[["idx"]] <- as.integer(rownames(rd))
    rd <- S4Vectors::merge(rd, match_summary, by = "idx", all.x = TRUE)
    # idx is a character so sorts like 1, 10, 100, ...
    rownames(rd) <- as.character(rd[["idx"]])
    rd[["idx"]] <- NULL
    SummarizedExperiment::rowData(object) <- rd

    if (length(rows_to_drop) > 0) {
      object <- object[-rows_to_drop, ]
    }

    object
  }
) # TODO progress bar in SE, in parallel?
