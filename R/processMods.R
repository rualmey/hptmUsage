#' Process Histone Modifications
#'
#' @description
#' hPTMs are rewritten in a common format ("aa|loc|mod;..."), for example
#' "K|27|Me3;K|36|Ac;K|37|Unmod". This involves mapping the hPTM location from
#' within a peptide to the matching histone variant(s). It also allows to add
#' "Unmod" pseudo-hPTMs, remove uninformative hPTMs (e.g., propionylation from
#' chemical derivatization), and rename hPTMs (e.g., from "Trimethyl" to "Me3").
#'
#' @param object The `QFeatures` or `SummarizedExperiment` object in which the
#'   mods should be processed.
#' @param msa A `list()` of `AAStringSetList` objects containing unaligned,
#'   aligned, and optionally reference histone sequences. Each `AAStringSetList`
#'   element contains sequences of one histone family. For example, the
#'   (default) result from [alignHistones()].
#' @param ... Additional arguments passed to specific methods.
#' @param mod_format The format of the modification string. One of:
#' * `"progenesis"` = "[loc] mod (aa)|..."
#' * `"progenesis_sw"` = "[loc] (aa) mod|..."
#' Defaults to `"progenesis"`.
#' @param unmods A `character()` of amino acids to which "Unmod" pseudo-hPTMs
#'   should be added. Optional, `NULL` disables adding unmods.
#' @param strip_mods A named `list(mod1 = character("aa1", ...), ...)` of
#'   modifications to strip. Optional, `NULL` disables stripping mods.
#' @param rename_mods A `list()` of two-sided formulas where the LHS is the old
#'   name and the RHS the new name (e.g., `"old" ~ "new"`) specifying how to
#'   rename hPTMs. Optional.
#' @returns
#' A `QFeatures` or `SummarizedExperiment` (same as supplied) with the rowData
#' containing four new columns, where each column differs in the locations of
#' hPTMs:
#' * "mods_pep" = location within peptide
#' * "mods_var" = location within histone variant(s)
#' * "mods_msa" = location within histone variant(s) after multiple sequence alignment (see [alignHistones()])
#' * "mods_ref" = corresponding location within the aligned reference sequence (see [alignHistones()]), only when `msa` contains "msa_ref"
#' One more column is added: "precursor", which is the ProForma combination of
#' sequence, mods, and charge.
#' @export
#' @name processMods
setGeneric(
  "processMods",
  function(
    object,
    msa,
    ...,
    mod_format = c("progenesis", "progenesis_sw"),
    unmods = c("K"),
    strip_mods = list(Propionyl = c("K", "N-term")),
    rename_mods = NULL
  ) {
    standardGeneric("processMods")
  },
  signature = c("object", "msa")
)

#' @rdname processMods
#' @examples
#' \dontrun{
#' se <- ncbtoy[[1]] |>
#'   matchHistones(aligned_histones$unaligned)
#'
#' # in this dataset, mods follow the "progenesis_sw" format
#' process_mods(se, aligned_histones, 1, mod_format = "progenesis_sw")
#'
#' # "Unmod" can be added to other amino acids and we can remove chemical artifacts
#' process_mods(
#'   se,
#'   aligned_histones,
#'   1,
#'   mod_format = "progenesis_sw",
#'   unmods = c("K", "R"),
#'   strip_mods = list(Fo = c("K"))
#' )
#'
#' # maybe we would like to write out hPTMs in full
#' process_mods(
#'   se,
#'   aligned_histones,
#'   1,
#'   mod_format = "progenesis_sw",
#'   rename_mods = list("Ac" ~ "Acetyl", "Me2" ~ "Dimethyl", "Me3" ~ "Trimethyl")
#' )
#'}
# TODO
setGeneric(
  "processMods2",
  function(
    object,
    msa,
    ...,
    mod_format = c("progenesis", "progenesis_sw"),
    unmods = c("K"),
    strip_mods = list(Propionyl = c("K", "N-term")),
    rename_mods = NULL
  ) {
    standardGeneric("processMods2")
  },
  signature = c("object", "msa")
)
setMethod(
  "processMods2",
  c("SummarizedExperiment", "list"),
  function(
    object,
    msa,
    mod_format = c("progenesis", "progenesis_sw"),
    unmods = c("K"),
    strip_mods = list(Propionyl = c("K", "N-term")),
    rename_mods = NULL
  ) {
    # argument checking
    # TODO

    # retrieve positions of all mods
    # mods <- rowData(object)[["mods"]] |>
    #   stringr::str_split(stringr::fixed("|"))
    # positions <- lapply(mods, function(x) {
    #   stringr::str_split_i(x, pattern = " ", i = 1) |>
    #     stringr::str_sub(start = 2L, end = -2L)
    # })

    # # extract proper names of modifications
    # if (mod_format) {
    #   # change to progenesis
    #   mods <- lapply(mods, stringr::str_split_i, pattern = " ", i = 3) # TODO regex "", see below
    # } else {
    #   mods <- lapply(mods, stringr::str_extract, pattern = r"(\] ([^(]+) \()", group = 1)
    # }
    # unique_mods <- sort(unique(unlist(mods)), na.last = NA)
    # print(paste("Unique modifications found:", paste0(unique_mods, collapse = ";")))
    # # change mod names using if requested
    # if (!is.null(rename_mods)) {
    #   print(paste("Renaming mods:", paste(names(rename_mods), rename_mods, sep = "=", collapse = "; ")))
    #   mods_not_in_lookup <- unique_mods[!unique_mods %in% names(rename_mods)]
    #   names(mods_not_in_lookup) <- mods_not_in_lookup
    #   rename_mods <- c(rename_mods, mods_not_in_lookup)
    #   mods <- lapply(mods, \(x) unname(rename_mods[x]))
    # }

    # # extract actual amino acid that each mod is located on
    # amino_acids <- aa_from_position(positions, rowData(object)[["sequence"]])

    # strip mods e.g. propionylation ~ fixed mod
    # if (!is.null(strip_mods)) {
    # print(paste(
    #   "Stripping mods:",
    #   paste0(
    #     paste(names(strip_mods), lapply(strip_mods, \(x) paste0("(", paste0(x, collapse = ", "), ")"))),
    #     collapse = "; "
    #   )
    # ))
    #   to_strip_mask <- mapply(
    #     \(mod, aa) mapply(\(x, y) !y %in% strip_mods[[x]], mod, aa),
    #     mods,
    #     amino_acids
    #   )
    #   mods <- mapply(\(x, y) x[y], mods, to_strip_mask) |>
    #     lapply(\(x) if (length(x) == 0) NA else x)
    #   positions <- mapply(\(x, y) x[y], positions, to_strip_mask) |>
    #     lapply(\(x) if (length(x) == 0) NA else x)
    #   amino_acids <- mapply(\(x, y) x[y], amino_acids, to_strip_mask) |>
    #     lapply(\(x) if (length(x) == 0) NA else x)
    # }

    # # add unmods where required
    # if (!is.null(unmods)) {
    #   print(paste("Adding unmods to:", paste0(unmods, collapse = "; ")))
    #   unmods <- rowData(object)[["sequence"]] |>
    #     stringr::str_split("") |>
    #     sapply(\(x) which(x %in% unmods))
    #   unmods <- mapply(
    #     function(x, y, z, a) {
    #       # ignore co-extracts
    #       if (!a) {
    #         return(list(pos = y, mods = z))
    #       }
    #       # sorted vector of positions: N-term, 1, 2, ..., C-term
    #       n_term <- if ("N-term" %in% y) "N-term" else character(0)
    #       c_term <- if ("C-term" %in% y) "C-term" else character(0)
    #       int_positions <- y[!y %in% c("N-term", "C-term")] |>
    #         as.integer() |>
    #         c(x) |>
    #         unique() |>
    #         sort()
    #       all_positions <- c(n_term, as.character(int_positions), c_term)
    #       # amino acid vector with unmod in new positions
    #       mods_unmods <- rep.int("Unmod", length(all_positions))
    #       mods_unmods[match(y, all_positions)] <- z
    #       # return both new positions and mods
    #       list(pos = all_positions, mods = mods_unmods)
    #     },
    #     unmods,
    #     positions,
    #     mods,
    #     rowData(object)[["histone"]],
    #     SIMPLIFY = FALSE
    #   )
    #   positions <- lapply(unmods, `[[`, "pos")
    #   mods <- lapply(unmods, `[[`, "mods")
    # }

    # # redo extract actual amino acid that each mod is located on due to added unmods/...
    # amino_acids <- aa_from_position(positions, rowData(object)[["sequence"]])

    # # create precursor annotation
    # precursors <- mapply(
    #   create_proforma,
    #   rowData(object)[["sequence"]],
    #   positions,
    #   mods
    # ) |>
    #   paste0(rowData(object)[["charge"]])

    # # add to rowdata
    # rowData(object)[["precursor"]] <- precursors

    # # return parsed within-peptide mod string
    # within_peptide_mod_string <- mod_string_from_parts(
    #   amino_acids,
    #   positions,
    #   mods
    # )
    # # add to rowdata
    # rowData(object)[["mods_pep"]] <- within_peptide_mod_string

    # # map mod positions from within-peptide to within-variant (locations can mismatch between variants)
    # positions <- mapply(
    #   \(x, y) stringr::str_replace_all(x, pattern = c("N-term" = "1", "C-term" = as.character(nchar(y)))),
    #   positions,
    #   rowData(object)[["sequence"]]
    # ) |>
    #   sapply(as.integer)
    # # for every seq returns a matrix with rows = variant and cols = mod
    # positions <- mapply(
    #   \(x, y) matrix(rep.int(x, length(y)), ncol = length(x), byrow = TRUE) + y - 1, # -1 because initiator M
    #   positions,
    #   rowData(object)[["start_index"]]
    # )
    # positions[is.na(mods) | is.na(rowData(object)[["start_index"]])] <- NA
    # # and collapse the matrix to get the unique PTM locs in the variant group
    # within_variant_position <- lapply(positions, function(x) {
    #   if (!all(is.na(x))) {
    #     unique(x) |> apply(MARGIN = 2, FUN = \(y) paste(sort(y), collapse = "/"))
    #   } else {
    #     NA
    #   }
    # })
    # # create a mod string from within-protein positions
    # within_variant_mod_string <- mod_string_from_parts(
    #   amino_acids,
    #   within_variant_position,
    #   mods,
    #   filter = !rowData(object)[["histone"]]
    # )
    # # add to rowdata
    # rowData(object)[["mods_var"]] <- within_variant_mod_string

    # map within-variant positions to within-msa-frame positions
    mapper <- lapply(msa$msa, mapper_from_msa)
    mapped_positions <- mapply(
      function(x, y, z) {
        if (is.na(x) || all(is.na(y))) {
          return(NA) # not a histone or not modified
        }
        # the row index of the variant (column) mapper that matches the within-variant index returns the within-frame idx
        idx_mask <- mapply(\(a, b) mapper[[x]][, b] == a, y, z, SIMPLIFY = FALSE)
        idx_mask <- sapply(idx_mask, which)
        matrix(idx_mask, ncol = ncol(y)) - 1 # -1 as initiator M gets index 1 in lookup table but should get position 0
      },
      rowData(object)[["histone_family"]],
      positions,
      stringr::str_split(rowData(object)[["histone_group"]], "/")
    ) |>
      lapply(unique) # collapse the matrix
    # check for ambiguous PTM location within the alignment frame
    ambiguous_ptm_loc <- lapply(mapped_positions, \(x) nrow(x) > 1) |>
      sapply(any)
    if (any(ambiguous_ptm_loc)) {
      print(paste(
        "Modified sequences with ambiguous location post-alignment:",
        rowData(object)[ambiguous_ptm_loc, c("histone_family", "sequence")] |>
          unique() |>
          apply(1, \(x) paste(x, collapse = "-")) |>
          paste0(collapse = "; "),
        "at rows:",
        paste(which(ambiguous_ptm_loc), collapse = "; ")
      ))
    }
    # and collapse the matrix to get the unique PTM locs in the variant group
    mapped_positions <- lapply(mapped_positions, function(x) {
      if (!all(is.na(x))) {
        apply(x, MARGIN = 2, FUN = \(y) paste(sort(as.integer(unique(y))), collapse = "/"))
      } else {
        NA
      }
    })
    # create a mod string from within-msa-frame positions
    within_frame_mod_string <- mod_string_from_parts(
      amino_acids,
      mapped_positions,
      mods,
      filter = !rowData(object)[["histone"]]
    )
    # add to rowdata
    rowData(object)[["mods_msaframe"]] <- within_frame_mod_string

    # translate within-msa-frame positions to refseq-relative locations if desired
    if ("msa_ref" %in% names(msa)) {
      refseq_mapper <- lapply(msa$msa_ref, mapper_from_msa) |>
        lapply(fill_msa_gaps)
      # use the refseq mapper to translate within-msa-frame to refseq-relative locations
      refseq_positions <- mapply(
        function(x, y) {
          if (is.na(x) || all(is.na(y))) {
            return(NA) # not a histone or not modified
          }
          # need to account for possible location ambiguity e.g. 145/170 -> 89/119
          pos <- stringr::str_split(y, "/")
          sapply(pos, \(z) paste(refseq_mapper[[x]][as.integer(z) + 1], collapse = "/")) # -1 due to initiator M
        },
        rowData(object)[["histone_family"]],
        mapped_positions
      )
      # create a mod string from within-protein positions
      within_refseq_mod_string <- mod_string_from_parts(
        amino_acids,
        refseq_positions,
        mods,
        filter = !rowData(object)[["histone"]]
      )
      # add to rowdata
      rowData(object)[["mods_ref"]] <- within_refseq_mod_string
    }

    object
  }
)

# # Create a mapper for aligned sequences (containing "-" symbols as gaps) back to
# # unaligned indexes: the aligned sequence is translated so that, at every index
# # in the msa frame, a gap becomes -1 and the initiator methionine gets 0 with
# # every subsequent amino acid increasing that index by 1. E.g. "--MTE-ST-"
# # becomes c(-1, -1, 0, 1, 2, -1, 3, 4, -1)
# mapper_from_msa <- function(msa) {
#   msa <- as.character(msa)
#   histone_names <- names(msa)
#   msa <- stringr::str_split(msa, "")
#   names(msa) <- histone_names
#   mapper <- sapply(msa, function(x) {
#     # initialize
#     unaligned_idx <- integer(length = length(x))
#     counter <- 0L
#     # - becomes -1 and amino acids increase the index by +1, starting from 0
#     for (i in seq_along(x)) {
#       if (x[i] == "-") {
#         unaligned_idx[i] <- -1L
#       } else {
#         unaligned_idx[i] <- counter
#         counter <- counter + 1L
#       }
#     }
#     unaligned_idx
#   })

#   mapper
# }

# helper function ... TODO
fill_msa_gaps <- function(msa_seq) {
  # gaps before 0 have a negative index that counts down from the 0 position
  seq_start <- which(msa_seq == 0)
  msa_seq[1:seq_start] <- -(seq_start - 1):0
  # find the -1 values that come after 0 (the start of the sequence)
  zero_index <- which(msa_seq == 0)
  after_zero_index <- vector(length = length(msa_seq))
  after_zero_index[zero_index:length(after_zero_index)] <- TRUE
  minus1_after_zero_index <- which(msa_seq == -1 & after_zero_index)
  # -1 after 0 becomes "preceding_position.position_in_-1_stretch"
  if (length(minus1_after_zero_index) > 0) {
    canonical_pos_str <- vector("character", length = length(msa_seq))
    for (i in seq_along(msa_seq)) {
      if (i %in% minus1_after_zero_index) {
        # no need to set the variables first as they will always be set
        minus1_stretch_pos <- minus1_stretch_pos + 1
        canonical_pos_str[[i]] <- paste(current_pos, as.character(minus1_stretch_pos), sep = ".")
      } else {
        minus1_stretch_pos <- 0
        current_pos <- as.character(msa_seq[[i]])
        canonical_pos_str[[i]] <- current_pos
      }
    }
  } else {
    canonical_pos_str <- as.character(msa_seq)
  }
  canonical_pos_str
}

# helper function ... TODO
mod_string_from_parts <- function(amino_acids, positions, mods, filter = NULL) {
  mod_string <- mapply(
    \(x, y, z) paste(x, y, z, sep = "|"),
    amino_acids,
    positions,
    mods
  ) |>
    sapply(paste, collapse = ";")
  # only change mod string where applicable e.g. histones only
  if (!is.null(filter)) {
    mod_string[which(filter)] <- NA
  }
  # no mods should return NA
  mod_string[which(sapply(mods, \(x) all(is.na(x))))] <- NA
  mod_string
}

# helper function ... TODO
aa_from_position <- function(positions, sequences) {
  mapply(
    function(x, y) {
      if (all(is.na(x))) {
        # no mods
        NA
      } else {
        aa_char <- character(length(x))
        for (j in seq_along(x)) {
          # special handling for N- and C-term mods
          aa_char[j] <- if (x[j] == "N-term") {
            "N-term"
          } else if (x[j] == "C-term") {
            "C-term"
          } else {
            # get amino acid character at position instead of relying on the amino acid from the mod name
            stringr::str_sub(y, start = as.integer(x[j]), end = as.integer(x[j]))
          }
        }
        aa_char
      }
    },
    positions,
    sequences
  )
}

# helper function ... TODO
create_proforma <- function(aa_seq, positions, mods) {
  # unchanged if no mods
  if (all(is.na(mods))) {
    return(aa_seq)
  }
  aa_seq <- stringr::str_split(aa_seq, "")[[1]]
  # N-term or C-term mod needs special handling
  if (positions[[1]] == "N-term") {
    aa_seq[1] <- paste0("[", mods[1], "]-", aa_seq[1])
  }
  if (tail(positions, n = 1) == "C-term") {
    aa_seq[length(aa_seq)] <- paste0(aa_seq[length(aa_seq)], "-[", mods[length(mods)], "]")
  }
  # add side chain mods
  side_chain_mods <- which(positions != "N-term" & positions != "C-term")
  for (i in side_chain_mods) {
    pos <- as.numeric(positions[i])
    aa_seq[pos] <- paste0(aa_seq[pos], "[", mods[i], "]")
  }

  # return collapsed string
  paste(aa_seq, collapse = "")
}

setMethod(
  "processMods",
  c("SummarizedExperiment", "list"),
  function(
    object,
    msa,
    mod_format = c("progenesis", "progenesis_sw"),
    unmods = c("K"),
    strip_mods = list(Propionyl = c("K", "N-term")),
    rename_mods = NULL
  ) {
    rd <- rowData(object)

    # parse mods
    mod_format <- match.arg(mod_format)
    parsed <- .parse_mods(rd$mods, mod_format)
    parsed$mods <- .rename_mods(parsed$mods, rename_mods)
    mod_info <- .process_modifications(
      parsed$locs,
      parsed$mods,
      rd$sequence,
      strip_mods,
      unmods,
      rd$histone
    )

    # within-peptide
    rd$mods_pep <- .create_mod_string(mod_info$aa, mod_info$loc, mod_info$mod)

    # within-peptide to within-variant
    pos_var_mat <- .map_pep_to_var(mod_info$loc, rd$start_index, rd$sequence)
    pos_var <- purrr::modify_if(
      pos_var_mat,
      \(p) !all(is.na(p)),
      \(p) apply(p, 2, \(x) paste(sort(unique(x)), collapse = "/"))
    )
    rd$mods_var <- .create_mod_string(mod_info$aa, pos_var, mod_info$mod, filter = !rd$histone)

    # within-variant to within-msa-frame
    # mapper translates aligned sequence index back to unaligned index
    msa_mappers <- purrr::map(msa$msa, .mapper_from_msa)
    pos_msa <- .map_var_to_msa(pos_var_mat, rd$histone_family, rd$histone_group, msa_mappers)
    rd$mods_msa <- .create_mod_string(mod_info$aa, pos_msa, mod_info$mod, filter = !rd$histone)
    return(rd)

    # within-msa-frame to within-reference for interpretability
    if ("msa_ref" %in% names(msa)) {
      ref_mappers <- purrr::map(msa$msa_ref, .mapper_from_msa) |>
        purrr::map(dplyr::first) |>
        purrr::map(.fill_msa_gaps)
      pos_ref <- .map_msa_to_ref(pos_msa, rd$histone_family, ref_mappers)
      rd$mods_ref <- .create_mod_string(mod_info$aa, pos_ref, mod_info$mod, filter = !rd$histone)
    }

    # add precursor column
    rd$precursor <- .create_proforma(rd$sequence, mod_info$loc, mod_info$mod, rd$charge)

    rowData(object) <- rd

    object
  }
)

.parse_mods <- function(mod_strings, format) {
  mods_split <- stringr::str_split(mod_strings, stringr::fixed("|"))
  pattern <- switch(
    format,
    "progenesis" = r"(^\[(?<loc>\d+|N-term|C-term)\] (?<mod>.+) \(.+\)$)",
    "progenesis_sw" = r"(^\[(?<loc>\d+|N-term|C-term)\] \(.+\) (?<mod>.+)$)"
  )
  mods <- purrr::map(mods_split, \(x) stringr::str_match(x, pattern))
  mods <- list(
    locs = purrr::map(mods, \(x) x[, "loc"]),
    mods = purrr::map(mods, \(x) x[, "mod"])
  ) |>
    # NA values have a name, which should be dropped for consistency
    purrr::map_depth(2, \(x) purrr::set_names(x, nm = NULL))
  # for easy renaming / sanity check
  message(
    "Unique modifications found: ",
    paste(sort(unique(unlist(mods$mods)), na.last = NA), collapse = "; ")
  )

  mods
}

.rename_mods <- function(mods, rename_map) {
  if (is.null(rename_map)) {
    return(mods)
  }
  message(
    "Renaming mods: ",
    paste(rename_map, collapse = "; ")
  )
  purrr::map(mods, \(x) dplyr::case_match(x, !!!rename_map, .default = x))
}

.process_modifications <- function(locs, mods, sequences, strip_mods, unmods, is_histone) {
  aa <- .aa_from_position(locs, sequences)

  # first strip mods, as "Unmod" may need to be added to the newly unmodified amino acids e.g. K-prop -> K-unmod
  if (!is.null(strip_mods)) {
    message(
      "Stripping mods: ",
      paste0(names(strip_mods), " on ", purrr::map_chr(strip_mods, \(x) paste0(x, collapse = ", ")), collapse = "; ")
    )
    # lookup tibble
    strip_df <- tibble::enframe(strip_mods, name = "mod", value = "aa") |> tidyr::unnest("aa")

    stripped <- purrr::pmap(list(locs, mods, aa), function(p, m, a) {
      # unchanged if no mods
      if (all(is.na(m))) {
        return(list(loc = p, mod = m, aa = a))
      }
      to_keep <- !paste(m, a) %in% paste(strip_df$mod, strip_df$aa)
      list(loc = p[to_keep], mod = m[to_keep], aa = a[to_keep])
    })
    locs <- purrr::map(stripped, "loc")
    mods <- purrr::map(stripped, "mod")
    aa <- purrr::map(stripped, "aa")
  }

  if (!is.null(unmods)) {
    message("Adding Unmods to: ", paste(unmods, collapse = ", "))
    unmod_info <- .add_unmods(locs, mods, sequences, unmods, is_histone)
    locs <- unmod_info$loc
    mods <- unmod_info$mod
    # redo extract actual amino acids due to newly added unmods
    aa <- .aa_from_position(locs, sequences)
  }

  list(loc = locs, mod = mods, aa = aa)
}

.aa_from_position <- function(locs, sequences) {
  purrr::map2(locs, sequences, function(loc, seq) {
    if (all(is.na(loc))) {
      return(NA_character_)
    }
    dplyr::case_when(
      loc == "N-term" ~ "N-term",
      loc == "C-term" ~ "C-term",
      .default = stringr::str_sub(seq, as.integer(loc), as.integer(loc))
    )
  })
}

.add_unmods <- function(locs, mods, sequences, unmods, is_histone) {
  pattern <- paste0("[", paste0(unmods, collapse = ""), "]")
  unmod_sites <- stringr::str_locate_all(sequences, pattern) |>
    purrr::map(\(x) as.character(x[, "start"]))
  new_info <- purrr::pmap(list(locs, mods, unmod_sites, is_histone), function(l, m, u, h) {
    # do not handle co-extracts (or if no unmods need to be added)
    if (!isTRUE(h) || length(u) == 0) {
      return(list(loc = l, mod = m))
    }

    unmod_to_add <- u[!u %in% l]
    all_pos <- c(l, unmod_to_add)
    # l can be NA, so need to remove this
    all_pos <- all_pos[!is.na(all_pos)]
    all_mods <- c(m, rep.int("Unmod", length(unmod_to_add)))
    # m can be NA, so need to remove this
    all_mods <- all_mods[!is.na(all_mods)]

    # sort mods ascending except N-term first and C-term last
    pos_order <- dplyr::case_when(
      all_pos == "N-term" ~ -Inf,
      all_pos == "C-term" ~ Inf,
      .default = as.integer(all_pos)
    )

    sorted_order <- order(pos_order)
    list(loc = all_pos[sorted_order], mod = all_mods[sorted_order])
  })

  list(
    loc = purrr::map(new_info, "loc"),
    mod = purrr::map(new_info, "mod")
  )
}

.create_mod_string <- function(amino_acids, locs, mods, filter = NULL) {
  mod_string <- purrr::pmap_chr(list(amino_acids, locs, mods), function(aa, loc, mod) {
    if (all(is.na(mod))) {
      return(NA_character_)
    }
    paste(aa, loc, mod, sep = "|", collapse = ";")
  })
  # can be used to return NA when not a histone, for example, as these mods do not need to be mapped
  if (!is.null(filter)) {
    mod_string[filter] <- NA_character_
  }
  mod_string
}

.map_pep_to_var <- function(locs, start_indices, sequences) {
  purrr::pmap(list(locs, start_indices, sequences), function(loc, starts, seq) {
    # NA if no mods or if not a histone (so no starts)
    if (all(is.na(loc)) || all(is.na(starts))) {
      return(NA_character_)
    }
    # get within-peptide index
    loc_num <- dplyr::case_when(
      loc == "N-term" ~ 1L,
      loc == "C-term" ~ nchar(seq),
      .default = as.integer(loc)
    )
    # -1 because initiator M
    outer(starts, loc_num, `+`) - 1L
  })
}

# returns a list where each element (histone variant) is an integer vector
# indexing the vector by within-variant location returns the index at the
# corresponding location within-msa-frame e.g. result$H31_HUMAN[[27]] -> idx_in_msa_frame
.mapper_from_msa <- function(msa) {
  aligned_lengths <- Biostrings::width(msa) |>
    purrr::set_names(names(msa))
  is_aa <- as.character(msa) |>
    stringr::str_split("") |>
    purrr::map(\(x) x != "-")
  purrr::map2(aligned_lengths, is_aa, \(seq_length, aa) seq_len(seq_length)[aa])
}

.map_var_to_msa <- function(locs, families, variant_groups, msa_mappers) {
  purrr::pmap(list(locs, families, variant_groups), function(loc, fam, grp) {
    # NA if no mods or if not a histone (so no assigned histone family)
    if (all(is.na(loc)) || is.na(fam)) {
      return(NA_character_)
    }

    mapper <- msa_mappers[[fam]]
    variants <- stringr::str_split_1(grp, "/")
    mapped_loc <- purrr::map2(
      # each element in the mapper corresponds to 1 variant
      mapper[variants],
      # each row in the loc matrix corresponds to 1 variant
      asplit(loc, MARGIN = 1),
      \(variant_map, mod_locs) variant_map[mod_locs]
    ) |>
      # keep only unique locations for each mod
      purrr::pmap(\(...) unique(c(...)))

    # notify if some sequences still give ptm loc mismatch after msa
    ambiguous_loc <- purrr::map_lgl(mapped_loc, \(x) length(x) > 1) |>
      any()
    if (ambiguous_loc) {
      message("yo", grp) # TODO
      mapped_loc <- purrr::map(mapped_loc, \(x) paste(sort(x), collapse = "/"))
    }

    mapped_loc

    # return(mapper[variants])
    # return(purrr::map2(
    #   # each column in the loc matrix corresponds to 1 mod, so this creates a list of locs per mod
    #   asplit(loc, MARGIN = 2),
    #   # then map each mod with its corresponding mapper
    #   mapper[variants],
    #   \(mod_locs, map) map[[mod_locs]]
    # ))
    # return(
    #   apply(loc, 2, \(x) purrr::map2(mapper[variants], x, \(y, z) y[[z]]), simplify = FALSE)
    # )
    # purrr::map_depth(mapper[variants], 2, )

    # purrr::map(loc, \(l) {
    #   var_loc <- as.integer(stringr::str_split_1(l, "/"))
    #   # mapper has variants as rows and msa-frame-idx as columns
    #   # +1 because mapper is 1-indexed, protein position is 0-based for M
    #   msa_loc <- mapper[variants, var_loc + 1, drop = FALSE] |>
    #     as.vector() |>
    #     unique()
    #   paste(sort(msa_loc[msa_loc != -1]), collapse = "/")
    # })
  })
}

.fill_msa_gaps <- function(ref_mapper_vec) {
  result <- as.character(ref_mapper_vec)
  is_gap <- ref_mapper_vec == -1
  if (!any(is_gap)) {
    return(result)
  }

  first_aa_idx <- min(which(ref_mapper_vec >= 0), na.rm = TRUE)
  is_gap[seq_len(first_aa_idx - 1)] <- FALSE
  if (!any(is_gap)) {
    return(result)
  }

  gap_id <- rle(is_gap)$lengths |> (function(x) rep(seq_along(x), x))()
  preceding_val <- (function(x) x[cummax(seq_along(x) * !is_gap)])(ref_mapper_vec)
  gap_count <- ave(is_gap, gap_id, FUN = cumsum)

  result[is_gap] <- paste(preceding_val[is_gap], gap_count[is_gap], sep = ".")
  result
}

.map_msa_to_ref <- function(locs, families, ref_mappers) {
  purrr::pmap(list(locs, families), function(loc, fam) {
    if (all(is.na(loc)) || is.na(fam)) {
      return(NA_character_)
    }
    mapper <- ref_mappers[[fam]]
    purrr::map_chr(loc, \(l) {
      if (l == "") {
        return("")
      }
      msa_loc <- as.integer(stringr::str_split_1(l, "/"))
      # +1 because mapper is 1-indexed, msa loc is 0-based
      paste(mapper[msa_loc + 1], collapse = "/")
    })
  })
}

.create_proforma <- function(sequences, locs, mods, charges) {
  purrr::pmap_chr(list(sequences, locs, mods, charges), function(seq, loc, mod, charge) {
    if (all(is.na(mod))) {
      return(paste0(seq, "/", charge))
    }

    seq_chars <- stringr::str_split(seq, "")[[1]]
    mods_df <- tibble::tibble(loc = loc, mod = mod)

    side_mods <- mods_df |> dplyr::filter(!.data$loc %in% c("N-term", "C-term"))
    if (nrow(side_mods) > 0) {
      side_mods <- side_mods |> dplyr::mutate(loc_idx = as.integer(.data$loc))
      for (i in seq_len(nrow(side_mods))) {
        p <- side_mods$loc_idx[i]
        m <- side_mods$mod[i]
        seq_chars[p] <- paste0(seq_chars[p], "[", m, "]")
      }
    }

    proforma_str <- paste(seq_chars, collapse = "")
    if ("N-term" %in% mods_df$loc) {
      proforma_str <- paste0("[", mods_df$mod[mods_df$loc == "N-term"], "]-", proforma_str)
    }
    if ("C-term" %in% mods_df$loc) {
      proforma_str <- paste0(proforma_str, "-[", mods_df$mod[mods_df$loc == "C-term"], "]")
    }

    paste0(proforma_str, "/", charge)
  })
}
