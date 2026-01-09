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
#' @param mod_format The format of the modification string. Defaults to
#'   `"progenesis"`. Must be one of:
#'   * `"progenesis"` = "\[loc\] mod (aa)|..."
#'   * `"progenesis_sw"` = "\[loc\] (aa) mod|..."
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
#'
#' One additional column is added: "precursor", which is the ProForma combination of
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
#' @param i The index (`integer()`) or name (`character()`) of the assay(s) to be processed.
#' @examples
#' \dontrun{
#' qf <- matchHistones(ncbtoy, aligned_histones$unaligned, 1)
#'
#' # in this dataset, mods follow the "progenesis_sw" format
#' processMods(qf, aligned_histones, 1, mod_format = "progenesis_sw")
#'
#' # "Unmod" can be added to other amino acids and we can remove chemical artifacts
#' processMods(
#'   qf,
#'   aligned_histones,
#'   1,
#'   mod_format = "progenesis_sw",
#'   unmods = c("K", "R"),
#'   strip_mods = list(Fo = c("K"))
#' )
#'
#' # maybe we would like to write out hPTMs in full
#' processMods(
#'   qf,
#'   aligned_histones,
#'   1,
#'   mod_format = "progenesis_sw",
#'   rename_mods = list("Ac" ~ "Acetyl", "Me2" ~ "Dimethyl", "Me3" ~ "Trimethyl")
#' )
#'}
setMethod(
  "processMods",
  c("QFeatures", "list"),
  function(
    object,
    msa,
    i,
    mod_format = c("progenesis", "progenesis_sw"),
    unmods = c("K"),
    strip_mods = list(Propionyl = c("K", "N-term")),
    rename_mods = NULL
  ) {
    i <- QFeatures:::.normIndex(object, i)
    for (j in i) {
      object <- QFeatures::replaceAssay(
        object,
        processMods(
          object[[j]],
          msa,
          mod_format = mod_format,
          unmods = unmods,
          strip_mods = strip_mods,
          rename_mods = rename_mods
        ),
        j
      )
    }
    object
  }
)

#' @rdname processMods
#' @examples
#' \dontrun{
#' se <- matchHistones(ncbtoy[[1]], aligned_histones$unaligned)
#'
#' # in this dataset, mods follow the "progenesis_sw" format
#' processMods(se, aligned_histones, mod_format = "progenesis_sw")
#'
#' # "Unmod" can be added to other amino acids and we can remove chemical artifacts
#' processMods(
#'   se,
#'   aligned_histones,
#'   mod_format = "progenesis_sw",
#'   unmods = c("K", "R"),
#'   strip_mods = list(Fo = c("K"))
#' )
#'
#' # maybe we prefer to write out hPTMs in full
#' processMods(
#'   se,
#'   aligned_histones,
#'   mod_format = "progenesis_sw",
#'   rename_mods = list("Ac" ~ "Acetyl", "Me2" ~ "Dimethyl", "Me3" ~ "Trimethyl")
#' )
#'}
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

    # argument checks
    #fmt: skip
    stopifnot(
      all(c("mods", "sequence", "histone", "start_index", "histone_family", "histone_group", "feature_number", "charge") %in% colnames(rd)),
      all(c("unaligned", "msa") %in% names(msa)),
      all(purrr::map_lgl(msa, inherits, "AAStringSetList")),
      is.null(unmods) || is.character(unmods),
      is.null(strip_mods) || (is.list(strip_mods) && rlang::is_named(strip_mods) && all(purrr::map_lgl(strip_mods, is.character))),
      is.null(rename_mods) || (is.list(rename_mods) && all(purrr::map_lgl(rename_mods, rlang::is_formula)))
    )

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

    # add precursor column and use for naming
    rd$precursor <- .create_proforma(rd$sequence, mod_info$loc, mod_info$mod, rd$charge)
    stopifnot("Duplicate precursors detected" = !any(duplicated(rd$precursor)))
    rownames(object) <- rd$precursor

    # within-peptide
    mod_info$loc <- purrr::map2(mod_info$loc, rd$sequence, function(locs, seq) {
      dplyr::case_when(
        locs == "N-term" ~ 1L,
        locs == "C-term" ~ nchar(seq),
        is.na(locs) ~ NA_integer_,
        # warning "NAs introduced by coercion" due possible N-term/C-term, works regardless
        .default = suppressWarnings(as.integer(locs))
      )
    })
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
    pos_msa <- .map_var_to_msa(pos_var_mat, rd, msa_mappers)
    rd$mods_msa <- .create_mod_string(mod_info$aa, pos_msa, mod_info$mod, filter = !rd$histone)

    # within-msa-frame to within-reference for interpretability
    if ("msa_ref" %in% names(msa)) {
      ref_mappers <- purrr::map(msa$msa_ref, .mapper_from_msa) |>
        purrr::map(1) |>
        purrr::map2(msa$msa_ref, .fill_msa_gaps)
      pos_ref <- .map_msa_to_ref(pos_msa, rd$histone_family, ref_mappers)
      rd$mods_ref <- .create_mod_string(mod_info$aa, pos_ref, mod_info$mod, filter = !rd$histone)
    }

    rowData(object) <- rd

    object
  }
)

.parse_mods <- function(mod_strings, format) {
  # workaround for historic bug in parsing script...
  if (format == "progenesis_sw") {
    mod_strings <- stringr::str_replace_all(mod_strings, r"(\(ST\)Ph)", "(ST) Ph")
  }

  pattern <- switch(
    format,
    "progenesis" = r"(\[(?<loc>\d+|N-term|C-term)\] (?<mod>.+?) \(.+?\))",
    "progenesis_sw" = r"(\[(?<loc>\d+|N-term|C-term)\] \(.+?\) (?<mod>.+?)(?:\||$))"
  )
  mods <- stringr::str_match_all(mod_strings, pattern) |>
    purrr::map(\(x) list(locs = x[, "loc"], mods = x[, "mod"]))
  mods <- list(
    locs = purrr::map(mods, "locs"),
    mods = purrr::map(mods, "mods")
  ) |>
    # NA values have a name, which should be dropped for consistency
    purrr::map_depth(2, \(x) purrr::set_names(x, nm = NULL))

  # empty mod strings result in a character(0) instead of NA
  mods$locs[lengths(mods$locs) == 0] <- NA_character_
  mods$mods[lengths(mods$mods) == 0] <- NA_character_

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
  stopifnot(is.list(rename_map) && length(rename_map) > 0)
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
      # warning "NAs introduced by coercion" due to seq and loc possibly having different length, works regardless
      .default = suppressWarnings(stringr::str_sub(seq, as.integer(loc), as.integer(loc)))
    )
  })
}

.add_unmods <- function(locs, mods, sequences, unmods, is_histone) {
  aa_unmods <- unmods[unmods != "N-term" & unmods != "C-term"]
  if (length(aa_unmods) > 0) {
    pattern <- paste0("[", paste0(aa_unmods, collapse = ""), "]")
    unmod_sites <- stringr::str_locate_all(sequences, pattern) |>
      purrr::map(\(x) as.character(x[, "start"]))
  } else {
    unmod_sites <- replicate(length(sequences), character(), simplify = FALSE)
  }
  # N- or C-term will not be found through the pattern match, so add manually
  if ("N-term" %in% unmods) {
    unmod_sites <- purrr::map(unmod_sites, \(x) c(x, "N-term"))
  }
  if ("C-term" %in% unmods) {
    unmod_sites <- purrr::map(unmod_sites, \(x) c(x, "C-term"))
  }

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
      # warning "NAs introduced by coercion" due to seq and loc possibly having different length, works regardless
      .default = suppressWarnings(as.integer(all_pos))
    )

    sorted_order <- order(pos_order)
    list(loc = all_pos[sorted_order], mod = all_mods[sorted_order])
  })

  list(
    loc = purrr::map(new_info, "loc"),
    mod = purrr::map(new_info, "mod")
  )
}

.create_proforma <- function(sequences, locs, mods, charges) {
  purrr::pmap_chr(list(sequences, locs, mods, charges), function(seq, loc, mod, charge) {
    if (all(is.na(mod))) {
      return(paste0(seq, "/", charge))
    }

    # turn loc into the amino acid index of the sequence
    mod_loc <- dplyr::case_when(
      loc == "N-term" ~ 1L,
      loc == "C-term" ~ nchar(seq),
      # warning "NAs introduced by coercion" due to seq and loc possibly having different length, works regardless
      .default = suppressWarnings(as.integer(loc))
    )
    # funky looking string for later evaluation with glue::glue()
    # {m} will be substituted by the mod and {aa} by the amino acid
    subs_str <- dplyr::case_when(
      loc == "N-term" ~ "[{m}]-{aa}",
      loc == "C-term" ~ "{aa}-[{m}]",
      .default = "{aa}[{m}]"
    )

    seq_chars <- stringr::str_split_1(seq, "")
    seq_chars[mod_loc] <- purrr::pmap(list(mod_loc, subs_str, mod), function(ml, ms, m) {
      # need to define aa for glue() to use, m is already defined inside this pmap call
      aa <- seq_chars[ml]
      (\(x) NULL)(aa) # just getting rid of warning that aa is unused...
      seq_chars <- stringr::str_glue(ms)
    })

    paste(seq_chars, collapse = "") |> paste0("/", charge)
  })
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
      return(NA_integer_)
    }
    # -1 because initiator M
    outer(starts, loc, `+`) - 1L
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
  # -1 as initiator M is deemed to start at 0
  purrr::map2(aligned_lengths, is_aa, \(seq_length, aa) seq_len(seq_length)[aa] - 1L)
}

.map_var_to_msa <- function(locs, rowData, msa_mappers) {
  res <- purrr::pmap(
    list(locs, rowData$histone_family, rowData$histone_group, rowData$feature_number, rowData$sequence),
    function(loc, fam, grp, feat, seq) {
      # NA if no mods or if not a histone (so no assigned histone family)
      if (all(is.na(loc)) || is.na(fam)) {
        return(list(NA_character_))
      }

      mapper <- msa_mappers[[fam]]
      variants <- stringr::str_split_1(grp, "/")
      stopifnot(all(variants %in% names(mapper)))

      purrr::map2(
        # each element in the mapper corresponds to 1 variant
        mapper[variants],
        # each row in the loc matrix corresponds to 1 variant
        asplit(loc, MARGIN = 1),
        # + 1 to correct for initiator M (i.e. this has index 0 at column 1 of the mapper)
        \(variant_map, mod_locs) as.character(variant_map[mod_locs + 1L])
      ) |>
        # keep only unique locations for each mod
        purrr::pmap(\(...) unique(c(...)))
    }
  )

  # some sequences may still have a ptm loc mismatch after msa
  ambiguous_loc <- purrr::map_depth(res, 2, \(x) length(x) > 1) |> purrr::map(as.logical) |> purrr::map_lgl(any)
  if (any(ambiguous_loc)) {
    message("Modified sequences with ambiguous location post-alignment:\nFeat\tSequence")
    message(paste(rowData$feature_number[ambiguous_loc], rowData$sequence[ambiguous_loc], sep = "\t", collapse = "\n"))
    res[ambiguous_loc] <- purrr::map_depth(res[ambiguous_loc], 2, \(x) paste(sort(x), collapse = "/"))
  }

  res
}

.fill_msa_gaps <- function(ref_mapper, msa_ref) {
  # start with the "unaligned" (1 to 1) index map, starting at 0
  res <- as.character(seq_along(ref_mapper) - 1L)
  # gaps before initiator M (idx = 0) will have a negative index that counts down from the 0 position
  gaps_before_seq <- as.character(seq(to = -1L, length.out = ref_mapper[[1]]))
  # number of gaps after sequence ends, to fill out the msa frame
  gaps_after_seq <- as.character(msa_ref) |>
    stringr::str_extract("(-*)$") |>
    nchar()
  # then add "gaps" after amino acids
  gaps_after_aa <- (diff(ref_mapper) - 1L) |>
    append(gaps_after_seq)
  interior_gaps <- gaps_after_aa[gaps_after_aa > 0] |>
    # "R - - - K -"" becomes "1 1.1 1.2 1.3 2 2.1"
    purrr::map(\(x) c("", paste0(".", seq_len(x)))) |>
    purrr::map2(res[gaps_after_aa > 0], \(x, y) paste0(y, x))

  res[gaps_after_aa > 0] <- interior_gaps
  c(gaps_before_seq, unlist(res))
}

.map_msa_to_ref <- function(locs, families, ref_mappers) {
  purrr::map2(locs, families, function(loc, fam) {
    if (all(is.na(loc)) || is.na(fam)) {
      return(list(NA_character_))
    }

    mapper <- ref_mappers[[fam]]
    purrr::map_chr(loc, \(l) {
      msa_loc <- as.integer(stringr::str_split_1(l, "/"))
      # +1 because mapper is 1-indexed, msa loc is 0-based
      paste(mapper[msa_loc + 1], collapse = "/")
    })
  })
}
