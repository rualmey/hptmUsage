#' Get histone sequences from UniProt
#'
#' Request histone families to the UniProt API and return an AAStringSetList of
#' all requested families.
#'
#' @param histone_families A `character()` of histone families to request. One
#'   or more of `"H1"`, `"H2A"`, `"H2B"`, `"H3"`, or `"H4"`. Requests all
#'   histone families by default.
#' @param query A `character()` containing all parts of the target query except
#'   the "family" field. For query fields see https://www.uniprot.org/help/query-fields.
#'   See also [UniProt.ws::queryUniProt()].
#' @param name_field A `character(1)` specifying the UniProt field to use for
#'   naming the sequences. For a full list of possible return fields see
#'   https://www.uniprot.org/help/return_fields.
#' @param collapse A `character(1)` indicating how to collapse multiple query
#'   clauses: either " OR " or " AND ". See also [UniProt.ws::queryUniProt()].
#' @returns An `AAStringSetList` where each element corresponds to a histone
#'   family and is an `AAStringSet` object of matching histone variant sequences.
#' @export
histonesFromUniprot <- function(
  histone_families = c("H1", "H2A", "H2B", "H3", "H4"),
  query = c("organism_id:9606", "reviewed:true"),
  name_field = "id",
  collapse = " AND "
) {
  # check arguments
  stopifnot(all(histone_families %in% c("H1", "H2A", "H2B", "H3", "H4")))
  stopifnot(collapse == " AND " || collapse == " OR ")

  # lookup table to translate family into the correct query term, see UniProt "SIMILARITY comments: index"
  # at https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/similar
  family_lookup <- c(
    "H1" = 'family:"histone H1/H5 family"',
    "H2A" = 'family:"histone H2A family"',
    "H2B" = 'family:"histone H2B family"',
    "H3" = 'family:"histone H3 family"',
    "H4" = 'family:"histone H4 family"'
  )

  # each histone family is individually requested to UniProt
  family_lookup[histone_families] |>
    purrr::map(\(family_query) {
      res <- UniProt.ws::queryUniProt(
        query = c(query, family_query),
        fields = c("sequence", name_field),
        collapse = collapse
      )
      # name_field "id" results in "Entry.Name", so use column index instead
      setNames(Biostrings::AAStringSet(res$Sequence), res[[2]])
    }) |>
    Biostrings::AAStringSetList()
}

#' Align histone sequences using MAFFT.
#'
#' @param unaligned_histones An `AAStringSetList` object containing sequences of
#'  each histone family to align. If left empty, a default call to
#'  `histonesFromUniprot()` is made, i.e., all five histone families from
#'  human and reviewed.
#' @param return_alignment Whether or not to return the aligned fasta.
#' @param output_path A `character()` of the path where to store the aligned
#'  fasta. Ignored if `return_alignment = FALSE`.
#' @param overwrite Overwrite aligned fasta if already present at
#'  `output_path`? Ignored if `return_alignment = FALSE`.
#' @param use_profiles Use the curated alignments from HistoneDB 2.0 (1) as a
#'  high-quality base MSA to which new sequences are added, ensuring consistent
#'  and accurate alignment.
#'
#'  Can be `TRUE` to use profiles of H1, H2A, H2B, H3, H4, respectively, or can
#'  be a `character()` specifying the histone families to use in order, e.g.,
#'  `c("H3", "H4")` of specifying paths to custom pre-aligned `.fasta` files.
#'
#'  (1) Draizen, E. J.; Shaytan, A. K.; Mariño-Ramírez, L.; Talbert, P. B.;
#'  Landsman, D.; Panchenko, A. R. HistoneDB 2.0: A Histone Database with
#'  Variants—an Integrated Resource to Explore Histones and Their Variants.
#'  Database (Oxford) 2016, 2016, baw014. https://doi.org/10.1093/database/baw014.
#' @param nondefault_refseq_names Names of the reference sequences for each histone family as
#'  present in the profiles, with `NULL` using H1.1 and canonical H2A/H2B/H3/H4.
#'  Be aware that changing this can change hPTM location definition! Ignored if
#'  not using profiles.
#' @param return_only_original ...
#' @param quiet ...
#' @param num_cores ...
#' @param ... Additional arguments passed to MAFFT, overriding defaults.
#' @returns A named list of "unaligned" and "msa", each a named `AAStringSetList`
#'  object, respectively.
#' @export
alignHistones <- function(
  unaligned_histones = NULL,
  return_alignment = TRUE,
  output_path = "./out/msa/",
  overwrite = FALSE,
  use_profiles = TRUE,
  nondefault_refseq_names = NULL,
  return_only_original = TRUE,
  quiet = TRUE,
  num_cores = -1,
  ...
) {
  # # check MAFFT installation
  # mafft_binary <- Sys.which("mafft")
  # if (all(mafft_binary == "")) {
  #   stop("MAFFT was not found, please ensure that MAFFT was properly installed and can be found on PATH")
  # }
  # # check arguments
  # stopifnot(class(unaligned_histones) == "AAStringSetList")
  # stopifnot(all(lapply(unaligned_histones, class) == "AAStringSet"))
  # output_path <- trimws(output_path, which = "right", whitespace = "/")
  # # if NULL input, retrieve histone fasta files for every family with default parameters
  # if (is.null(unaligned_histones)) {
  #   unaligned_histones <- histonesFromUniprot()
  # }
  # # if NULL as refseqs, use defaults = H1.1 and canonical H2A/H2B/H3/H4 from HistoneDB 2.0
  # if (is.null(nondefault_refseq_names)) {
  #   ref_seqs <- c(
  #     H1 = "generic_H1|Homo|NP_005316.1 Homo|NP_005316.1|generic_H1 Homo_generic_H1_4885373", # H11_HUMAN
  #     H2A = "canonical_H2A|Homo|NP_066390.1 Homo|NP_066390.1|canonical_H2A Homo_canonical_H2A_10645195", # H2A1B_HUMAN
  #     H2B = "canonical_H2B|Homo|NP_066402.2 Homo|NP_066402.2|canonical_H2B Homo_canonical_H2B_20336754", # H2B1J_HUMAN
  #     H3 = "canonical_H3|Homo|NP_003520.1 Homo|NP_003520.1|canonical_H3 Homo_canonical_H3_4504281", # H31_HUMAN
  #     H4 = "canonical_H4|Homo|NP_001029249.1 Homo|NP_001029249.1|canonical_H4 Homo_canonical_H4_77539758" # H4_HUMAN
  #   )
  # } else {
  #   stopifnot(is.character(nondefault_refseq_names))
  #   ref_seqs <- nondefault_refseq_names
  # }

  # # pre-aligned profiles
  # if (isTRUE(use_profiles)) {
  #   msa_profile <- hptmUsageData(paste0(c("H1", "H2A", "H2B", "H3", "H4"), ".fasta"))
  #   stopifnot(length(unaligned_histones) == length(msa_profile)) # all five histone families
  # } else if (is.character(use_profiles)) {
  #   msa_profile <- readRDS("./data/msa/histone_msa.rds")
  #   if (all(use_profiles %in% c("H1", "H2A", "H2B", "H3", "H4"))) {
  #     # user specified which histone families to use, e.g., through names(unaligned_histones)
  #     msa_profile <- hptmUsageData(paste0(use_profiles, ".fasta"))
  #     ref_seqs <- ref_seqs[use_profiles] # TODO?
  #   } else {
  #     # user directly supplied files
  #     stopifnot(
  #       all(file.exists(use_profiles)) &&
  #         length(use_profiles) == length(unaligned_histones)
  #     )
  #     msa_profile <- use_profiles
  #     ref_seqs <- ref_seqs[names(unaligned_histones)] # TODO?
  #   }
  # } else {
  #   msa_profile <- rep.int(FALSE, length(unaligned_histones))
  # }

  # # prepare params
  # extra_params <- list(
  #   # general parameters
  #   amino = TRUE,
  #   quiet = quiet,
  #   thread = num_cores,
  #   threadit = 0, # --threadit 0 for consistency
  #   # return
  #   treeout = FALSE,
  #   reorder = TRUE,
  #   # performance parameters, use E-INS-i https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html#GLE
  #   genafpair = TRUE,
  #   maxiterate = 1000
  # )
  # # add in user-supplied params and turn into a usable CLI string
  # extra_params <- utils::modifyList(extra_params, list(...))
  # extra_params <- .params_to_mafft_str(extra_params)

  # # prepare .fasta files for MAFFT to read, cannot use stdin?
  # in_fasta <- tempfile(pattern = names(unaligned_histones), fileext = ".fasta")
  # out_fasta <- tempfile(pattern = names(unaligned_histones), fileext = "_msa.fasta")
  # invisible(mapply(
  #   Biostrings::writeXStringSet,
  #   unaligned_histones,
  #   in_fasta
  # ))
  # # align
  # mafft_error_code <- mapply(
  #   function(x, y, z) {
  #     system2(
  #       mafft_binary,
  #       paste0(
  #         extra_params,
  #         if (is.character(y)) paste(" --add", x, y, ">", z) else paste(" ", x, ">", z)
  #       )
  #     )
  #   },
  #   in_fasta,
  #   msa_profile,
  #   out_fasta
  # )
  # stopifnot(all(mafft_error_code == 0))
  # # read in the resulting alignments
  # aligned_histones <- lapply(out_fasta, Biostrings::readAAMultipleAlignment)
  # names(aligned_histones) <- names(unaligned_histones)
  # # clean up
  # invisible(file.remove(in_fasta, out_fasta))

  # remove sequences that were not user supplied or the ref seq in profile mode
  if (is.character(msa_profile)) {
    if (return_only_original) {
      aligned_histones <- mapply(
        \(x, y, z) {
          # mask sequences that were not user supplied or the ref seq
          # match ensures that no duplicate sequences are kept (possibility in profile mode)
          Biostrings::rowmask(x, "replace", TRUE) <- match(
            c(names(y), z),
            names(Biostrings::unmasked(x))
          )
          # mask positions that are now all gaps
          x <- Biostrings::maskGaps(x, 1, 1)
          # export the aligned refseq for later mapping
          ref_seq_aligned <- as(x, "AAStringSet")[z]
          # additionally mask the ref seq if not part of original sequences
          if (!z %in% names(y)) {
            Biostrings::rowmask(x, "union", FALSE) <- match(
              c(z),
              names(Biostrings::unmasked(x))
            )
          }
          list(msa = as(x, "AAStringSet"), ref = ref_seq_aligned)
        },
        aligned_histones,
        unaligned_histones,
        ref_seqs,
        SIMPLIFY = FALSE
      )
      msa_ref <- lapply(aligned_histones, \(x) x$ref) |>
        do.call(Biostrings::AAStringSetList, args = _)
      aligned_histones <- lapply(aligned_histones, \(x) x$msa) |>
        do.call(Biostrings::AAStringSetList, args = _)
    } # else {
    # aligned_histones <- do.call(Biostrings::AAStringSetList, aligned_histones)
    #   msa_ref <- mapply(
    #     \(x, y) x[y],
    #     aligned_histones,
    #     ref_seqs
    #   ) |>
    #     do.call(Biostrings::AAStringSetList, args = _)
    # }
    # # clean up default names from HistoneDB 2.0
    # if (is.null(nondefault_refseq_names)) {
    #   clean_ref_names <- c(
    #     H1 = "ref_H11_HUMAN",
    #     H2A = "ref_H2A1B_HUMAN",
    #     H2B = "ref_H2B1J_HUMAN",
    #     H3 = "ref_H31_HUMAN",
    #     H4 = "ref_H4_HUMAN"
    #   )
    #   for (i in seq_along(msa_ref)) {
    #     names(msa_ref[[i]]) <- clean_ref_names[[names(unaligned_histones)[[i]]]]
    #   }
    # }
    # # to return: unaligned sequences, aligned sequences, aligned ref sequence
    # result <- list(
    #   unaligned = unaligned_histones,
    #   msa = aligned_histones,
    #   msa_ref = msa_ref
    # )
    # } #else {
    #   # no need to modify the msa fasta if not in profile mode
    #   result <- list(
    #     unaligned = unaligned_histones,
    #     msa = do.call(Biostrings::AAStringSetList, aligned_histones)
    #   )
  }

  # write (cleaned up if profile mode) fasta
  if (return_alignment) {
    # prepare paths
    dir.create(output_path, showWarnings = FALSE)
    output_paths <- sapply(names(unaligned_histones), \(x) file.path(output_path, paste0(x, "_msa.fasta")))
    # write
    if (any(file.exists(output_paths)) && !overwrite) {
      warning(
        "Output files already exist at ",
        output_paths[file.exists(output_paths)],
        ", skipping, consider setting `overwrite = TRUE`"
      )
    } else {
      mapply(\(x, y) Biostrings::writeXStringSet(x, y), result$msa, output_paths)
    }
  }

  result
}

# OLD --------------------------------------------------------------------------
#
#

alignHistones <- function(
  unaligned_histones = NULL,
  return_alignment = TRUE,
  output_path = "./out/msa/",
  overwrite = FALSE,
  use_profiles = TRUE,
  nondefault_refseq_names = NULL,
  return_only_original = TRUE,
  quiet = TRUE,
  num_cores = -1,
  ...
) {
  # check MAFFT installation
  mafft_binary <- .check_mafft()
  # retrieve sequences for all histone families if NULL
  if (is.null(unaligned_histones)) {
    unaligned_histones <- histonesFromUniprot()
  }
  stopifnot(
    is(unaligned_histones, "AAStringSetList"),
    all(vapply(unaligned_histones, is, logical(1), "AAStringSet"))
  )
  family_names <- names(unaligned_histones)
  # default refseqs or custom?
  ref_seqs <- nondefault_refseq_names %||% .get_default_ref_seqs()
  stopifnot(is.character(ref_seqs))

  # align against existing profiles?
  msa_profiles <- .prepare_msa_profiles(use_profiles, family_names)

  # prepare params for mafft
  mafft_params <- list(
    # general parameters
    amino = TRUE,
    quiet = quiet,
    thread = num_cores,
    # set --threadit 0 for consistency
    threadit = 0,
    # return
    treeout = FALSE,
    reorder = TRUE,
    # performance parameters, use E-INS-i https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html#GLE
    genafpair = TRUE,
    maxiterate = 1000
  ) |>
    # add in user-supplied params and turn into a usable CLI string
    utils::modifyList(list(...)) |>
    .params_to_mafft_str()

  # align using mafft
  alignment_results <- purrr::map2(
    unaligned_histones,
    msa_profiles,
    \(x, y) .run_mafft_alignment(x, y, mafft_binary, mafft_params)
  )

  # check if MSA was successful
  statuses <- purrr::map_int(alignment_results, "status")
  if (any(statuses != 0)) {
    failed <- names(statuses[statuses != 0])
    stop("MAFFT failed for: ", paste(failed, collapse = ", "), call. = FALSE)
  }
  aligned_histones_raw <- purrr::map(alignment_results, "result")

  # no need to modify the msa result if not in profile mode
  if (isFALSE(use_profiles)) {
    result <- list(unaligned = unaligned_histones, msa = Biostrings::AAStringSetList(aligned_histones_raw))
  } else {
    # remove sequences that were not user supplied or the ref seq?
    processed <- if (return_only_original) {
      purrr::pmap(list(aligned_histones_raw, unaligned_histones, ref_seqs[family_names]), .filter_alignment)
    } else {
      list(
        msa = Biostrings::AAStringSetList(aligned_histones_raw),
        ref = purrr::map2(aligned_histones_raw, ref_seqs[family_names], \(x, y) as(x, "AAStringSet")[y])
      )
    }

    # clean up default names from HistoneDB 2.0
    if (is.null(nondefault_refseq_names)) {
      clean_names <- c(
        H1 = "ref_H11_HUMAN",
        H2A = "ref_H2A1B_HUMAN",
        H2B = "ref_H2B1J_HUMAN",
        H3 = "ref_H31_HUMAN",
        H4 = "ref_H4_HUMAN"
      )
      processed$ref <- purrr::map2(
        processed$ref,
        clean_names[names(processed$ref)],
        \(x, y) `names<-`(x, y)
      ) |>
        Biostrings::AAStringSetList()
    } else {
      processed$ref <- Biostrings::AAStringSetList(processed$ref)
    }

    result <- list(unaligned = unaligned_histones, msa = processed$msa, msa_ref = processed$ref)
  }

  return(result) # TODO

  if (return_alignment) {
    .write_msa_files(result$msa, trimws(output_path, "right", "/"), overwrite)
  }

  result
}


.check_mafft <- function() {
  mafft_binary <- Sys.which("mafft")
  if (all(mafft_binary == "")) {
    stop(
      "MAFFT was not found, please ensure that MAFFT was properly installed and can be found on PATH.",
      call. = FALSE
    )
  }
  mafft_binary
}

.get_default_ref_seqs <- function() {
  c(
    H1 = "generic_H1|Homo|NP_005316.1 Homo|NP_005316.1|generic_H1 Homo_generic_H1_4885373", # H11_HUMAN
    H2A = "canonical_H2A|Homo|NP_066390.1 Homo|NP_066390.1|canonical_H2A Homo_canonical_H2A_10645195", # H2A1B_HUMAN
    H2B = "canonical_H2B|Homo|NP_066402.2 Homo|NP_066402.2|canonical_H2B Homo_canonical_H2B_20336754", # H2B1J_HUMAN
    H3 = "canonical_H3|Homo|NP_003520.1 Homo|NP_003520.1|canonical_H3 Homo_canonical_H3_4504281", # H31_HUMAN
    H4 = "canonical_H4|Homo|NP_001029249.1 Homo|NP_001029249.1|canonical_H4 Homo_canonical_H4_77539758" # H4_HUMAN
  )
}

.prepare_msa_profiles <- function(use_profiles, family_names) {
  if (isFALSE(use_profiles)) {
    return(rep(NA_character_, length(family_names)))
  }
  if (isTRUE(use_profiles)) {
    return(hptmUsageData(paste0(family_names, ".fasta")))
  }
  if (is.character(use_profiles)) {
    # user specified which histone families to use
    if (all(use_profiles %in% c("H1", "H2A", "H2B", "H3", "H4"))) {
      return(hptmUsageData(paste0(use_profiles, ".fasta")))
    }
    # user directly supplied files
    stopifnot(
      all(file.exists(use_profiles)),
      length(use_profiles) == length(family_names)
    )
    return(use_profiles)
  }
  stop("`use_profiles` must be logical or a character vector of paths/families.", call. = FALSE)
}

.params_to_mafft_str <- function(params) {
  purrr::imap(params, \(value, arg) {
    if (isTRUE(value)) {
      return(paste0("--", arg))
    }
    if (isFALSE(value)) {
      return(NULL)
    }
    if (is.numeric(value) || is.character(value)) {
      return(paste0("--", arg, " ", value))
    }
    # only character/integer/boolean arguments are possible
    warning("Skipping argument `", arg, "` with unsupported type `", typeof(value), "`.", call. = FALSE)
    NULL
  }) |>
    purrr::compact() |>
    paste(collapse = " ")
}

.run_mafft_alignment <- function(unaligned_set, profile_path, mafft_binary, mafft_params) {
  in_fasta <- tempfile(fileext = ".fasta")
  out_fasta <- tempfile(fileext = "_msa.fasta")
  on.exit(file.remove(in_fasta, out_fasta), add = TRUE)

  Biostrings::writeXStringSet(unaligned_set, in_fasta)

  cmd_args <- if (!is.na(profile_path)) {
    glue::glue("{mafft_params} --add {shQuote(in_fasta)} {shQuote(profile_path)} > {shQuote(out_fasta)}")
  } else {
    glue::glue("{mafft_params} {shQuote(in_fasta)} > {shQuote(out_fasta)}")
  }

  status <- system2(mafft_binary, args = cmd_args)
  list(status = status, result = if (status == 0) Biostrings::readAAMultipleAlignment(out_fasta) else NULL)
}

.filter_alignment <- function(aligned_set, original_set, ref_name) {
  seqs_to_keep <- c(names(original_set), ref_name) |>
    unique()
  keep_indices <- match(seqs_to_keep, names(as(aligned_set, "AAStringSet")))

  Biostrings::rowmask(aligned_set, "replace", TRUE) <- keep_indices
  aligned_set <- Biostrings::maskGaps(aligned_set, min.occupancy = 1, min.fraction = 1)
  ref_seq_aligned <- as(aligned_set, "AAStringSet")[ref_name]

  if (!ref_name %in% names(original_set)) {
    ref_idx <- match(ref_name, names(Biostrings::unmasked(aligned_set)))
    Biostrings::rowmask(aligned_set, "union", FALSE) <- ref_idx
  }

  list(msa = as(aligned_set, "AAStringSet"), ref = ref_seq_aligned)
}

.write_msa_files <- function(msa_list, output_path, overwrite) {
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  output_paths <- file.path(output_path, paste0(names(msa_list), "_msa.fasta"))
  names(output_paths) <- names(msa_list)

  to_write_paths <- if (overwrite) output_paths else output_paths[!file.exists(output_paths)]
  to_skip_paths <- setdiff(output_paths, to_write_paths)

  if (length(to_skip_paths) > 0) {
    warning("Output files already exist, skipping: \n", paste(to_skip_paths, collapse = "\n"), call. = FALSE)
  }
  if (length(to_write_paths) > 0) {
    purrr::pwalk(list(msa_list[names(to_write_paths)], to_write_paths), Biostrings::writeXStringSet)
  }
}
