#' Retrieve Histone Sequences From UniProt
#'
#' @description
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
#' @examples
#' \dontrun{
#' # Retrieve all human histone sequences from UniProt (reviewed entries only)
#' human_histones <- histonesFromUniprot()
#'
#' # Retrieve only H3 and H4 sequences
#' h3_h4_histones <- histonesFromUniprot(histone_families = c("H3", "H4"))
#'
#' # Retrieve mouse histones instead of human (Taxon ID = 10090)
#' mouse_histones <- histonesFromUniprot(query = c("organism_id:10090", "reviewed:true"))
#'
#' # Use accession numbers for sequence names instead of entry names
#' histones_by_accession <- histonesFromUniprot(name_field = "accession")
#'
#' # Retrieve human histones (reviewed only) and mouse histones (including unreviewed)
#' human_mouse_histones <- histonesFromUniprot(
#'   query = c("organism_id:9606 AND reviewed:true", "organism_id:10090"),
#'   collapse = " OR "
#' )
#'}

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

  # the query should be self-contained, with the histone family being added on top, i.e. "(query) AND histone_family"
  query <- paste0("(", paste0(query, collapse = collapse), ")")

  # each histone family is individually requested to UniProt
  family_lookup[histone_families] |>
    purrr::map(\(family_query) {
      res <- UniProt.ws::queryUniProt(
        query = c(query, family_query),
        fields = c("sequence", name_field),
        collapse = " AND "
      )
      # name_field "id" results in "Entry.Name", so use column index instead
      stats::setNames(Biostrings::AAStringSet(res$Sequence), res[[2]])
    }) |>
    Biostrings::AAStringSetList()
}

#' Align Histone Sequences Using MAFFT
#'
#' @description
#' Histone sequences can be aligned to prevent location mismatches between
#' corresponding hPTMs, e.g, H31S57 vs. H3CS56.
#'
#' @details
#' To ensure consistent and accurate alignment regardless of the underlying
#' sequences, this alignment is by default added into a predefined MSA profile
#' using the MAFFT `--add` functionality (1). This predefined MSA profile is by
#' default the curated alignment of the corresponding histone family retrieved
#' from HistoneDB 2.0 (2).
#'
#' It is therefore highly recommended to leave these settings at default, so
#' that hPTM definition will be consistent, even when the sequences to align
#' change.
#'
#' @param unaligned_histones An optional `AAStringSetList` of histone sequences
#'   to align, by family. If `NULL`, a default call to [histonesFromUniprot()]
#'   is made, i.e., sequences from all five histone families (human and
#'   reviewed) are retrieved.
#' @param return_alignment,output_path,overwrite Return the alignment? In order:
#'   * A `logical(1)` indicating whether to save the alignment to `.fasta` files.
#'   * A `character(1)` specifying the directory to save alignment files.
#'   * A `logical(1)` indicating whether to overwrite existing files at `output_path`
#' @param use_profiles A `logical(1)` or `character()` indicating whether to
#'   add sequences to existing histone MSA profiles. If `TRUE`, uses the curated
#'   alignments from HistoneDB 2.0 (2). Can also be a `character()` specifying
#'   paths to custom, pre-aligned `.fasta` files.
#' @param nondefault_refseq_names An optional `character()` of non-default
#'   reference sequence names. If `NULL`, H1.1 and canonical H2A/H2B/H3/H4 are
#'   used. Be aware that changing this can change hPTM location definition!
#' @param return_only_original A `logical(1)` indicating whether to return only
#'   the sequences present in `unaligned_histones` in the MSA. Optional, only
#'   used in profile mode.
#' @param quiet A `logical(1)` indicating whether to suppress MAFFT output.
#' @param num_cores An `integer(1)` specifying the number of CPU cores to use
#'   during alignment. Optional, `-1` uses all cores.
#' @param ... Additional arguments passed to MAFFT, overriding any defaults.
#' @returns A list containing the unaligned input sequences (`unaligned`), the
#'   multiple sequence alignment (`msa`), and optionally the aligned reference
#'   sequences (`msa_ref`) if `use_profiles` is not `FALSE`. Each list element
#'   is a named `AAStringSetList`. An error is raised if MAFFT fails alignment
#'   for any histone family.
#' @references
#' (1) Katoh, K.; Standley, D. M. MAFFT Multiple Sequence Alignment Software
#' Version 7: Improvements in Performance and Usability. Molecular Biology and
#' Evolution 2013, 30 (4), 772–780. <https://doi.org/10.1093/molbev/mst010>.
#'
#' (2) Draizen, E. J.; Shaytan, A. K.; Mariño-Ramírez, L.; Talbert, P. B.;
#' Landsman, D.; Panchenko, A. R. HistoneDB 2.0: A Histone Database with
#' Variants—an Integrated Resource to Explore Histones and Their Variants.
#' Database (Oxford) 2016, 2016, baw014. <https://doi.org/10.1093/database/baw014>.
#' @export
#' @examples
#' \dontrun{
#' # Align default human histone sequences from UniProt (reviewed entries only)
#' aligned_human_histones <- alignHistones()
#'
#' # Align a custom set of histone sequences
#' library(Biostrings)
#' h3_seqs <- AAStringSet(c(H3.1-Ntail = "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPH",
#'                          H3.3-Ntail = "ARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVKKPH"))
#' h4_seqs <- AAStringSet(c(H4-Ntail = "SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRR"))
#' unaligned_seqs <- AAStringSetList(H3 = h3_seqs, H4 = h4_seqs)
#' aligned_custom_seqs <- alignHistones(unaligned_histones = unaligned_seqs)
#'
#' # Align mouse histones instead of human
#' mouse_unaligned <- histonesFromUniprot(query = c("organism_id:10090", "reviewed:true"))
#' aligned_mouse <- alignHistones(unaligned_histones = mouse_unaligned)
#'
#' # Perform a de novo alignment instead of adding to a profile
#' denovo_aligned <- alignHistones(use_profiles = FALSE)
#'
#' # Pass additional arguments to MAFFT, e.g., to change the gap opening penalty
#' alignHistones(op = 3.0)
#'
#' # Use non-default reference sequences, note that this is not recommended!
#' alignHistones(
#'   nondefault_refseq_names = c(
#'     H1 = "H1.0|Bos|NP_001069955.1 Bos|NP_001069955.1|H1.0 Bos_H1.0_115496898",
#'     H2A = "H2A.1|Homo|NP_734466.1 Homo|NP_734466.1|H2A.1 Homo_H2A.1_25092737",
#'     H2B = "H2B.1|Canis|XP_005640164.1 Canis|XP_005640164.1|H2B.1 Canis_H2B.1_545554624",
#'     H3 = "H3.3|Arabidopsis|NP_195713.1 Arabidopsis|NP_195713.1|H3.3 Arabidopsis_H3.3_15236103",
#'     H4 = "canonical_H4|Mus|NP_291074.1 Mus|NP_291074.1|canonical_H4 Mus_canonical_H4_21361209"
#'   )
#' )
#'}
alignHistones <- function(
  unaligned_histones = NULL,
  return_alignment = FALSE,
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
    result <- list(
      unaligned = unaligned_histones,
      msa = Biostrings::AAStringSetList(aligned_histones_raw)
    )
  } else {
    # remove sequences that were not user supplied or the ref seq?
    processed <- if (return_only_original) {
      filtered_msa <- purrr::pmap(
        list(aligned_histones_raw, unaligned_histones, ref_seqs[family_names]),
        .filter_alignment
      )
      list(
        msa = purrr::map(filtered_msa, "msa") |> Biostrings::AAStringSetList(),
        ref = purrr::map(filtered_msa, "ref")
      )
    } else {
      list(
        msa = Biostrings::AAStringSetList(aligned_histones_raw),
        # retrieve aligned ref seq from MSA frame
        ref = purrr::map2(
          aligned_histones_raw,
          ref_seqs[family_names],
          \(x, y) as(x, "AAStringSet")[y]
        )
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
      )
    }

    result <- list(
      unaligned = unaligned_histones,
      msa = processed$msa,
      msa_ref = Biostrings::AAStringSetList(processed$ref)
    )
  }

  if (return_alignment) {
    .write_msa_files(result$msa, trimws(output_path, "right", "/"), overwrite)
  }

  result
}

.check_mafft <- function() {
  mafft_binary <- Sys.which("mafft")
  if (!nzchar(mafft_binary)) {
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
    stopifnot(
      all(file.exists(use_profiles)),
      length(use_profiles) == length(family_names)
    )
    return(use_profiles)
  }
  stop("`use_profiles` must be logical or a character vector of paths", call. = FALSE)
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
    stringr::str_glue("{mafft_params} --add {shQuote(in_fasta)} {shQuote(profile_path)} > {shQuote(out_fasta)}")
  } else {
    stringr::str_glue("{mafft_params} {shQuote(in_fasta)} > {shQuote(out_fasta)}")
  }

  status <- system2(mafft_binary, args = cmd_args)
  list(
    status = status,
    result = if (status == 0) Biostrings::readAAMultipleAlignment(out_fasta) else NULL
  )
}

.filter_alignment <- function(aligned_set, original_set, ref_name) {
  seqs_to_keep <- c(names(original_set), ref_name) |>
    unique()
  keep_indices <- match(seqs_to_keep, names(Biostrings::unmasked(aligned_set)))

  # mask sequences that were not user supplied or the ref seq
  Biostrings::rowmask(aligned_set, "replace", TRUE) <- keep_indices
  # mask positions that are now all gaps across the remaining sequences
  aligned_set <- Biostrings::maskGaps(aligned_set, min.fraction = 1, min.block.width = 1)
  # export the aligned refseq for later mapping functions
  ref_seq_aligned <- as(aligned_set, "AAStringSet")[ref_name]
  # now additionally mask the ref seq if not part of original sequences
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
    warning(
      "Output files already exist, skipping: \n",
      paste(to_skip_paths, collapse = "\n"),
      call. = FALSE
    )
  }
  if (length(to_write_paths) > 0) {
    purrr::walk2(
      msa_list[names(to_write_paths)],
      to_write_paths,
      Biostrings::writeXStringSet
    )
  }
}
