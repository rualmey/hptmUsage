#' Read a Progenesis QIP Peptide Ion Data Export
#'
#' A Progenesis QIP peptide ion data `.csv` file is cleaned up to make it
#' compatible with the `QFeatures` framework. This involves subsetting the
#' quantitative data to one type, e.g., raw abundances, sanitizing sample names,
#' and generating initial metadata from the Progenesis experimental design.
#'
#' @param file The path (`character(1)`) to the Progenesis QIP peptide ion data
#'   `.csv` file. This file should contain all identified features, including
#'   co-extracts, and at least the properties of feature number, charge,
#'   protein, sequence, variable modifications, and a quantitative value (see
#'   `quant`).
#' @param quant The quantitative value to use: `"Raw abundance"` (default),
#'   `"Normalized abundance"`, or `"Intensity"`.
#' @param generate_metadata Generate a metadata `.csv` file that can be manually
#'   edited (default `FALSE`, see also `overwrite_metadata`)? Can also be a path
#'   to the location where the metadata file will be written (e.g.,
#'   `"./data/metadata.csv"`). Updated metadata can later be added to a
#'   `QFeatures` object using [replace_metadata()].
#' @param overwrite_metadata Overwrite the metadata file if it already exists
#'   (default `FALSE`)?
#' @param simplify_column_names Remove pre-/suffixes from sample names (default
#'   `TRUE`)?
#' @returns A `QFeatures` containing the Progenesis QIP dataset.
#' @export
#' @seealso The [QFeatures::QFeatures()] class to read about how to manipulate
#'   the resulting `QFeatures` object (if `return_qfeatures` was `TRUE`).
#' @examples
#' \dontshow{
#' .old_wd <- setwd(tempdir())
#' file.copy(hptmUsageData("all_ion_export.csv"), ".")
#' }
#' # Simplest case
#' readProgenesis("./all_ion_export.csv")
#'
#' # Generate a metadata `.csv` file that can be manually edited
#' readProgenesis("./all_ion_export.csv", generate_metadata = TRUE)
#' # A path can also be supplied, instead of having the file generated next to
#' # the `.csv` source.
#' readProgenesis("./all_ion_export.csv", generate_metadata = "./metadata.csv")
#' # Note, if the metadata file already exists, nothing will be written unless
#' # `overwrite_metadata` is set to `TRUE`
#'
#' # Use feature intensities instead of raw abundances
#' readProgenesis("./all_ion_export.csv", quant = "Intensity")
#'
#' # Do not simplify (remove pre- and suffixes) sample names
#' readProgenesis("./all_ion_export.csv", simplify_column_names = FALSE)
#' \dontshow{
#' setwd(.old_wd)
#' }
readProgenesis <- function(
  file,
  quant = "Raw abundance",
  generate_metadata = FALSE,
  overwrite_metadata = FALSE,
  simplify_column_names = TRUE
) {
  # check arguments
  stopifnot(file.exists(file))
  stopifnot(quant %in% c("Raw abundance", "Normalized abundance", "Intensity") && length(quant) == 1)
  stopifnot(
    isFALSE(generate_metadata) ||
      isTRUE(generate_metadata) ||
      (is.character(generate_metadata) && length(generate_metadata) == 1)
  )
  # check header
  header <- utils::read.csv(file, header = FALSE, nrows = 3)
  stopifnot(quant %in% header[1, ])
  stopifnot(all(
    c("#", "Charge", "Protein", "Sequence", "Variable modifications ([position] description)") %in% header[3, ]
  ))

  # read progenesis csv, need to ignore header as readr cannot handle multi-level
  df <- file |>
    readr::read_csv(
      col_names = FALSE,
      col_types = readr::cols(.default = readr::col_character()),
      progress = FALSE
    )

  # subselect quant data
  quant_idx <- get_quant_cols(df, quant)
  sample_names <- as.character(df[3, quant_idx])
  other_quant_idx <- setdiff(
    which(df[3, ] %in% sample_names),
    quant_idx
  )
  df <- df |> dplyr::select(!tidyselect::all_of(other_quant_idx))
  quant_idx <- which(df[3, ] %in% sample_names)

  # store group assignment from Progenesis QIP for later metadata, will be dropped soon
  metadata <- tibble::tibble(
    original_name = sample_names,
    group = as.character(df[2, quant_idx])
  ) |>
    tidyr::fill("group") |>
    mutate(group = forcats::as_factor(.data[["group"]]))

  # simplify sample names if requested and clean up header
  if (simplify_column_names) {
    sample_names <- remove_pre_suffix(sample_names)
    df[3, quant_idx] <- as.list(sample_names)
  }
  # first two rows contain no useful info, ignore warning as we fix it below
  df <- suppressWarnings(janitor::row_to_names(df, 3), classes = "janitor_warn_row_to_names_not_unique") |>
    # remove duplicate names
    janitor::clean_names() |>
    dplyr::rename(
      # no need to use `any_of(lookup)` as both columns have to be present, see argument checking
      feature_number = "number",
      mods = "variable_modifications_position_description"
    )
  cleaned_names <- names(df)[quant_idx]

  # clean up the dataframe and fix dtypes
  df <- df |>
    # ID transfer from different charge states in ProgenesisQIP gives NA score
    mutate(score = dplyr::na_if(.data[["score"]], "---")) |>
    fix_progenesis_dtypes(cleaned_names)

  # check for notes, if present
  if ("notes" %in% colnames(df)) {
    notes <- df |>
      tidyr::drop_na("notes") |>
      dplyr::select(tidyselect::all_of(c("feature_number", "notes")))
    if (nrow(notes) >= 1) {
      message(
        "Some features had a note:\n",
        paste0("  Feature ", notes[["feature_number"]], ": ", notes[["notes"]], "\n")
      )
    }
  }

  # check for empty sequences
  missing_seq <- df[["sequence"]] |>
    vctrs::vec_detect_missing() |>
    which()
  missing_seq <- df |>
    dplyr::slice(missing_seq) |>
    _[["feature_number"]]
  if (length(missing_seq) >= 1) {
    warning(
      "Some features have no assigned sequence, please verify. These will be dropped: ",
      paste(missing_seq, collapse = ", ")
    )
  }
  # drop the features with no assigned sequence
  df <- df |>
    tidyr::drop_na("sequence")

  # generate metadata if requested
  metadata <- metadata |>
    mutate(quantCols = cleaned_names, .before = 1) |>
    mutate(
      # "include" = include the sample for statistics or only visualization, so QC samples should be FALSE
      include = replace(rep.int(TRUE, nrow(metadata)), .data[["group"]] == "QC", FALSE),
      outlier = rep.int(FALSE, nrow(metadata))
    )
  if (!isFALSE(generate_metadata)) {
    if (is.character(generate_metadata)) {
      outfile <- generate_metadata
    } else {
      outfile <- paste0(tools::file_path_sans_ext(file), "_metadata.csv")
    }
    if (file.exists(outfile) && !overwrite_metadata) {
      warning(paste0("File already exists at \"", outfile, "\", see argument `overwrite_metadata`."))
    } else {
      message(
        if (file.exists(outfile)) "Overwriting metadata... " else "",
        "Metadata written at ",
        outfile,
        "\n"
      )
      readr::write_csv(metadata, file = outfile)
    }
  }

  # return as a QFeatures object
  QFeatures::readQFeatures(
    df,
    colData = metadata,
    name = paste0("precursor", stringr::str_split_i(quant, " ", 1)),
    verbose = FALSE
  )
}

# helper to retrieve the column indexes belonging to the quant type of choice
get_quant_cols <- function(df, quant) {
  quant_start_index <- which(df[1, ] == quant)
  first_sample_name <- df[[3, quant_start_index]]
  # the quant columns stop either at the next type of quant (e.g. Intensity), so when the sample names repeat
  repeat_first_sample_idx <- which(df[3, ] == first_sample_name)
  repeat_first_sample_idx <- repeat_first_sample_idx[repeat_first_sample_idx > quant_start_index]
  # or when one of the following columns is found (note that "Protein" must always be present)
  first_property_after_samples <- which(
    df[3, ] %in% c("Notes", "Score", "Mass error (u)", "Mass error (ppm)", "Protein")
  )
  quant_end_index <- min(c(repeat_first_sample_idx, first_property_after_samples))
  if (quant_end_index <= quant_start_index) {
    stop("Error in selecting quant columns")
  }
  quant_start_index:(quant_end_index - 1)
}

# helper to remove pre- and/or suffixes from a character vector
remove_pre_suffix <- function(sample_names, prefix = TRUE, suffix = TRUE) {
  if (prefix) {
    prefix_len <- nchar(Biobase::lcPrefix(sample_names))
    sample_names <- stringr::str_sub(sample_names, start = prefix_len + 1)
  }
  if (suffix) {
    suffix_len <- nchar(Biobase::lcSuffix(sample_names))
    sample_names <- stringr::str_sub(sample_names, end = -suffix_len - 1)
  }
  sample_names
}


fix_progenesis_dtypes <- function(df, cleaned_names) {
  df |>
    mutate(
      # fix quant column types
      across(tidyselect::all_of(cleaned_names), as.double),
      # fix integer column types
      across(c("feature_number", "charge"), as.integer),
      # fix double column types
      across(
        tidyselect::any_of(c(
          "m_z",
          "retention_time_min",
          "retention_time_window_min",
          "mass",
          "max_fold_change",
          "anova",
          "maximum_cv",
          "score",
          "mass_error_u",
          "mass_error_ppm"
        )),
        as.double
      ),
      # all except below as factor, this due to unknown number of tag columns
      across(
        -tidyselect::any_of(c(
          cleaned_names,
          "feature_number",
          "charge",
          "m_z",
          "retention_time_min",
          "retention_time_window_min",
          "mass",
          "max_fold_change",
          "anova",
          "maximum_cv",
          "score",
          "mass_error_u",
          "mass_error_ppm",
          "notes",
          "protein",
          "sequence",
          "mods",
          "description"
        )),
        as.factor
      )
    )
}
