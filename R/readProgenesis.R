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
#' @param return_qfeatures Whether to return a `QFeatures` (`TRUE`, default)
#'   object or a `tibble` (`FALSE`).
#' @param simplify_column_names Remove pre-/suffixes from sample names (default
#'   `TRUE`)?
#' @param generate_metadata Generate a metadata `.csv` file that can be manually
#'   edited (default `TRUE`, see also `overwrite_metadata`)? Can also be a path
#'   to the location where the metadata file will be written. Updated metadata
#'   can later be added to a `QFeatures` object using [replace_metadata()].
#' @param overwrite_metadata Overwrite the metadata file if it already exists
#'   (default `FALSE`)?
#' @returns A `QFeatures` or `tibble` (see `return_qfeatures`) containing the
#'   Progenesis QIP dataset.
#' @export
#' @seealso The [QFeatures::QFeatures()] class to read about how to manipulate
#'   the resulting `QFeatures` object (if `return_qfeatures` was `TRUE`).
#' @examples
#' # Simplest case
readProgenesis <- function(
  file,
  quant = "Raw abundance",
  return_qfeatures = TRUE,
  simplify_column_names = TRUE,
  generate_metadata = TRUE,
  overwrite_metadata = FALSE
) {
  # check arguments
  stopifnot(file.exists(file))
  stopifnot(quant %in% c("Raw abundance", "Normalized abundance", "Intensity"))
  # check header
  header <- read.csv(file, header = FALSE, nrows = 3)
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

  return(df)

  # check for empty sequences
  sequence_col <- which(df[3, ] == "Sequence")
  sequence_df <- df[-(1:2), ][is.na(df[-(1:2), sequence_col]), c(1, sequence_col)]
  if (nrow(sequence_df) > 2) {
    print("Following features have no assigned sequence, please verify in Progenesis QIP. These will be removed...")
    print(sequence_df)
  }
  # check for notes
  notes_col <- which(df[3, ] == "Notes")
  if (length(notes_col)) {
    notes_df <- df[-(1:2), ][!is.na(df[-(1:2), notes_col]), c(1, notes_col)]
    if (nrow(notes_df) > 1) {
      print("Following features had a note, please make sure these do not need to be addressed:")
      print(notes_df)
    }
  }

  # to slice quant columns
  quant_start_index <- which(df[1, ] == quant) # can only occur once so this is fine
  first_sample_name <- df[[3, quant_start_index]]
  repeat_first_sample_idx <- which(df[3, ] == first_sample_name)
  repeat_first_sample_idx <- repeat_first_sample_idx[repeat_first_sample_idx > quant_start_index]
  first_property_after_samples <- which(
    df[3, ] == "Notes" |
      df[3, ] == "Score" |
      df[3, ] == "Mass error (u)" |
      df[3, ] == "Mass error (ppm)" |
      df[3, ] == "Protein" # This column must be present or it will have errored out during argument checking
  )
  quant_end_index <- min(c(repeat_first_sample_idx, first_property_after_samples))
  stopifnot(quant_start_index < quant_end_index)
  quant_idx <- quant_start_index:(quant_end_index - 1)

  # store group assignment from Progenesis QIP for metadata
  sample_names <- df[3, quant_idx]
  if (generate_metadata) {
    metadata <- tibble(
      original_name = as.character(sample_names),
      group = as.character(df[2, quant_idx])
    ) |>
      tidyr::fill(group) |>
      mutate(group = as.factor(group))
  }

  # drop other quant columns
  other_quant_idx <- which(df[3, ] %in% sample_names)
  other_quant_idx <- other_quant_idx[!other_quant_idx %in% quant_idx]
  df <- df[, -other_quant_idx]
  new_quant_idx <- which(df[3, ] %in% sample_names)
  # simplify sample names if requested
  if (simplify_column_names) {
    prefix <- Biobase::lcPrefix(sample_names)
    if (prefix != "") {
      sample_names_clean <- stringr::str_remove(sample_names, stringr::fixed(prefix))
    }
    suffix <- Biobase::lcSuffix(sample_names_clean)
    if (suffix != "") {
      sample_names_clean <- stringr::str_remove(sample_names_clean, stringr::fixed(suffix))
    }
    df[3, new_quant_idx] <- as.list(sample_names_clean)
  }

  # clean up the dataframe
  df <- df |>
    # drop first two header rows and set the third as (clean) column names
    # ignore warning as we fix it below
    (function(x) {
      suppressWarnings(janitor::row_to_names(x, 3), classes = "janitor_warn_row_to_names_not_unique")
    })() |>
    # remove duplicate names
    janitor::clean_names() |>
    dplyr::rename(
      # no need to use `any_of(lookup)` as both columns have to be present, see argument checking
      feature_number = number,
      mods = variable_modifications_position_description
    )
  cleaned_names <- names(df)[new_quant_idx]
  df <- df |>
    dplyr::mutate(
      # ID transfer from different charge states in ProgenesisQIP gives NA score
      score = na_if(score, "---"),
      # fix quant column types
      across(all_of(cleaned_names), as.double),
      # fix integer column types
      across(any_of(c("feature_number", "charge")), as.integer),
      # fix double column types
      across(
        any_of(c(
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
        -any_of(c(
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
    ) |>
    drop_na(sequence) # drop features with no assigned sequence

  # generate metadata if requested TODO allow path
  metadata <- dplyr::mutate(metadata, quantCols = names(df)[new_quant_idx], .before = 1) # original sample names
  if (generate_metadata) {
    metadata <- dplyr::mutate(
      metadata,
      # "include" = include the sample for statistics or only visualization, so QC samples should be FALSE
      include = replace(rep.int(TRUE, nrow(metadata)), stringr::str_detect(quantCols, "qc"), FALSE),
      outlier = rep.int(FALSE, nrow(metadata))
    )
    outfile <- paste0(tools::file_path_sans_ext(file), "_metadata.csv")
    if (file.exists(outfile) && !overwrite_metadata) {
      warning(paste0("File already exists at \"", outfile, "\", see argument `overwrite_metadata`."))
    } else {
      if (file.exists(outfile)) {
        print("Overwriting metadata...")
      }
      readr::write_csv(metadata, file = outfile)
      print(paste("Metadata written at", outfile))
    }
  }

  # transfer to a QFeatures object
  if (return_qfeatures) {
    df <- QFeatures::readQFeatures(
      df,
      colData = metadata,
      name = paste0("precursor", stringr::str_split_i(quant, " ", 1)),
      verbose = FALSE
    )
  }

  df
}
