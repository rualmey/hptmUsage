#' Generate Differential Usage Analysis Report
#'
#' @description
#' TODO
#'
#' @details
#' TODO
#'
#' @param dataset The hPTM dataset generated using [hptmUsage::readProgenesis()].
#' @param output_dir The path where the final report should be saved (e.g.,
#'   "./output").
#' @param overwrite If a report with the same name already exists at
#'   `output_dir`, should it be overwritten (`logical(1)`?
#' @param histone_params A `list()` of arguments passed to
#'   [histonesFromUniprot()].
#' @param msa_params A `list()` of arguments passed to [alignHistones()] (except
#'   for "unaligned_histones").
#' @param mod_params A `list()` of arguments passed to [processMods()] (except
#'   for "object", "msa", or "i").
#' @param contaminant_params A `list()` of arguments passed to
#'   [tagContaminants()] (except for "object" or "i").
#' @param usage_params A `list()` of arguments passed to [calculateUsage()]
#'   (except for "object", "i", "name", or "target").
#' @param variant_usage_level A `character(1)` specifying the `usage_level`
#'   (see [calculateUsage()]) for variant usage calculation, i.e., one of
#'   "histone" or "histone_family".
#' @param ... Additional parameters passed to [quarto::quarto_render()].
#' @returns Invisibly returns the final dataset, including aggregated hPTM assays, for example.
#' @export
#' @examples
#' \dontrun{
#' # TODO
#' # Usage level
#' ...
#'}
generateReport <- function(
  dataset,
  output_dir,
  overwrite = FALSE,
  # TODO
  # preprocessing
  histone_params = list(),
  msa_params = list(),
  mod_params = list(),
  contaminant_params = list(),
  usage_params = list(),
  variant_usage_level = c("histone", "histone_family"),
  ...
) {
  # TODO check params
  .paramcheck(histone_params, hptmUsage::histonesFromUniprot)
  .paramcheck(msa_params, hptmUsage::alignHistones, exclusions = "unaligned_histones")
  .paramcheck(mod_params, hptmUsage::processMods, exclusions = c("object", "msa", "i"))
  .paramcheck(contaminant_params, hptmUsage::tagContaminants, exclusions = c("object", "i"))
  .paramcheck(usage_params, hptmUsage::calculateUsage, exclusions = c("object", "i", "name", "target"))
  variant_usage_level <- match.arg(variant_usage_level)

  # quarto needs to be available
  quarto::quarto_available(min = "1.8", max = NULL, error = TRUE)

  # quarto will render in a temporary project directory and the report exported
  resource_path <- fs::path_package("quarto", package = "hptmUsage")
  stopifnot("Report templates not found in package. Try reinstalling the package." = resource_path != "")
  temp_dir <- withr::local_tempdir(pattern = "hptm_report_")
  fs::dir_copy(path = resource_path, new_path = temp_dir, overwrite = TRUE)

  # serialize dataset first instead of passing to quarto::quarto_render() as this should be safer
  ds_path <- fs::path(temp_dir, "input_data", ext = "rds")
  saveRDS(dataset, ds_path)

  # render
  template_qmd <- fs::path(temp_dir, "report_template", ext = "qmd")
  report_basename <- paste0(format(Sys.time(), format = "%y%m%d_%H%M"), "_hptmUsage.html")
  stopifnot("Template not found in package resources. Try reinstalling the package." = fs::file_exists(template_qmd))
  message("Rendering report...")
  tryCatch(
    {
      quarto::quarto_render(
        input = template_qmd,
        output_format = "html",
        output_file = report_basename,
        execute_params = list(
          # passthrough of function parameters, see template YAML header
          ds_path = ds_path,
          histone_params = histone_params,
          msa_params = msa_params,
          mod_params = mod_params,
          contaminant_params = contaminant_params,
          usage_params = usage_params,
          variant_usage_level = variant_usage_level
          # TODO,
        ),
        ...
      )
    },
    error = function(e) stop("Quarto rendering failed: ", e$message)
  )

  # move results to output directory
  generated_html <- fs::path(temp_dir, report_basename)
  stopifnot("Report generated but output file could not be located." = fs::file_exists(generated_html))
  html_output <- fs::path(output_dir, report_basename)
  if (fs::file_exists(html_output) && !overwrite) {
    stop("Output file already exists: ", html_output, ". Set overwrite = TRUE if so desired.")
  }
  fs::dir_create(output_dir)
  fs::file_move(generated_html, output_dir)
  message("Report saved to: ", normalizePath(output_dir))

  # the processed QFeatures object was serialized at the end of the quarto_render call
  dataset <- readRDS(fs::path(temp_dir, "ds_processed", ext = "rds"))
  return(invisible(dataset))
}

.paramcheck <- function(param, fun, exclusions = NULL) {
  wrong_names <- names(param)[!names(param) %in% methods::formalArgs(fun)]
  if (length(wrong_names) > 0) {
    stop(
      paste0(
        "Invalid argument",
        if (length(wrong_names) > 1) "s" else "",
        " in '",
        deparse(substitute(param)),
        "': ",
        paste(wrong_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  excluded_names <- names(param)[names(param) %in% exclusions]
  if (length(excluded_names) > 0) {
    stop(
      paste0(
        "Argument",
        if (length(excluded_names) > 1) "s" else "",
        " not allowed: ",
        paste(excluded_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}
