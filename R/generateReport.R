#' Generate Differential Usage Analysis Report
#'
#' @description
#' TODO
#'
#' @details
#' TODO
#'
#' @param dataset A dataset object produced by your package (e.g., from `deconvolute()`).
#' @param output_file The path where the final report should be saved (e.g., "analysis_report.html").
#' @param contrast A string describing the contrast (passed to report parameters).
#' @param ... Additional parameters passed to `quarto::quarto_render`.
#' @returns The dataset, invisible
#' @export
#' @examples
#' \dontrun{
#' # TODO
#'}
generateReport <- function(
  dataset,
  output_dir,
  overwrite = TRUE,
  # TODO = params,
  ...
) {
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
          ds_path = ds_path
          # TODO = params,
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
