#' Generate Differential Usage Analysis Report
#'
#' @description
#' Dynamically generates a differential usage report using Quarto. All that
#' needs to be supplied is the histone dataset (see
#' [hptmUsage::readProgenesis()] and [hptmUsage::replaceColData()]), an output
#' directory, and optional settings that modify the analysis. Note that the
#' final QFeatures dataset is invisibly returned so that subsequent custom
#' analysis/plotting is still possible.
#'
#' @details
#' Under the hood, this function creates a temporary Quarto project directory
#' where a template `.qmd` is then rendered. This template file is dynamically
#' rendered using knitr YAML params that are defined by this function's
#' parameters. The document is rendered as a self-contained HTML, which includes
#' a download button to static assets generated during rendering, i.e., figures
#' and `.csv` assays of relevant usage values.
#'
#' @param dataset The hPTM dataset generated using [hptmUsage::readProgenesis()].
#' @param output_dir The path where the final report should be saved (e.g.,
#'   "./output/").
#' @param overwrite If a report with the same name already exists at
#'   `output_dir`, should it be overwritten (`logical(1)`?
#' @param design_formula A named `character()` of design formula(e), e.g., the
#'   default `c(factor = "~ 0 + group")`. The experimental design matrix will
#'   only be shown for models that contain "factor" in the name, i.e., the
#'   default formula name.
#' @param random_effects Optional, a `character()` of categorical variables in
#'   the colData that specifies blocking factors for which a random intercept
#'   should be fit, e.g., sample preparation batch or study individual.
#' @param contrasts A named `list()` of named `character()`s containing relevant
#'   contrasts for differential usage testing; for example,
#'   `list(factor = c("A vs B" = "groupB - groupA"))` for design_formula
#'   `c(factor = "~ 0 + group")` will return a positive LogFC if usage is higher
#'   in "conditionB" according to the "factor" model. If no contrasts are
#'   provided, differential usage testing will be skipped but the design matrix
#'   and QC plots will still be generated for quick finetuning of the analysis.
#' @param reference_levels Optional, a named collection of factors (as present
#'   in the colData) and their reference level, e.g., `c(group = "wild_type")`.
#'   This allows setting the reference group (i.e., intercept) in a
#'   "mean-reference" model.
#' @param generate_lineplots For which hPTMs/peptidoforms/variants should
#'   lineplots be generated? One of "none", "significant", or "all".
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
#' # By default, a "variant-agnostic" usage analysis of human histone samples is performed
#' # The "ncbtoy" dataset was created by supplying a ".csv" to readProgenesis() and then adding metadata using replaceColData()
#' generateReport(ncbtoy, "./out/")
#'
#' # Inspecting this HTML report, we can see that hPTMs are not properly parsed because they follow a non-standard format
#' # We can use the "*_params" parameters to supply this sort of information, see the relevant parameter descriptions for more explanation
#' generateReport(ncbtoy, "./out/", mod_params = list(mod_format = "progenesis_sw"))
#'
#' # Note that this analysis ends at the design matrix of a simple means model that captures the experimental groups
#' # This matrix can help us set up the relevant contrasts
#' generateReport(
#'   dataset = ncbtoy,
#'   output_dir = "./out/",
#'   mod_params = list(mod_format = "progenesis_sw"),
#'   contrasts = list(factor = c("A vs B" = "groupcondition_B - groupcondition_A"))
#' )
#'
#' # We can of course supply our own design formula(e)
#' generateReport(
#'   dataset = ncbtoy,
#'   output_dir = "./out/",
#'   mod_params = list(mod_format = "progenesis_sw"),
#'   design_formula = c(factor = "~ 0 + group", factor_batch = "~ 0 + group + prep_batch"),
#'   contrasts = list(
#'     factor = c("A vs B" = "groupcondition_B - groupcondition_A"),
#'     factor_batch = c(
#'       "A vs B" = "groupcondition_B - groupcondition_A",
#'       "batch_B vs batch_A" = "prep_batchB"
#'     )
#'   )
#' )
#' # or instead supply such variables as a random effect if this makes more sense
#' generateReport(
#'   dataset = ncbtoy,
#'   output_dir = "./out/",
#'   mod_params = list(mod_format = "progenesis_sw"),
#'   random_effects = "prep_batch",
#'   contrasts = list(
#'     factor = c("A vs B" = "groupcondition_B - groupcondition_A")
#'   )
#' )
#'
#' # Lineplots can be generated if so desired, for example of the significant hPTMs/peptidoforms/variants
#' # These will be available through a download button at the bottom of the HTML report
#' # Do note that this can take quite some time/resources to generate all plots
#' generateReport(
#'   dataset = ncbtoy,
#'   output_dir = "./out/",
#'   mod_params = list(mod_format = "progenesis_sw"),
#'   contrasts = list(factor = c("A vs B" = "groupcondition_B - groupcondition_A")),
#'   generate_lineplots = "significant"
#' )
#'
#' # Finally, we can easily define the level at which usage is defined
#' # By default, histone precursor usage and therefore hPTM/histone variant usage are defined against the entire chromatosome (i.e., all histone proteins)
#' # Instead, we can choose a "variant-corrected" usage workflow for histone precursors/hPTMs or calculate usage of histone variants against their corresponding histone families
#' generateReport(
#'   dataset = ncbtoy,
#'   output_dir = "./out/",
#'   mod_params = list(mod_format = "progenesis_sw"),
#'   contrasts = list(factor = c("A vs B" = "groupcondition_B - groupcondition_A")),
#'   usage_params = list(usage_level = "histone_group"),
#'   variant_usage_level = "histone_family"
#' )
#'}
generateReport <- function(
  dataset,
  output_dir,
  overwrite = FALSE,
  design_formula = c(factor = "~ 0 + group"),
  random_effects = NULL,
  contrasts = NULL,
  reference_levels = NULL,
  generate_lineplots = c("none", "significant", "all"),
  histone_params = list(),
  msa_params = list(),
  mod_params = list(),
  contaminant_params = list(),
  usage_params = list(),
  variant_usage_level = c("histone", "histone_family"),
  ...
) {
  .paramcheck(histone_params, hptmUsage::histonesFromUniprot)
  .paramcheck(msa_params, hptmUsage::alignHistones, exclusions = "unaligned_histones")
  .paramcheck(mod_params, hptmUsage::processMods, exclusions = c("object", "msa", "i"))
  .paramcheck(contaminant_params, hptmUsage::tagContaminants, exclusions = c("object", "i"))
  .paramcheck(usage_params, hptmUsage::calculateUsage, exclusions = c("object", "i", "name", "target"))
  variant_usage_level <- match.arg(variant_usage_level)
  generate_lineplots <- match.arg(variant_usage_level)

  # quarto needs to be available
  quarto::quarto_available(min = "1.8", max = NULL, error = TRUE)

  # quarto will render in a temporary project directory from which the report will be exported
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
          variant_usage_level = variant_usage_level,
          reference_levels = reference_levels,
          reference_levels_names = names(reference_levels),
          design_formula = design_formula,
          design_formula_names = names(design_formula),
          random_effects = random_effects,
          contrasts = contrasts,
          contrasts_names = purrr::map(contrasts, names),
          generate_lineplots = generate_lineplots
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
