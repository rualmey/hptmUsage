#' Get Path to Data File(s)
#'
#' Convenience function for accessing data files (located in source at
#' `inst/extdata`) included with the hptmUsage package.
#'
#' @param file Package data file name. If `NULL` (default), all data files are
#'   listed.
#' @export
#' @examples
#' hptmUsageData()
#' hptmUsageData("H3.fasta")
hptmUsageData <- function(file = NULL) {
  if (is.null(file)) {
    dir(fs::path_package("extdata", package = "hptmUsage"))
  } else {
    fs::path_package("extdata", file, package = "hptmUsage")
  }
}

#' hPTM Benchmark Data
#'
#' A preprocessed `QFeatures` of hPTM data, functioning as a benchmark for
#' downstream analysis and statistics. The samples consist of commercial histone
#' extracts treated with a histone deacetylase (HDAC1) in a time-lapse design as
#' described in the reference below.
#'
#' @format TODO A QFeatures with 7,240 rows and 60 columns:
#' \describe{
#'   \item{assay}{contains the raw peptide intensities}
#'   \item{rowData}{contains a variable "Proteins" with the protein accession and an variable ecoli to indicate if the protein is a spikin}
#'   \item{colData}{contains a factor condition indicating the spike-in condition}
#'   ...
#' }
#' @source <https://www.ebi.ac.uk/pride/archive/projects/PXD009910>
#' @references Clerck, L. D.; Willems, S.; Daled, S.; Puyvelde, B. V.;
#'   Verhelst, S.; Corveleyn, L.; Deforce, D.; Dhaenens, M. An Experimental
#'   Design to Extract More Information from MS-Based Histone Studies. Molecular
#'   Omics 2021, 17 (6), 929–938. <https://doi.org/10.1039/D1MO00201E>.
#' @examples
#' hPTM_benchmark |>
#'   show()
"hPTM_benchmark"
