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

#' Example hPTM Data
#'
#' A minimally preprocessed `QFeatures` of hPTM data for testing and
#' demonstration purposes. The samples consist of bottom-up histone extracts
#' from naive (condition_A) and primed (condition_B) stem cell cultures as
#' described in the reference below.
#'
#' @format A QFeatures object with 1 assay named `precursorRaw`. This assay
#' contains 5663 features (peptide ions) and 10 samples from two conditions
#' ('condition_A' and 'condition_B').
#' \describe{
#'   \item{assay}{contains the raw peptide ion abundances}
#'   \item{rowData}{contains 21 variables, including (most relevant) "protein", "sequence", and "mods"}
#'   \item{colData}{contains 8 columns of sample metadata, including the "group" (biological factor), "ms_run"/"date_collected"/"prep_batch" (technical factors), and "include"/"outlier" (tags for analysis)}
#' }
#' @source <https://www.ebi.ac.uk/pride/archive/projects/PXD028162>
#' @references Zijlmans, D. W.; Talon, I.; Verhelst, S.; Bendall, A.;
#'   Van Nerum, K.; Javali, A.; Malcolm, A. A.; Van Knippenberg, S. S. F. A.;
#'   Biggins, L.; To, S. K.; Janiszewski, A.; Admiraal, D.; Knops, R.;
#'   Corthout, N.; Balaton, B. P.; Georgolopoulos, G.; Panda, A.; Bhanu, N. V.;
#'   Collier, A. J.; Fabian, C.; Allsop, R. N.; Chappell, J.; Pham, T. X. A.;
#'   Oberhuemer, M.; Ertekin, C.; Vanheer, L.; Athanasouli, P.; Lluis, F.;
#'   Deforce, D.; Jansen, J. H.; Garcia, B. A.; Vermeulen, M.; Rivron, N.;
#'   Dhaenens, M.; Marks, H.; Rugg-Gunn, P. J.; Pasque, V. Integrated
#'   Multi-Omics Reveal Polycomb Repressive Complex 2 Restricts Human
#'   Trophoblast Induction. Nat Cell Biol 2022, 24 (6), 858–871.
#'   <https://doi.org/10.1038/s41556-022-00932-w>.
#' @examples
#' ncbtoy |>
#'   show()
"ncbtoy"

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
#' benchmark |>
#'   show()
"benchmark"
