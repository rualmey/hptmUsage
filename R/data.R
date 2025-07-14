#' Get Path to Data File(s)
#'
#' Convenience function for accessing data files (located in source at
#' `inst/extdata`) included with the hptmUsage package.
#'
#' @param file Package data file name. If `NULL` (default), all data files are
#'   listed.
#' @export
#' @examples
#' hptmUsage_data()
#' hptmUsage_data("H3.fasta")
hptmUsage_data <- function(file = NULL) {
  if (is.null(file)) {
    dir(fs::path_package("extdata", package = "hptmUsage"))
  } else {
    fs::path_package("extdata", file, package = "hptmUsage")
  }
}
