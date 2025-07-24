skip_if_no_mafft <- function() {
  if (!nzchar(Sys.which("mafft"))) {
    testthat::skip("MAFFT not found")
  }
}
