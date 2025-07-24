# Default histones are retrieved (reuse the same request as during testing)
# The alignment is manually verified and is considered a ground truth
vcr::use_cassette("histones_default", {
  aligned_histones <- alignHistones()
  # manually check output
  purrr::walk(aligned_histones$unaligned, print)
  purrr::walk(aligned_histones$msa, print)
  purrr::walk(aligned_histones$msa_ref, print)
  # looking good
  saveRDS(aligned_histones, file = testthat::test_path("fixtures", "aligned_histones.rds"))
})
