# QFeatures object of correctly read minimal example
readProgenesis(
  testthat::test_path("fixtures", "all_ion_export.csv"),
  generate_metadata = testthat::test_path("fixtures", "all_ion_export_metadata.csv")
) |>
  saveRDS(file = testthat::test_path("fixtures", "all_ion_export_read_correct.rds"))
