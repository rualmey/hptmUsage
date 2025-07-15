test_that("correct quant data is subselected", {
  df <- hptmUsageData("all_ion_export.csv") |>
    readr::read_csv(
      col_names = FALSE,
      col_types = readr::cols(.default = readr::col_character()),
      progress = FALSE
    )
  expect_equal(get_quant_cols(df, "Raw abundance"), 22:31)
  expect_equal(get_quant_cols(df, "Normalized abundance"), 12:21)
  expect_equal(get_quant_cols(df, "Intensity"), 32:41)
})

test_that("correct quant data is subselected", {
  expect_equal(remove_pre_suffix(paste("pre", c("a", "b", "c"))), c("a", "b", "c"))
  expect_equal(remove_pre_suffix(paste(c("a", "b", "c"), "suf")), c("a", "b", "c"))
  expect_equal(remove_pre_suffix(paste("pre", c("a", "b", "c"), "suf")), c("a", "b", "c"))
})
