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

test_that("names are simplified", {
  expect_equal(.remove_pre_suffix(paste("pre", c("a", "b", "c"))), c("a", "b", "c"))
  expect_equal(.remove_pre_suffix(paste(c("a", "b", "c"), "suf")), c("a", "b", "c"))
  expect_equal(.remove_pre_suffix(paste("pre", c("a", "b", "c"), "suf")), c("a", "b", "c"))
})

test_that("non-existent file throws error", {
  expect_snapshot(
    error = TRUE,
    readProgenesis("foo")
  )
})

test_that("readProgenesis catches invalid quant", {
  expect_snapshot(
    error = TRUE,
    readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"), quant = NA)
  )
  expect_snapshot(
    error = TRUE,
    readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"), quant = -1)
  )
  expect_snapshot(
    error = TRUE,
    readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"), quant = Inf)
  )
  expect_snapshot(
    error = TRUE,
    readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"), quant = NULL)
  )
  expect_snapshot(
    error = TRUE,
    readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"), quant = "Spectral count")
  )
  expect_snapshot(
    error = TRUE,
    readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"), quant = c("Raw abundance", "Intensity"))
  )
})

test_that("readProgenesis catches invalid generate_metadata", {
  expect_snapshot(
    error = TRUE,
    readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"), generate_metadata = NA)
  )
  expect_snapshot(
    error = TRUE,
    readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"), generate_metadata = -1)
  )
  expect_snapshot(
    error = TRUE,
    readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"), generate_metadata = Inf)
  )
  expect_snapshot(
    error = TRUE,
    readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"), generate_metadata = NULL)
  )
  expect_snapshot(
    error = TRUE,
    readProgenesis(
      testthat::test_path("fixtures", "all_ion_export.csv"),
      generate_metadata = c("./path_foo", "./path_bar")
    )
  )
})

test_that("readProgenesis catches invalid headers", {
  expect_snapshot(
    error = TRUE,
    readProgenesis(testthat::test_path("fixtures", "all_ion_export_no_quant.csv"))
  )
  expect_snapshot(
    error = TRUE,
    readProgenesis(testthat::test_path("fixtures", "all_ion_export_no_quant.csv"), quant = "Intensity")
  )
  expect_snapshot(
    error = TRUE,
    readProgenesis(testthat::test_path("fixtures", "all_ion_export_no_feat.csv"))
  )
})

test_that("messages and warnings work as expected", {
  expect_snapshot(
    invisible(readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv")))
  )
})

test_that("reading works as expected", {
  # suppress output for clean test() result
  out <- suppressMessages(suppressWarnings(
    readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"))
  ))
  # compare against manually checked QFeatures, created in "./fixtures/all_ion_export_make_rds.R"
  expect_identical(out, readRDS(test_path("fixtures", "all_ion_export_read_correct.rds")))
})

test_that("changing quant column works", {
  # suppress output for clean test() result
  out <- suppressMessages(suppressWarnings(
    readProgenesis(
      testthat::test_path("fixtures", "all_ion_export.csv"),
      quant = "Intensity"
    )
  ))
  reference <- readRDS(test_path("fixtures", "all_ion_export_read_correct.rds"))
  # everything stays the same
  expect_identical(colData(out), colData(reference))
  expect_identical(rowData(out[[1]]), rowData(reference[[1]]))
  # except the assay values
  # fmt: skip
  expect_identical(
    assay(out),
    matrix(
      data = c(
        6917960, 7119882, 6482082, 8740792, 6731577, 3733435, 12306793, 5283642, 11375076, 7993145,
        4461476, 4583903, 4320542, 5268651, 4067003, 2860968, 6746993, 4771641, 6541611, 3972114,
        2224393, 2194310, 2009521, 2904486, 2205489, 1467932, 6943522, 2117334, 5315106, 3128014,
        3315504, 3544847, 4728541, 4205625, 3274177, 3523505, 6137092, 4061494, 5462489, 4516166
      ),
      nrow = 4,
      byrow = TRUE,
      dimnames = list(c("1", "2", "3", "4"), c("a_1", "a_2", "a_3", "a_4", "a_5", "b_1", "b_2", "b_3", "b_4", "b_5"))
    )
  )
})

test_that("unsimplified names work", {
  # suppress output for clean test() result
  out <- suppressMessages(suppressWarnings(
    readProgenesis(
      testthat::test_path("fixtures", "all_ion_export.csv"),
      simplify_column_names = FALSE
    )
  ))
  expect_identical(
    rownames(colData(out)),
    c(
      "x250715_cond_a_1",
      "x250715_cond_a_2",
      "x250715_cond_a_3",
      "x250715_cond_a_4",
      "x250715_cond_a_5",
      "x250715_cond_b_1",
      "x250715_cond_b_2",
      "x250715_cond_b_3",
      "x250715_cond_b_4",
      "x250715_cond_b_5"
    )
  )
})

test_that("reading works as expected", {
  # suppress output for clean test() result
  out <- suppressMessages(suppressWarnings(
    readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"))
  ))
  expect_identical(out, readRDS(test_path("fixtures", "all_ion_export_read_correct.rds")))
})

test_that("metadata is generated as expected", {
  temp_path <- withr::local_tempfile()
  expect_snapshot(
    out <- readProgenesis(
      testthat::test_path("fixtures", "all_ion_export.csv"),
      generate_metadata = temp_path
    ),
    transform = \(x) sub(temp_path, "TEMP_FILEPATH", x, fixed = TRUE)
  )
  expect_identical(read.csv(temp_path), read.csv(test_path("fixtures", "all_ion_export_metadata.csv")))
})

test_that("overwrite works as expected", {
  # create a temp file and make sure it exists
  temp_path <- withr::local_tempfile()
  file.create(temp_path)
  # do not overwrite, a warning should reflect this
  expect_snapshot(
    out <- readProgenesis(
      testthat::test_path("fixtures", "all_ion_export.csv"),
      generate_metadata = temp_path,
      overwrite_metadata = FALSE
    ),
    transform = \(x) sub(temp_path, "TEMP_FILEPATH", x, fixed = TRUE)
  )
  # and nothing should have been written
  expect_equal(file.size(temp_path), 0)
  # with "overwrite" it should change the message
  expect_snapshot(
    out <- readProgenesis(
      testthat::test_path("fixtures", "all_ion_export.csv"),
      generate_metadata = temp_path,
      overwrite_metadata = TRUE
    ),
    transform = \(x) sub(temp_path, "TEMP_FILEPATH", x, fixed = TRUE)
  )
  # and the correct metadata should have been written to the temp file
  expect_identical(read.csv(temp_path), read.csv(test_path("fixtures", "all_ion_export_metadata.csv")))
})
