test_that("replaceColData works with tbl_df", {
  # objects to replace colData in
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = matrix(1:12, 3, 4, dimnames = list(NULL, paste0("S", 1:4))),
    colData = S4Vectors::DataFrame(row.names = paste0("S", 1:4))
  )
  qf <- QFeatures::QFeatures(list(assay1 = se))
  # new colData
  new_cd_tbl <- tibble::tibble(
    quantCols = paste0("S", 1:4),
    original_name = paste0("Sample_", 1:4),
    group = factor(c("A", "A", "B", "B")),
    include = TRUE,
    outlier = FALSE
  )
  expected_cd <- as(new_cd_tbl, "DataFrame")
  rownames(expected_cd) <- colnames(se)
  # test for SummarizedExperiment
  se_new <- replaceColData(se, new_cd_tbl)
  expect_s4_class(se_new, "SummarizedExperiment")
  expect_identical(
    colData(se_new),
    expected_cd
  )
  # test for QFeatures
  qf_new <- replaceColData(qf, new_cd_tbl)
  expect_s4_class(qf_new, "QFeatures")
  expect_identical(
    colData(qf_new),
    expected_cd
  )
})

test_that("replaceColData works with file path", {
  # objects to replace colData in
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = matrix(1:12, 3, 4, dimnames = list(NULL, paste0("S", 1:4))),
    colData = S4Vectors::DataFrame(row.names = paste0("S", 1:4))
  )
  qf <- QFeatures::QFeatures(list(assay1 = se))
  # new colData
  new_cd_tbl <- tibble::tibble(
    quantCols = paste0("S", 1:4),
    original_name = paste0("Sample_", 1:4),
    group = factor(c("A", "A", "B", "B")),
    include = TRUE,
    outlier = FALSE,
    extra_num = 1:4
  )
  temp_csv <- withr::local_tempfile(fileext = ".csv")
  readr::write_csv(new_cd_tbl, temp_csv)
  # test for SummarizedExperiment, do not test for custom coltype
  se_new <- replaceColData(se, temp_csv)
  expect_s4_class(se_new, "SummarizedExperiment")
  expect_s3_class(colData(se_new)$extra_num, "factor")
  # test for QFeatures, do test for custom coltype
  qf_new <- replaceColData(qf, temp_csv, custom_coltypes = list(extra_num = readr::col_integer()))
  expect_s4_class(qf_new, "QFeatures")
  expect_type(colData(qf_new)$extra_num, "integer")
})

test_that("replaceColData handles errors correctly", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = matrix(1:12, 3, 4, dimnames = list(NULL, paste0("S", 1:4))),
    colData = DataFrame(row.names = paste0("S", 1:4))
  )
  # non-existent file
  expect_error(replaceColData(se, "nonexistent.csv"))
  # not all expected columns in tbl
  bad_cd_tbl <- tibble::tibble(quantCols = paste0("S", 1:4))
  expect_error(replaceColData(se, bad_cd_tbl))
  # not all expected columns in csv
  temp_csv <- withr::local_tempfile(fileext = ".csv")
  readr::write_csv(bad_cd_tbl, temp_csv)
  expect_error(replaceColData(se, temp_csv))
  # mismatch between samples
  mismatch_cd_tbl_1 <- tibble::tibble(
    quantCols = paste0("S", c(1:3, 5)),
    original_name = "n",
    group = factor("A"),
    include = TRUE,
    outlier = FALSE
  )
  expect_error(replaceColData(se, mismatch_cd_tbl_1), "Samples do not match")
  # extra sample in colData
  mismatch_cd_tbl_1 <- tibble::tibble(
    quantCols = paste0("S", c(1:5)),
    original_name = "n",
    group = factor("A"),
    include = TRUE,
    outlier = FALSE
  )
  expect_error(replaceColData(se, mismatch_cd_tbl_1), "Samples do not match")
  # sample missing in colData
  mismatch_cd_tbl_2 <- tibble::tibble(
    quantCols = paste0("S", 1:3),
    original_name = "n",
    group = factor("A"),
    include = TRUE,
    outlier = FALSE
  )
  expect_error(replaceColData(se, mismatch_cd_tbl_2), "Samples do not match")
})

test_that("replaceColData warns for incorrect column types", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = matrix(1:12, 3, 4, dimnames = list(NULL, paste0("S", 1:4))),
    colData = DataFrame(row.names = paste0("S", 1:4))
  )
  bad_types_tbl <- tibble::tibble(
    quantCols = paste0("S", 1:4),
    original_name = paste0("Sample_", 1:4),
    group = "A",
    include = 1,
    outlier = FALSE
  )
  expect_warning(replaceColData(se, bad_types_tbl))
})

test_that("replaceColData preserves original sample order", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = matrix(1:12, 3, 4, dimnames = list(NULL, paste0("S", c(3, 1, 4, 2)))),
    colData = DataFrame(row.names = paste0("S", c(3, 1, 4, 2)))
  )
  shuffled_cd_tbl <- tibble::tibble(
    quantCols = paste0("S", 4:1),
    original_name = paste0("Sample_", 4:1),
    group = factor(c("A", "A", "B", "B")),
    include = TRUE,
    outlier = FALSE
  )

  se_new <- replaceColData(se, shuffled_cd_tbl)
  expect_identical(rownames(colData(se)), colData(se_new)$quantCols)
})
