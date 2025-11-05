# helper functions -------------------------------------------------------------

## .deconv_grouper -------------------------------------------------------------

test_that(".deconv_grouper works with standard input", {
  mock_se <- new_mock_se()
  rowData(mock_se)$mods_ref <- c(
    "K|4|Ac",
    "K|27|Ac;K|36|Unmod",
    "K|56|Unmod;S|57|Ph",
    "S|57|Ph",
    "K|5|Unmod;K|8|Unmod;K|12|Unmod;K|16|Unmod",
    NA_character_
  )
  result <- .deconv_grouper(mock_se, "mods_ref", ";", "histone_group")

  expect_s3_class(result, "tbl_df")
  expect_named(result, c("rows_same_ptm", "ptms"))
  expect_equal(nrow(result), 6)
  expect_equal(sort(unique(result$rows_same_ptm)), as.character(1:6))
  expect_equal(
    sort(result$ptms),
    sort(c(
      "H31_HUMAN#K|4|Ac",
      "H31_HUMAN#K|27|Ac;K|36|Unmod",
      "H31_HUMAN/H3C_HUMAN#K|56|Unmod;S|57|Ph",
      "H31_HUMAN/H3C_MISALIGNED#S|57|Ph",
      "H4_HUMAN#K|5|Unmod;K|8|Unmod;K|12|Unmod;K|16|Unmod"
    ))
  )
})

test_that(".deconv_grouper handles degenerate hPTMs", {
  # rows 1 and 2 share PTMs K27Ac and K36Me3, making them degenerate
  se_degen <- new_mock_se(2)
  rowData(se_degen)$mods_ref <- c("K|27|Ac;K|36|Me3;K|37|Me2", "K|27|Ac;K|36|Me3")
  rowData(se_degen)$histone_group <- c("H31", "H31")
  result <- .deconv_grouper(se_degen, "mods_ref", ";", "histone_group")

  expect_equal(nrow(result), 3)
  # The PTMs are combined because they are indistinguishable
  expect_equal(result$ptms, c("H31#K|27|Ac;K|36|Me3", "H31#K|27|Ac;K|36|Me3", "H31#K|37|Me2"))
  expect_equal(result$rows_same_ptm, c("1", "2", "1"))
})

test_that(".deconv_grouper works with group = NULL", {
  mock_se <- new_mock_se(1)
  rowData(mock_se)$mods_ref <- "K|4|Ac"
  result <- .deconv_grouper(mock_se, "mods_ref", ";", group = NULL)

  expect_named(result, c("rows_same_ptm", "ptms"))
  expect_equal(result$ptms, "K|4|Ac")
  expect_false(any(grepl("#", result$ptms)))
})

test_that(".deconv_grouper works with a different separator", {
  mock_se <- new_mock_se(1)
  rowData(mock_se)$mods_ref <- "K4Ac|K9Me"
  result <- .deconv_grouper(mock_se, "mods_ref", "|", "histone_group")
  expect_equal(result$ptms, "H31_HUMAN#K4Ac|K9Me")
})

test_that(".deconv_grouper throws an error for missing columns", {
  mock_se <- new_mock_se(1)
  expect_error(.deconv_grouper(mock_se, "missing_col", ";", "histone_group"))
  expect_error(.deconv_grouper(mock_se, "mods_ref", ";", "missing_group"))
})

# S4 Methods -------------------------------------------------------------------

## deconvolute.SummarizedExperiment --------------------------------------------

test_that("deconvolute,SE-method works with standard inputs", {
  mock_se <- new_mock_se(5)
  rowData(mock_se)$mods_ref <- c("K|4|Ac", "K|27|Ac;K|36|Me2", "K|56|Unmod;S|57|Ph", "S|57|Ph", "K|27|Ac;K|36|Me2")
  rowData(mock_se)$histone_group <- "H3"
  result <- deconvolute(mock_se)

  expect_s4_class(result, "SummarizedExperiment")
  expect_true("hptm" %in% names(rowData(result)))
  expect_equal(nrow(result), 6)
  expect_equal(
    rowData(result)$hptm,
    c("H3#K|4|Ac", "H3#K|27|Ac;K|36|Me2", "H3#K|27|Ac;K|36|Me2", "H3#K|56|Unmod", "H3#S|57|Ph", "H3#S|57|Ph")
  )
  # Check that assay data is correctly replicated
  expect_identical(assay(result)[4, ], assay(mock_se)[3, ])
  expect_identical(assay(result)[5, ], assay(mock_se)[3, ])
  # Check rownames are new and unique
  expect_equal(rownames(result), as.character(1:6))
})

test_that("deconvolute,SE-method correctly filters out NA features", {
  mock_se <- new_mock_se(3)
  rowData(mock_se)$mods_ref <- c("K|4|Ac", NA, "S|57|Ph")
  result <- deconvolute(mock_se)

  expect_equal(nrow(result), 2)
  expect_false(any(is.na(rowData(result)$mods_ref)))
  expect_equal(rowData(result)$feature_number, c(1, 3))
})

test_that("deconvolute,SE-method forwards arguments correctly", {
  mock_se <- new_mock_se(2)
  rowData(mock_se)$mods_alt <- c("K|4|Ac", "K|27|Ac;K|36|Unmod")
  rowData(mock_se)$group_alt <- "H3_alt"

  # deconv and sep
  res_deconv <- deconvolute(mock_se, deconv = "mods_alt", group = "group_alt")
  expect_equal(nrow(res_deconv), 2)
  expect_equal(
    rowData(res_deconv)$hptm,
    c("H3_alt#K|4|Ac", "H3_alt#K|27|Ac;K|36|Unmod")
  )

  # group = NULL
  res_no_group <- deconvolute(mock_se, deconv = "mods_alt", group = NULL)
  expect_equal(rowData(res_no_group)$hptm, c("K|4|Ac", "K|27|Ac;K|36|Unmod"))
  expect_false(any(grepl("#", rowData(res_no_group)$hptm)))
})

test_that("deconvolute,SE-method fails with missing columns", {
  mock_se <- new_mock_se()
  expect_error(deconvolute(mock_se, deconv = "nonexistent_column"))
  expect_error(deconvolute(mock_se, group = "nonexistent_group"))
})

## deconvolute.QFeatures -------------------------------------------------------

test_that("deconvolute,QFeatures-method works with a single assay", {
  mock_qf <- new_mock_qf(5)
  rowData(mock_qf[["assay1"]])$mods_ref <- c(
    "K|4|Ac",
    "K|27|Ac;K|36|Me2",
    "K|56|Unmod;S|57|Ph",
    "S|57|Ph",
    "K|27|Ac;K|36|Me2"
  )
  rowData(mock_qf[["assay1"]])$histone_group <- "H3"

  # By name
  res_qf_name <- deconvolute(mock_qf, i = "assay1")
  expect_s4_class(res_qf_name, "QFeatures")
  expect_true("precursorDeconv" %in% names(res_qf_name))
  expect_equal(length(res_qf_name), 3)
  expect_equal(nrow(res_qf_name[["precursorDeconv"]]), 6)
  expect_true("hptm" %in% names(rowData(res_qf_name[["precursorDeconv"]])))
  expect_identical(rowData(mock_qf[["assay2"]]), rowData(res_qf_name[["assay2"]])) # other assay untouched
  link <- assayLink(res_qf_name, "precursorDeconv")
  expect_equal(link@from, "assay1")
  expect_equal(link@fcol, "feature_number")

  # By index
  res_qf_idx <- deconvolute(mock_qf, i = 1)
  expect_identical(res_qf_name, res_qf_idx)
})

test_that("deconvolute,QFeatures-method works with multiple assays", {
  mock_qf <- new_mock_qf()
  rowData(mock_qf[["assay1"]])$mods_ref <- c("K|4|Ac", "K|27|Ac", rep(NA, 4))
  rowData(mock_qf[["assay2"]])$mods_ref <- c("K|4|Ac", rep(NA, 5))
  new_names <- c("assay1_d", "assay2_d")

  res_qf <- deconvolute(mock_qf, i = c("assay1", "assay2"), name = new_names)
  expect_true(all(new_names %in% names(res_qf)))
  expect_equal(length(res_qf), 4)
  expect_equal(nrow(res_qf[["assay1_d"]]), 2)
  expect_equal(nrow(res_qf[["assay2_d"]]), 1)
  expect_true(assayLink(res_qf, "assay1_d")@from == "assay1")
  expect_true(assayLink(res_qf, "assay2_d")@from == "assay2")
})

test_that("deconvolute,QFeatures-method passes arguments to SE method", {
  mock_qf <- new_mock_qf()
  rowData(mock_qf[["assay1"]])$mods_alt <- c("K|4|Ac", NA, NA, NA, NA, NA)

  res_qf <- deconvolute(
    mock_qf,
    i = "assay1",
    name = "deconv_test",
    deconv = "mods_alt",
    group = NULL
  )
  se_res <- res_qf[["deconv_test"]]
  expect_false(any(grepl("#", rowData(se_res)$hptm)))
  expect_equal(rowData(se_res)$hptm, "K|4|Ac")
})

test_that("deconvolute,QFeatures-method throws errors for invalid arguments", {
  mock_qf <- new_mock_qf()
  # Mismatched lengths of i and name
  expect_error(deconvolute(mock_qf, i = "assay1", name = c("a", "b")))
  expect_error(deconvolute(mock_qf, i = c("assay1", "assay2"), name = "a"))
  # Non-existent assay
  expect_error(deconvolute(mock_qf, i = "nonexistent_assay", name = "a"))
  # Missing identifier for assay link
  mock_qf_no_id <- mock_qf
  rowData(mock_qf_no_id[[1]])$feature_number <- NULL
  expect_error(deconvolute(mock_qf_no_id, i = 1, name = "a"))
})
