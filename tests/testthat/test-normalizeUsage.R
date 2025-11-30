# S4 Methods -------------------------------------------------------------------

## normalizeUsage.SummarizedExperiment -----------------------------------------

test_that("normalizeUsage,SE works with usage_level = 'histone'", {
  mock_se <- new_mock_se()
  # normalizeUsage subtracts the aggregated abundance from the features
  # usage_level = "histone" uses the boolean column.
  # features 1-5 are histone=TRUE, feature 6 is histone=FALSE.

  # Calculate expected manually using robustSummary (default)
  # We use colMeans here for simplicity by overriding 'fun' to verify the logic
  res_se <- suppressMessages(normalizeUsage(mock_se, usage_level = "histone", fun = colMeans, na.rm = TRUE))

  scaling_factors <- matrix(NA_real_, nrow = nrow(mock_se), ncol = ncol(mock_se))

  # Group TRUE
  idx_true <- which(rowData(mock_se)$histone)
  agg_true <- colMeans(assay(mock_se)[idx_true, ], na.rm = TRUE)
  scaling_factors[idx_true, ] <- rep(agg_true, each = length(idx_true))

  # Group FALSE
  idx_false <- which(!rowData(mock_se)$histone)
  agg_false <- colMeans(assay(mock_se)[idx_false, , drop = FALSE], na.rm = TRUE)
  scaling_factors[idx_false, ] <- rep(agg_false, each = length(idx_false))

  expected_assay <- assay(mock_se) - scaling_factors

  expect_s4_class(res_se, "SummarizedExperiment")
  expect_equal(assay(res_se), expected_assay)
  # rowData should remain unchanged
  expect_identical(rowData(res_se), rowData(mock_se))
})

test_that("normalizeUsage,SE works with usage_level = 'histone_family'", {
  mock_se <- new_mock_se()
  # H3 family: rows 1, 2, 3, 4
  # H4 family: row 5
  # NA family: row 6
  mock_se <- mock_se[1:5, ]

  res_se <- suppressMessages(normalizeUsage(mock_se, usage_level = "histone_family", fun = colMeans, na.rm = TRUE))

  # Calculate expected for H3
  idx_h3 <- which(rowData(mock_se)$histone_family == "H3")
  agg_h3 <- colMeans(assay(mock_se)[idx_h3, ], na.rm = TRUE)

  # Check row 1 (H3)
  expect_equal(
    assay(res_se)[1, ],
    assay(mock_se)[1, ] - agg_h3
  )

  # Calculate expected for H4 (single row)
  idx_h4 <- which(rowData(mock_se)$histone_family == "H4")
  agg_h4 <- assay(mock_se)[idx_h4, ] # mean of 1 row is the row itself

  # Check row 5 (H4) -> Should be 0 (value - itself)
  expect_equal(
    assay(res_se)[5, ],
    setNames(rep(0, 7), LETTERS[1:7])
  )
})

test_that("normalizeUsage,SE works with usage_level = 'histone_group'", {
  mock_se <- new_mock_se()[1:5, ]
  # H31_HUMAN: rows 1, 2
  # H31_HUMAN/H3C_HUMAN: row 3
  # ...

  res_se <- suppressMessages(normalizeUsage(mock_se, usage_level = "histone_group", fun = colMeans, na.rm = TRUE))

  # Check group H31_HUMAN
  idx_grp <- which(rowData(mock_se)$histone_group == "H31_HUMAN")
  agg_grp <- colMeans(assay(mock_se)[idx_grp, ], na.rm = TRUE)

  expect_equal(
    assay(res_se)[1, ],
    assay(mock_se)[1, ] - agg_grp
  )
  expect_equal(
    assay(res_se)[2, ],
    assay(mock_se)[2, ] - agg_grp
  )
})

test_that("normalizeUsage,SE handles NAs in data correctly", {
  mock_se <- new_mock_se()
  # Introduce NAs in assay
  assay(mock_se)[1, 1:3] <- NA
  assay(mock_se)[2, 1] <- NA

  # Without na.rm = TRUE, this would result in NAs in the aggregation
  res_se <- suppressMessages(normalizeUsage(mock_se, usage_level = "histone", fun = colMeans, na.rm = TRUE))

  # Verify aggregation ignored NAs
  idx_true <- which(rowData(mock_se)$histone)
  agg_col1 <- mean(assay(mock_se)[idx_true, 1], na.rm = TRUE)

  # Expected result for a non-NA value in col 1 (e.g., row 3)
  expect_equal(
    assay(res_se)[3, 1],
    assay(mock_se)[3, 1] - agg_col1
  )

  # Expected result for an NA value remains NA
  expect_true(is.na(assay(res_se)[1, 1]))
})

test_that("normalizeUsage,SE handles custom aggregation functions via ...", {
  mock_se <- new_mock_se()[1:5, ]
  # Use median instead of default robustSummary
  res_se <- suppressMessages(normalizeUsage(
    mock_se,
    usage_level = "histone_family",
    fun = MatrixGenerics::colMedians,
    na.rm = TRUE
  ))

  idx_h3 <- which(rowData(mock_se)$histone_family == "H3")
  agg_h3_median <- MatrixGenerics::colMedians(assay(mock_se)[idx_h3, ], na.rm = TRUE)

  expect_equal(
    assay(res_se)[1, ],
    assay(mock_se)[1, ] - agg_h3_median
  )
})

test_that("normalizeUsage,SE throws error for invalid usage_level argument", {
  mock_se <- new_mock_se()
  expect_error(
    normalizeUsage(mock_se, usage_level = "invalid_level"),
    "'arg' should be one of"
  )
})

test_that("normalizeUsage,SE throws error if usage_level column is missing", {
  mock_se <- new_mock_se()
  rowData(mock_se)$histone_family <- NULL
  expect_error(
    normalizeUsage(mock_se, usage_level = "histone_family"),
    "'fcol' not found"
  )
})

test_that("normalizeUsage,SE handles NA in usage_level column", {
  mock_se <- new_mock_se()
  # Row 6 has histone_family = NA
  res_se <- expect_error(
    normalizeUsage(mock_se, usage_level = "histone_family", fun = colMeans, na.rm = TRUE),
    "'usage_level' should not contain"
  )
})

## normalizeUsage.QFeatures ----------------------------------------------------

test_that("normalizeUsage,QFeatures works with default parameters", {
  mock_qf <- new_mock_qf()

  res_qf <- suppressWarnings(suppressMessages(normalizeUsage(mock_qf, i = "assay1")))

  expect_s4_class(res_qf, "QFeatures")
  expect_true("precursorNorm" %in% names(res_qf))
  # Should add 1 new assay
  expect_length(res_qf, length(mock_qf) + 1)

  # Check assay link (should be OneToOne)
  link <- assayLink(res_qf, "precursorNorm")
  expect_equal(link@from, "assay1")
  expect_equal(link@fcol, "._oneToOne")
})

test_that("normalizeUsage,QFeatures works with custom name and usage_level", {
  mock_qf <- new_mock_qf()[1:5, ]

  res_qf <- suppressMessages(normalizeUsage(
    mock_qf,
    i = "assay1",
    name = "norm_h_fam",
    usage_level = "histone_family",
    fun = colMeans,
    na.rm = TRUE
  ))

  expect_true("norm_h_fam" %in% names(res_qf))

  # Verify calculation matches SE equivalent
  se_expected <- suppressMessages(normalizeUsage(
    mock_qf[["assay1"]],
    usage_level = "histone_family",
    fun = colMeans,
    na.rm = TRUE
  ))
  expect_equal(assay(res_qf[["norm_h_fam"]]), assay(se_expected))
})

test_that("normalizeUsage,QFeatures processes multiple assays", {
  mock_qf <- new_mock_qf()
  # mock_qf has assay1 and assay2

  res_qf <- suppressMessages(normalizeUsage(
    mock_qf,
    i = c("assay1", "assay2"),
    name = c("norm1", "norm2")
  ))

  expect_true(all(c("norm1", "norm2") %in% names(res_qf)))
  expect_length(res_qf, length(mock_qf) + 2)

  expect_equal(assayLink(res_qf, "norm1")@from, "assay1")
  expect_equal(assayLink(res_qf, "norm2")@from, "assay2")
})

test_that("normalizeUsage,QFeatures errors with mismatched i and name lengths", {
  mock_qf <- new_mock_qf()
  expect_error(
    normalizeUsage(mock_qf, i = "assay1", name = c("n1", "n2")),
    "length\\(i\\) == length\\(name\\) is not TRUE"
  )
})

test_that("normalizeUsage,QFeatures errors with invalid assay index", {
  mock_qf <- new_mock_qf()
  expect_error(
    normalizeUsage(mock_qf, i = "non_existent"),
    "The following assay" # The following assay(s) is/are not found:non_existent
  )
})

test_that("normalizeUsage,QFeatures passes ... arguments correctly", {
  mock_qf <- new_mock_qf()

  # Pass a dummy function that returns 0s to verify it's being used
  dummy_agg <- function(x, ...) {
    rep(0, ncol(x))
  }

  res_qf <- suppressMessages(normalizeUsage(mock_qf, i = "assay1", fun = dummy_agg))

  # If aggregation is 0, result = original - 0 = original
  expect_equal(
    assay(res_qf[["precursorNorm"]]),
    assay(mock_qf[["assay1"]])
  )
})
