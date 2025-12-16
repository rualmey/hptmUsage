# S4 Methods -------------------------------------------------------------------

## calculateUsage.SummarizedExperiment -----------------------------------------

test_that("calculateUsage,SE works with standard workflow", {
  # Setup: Create an SE where two features will map to the same hPTM
  mock_se <- new_mock_se(nrow = 4)
  # Rows 1 and 2 will share the same hPTM
  rowData(mock_se)$mods_ref <- c("K27Me3", "K27Me3", "K9Ac", "K36Me2")
  rowData(mock_se)$histone_group <- c("H3", "H3", "H3", "H3")
  rowData(mock_se)$histone <- TRUE

  # Execution
  # Using colMeans for aggregation to make verification easy
  res_se <- calculateUsage(
    mock_se,
    usage_level = "histone",
    fun = colMeans,
    na.rm = TRUE
  )

  expect_s4_class(res_se, "SummarizedExperiment")

  # Check Dimensions:
  # Input 4 rows -> Normalization (4) -> Deconvolution (4) -> Aggregation
  # hPTMs: "H3#K27Me3" (from rows 1&2), "H3#K9Ac" (row 3), "H3#K36Me2" (row 4)
  # Result should have 3 rows
  expect_equal(nrow(res_se), 3)

  # Check Rownames
  expect_true("H3#K27Me3" %in% rownames(res_se))

  # Check Calculation Logic (Simplified):
  # 1. Normalization: usage_level="histone".
  #    Subtract mean of all histones (rows 1-4) from each value.
  # 2. Aggregation: Group rows 1 and 2 (H3#K27Me3).
  #    Value = Mean(NormValue1, NormValue2).

  vals <- assay(mock_se)
  # 1. Normalization factor (mean of all 4 rows per sample)
  norm_factors <- colMeans(vals)
  norm_vals <- sweep(vals, 2, norm_factors, "-")

  # 2. Aggregation for H3#K27Me3 (rows 1 and 2)
  expected_k27 <- colMeans(norm_vals[1:2, ])

  expect_equal(assay(res_se)["H3#K27Me3", ], expected_k27)
})

test_that("calculateUsage,SE works with different usage_levels", {
  mock_se <- new_mock_se(nrow = 3)
  rowData(mock_se)$mods_ref <- c("K27Me3", "K27Me3", "K20Ac")
  rowData(mock_se)$histone_family <- c("H3", "H3", "H4")
  rowData(mock_se)$histone_group <- c("H3_V1", "H3_V1", "H4_V1")

  # usage_level = "histone_family"
  res_se <- calculateUsage(
    mock_se,
    usage_level = "histone_family",
    fun = colMeans,
    na.rm = TRUE
  )

  # Logic:
  # H3 rows (1&2) normalized by mean of H3 rows.
  # H4 row (3) normalized by itself (mean of 1 row).
  # Aggregation: Rows 1&2 combine to H3_V1#K27Me3.

  expect_equal(nrow(res_se), 2) # H3_V1#K27Me3 and H4_V1#K20Ac

  # H4 check: Value - Mean(Value) = 0
  expect_equal(as.vector(assay(res_se)["H4_V1#K20Ac", ]), rep(0, 7))
})

test_that("calculateUsage,SE works for target = 'variant'", {
  mock_se <- new_mock_se(nrow = 3)
  rowData(mock_se)$histone_group <- c("G1", "G1", "G2")
  rowData(mock_se)$histone <- TRUE # for usage_level="histone" normalization

  # Execution
  res_se <- calculateUsage(
    mock_se,
    target = "variant",
    usage_level = "histone",
    fun = colMeans, # Simplify aggregation for verification
    na.rm = TRUE
  )

  expect_s4_class(res_se, "SummarizedExperiment")

  # Check Dimensions: 2 groups -> 2 rows
  expect_equal(nrow(res_se), 2)
  expect_equal(sort(rownames(res_se)), c("G1", "G2"))

  # Check Logic:
  # 1. Normalization happens first (subtract mean of all)
  # 2. Deconvolution is SKIPPED
  # 3. Aggregation happens by histone_group

  vals <- assay(mock_se)
  norm_factors <- colMeans(vals)
  norm_vals <- sweep(vals, 2, norm_factors, "-")

  # Expected value for G1 (mean of first two rows)
  expected_g1 <- colMeans(norm_vals[1:2, , drop = FALSE])
  # Expected value for G2 (value of third row)
  expected_g2 <- norm_vals[3, ]

  expect_equal(assay(res_se)["G1", ], expected_g1)
  expect_equal(assay(res_se)["G2", ], expected_g2)
})

test_that("calculateUsage,SE passes arguments to internal functions", {
  mock_se <- new_mock_se(nrow = 2)
  rowData(mock_se)$alt_mods <- c("ModA|ModB", "ModA|ModB")
  rowData(mock_se)$alt_group <- c("G1", "G1")

  # Test custom deconv column, separator, grouping, and aggregation function
  res_se <- calculateUsage(
    mock_se,
    deconv = "alt_mods",
    sep = "|",
    group = "alt_group",
    usage_level = "histone"
  )

  expect_equal(nrow(res_se), 1)
  expect_equal(rownames(res_se), "G1#ModA|ModB")
})

test_that("calculateUsage,SE handles empty inputs/outputs", {
  # Case: Input has no valid mods for deconvolution
  mock_se <- new_mock_se(nrow = 2)
  rowData(mock_se)$mods_ref <- c("K27Me3", NA_character_)

  # deconvolute filters out NAs, resulting in 1 row
  res_se <- suppressWarnings(suppressMessages(calculateUsage(mock_se)))

  expect_equal(nrow(res_se), 1)
})

## calculateUsage.QFeatures ----------------------------------------------------

test_that("calculateUsage,QFeatures works with default parameters", {
  mock_qf <- new_mock_qf(nrow = 4)
  rowData(mock_qf[["assay1"]])$mods_ref <- c("K1", "K1", "K2", "K3")
  rowData(mock_qf[["assay1"]])$histone_group <- "H3"
  rowData(mock_qf[["assay1"]])$histone <- TRUE

  res_qf <- calculateUsage(mock_qf, i = "assay1")

  # Check if expected assays are added
  # 1. Normalization: "assay1NormHistone"
  # 2. Deconvolution: "assay1NormHistoneDeconv"
  # 3. Aggregation (Result): "ptm" (default name)
  expected_names <- c("assay1", "assay1NormHistone", "assay1NormHistoneDeconv", "ptm")
  expect_true(all(expected_names %in% names(res_qf)))

  # Check final assay dimensions
  # "K1" rows combine -> 3 unique hPTMs
  expect_equal(nrow(res_qf[["ptm"]]), 3)

  # Check links chain
  # assay1 -> assay1NormHistone
  link1 <- assayLink(res_qf, "assay1NormHistone")
  expect_equal(link1@from, "assay1")

  # assay1NormHistone -> assay1NormHistoneDeconv
  link2 <- assayLink(res_qf, "assay1NormHistoneDeconv")
  expect_equal(link2@from, "assay1NormHistone")

  # assay1NormHistoneDeconv -> ptm
  link3 <- assayLink(res_qf, "ptm")
  expect_equal(link3@from, "assay1NormHistoneDeconv")
})

test_that("calculateUsage,QFeatures works with custom names and multiple assays", {
  mock_qf <- new_mock_qf(nrow = 2)
  rowData(mock_qf[["assay1"]])$mods_ref <- c("K1", "K2")
  rowData(mock_qf[["assay1"]])$histone_group <- "H3"
  rowData(mock_qf[["assay1"]])$histone_family <- "H3"
  rowData(mock_qf[["assay2"]])$mods_ref <- c("K1", "K2")
  rowData(mock_qf[["assay2"]])$histone_group <- "H3"
  rowData(mock_qf[["assay2"]])$histone_family <- "H3"

  res_qf <- calculateUsage(
    mock_qf,
    i = c("assay1", "assay2"),
    name = c("ptm1", "ptm2"),
    usage_level = "histone_family"
  )

  # Check for intermediate names based on usage_level
  # "assay1NormHistone_family", "assay1NormHistone_familyDeconv", "ptm1"
  expect_true("assay1NormHistone_family" %in% names(res_qf))
  expect_true("assay1NormHistone_familyDeconv" %in% names(res_qf))
  expect_true("ptm1" %in% names(res_qf))
  expect_true("ptm2" %in% names(res_qf))

  expect_equal(nrow(res_qf[["ptm1"]]), 2)
})

test_that("calculateUsage,QFeatures works for target = 'variant'", {
  mock_qf <- new_mock_qf()
  # Ensure assay1 has necessary columns
  rowData(mock_qf[["assay1"]])$histone <- TRUE
  rowData(mock_qf[["assay1"]])$histone_group <- c(rep("G1", 3), rep("G2", 3))

  res_qf <- suppressWarnings(suppressMessages(calculateUsage(
    mock_qf,
    i = "assay1",
    target = "variant",
    usage_level = "histone"
  )))

  # Check Assays
  # Should have:
  # 1. "assay1" (Original)
  # 2. "assay1NormHistone" (Normalized)
  # 3. "variant" (Result, default name)
  # Should NOT have: "...Deconv" assay

  expect_true("assay1NormHistone" %in% names(res_qf))
  expect_true("variant" %in% names(res_qf))
  expect_false("assay1NormHistoneDeconv" %in% names(res_qf))

  # Check Dimensions
  # 2 groups (G1, G2) -> 2 rows in final assay
  expect_equal(nrow(res_qf[["variant"]]), 2)

  # Check Linkage
  # The 'variant' assay should be linked directly to 'assay1NormHistone'
  # skipping the deconvolution step
  link <- assayLink(res_qf, "variant")
  expect_equal(link@from, "assay1NormHistone")
})

test_that("calculateUsage,QFeatures passes arguments to inner methods", {
  mock_qf <- new_mock_qf(nrow = 2)
  rowData(mock_qf[["assay1"]])$my_mods <- c("A", "A")
  rowData(mock_qf[["assay1"]])$my_grp <- c("G", "G")
  rowData(mock_qf[["assay1"]])$my_id <- c("ID1", "ID2")

  # Testing identifier passing and custom deconvolution args
  res_qf <- calculateUsage(
    mock_qf,
    i = "assay1",
    name = "result",
    identifier = "my_id",
    deconv = "my_mods",
    group = "my_grp",
    fun = colSums # Aggregation arg
  )

  # Check linkage uses correct identifier
  # Link: assay1NormHistone -> assay1NormHistoneDeconv
  # This link is created by deconvolute() which uses 'identifier'
  link <- assayLink(res_qf, "assay1NormHistoneDeconv")
  expect_equal(link@from, "assay1NormHistone")
  expect_equal(link@fcol, "my_id")

  # Check aggregation occurred (2 rows -> 1 row)
  expect_equal(nrow(res_qf[["result"]]), 1)
  expect_equal(rownames(res_qf[["result"]]), "G#A")
})

test_that("calculateUsage,QFeatures errors with mismatched arguments", {
  mock_qf <- new_mock_qf()
  expect_error(
    calculateUsage(mock_qf, i = "assay1", name = c("n1", "n2")),
    "length\\(i\\) == length\\(name\\) is not TRUE"
  )
})

test_that("calculateUsage,QFeatures errors with invalid usage_level", {
  mock_qf <- new_mock_qf()
  expect_error(
    calculateUsage(mock_qf, i = "assay1", usage_level = "invalid"),
    "'arg' should be one of"
  )
})
