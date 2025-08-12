## tagContaminants.SummarizedExperiment ----------------------------------------

test_that("SummarizedExperiment tags with default 'universal' library", {
  mock_se <- new_mock_se()
  rowData(mock_se)$sequence[[6]] <- "VATVSLPR" # TRYP_PIG

  mock_se <- tagContaminants(mock_se)
  expect_s4_class(mock_se, "SummarizedExperiment")
  expect_true("contaminant" %in% colnames(rowData(mock_se)))
  expect_equal(rowData(mock_se)$contaminant, c(rep.int(FALSE, 5), TRUE))
})

test_that("SummarizedExperiment tags with a specified library", {
  mock_se <- new_mock_se()
  rowData(mock_se)$sequence[[6]] <- "KYTSWYVALKR" # FGF2_HUMAN from iPSC Medium not in universal library

  mock_se <- tagContaminants(mock_se)
  expect_equal(rowData(mock_se)$contaminant, rep.int(FALSE, 6))
  mock_se <- tagContaminants(mock_se, library = "stem_cell_culture")
  expect_equal(rowData(mock_se)$contaminant, c(rep.int(FALSE, 5), TRUE))
})

test_that("SummarizedExperiment works with custom sequence column name", {
  mock_se <- new_mock_se()
  rowData(mock_se)$sequence[[6]] <- "VATVSLPR" # TRYP_PIG
  colnames(rowData(mock_se))[[2]] <- "pep_seq"

  mock_se <- tagContaminants(mock_se, sequence = "pep_seq")
  expect_equal(rowData(mock_se)$contaminant, c(rep.int(FALSE, 5), TRUE))
})

test_that("SummarizedExperiment handles no matching contaminants", {
  mock_se <- new_mock_se()
  mock_se <- tagContaminants(mock_se)
  expect_equal(rowData(mock_se)$contaminant, rep.int(FALSE, 6))
})

test_that("SummarizedExperiment handles all features being contaminants", {
  mock_se <- new_mock_se()
  rowData(mock_se)$sequence <- "VATVSLPR" # TRYP_PIG
  mock_se <- tagContaminants(mock_se)
  expect_equal(rowData(mock_se)$contaminant, rep.int(TRUE, 6))
})

test_that("SummarizedExperiment throws error for invalid library name", {
  mock_se <- new_mock_se()
  expect_error(tagContaminants(mock_se, library = "invalid_library"))
  expect_error(tagContaminants(mock_se, library = 1))
  expect_error(tagContaminants(mock_se, library = c("universal", "stem_cell_culture")))
})

test_that("SummarizedExperiment throws error for invalid sequence column", {
  mock_se <- new_mock_se()
  expect_error(tagContaminants(mock_se, sequence = "invalid_sequence_col"))
  expect_error(tagContaminants(mock_se, sequence = 1))
  expect_error(tagContaminants(mock_se, sequence = c("sequence", "other")))
})

## tagContaminants.QFeatures ---------------------------------------------------

test_that("QFeatures tags a single assay by name", {
  mock_qf <- new_mock_qf()
  rowData(mock_qf)$assay1$sequence[[6]] <- "VATVSLPR" # TRYP_PIG

  res_qf <- tagContaminants(mock_qf, i = "assay1")
  # unchanged
  expect_s4_class(res_qf, "QFeatures")
  expect_length(assayLinks(res_qf, 2), 2)
  expect_equal(rowData(res_qf[[2]]), rowData(mock_qf[[2]]))
  # changed
  expect_true("contaminant" %in% colnames(rowData(res_qf[["assay1"]])))
  expect_equal(rowData(res_qf)$assay1$contaminant, c(rep.int(FALSE, 5), TRUE))
})

test_that("QFeatures tags a single assay by index", {
  mock_qf <- new_mock_qf()
  rowData(mock_qf)$assay1$sequence[[6]] <- "VATVSLPR" # TRYP_PIG

  res_qf <- tagContaminants(mock_qf, i = 1)
  # unchanged
  expect_s4_class(res_qf, "QFeatures")
  expect_length(assayLinks(res_qf, 2), 2)
  expect_equal(rowData(res_qf[[2]]), rowData(mock_qf[[2]]))
  # changed
  expect_true("contaminant" %in% colnames(rowData(res_qf[["assay1"]])))
  expect_equal(rowData(res_qf)$assay1$contaminant, c(rep.int(FALSE, 5), TRUE))
})

test_that("QFeatures tags multiple assays", {
  mock_qf <- new_mock_qf()
  rowData(mock_qf)$assay1$sequence[[6]] <- "VATVSLPR" # TRYP_PIG

  res_qf <- tagContaminants(mock_qf, i = c("assay1", "assay2"))
  # unchanged
  expect_s4_class(res_qf, "QFeatures")
  expect_length(assayLinks(res_qf, 2), 2)
  # changed
  expect_true("contaminant" %in% colnames(rowData(res_qf[["assay1"]])))
  expect_true("contaminant" %in% colnames(rowData(res_qf[["assay2"]])))
  expect_equal(rowData(res_qf)$assay1$contaminant, c(rep.int(FALSE, 5), TRUE))
  expect_equal(rowData(res_qf)$assay2$contaminant, rep.int(FALSE, 6))
})

test_that("QFeatures forwards arguments correctly", {
  mock_qf <- new_mock_qf()

  rowData(mock_qf)$assay1$sequence[[6]] <- "KYTSWYVALKR" # FGF2_HUMAN from iPSC Medium not in universal library
  res_qf <- tagContaminants(mock_qf, i = "assay1")
  expect_equal(rowData(res_qf)$assay1$contaminant, rep.int(FALSE, 6))
  res_qf <- tagContaminants(mock_qf, i = "assay1", library = "stem_cell_culture")
  expect_equal(rowData(res_qf)$assay1$contaminant, c(rep.int(FALSE, 5), TRUE))

  rowData(mock_qf)$assay2$sequence[[6]] <- "VATVSLPR" # FGF2_HUMAN from iPSC Medium not in universal library
  colnames(rowData(mock_qf)$assay2)[[2]] <- "pep_seq"
  res_qf <- tagContaminants(mock_qf, i = "assay2", sequence = "pep_seq")
  expect_equal(rowData(res_qf)$assay2$contaminant, c(rep.int(FALSE, 5), TRUE))
})

test_that("QFeatures propagates errors from SE method", {
  mock_qf <- new_mock_qf()
  expect_error(tagContaminants(mock_qf, i = "assay1", library = "invalid"))
  expect_error(tagContaminants(mock_qf, i = "assay1", sequence = "invalid"))
})

test_that("QFeatures errors on invalid assay index", {
  mock_qf <- new_mock_qf()
  expect_error(tagContaminants(mock_qf, i = "invalid_assay_name"))
  expect_error(tagContaminants(mock_qf, i = 99))
})
