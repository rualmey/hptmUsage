test_that("matchHistones works for SummarizedExperiment", {
  histone_seqs <- Biostrings::AAStringSetList(
    H1 = Biostrings::AAStringSet(c(H1.1 = "MSET", AMBIG = "MKQT")),
    H3 = Biostrings::AAStringSet(c(H3.1 = "MARTKQTAR", H3.ALT = "MTKQTAR")),
    H4 = Biostrings::AAStringSet(c(H4 = "MSGRGSGRG"))
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 7, ncol = 2, dimnames = list(c(as.character(1:7)), c("S1", "S2")))),
    rowData = S4Vectors::DataFrame(
      sequence = c(
        # H1.1
        "SET",
        # H3.1/H3.ALT
        "TKQTAR",
        # H3.1
        "RTKQTAR",
        # family ambiguous
        "KQT",
        # one vs two matches per variant
        "MSGRG",
        "SGRG",
        # no match
        "NOTAHISTONE"
      ),
      row.names = as.character(1:7)
    )
  )

  suppressMessages(expect_message(
    matchHistones(se, histone_seqs, progress = FALSE),
    "Sequences with ambiguous family assignment: KQT"
  ))
  suppressMessages(expect_message(
    matchHistones(se, histone_seqs, progress = FALSE),
    "Sequences with multiple possible positions per variant: SGRG"
  ))
  res_se <- suppressMessages(matchHistones(se, histone_seqs, progress = FALSE))

  expect_s4_class(res_se, "SummarizedExperiment")
  expect_equal(nrow(res_se), 5)
  expect_equal(rownames(res_se), c("1", "2", "3", "5", "7"))

  rd_res <- rowData(res_se)
  expect_true(all(
    c("histone", "histone_family", "core_histone", "histone_group", "start_index", "end_index") %in% colnames(rd_res)
  ))

  expect_true(rd_res["1", "histone"])
  expect_equal(rd_res["1", "histone_family"], "H1")
  expect_false(rd_res["1", "core_histone"])
  expect_equal(rd_res["1", "histone_group"], "H1.1")
  expect_equal(unlist(rd_res["1", "start_index"]), 1)
  expect_equal(unlist(rd_res["1", "end_index"]), 3)

  expect_true(rd_res["2", "histone"])
  expect_equal(rd_res["2", "histone_family"], "H3")
  expect_true(rd_res["2", "core_histone"])
  expect_equal(rd_res["2", "histone_group"], "H3.1/H3.ALT")
  expect_equal(unlist(rd_res["2", "start_index"]), c(3, 1))
  expect_equal(unlist(rd_res["2", "end_index"]), c(8, 6))

  expect_true(rd_res["3", "histone"])
  expect_equal(rd_res["3", "histone_family"], "H3")
  expect_true(rd_res["3", "core_histone"])
  expect_equal(rd_res["3", "histone_group"], "H3.1")
  expect_equal(unlist(rd_res["3", "start_index"]), 2)
  expect_equal(unlist(rd_res["3", "end_index"]), 8)

  expect_true(rd_res["5", "histone"])
  expect_equal(rd_res["5", "histone_family"], "H4")
  expect_true(rd_res["5", "core_histone"])
  expect_equal(rd_res["5", "histone_group"], "H4")
  expect_equal(unlist(rd_res["5", "start_index"]), 0)
  expect_equal(unlist(rd_res["5", "end_index"]), 4)

  expect_true(is.na(rd_res["7", "histone"]))
  expect_true(is.na(rd_res["7", "histone_family"]))
  expect_true(is.na(rd_res["7", "core_histone"]))
  expect_true(is.na(rd_res["7", "histone_group"]))
  expect_true(is.na(rd_res["7", "start_index"]))
  expect_true(is.na(rd_res["7", "end_index"]))
})

test_that("matchHistones handles other arguments for SE", {
  histone_seqs <- Biostrings::AAStringSetList(
    H1 = Biostrings::AAStringSet(c(H1.1 = "MSET", AMBIG = "MKQT")),
    H3 = Biostrings::AAStringSet(c(H3.1 = "MARTKQTAR"))
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 4, ncol = 2, dimnames = list(c(as.character(1:4)), c("S1", "S2")))),
    rowData = S4Vectors::DataFrame(
      peptide_sequence = c(
        # H1.1
        "SET",
        # family ambiguous
        "KQT",
        # no match
        "NOTAHISTONE",
        "ALSONOTAHISTONE"
      ),
      row.names = as.character(1:4)
    )
  )

  res_alt_col <- suppressMessages(
    matchHistones(se, histone_seqs, sequence_col = "peptide_sequence", progress = FALSE)
  )
  expect_equal(nrow(res_alt_col), 3)

  se_no_match <- se[3:4, ]
  expect_error(
    matchHistones(se_no_match, histone_seqs, sequence_col = "peptide_sequence", progress = FALSE),
    "No histone sequences could be matched"
  )
})

test_that("matchHistones works for QFeatures", {
  histone_seqs <- Biostrings::AAStringSetList(
    H1 = Biostrings::AAStringSet(c(H1.1 = "MSET", AMBIG = "MKQT")),
    H3 = Biostrings::AAStringSet(c(H3.1 = "MARTKQTAR", H3.ALT = "MAKTKQTAR")),
    H4 = Biostrings::AAStringSet(c(H4 = "MSGRGSGRG"))
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 7, ncol = 2, dimnames = list(c(as.character(1:7)), c("S1", "S2")))),
    rowData = S4Vectors::DataFrame(
      sequence = c(
        # H1.1
        "SET",
        # H3.1/H3.ALT
        "TKQTAR",
        # H3.1
        "RTKQTAR",
        # family ambiguous
        "KQT",
        # one vs two matches per variant
        "MSGRG",
        "SGRG",
        # no match
        "NOTAHISTONE"
      ),
      row.names = as.character(1:7)
    )
  )
  qf <- QFeatures::QFeatures(
    list(assay1 = se, assay2 = se),
    colData = S4Vectors::DataFrame(Var1 = rnorm(2), Var2 = LETTERS[1:2], row.names = c("S1", "S2"))
  ) |>
    QFeatures::addAssayLinkOneToOne("assay1", "assay2")

  res_qf <- suppressMessages(matchHistones(qf, histone_seqs, i = "assay1", progress = FALSE))

  expect_s4_class(res_qf, "QFeatures")
  expect_length(assayLinks(res_qf, 2), 2)
  expect_equal(nrow(res_qf[[1]]), 5)
  expect_equal(nrow(res_qf[[2]]), 7) # unchanged
  expect_true("histone_family" %in% names(rowData(res_qf[[1]])))
  expect_false("histone_family" %in% names(rowData(res_qf[[2]])))

  res_qf_idx <- suppressMessages(matchHistones(qf, histone_seqs, i = 1, progress = FALSE))
  expect_identical(res_qf, res_qf_idx)

  res_qf_multi <- suppressMessages(matchHistones(qf, histone_seqs, i = 1:2, progress = FALSE))
  expect_equal(nrow(res_qf_multi[[1]]), 5)
  expect_equal(nrow(res_qf_multi[[2]]), 5)
  expect_true("histone_family" %in% names(rowData(res_qf_multi[[2]])))
})
