# histonesFromUniprot ----------------------------------------------------------

# vcr makes it so that http requests are not repeated on every test
# this does mean that if the request were to change, this will not be tested
# see https://docs.ropensci.org/vcr/articles/vcr.html#what-happens-if-the-api-changes
test_that("histonesFromUniprot works with default arguments", {
  vcr::local_cassette("histones_default")
  res <- histonesFromUniprot()
  expect_s4_class(res, "AAStringSetList")
  expect_named(res, c("H1", "H2A", "H2B", "H3", "H4"))
  expect_true(all(vapply(res, methods::is, logical(1), "AAStringSet")))
  expect_true(all(lengths(res) > 0))
  # all sequences are the same as in a manually verified retrieval from UniProt
  verified_alignment <- readRDS(testthat::test_path("fixtures", "aligned_histones.rds"))
  expect_identical(
    purrr::map(res, as.character),
    purrr::map(verified_alignment$unaligned, as.character)
  )
})

test_that("histonesFromUniprot works with a subset of families", {
  vcr::local_cassette("histones_subset")
  res <- histonesFromUniprot(histone_families = c("H3", "H4"))
  expect_s4_class(res, "AAStringSetList")
  expect_named(res, c("H3", "H4"))
})

test_that("histonesFromUniprot works with a custom query", {
  vcr::local_cassette("histones_custom_query")
  res <- histonesFromUniprot(
    histone_families = "H3",
    query = c("organism_id:10090", "reviewed:true")
  )
  expect_s4_class(res, "AAStringSetList")
  expect_named(res, "H3")
  expect_gt(length(res[[1]]), 0)
})

test_that("histonesFromUniprot works with custom name_field", {
  vcr::local_cassette("histones_custom_name_field")
  res <- histonesFromUniprot(
    histone_families = "H3",
    name_field = "gene_names"
  )
  expect_s4_class(res, "AAStringSetList")
  expect_named(res, "H3")
  # H31_HUMAN has gene name H3C1, this has to be in the result among others
  expect_true(any(grepl("H3C1", names(res$H3))))
})

test_that("histonesFromUniprot works with 'OR' collapse", {
  vcr::local_cassette("histones_or_collapse")
  res <- histonesFromUniprot(
    histone_families = "H3",
    # human OR mouse
    query = c("organism_id:9606", "organism_id:10090"),
    collapse = " OR "
  )
  expect_s4_class(res, "AAStringSetList")
  expect_named(res, "H3")
  expect_gt(length(res[[1]]), 0)
  expect_true(any(grepl("_HUMAN", names(res$H3))))
  expect_true(any(grepl("_MOUSE", names(res$H3))))
})

test_that("histonesFromUniprot errors on invalid arguments", {
  expect_error(histonesFromUniprot(histone_families = "invalid_family"))
  expect_error(histonesFromUniprot(collapse = "INVALID"))
})

# alignHistones ----------------------------------------------------------------

test_that("alignHistones works with default Uniprot call", {
  skip_if_no_mafft()
  # reuse casette, as this is the same call as test "histonesFromUniprot works with default arguments"
  vcr::local_cassette("histones_default")
  res <- alignHistones()
  expect_true(is.list(res))
  expect_named(res, c("unaligned", "msa", "msa_ref"))
  expect_s4_class(res$unaligned, "AAStringSetList")
  expect_s4_class(res$msa, "AAStringSetList")
  expect_s4_class(res$msa_ref, "AAStringSetList")
  expect_equal(names(res$unaligned), c("H1", "H2A", "H2B", "H3", "H4"))
  expect_equal(names(res$msa), c("H1", "H2A", "H2B", "H3", "H4"))
  expect_equal(names(res$msa_ref), c("H1", "H2A", "H2B", "H3", "H4"))
  expect_named(res$msa_ref$H1, "ref_H11_HUMAN")
  expect_named(res$msa_ref$H2A, "ref_H2A1B_HUMAN")
  expect_named(res$msa_ref$H2B, "ref_H2B1J_HUMAN")
  expect_named(res$msa_ref$H3, "ref_H31_HUMAN")
  expect_named(res$msa_ref$H4, "ref_H4_HUMAN")
  # check against manually verified alignment to see if MAFFT performs as expected
  # unpack with as.character to prevent changes in the Biostrings structure to cause failure
  verified_alignment <- readRDS(testthat::test_path("fixtures", "aligned_histones.rds"))
  expect_identical(
    purrr::map(res$unaligned, as.character),
    purrr::map(verified_alignment$unaligned, as.character)
  )
  expect_identical(
    purrr::map(res$msa, as.character),
    purrr::map(verified_alignment$msa, as.character)
  )
  expect_identical(
    purrr::map(res$msa_ref, as.character),
    purrr::map(verified_alignment$msa_ref, as.character)
  )
})

test_that("alignHistones works in de novo mode", {
  skip_if_no_mafft()
  # reuse casette, as this is the same call as test "histonesFromUniprot works with default arguments"
  vcr::local_cassette("histones_default")
  res <- alignHistones(use_profiles = FALSE)

  expect_named(res, c("unaligned", "msa"))
  expect_null(res$msa_ref)
  expect_s4_class(res$msa, "AAStringSetList")
  expect_equal(names(res$msa), c("H1", "H2A", "H2B", "H3", "H4"))
  expect_vector(
    purrr::map2(res$unaligned, res$msa, \(x, y) setdiff(names(x), names(y))),
    ptype = list(character(0)),
    size = 5
  )
})

test_that("alignHistones works with custom AAStringSetList input", {
  skip_if_no_mafft()
  h3_seqs <- Biostrings::AAStringSet(c(
    H31_HUMAN = "MARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIR",
    H3C_HUMAN = "MARTKQTARKSTGGKAPRKQLATKAARKSTPSTCGVKPHRYRPGTVALREIRRYQKSTELLIR"
  ))
  unaligned_seqs <- Biostrings::AAStringSetList(H3 = h3_seqs)

  res <- alignHistones(unaligned_histones = unaligned_seqs, use_profiles = FALSE)

  expect_equal(as.character(res$msa$H3$H31_HUMAN), "MARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIR")
  expect_equal(as.character(res$msa$H3$H3C_HUMAN), "MARTKQTARKSTGGKAPRKQLATKAARKSTPSTCGV-KPHRYRPGTVALREIRRYQKSTELLIR")
})

test_that("alignHistones can write alignments to file", {
  skip_if_no_mafft()
  tmp_dir <- withr::local_tempdir()
  # reuse casette, as this is the same call as test "histonesFromUniprot works with a subset of families"
  vcr::local_cassette("histones_subset")
  # can only call once within a test, so save this result
  default_histones <- histonesFromUniprot(c("H3", "H4"))
  alignHistones(default_histones, return_alignment = TRUE, output_path = tmp_dir, overwrite = FALSE)
  expect_true(
    file.path(tmp_dir, paste0(c("H3", "H4"), "_msa.fasta")) |>
      file.exists() |>
      all()
  )

  expect_warning(
    alignHistones(default_histones, return_alignment = TRUE, output_path = tmp_dir, overwrite = FALSE),
    "Output files already exist, skipping"
  )
  expect_no_warning(
    alignHistones(default_histones, return_alignment = TRUE, output_path = tmp_dir, overwrite = TRUE)
  )
})

test_that("alignHistones handles additional MAFFT arguments", {
  skip_if_no_mafft()
  # reuse casette, as this is the same call as test "histonesFromUniprot works with a subset of families"
  vcr::local_cassette("histones_subset")
  expect_no_error(
    alignHistones(histonesFromUniprot(c("H3", "H4")), op = 3.0, maxiterate = 10)
  )
})


test_that("alignHistones handles additional MAFFT arguments", {
  skip_if_no_mafft()
  # reuse casette, as this is the same call as test "histonesFromUniprot works with a subset of families"
  vcr::local_cassette("histones_subset")
  expect_no_error(
    alignHistones(histonesFromUniprot(c("H3", "H4")), op = 3.0, maxiterate = 10)
  )
})

test_that("alignHistones returns only original sequences when requested", {
  skip_if_no_mafft()
  # reuse casette, as this is the same call as test "histonesFromUniprot works with a subset of families"
  vcr::local_cassette("histones_subset")
  # can only call once within a test, so save this result
  default_histones <- histonesFromUniprot(c("H3", "H4"))
  res <- alignHistones(default_histones, return_only_original = TRUE)
  res_all <- alignHistones(default_histones, return_only_original = FALSE)
  expect_true(all(lengths(res$msa) < lengths(res_all$msa)))
})

test_that("alignHistones errors with invalid input", {
  expect_error(alignHistones(unaligned_histones = list(a = 1)))
  expect_error(alignHistones(use_profiles = "not_a_path.fasta"))
})
