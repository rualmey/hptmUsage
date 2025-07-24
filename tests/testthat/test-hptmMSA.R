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
})


test_that("histonesFromUniprot works with a subset of families", {
  vcr::local_cassette("histones_subset")
  res <- histonesFromUniprot(histone_families = c("H2B", "H4"))
  expect_s4_class(res, "AAStringSetList")
  expect_named(res, c("H2B", "H4"))
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
