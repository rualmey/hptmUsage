# Test setup
mock_msa <- list(
  unaligned = Biostrings::AAStringSetList(
    H3 = Biostrings::AAStringSet(c(H31_HUMAN = "MARTKQTARKSTGGKAPRKQ")),
    H4 = Biostrings::AAStringSet(c(H4_HUMAN = "SGRGKGGKGLGKGGAKRHRK"))
  ),
  msa = Biostrings::AAStringSetList(
    H3 = Biostrings::AAStringSet(c(H31_HUMAN = "M-ARTKQTARKSTGGKAPRKQ", H3C_HUMAN = "")),
    H4 = Biostrings::AAStringSet(c(H4_HUMAN = "SGRGK-GGKGLGKGGAKRHRK"))
  ),
  msa_ref = Biostrings::AAStringSetList(
    H3 = Biostrings::AAStringSet(c(H3_REF = "M-ARKSTGGKAPRTKQTARKQ")),
    H4 = Biostrings::AAStringSet(c(H4_REF = "SGRHRGK-GGKGLGKGGAKRK"))
  )
)
#fmt: skip
mock_se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(intensity = matrix(runif(49), nrow = 7, dimnames = list(c(as.character(1:7)), LETTERS[1:7]))),
  rowData = S4Vectors::DataFrame(
    feature_number = 1:7,
    sequence = c("KQTARK", "KSTGGK", "KQTARK", "PEPTIDE", "GKGGK", "KQTARK", "YQKSTELLIR"),
    mods = c("[2] (Q) Ac|[5] (R) Me2", "[1] (K) Propionyl|[6] (K) Me3", "[2] (Q) Me3", NA_character_, "[2] (K) Ac", "[2] (Q) Dimethyl", NA_character_),
    charge = c(2L, 3L, 2L, 2L, 2L, 3L, 2L),
    histone = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
    histone_family = c("H3", "H3", "H3", NA_character_, "H4", "H3", "H3"),
    histone_group = c("H31_HUMAN", "H31_HUMAN", "H31_HUMAN", NA_character_, "H4_HUMAN", "H31_HUMAN", "H31_HUMAN/H3C_HUMAN"),
    start_index = I(list(4L, 9L, 4L, NA_integer_, 4L, 4L, c(54L, 53L))),
    row.names = as.character(1:7)
  )
)
mock_qf <- QFeatures::QFeatures(
  list(assay1 = mock_se, assay2 = mock_se),
  colData = S4Vectors::DataFrame(Var1 = rnorm(7), Var2 = LETTERS[1:7], row.names = LETTERS[1:7])
) |>
  QFeatures::addAssayLinkOneToOne("assay1", "assay2")

# helper functions
test_that(".parse_mods works with progenesis format", {
  mods <- c("[2] Ac (Q)|[5] Me2 (R)", "[1] Propionyl (K)")
  expect_message(res <- .parse_mods(mods, "progenesis"), "Unique modifications found: Ac; Me2; Propionyl")
  expect_equal(res$locs, list(c("2", "5"), "1"))
  expect_equal(res$mods, list(c("Ac", "Me2"), "Propionyl"))
})

test_that(".parse_mods works with progenesis_sw format", {
  mods <- c("[2] (Q) Ac|[5] (R) Me2", "[1] (K) Propionyl")
  expect_message(res <- .parse_mods(mods, "progenesis_sw"), "Unique modifications found: Ac; Me2; Propionyl")
  expect_equal(res$locs, list(c("2", "5"), "1"))
  expect_equal(res$mods, list(c("Ac", "Me2"), "Propionyl"))
})

test_that(".rename_mods works", {
  mods <- list(c("Ac", "Me2"), NA_character_, "Propionyl")
  rename_map <- list("Ac" ~ "Acetyl", "Me2" ~ "Dimethyl")
  expect_message(res <- .rename_mods(mods, rename_map), "Renaming mods: \"Ac\" ~ \"Acetyl\"; \"Me2\" ~ \"Dimethyl\"")
  expect_equal(res, list(c("Acetyl", "Dimethyl"), NA_character_, "Propionyl"))
  expect_identical(.rename_mods(mods, NULL), mods)
})

test_that(".aa_from_position works", {
  locs <- list(c("N-term", "5", "C-term"), NA_character_)
  seqs <- c("KQTARK", "KQTARK")
  res <- .aa_from_position(locs, seqs)
  expect_equal(res, list(c("N-term", "R", "C-term"), NA_character_))
})

test_that(".add_unmods works", {
  locs <- list(c("2"), c("2"), NA_character_, c("2"))
  mods <- list(c("Ac"), c("Ac"), NA_character_, c("Ac"))
  seqs <- c("KQTARK", "AKAK", "AKAK", "PEPTIDE")
  is_histone <- c(TRUE, TRUE, TRUE, FALSE)
  res <- .add_unmods(locs, mods, seqs, "K", is_histone)
  expect_equal(res$loc, list(c("1", "2", "6"), c("2", "4"), c("2", "4"), c("2")))
  expect_equal(res$mod, list(c("Unmod", "Ac", "Unmod"), c("Ac", "Unmod"), c("Unmod", "Unmod"), c("Ac")))
})

test_that(".process_modifications works", {
  locs <- list(c("2"), c("2", "4"), NA_character_, c("2"))
  mods <- list(c("Ac"), c("Propionyl", "Propionyl"), NA_character_, c("Ac"))
  seqs <- c("KQTARK", "AKAK", "AKAK", "PEPTIDE")
  strip_mods <- list(Propionyl = c("K", "N-term"))
  is_histone <- c(TRUE, TRUE, TRUE, FALSE)
  expect_message(
    expect_message(
      res <- .process_modifications(locs, mods, seqs, strip_mods, "K", is_histone),
      "Stripping mods: Propionyl on K, N-term"
    ),
    "Adding Unmods to: K"
  )
  expect_equal(res$loc, list(c("1", "2", "6"), c("2", "4"), c("2", "4"), c("2")))
  expect_equal(res$mod, list(c("Unmod", "Ac", "Unmod"), c("Unmod", "Unmod"), c("Unmod", "Unmod"), c("Ac")))
  expect_equal(res$aa, list(c("K", "Q", "K"), c("K", "K"), c("K", "K"), c("E")))
})

test_that(".create_mod_string works", {
  aa <- list(c("K", "T", "R"), NA_character_, c("K", "T", "R"))
  loc <- list(c("4", "5", "8"), NA_character_, c("4", "5", "8"))
  mod <- list(c("Unmod", "Ac", "Me2"), NA_character_, c("Unmod", "Ac", "Me2"))
  filter <- c(FALSE, FALSE, TRUE)
  expect_equal(.create_mod_string(aa, loc, mod, filter), c("K|4|Unmod;T|5|Ac;R|8|Me2", NA_character_, NA_character_))
})

# .map_pep_to_var --------------------------------------------------------------
# note the -1 offset due to initiator M getting index 0

test_that(".map_pep_to_var correctly maps a single numeric mod to one variant", {
  locs <- "5"
  starts <- 27L
  seqs <- "PEPTIDEK"

  expected <- list(matrix(31L))
  expect_equal(.map_pep_to_var(locs, starts, seqs), expected)
})

test_that(".map_pep_to_var handles multiple numeric mods on one variant", {
  locs <- list(c("2", "8"))
  starts <- 10L
  seqs <- "PEPTIDEKGG"

  expected <- list(matrix(c(11L, 17L), nrow = 1))
  expect_equal(.map_pep_to_var(locs, starts, seqs), expected)
})

test_that(".map_pep_to_var handles mapping to multiple variants", {
  locs <- "3"
  starts <- list(c(27L, 10L))
  seqs <- "PEPTIDEK"

  expected <- list(matrix(c(29L, 12L), ncol = 1))
  expect_equal(.map_pep_to_var(locs, starts, seqs), expected)
})

test_that(".map_pep_to_var handles multiple mods on multiple variants", {
  locs <- list(c("3", "6"))
  starts <- list(c(27L, 10L))
  seqs <- "PEPTIDEK"

  expected <- list(matrix(c(29L, 12L, 32L, 15L), nrow = 2))
  expect_equal(.map_pep_to_var(locs, starts, seqs), expected)
})

test_that(".map_pep_to_var correctly maps an N-terminal mod", {
  locs <- "N-term"
  starts <- list(c(27L, 10L))
  seqs <- "PEPTIDEK"

  expected <- list(matrix(c(27L, 10L), ncol = 1))
  expect_equal(.map_pep_to_var(locs, starts, seqs), expected)
})

test_that(".map_pep_to_var correctly maps a C-terminal mod", {
  locs <- "C-term"
  starts <- list(c(27L, 10L))
  seqs <- "PEPTIDEK"

  expected <- list(matrix(c(34L, 17L), ncol = 1))
  expect_equal(.map_pep_to_var(locs, starts, seqs), expected)
})

test_that(".map_pep_to_var handles mixed (N-term, numeric, C-term) mods", {
  locs <- list(c("N-term", "5", "C-term"))
  starts <- list(c(27L, 10L))
  seqs <- "PEPTIDEK"

  expected <- list(matrix(c(27L, 10L, 31L, 14L, 34L, 17L), nrow = 2))
  expect_equal(.map_pep_to_var(locs, starts, seqs), expected)
})

test_that(".map_pep_to_var returns NA when mods are NA", {
  locs <- NA_character_
  starts <- 27L
  seqs <- "PEPTIDEK"

  expect_equal(.map_pep_to_var(locs, starts, seqs), list(NA_integer_))
})

test_that(".map_pep_to_var returns NA when start indices are NA (i.e., not a histone)", {
  locs <- "5"
  starts <- NA_integer_
  seqs <- "PEPTIDEK"

  expect_equal(.map_pep_to_var(locs, starts, seqs), list(NA_integer_))
})

test_that(".map_pep_to_var works vectorized", {
  locs <- list("5", "5", NA_character_, c("N-term", "3"))
  starts <- list(27L, NA_integer_, 27L, c(10L, 15L))
  seqs <- list("PEPTIDEK", "AAAA", "BBBB", "GGG")

  expected <- list(
    matrix(31L),
    NA_integer_,
    NA_integer_,
    matrix(c(10L, 15L, 12L, 17L), nrow = 2)
  )
  expect_equal(.map_pep_to_var(locs, starts, seqs), expected)
})

test_that(".map_pep_to_var works with empty inputs", {
  expect_equal(.map_pep_to_var(list(), list(), list()), list())
})

test_that(".map_pep_to_var handles mod at position 1 equivalent to protein N-term", {
  locs <- "1"
  starts <- 0L
  seqs <- "PEPTIDEK"

  expect_equal(
    .map_pep_to_var(locs, starts, seqs),
    list(matrix(0))
  )
})

# .mapper_from_msa -------------------------------------------------------------

test_that(".mapper_from_msa works with a typical MSA", {
  msa <- Biostrings::AAStringSet(c(H1 = "A-C-E", H2 = "-B-DF"))
  expected <- list(
    H1 = c(0L, 2L, 4L),
    H2 = c(1L, 3L, 4L)
  )
  expect_equal(.mapper_from_msa(msa), expected)
})

test_that(".mapper_from_msa handles sequences without gaps", {
  msa <- Biostrings::AAStringSet(c(H1 = "ABCDE"))
  expected <- list(H1 = 0:4)
  expect_equal(.mapper_from_msa(msa), expected)
})

test_that(".mapper_from_msa handles sequences with only gaps", {
  msa <- Biostrings::AAStringSet(c(H1 = "ABC", H2 = "---"))
  expected <- list(H1 = 0:2, H2 = integer(0))
  expect_equal(.mapper_from_msa(msa), expected)
})

test_that(".mapper_from_msa handles leading and trailing gaps", {
  msa <- Biostrings::AAStringSet(c(H1 = "--ABC--", H2 = "---DEF", H3 = "DEF----"))
  expected <- list(
    H1 = c(2L, 3L, 4L),
    H2 = c(3L, 4L, 5L),
    H3 = c(0L, 1L, 2L)
  )
  expect_equal(.mapper_from_msa(msa), expected)
})

test_that(".mapper_from_msa preserves names", {
  msa <- Biostrings::AAStringSet(c(ProtX = "A-C", ProtY = "-BC"))
  result <- .mapper_from_msa(msa)
  expect_named(result, c("ProtX", "ProtY"))
})

test_that(".mapper_from_msa handles an empty AAStringSet", {
  msa <- Biostrings::AAStringSet(character(0))
  expect_equal(.mapper_from_msa(msa), list())
})

# .map_var_to_msa --------------------------------------------------------------

test_that(".map_var_to_msa single variant, single mod mapping works", {
  rd <- S4Vectors::DataFrame(
    feature_number = 1L,
    sequence = "PEPTIDEK",
    histone_family = "H3",
    histone_group = "H3_V1"
  )
  locs <- list(matrix(2))
  msa_mappers <- list(
    H3 = list(H3_V1 = c(0L, 1L, 2L, 4L, 5L, 6L, 7L, 8L, 9L, 10L), H3_V2 = c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L)),
    H4 = list(H4_V1 = c(0L, 1L, 2L, 3L))
  )

  expected <- list(list("2"))
  expect_equal(.map_var_to_msa(locs, rd, msa_mappers), expected)
})

test_that(".map_var_to_msa single variant, multiple mods mapping works", {
  rd <- S4Vectors::DataFrame(
    feature_number = 1L,
    sequence = "PEPTIDEK",
    histone_family = "H3",
    histone_group = "H3_V1"
  )
  locs <- list(matrix(c(2, 5), ncol = 2))
  msa_mappers <- list(
    H3 = list(H3_V1 = c(0L, 1L, 2L, 4L, 5L, 6L, 7L, 8L, 9L, 10L), H3_V2 = c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L)),
    H4 = list(H4_V1 = c(0L, 1L, 2L, 3L))
  )

  expected <- list(list("2", "6"))
  expect_equal(.map_var_to_msa(locs, rd, msa_mappers), expected)
})

test_that(".map_var_to_msa concordant mapping across multiple variants works", {
  rd <- S4Vectors::DataFrame(
    feature_number = 2L,
    sequence = "PEPTIDEK",
    histone_family = "H3",
    histone_group = "H3_V1/H3_V2"
  )
  locs <- list(matrix(c(1, 1, 2, 2), nrow = 2))
  msa_mappers <- list(
    H3 = list(H3_V1 = c(0L, 1L, 2L, 4L, 5L, 6L, 7L, 8L, 9L, 10L), H3_V2 = c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L)),
    H4 = list(H4_V1 = c(0L, 1L, 2L, 3L))
  )

  expected <- list(list("1", "2"))
  expect_equal(.map_var_to_msa(locs, rd, msa_mappers), expected)
  expect_no_message(.map_var_to_msa(locs, rd, msa_mappers))
})

test_that(".map_var_to_msa discordant mapping is handled and messaged", {
  rd <- S4Vectors::DataFrame(
    feature_number = 3L,
    sequence = "PEPTIDEK",
    histone_family = "H3",
    histone_group = "H3_V1/H3_V2"
  )
  locs <- list(matrix(c(3, 3, 5, 5), nrow = 2))
  msa_mappers <- list(
    H3 = list(H3_V1 = c(0L, 1L, 2L, 4L, 5L, 6L, 7L, 8L, 9L, 10L), H3_V2 = c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L)),
    H4 = list(H4_V1 = c(0L, 1L, 2L, 3L))
  )

  expected <- list(list("3/4", "5/6"))
  expect_message(
    expect_message(
      result <- .map_var_to_msa(locs, rd, msa_mappers),
      "Modified sequences with ambiguous location post-alignment:\nFeat\tSequence"
    ),
    "3\tPEPTIDEK"
  )
  expect_equal(result, expected)
})

test_that(".map_var_to_msa works vectorized", {
  rd <- S4Vectors::DataFrame(
    feature_number = c(1L, 2L),
    sequence = c("PEPTIDEK", "PEPTIDEL"),
    histone_family = c("H3", "H3"),
    histone_group = c("H3_V1", "H3_V2")
  )
  locs <- list(matrix(2), matrix(2))
  msa_mappers <- list(
    H3 = list(H3_V1 = c(0L, 1L, 2L, 4L, 5L, 6L, 7L, 8L, 9L, 10L), H3_V2 = c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L)),
    H4 = list(H4_V1 = c(0L, 1L, 2L, 3L))
  )

  expected <- list(list("2"), list("2"))
  expect_equal(.map_var_to_msa(locs, rd, msa_mappers), expected)
})

test_that(".map_var_to_msa NA inputs are handled gracefully", {
  rd_na_locs <- S4Vectors::DataFrame(
    feature_number = 4L,
    sequence = "NONHIST",
    histone_family = NA_character_,
    histone_group = NA_character_
  )
  locs_na <- list(NA_integer_)
  msa_mappers <- list(
    H3 = list(H3_V1 = c(0L, 1L, 2L, 4L, 5L, 6L, 7L, 8L, 9L, 10L), H3_V2 = c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L)),
    H4 = list(H4_V1 = c(0L, 1L, 2L, 3L))
  )

  expected_na <- list(list(NA_character_))
  expect_equal(.map_var_to_msa(locs_na, rd_na_locs, msa_mappers), expected_na)
})

test_that(".map_var_to_msa function works with empty inputs", {
  rd <- S4Vectors::DataFrame(
    feature_number = integer(),
    sequence = character(),
    histone_family = character(),
    histone_group = character()
  )
  locs <- list()
  msa_mappers <- list()

  expect_equal(.map_var_to_msa(locs, rd, msa_mappers), list())
})

test_that(".map_var_to_msa error on missing family or variant in mapper", {
  locs <- list(matrix(2))
  msa_mappers <- list(
    H3 = list(H3_V1 = c(0L, 1L, 2L, 4L, 5L, 6L, 7L, 8L, 9L, 10L), H3_V2 = c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L)),
    H4 = list(H4_V1 = c(0L, 1L, 2L, 3L))
  )

  rd_bad_fam <- S4Vectors::DataFrame(
    feature_number = 1L,
    sequence = "PEPTIDEK",
    histone_family = "H1",
    histone_group = "H1_V1"
  )
  expect_error(.map_var_to_msa(locs, rd_bad_fam, msa_mappers))

  rd_bad_var <- S4Vectors::DataFrame(
    feature_number = 1L,
    sequence = "PEPTIDEK",
    histone_family = "H3",
    histone_group = "H3_V3"
  )
  expect_error(.map_var_to_msa(locs, rd_bad_var, msa_mappers))
})

# .fill_msa_gaps ---------------------------------------------------------------

test_that(".fill_msa_gaps creates correct reference map", {
  ref_mapper <- 0:4
  msa_ref <- Biostrings::AAString("ARTKL")
  expect_equal(.fill_msa_gaps(ref_mapper, msa_ref), as.character(0:4))
})

test_that(".fill_msa_gaps handles leading gaps", {
  ref_mapper <- 2:6
  msa_ref <- Biostrings::AAString("--ARTKL")
  expect_equal(.fill_msa_gaps(ref_mapper, msa_ref), c("-2", "-1", as.character(0:4)))
})

test_that(".fill_msa_gaps handles trailing gaps", {
  ref_mapper <- 0:4
  msa_ref <- Biostrings::AAString("ARTKL---")
  expect_equal(.fill_msa_gaps(ref_mapper, msa_ref), c(as.character(0:3), "4", "4.1", "4.2", "4.3"))
})

test_that(".fill_msa_gaps handles internal gaps", {
  ref_mapper <- c(0L, 2L, 3L, 6L)
  msa_ref <- Biostrings::AAString("A-RT--K")
  expect_equal(
    .fill_msa_gaps(ref_mapper, msa_ref),
    c("0", "0.1", "1", "2", "2.1", "2.2", "3")
  )
})

test_that(".fill_msa_gaps handles all gap types simultaneously", {
  ref_mapper <- c(2L, 4L, 7L)
  msa_ref <- Biostrings::AAString("--A-R--T--")
  expect_equal(
    .fill_msa_gaps(ref_mapper, msa_ref),
    c("-2", "-1", "0", "0.1", "1", "1.1", "1.2", "2", "2.1", "2.2")
  )
})

test_that(".fill_msa_gaps handles single amino acid with gaps", {
  ref_mapper <- c(1L)
  msa_ref <- Biostrings::AAString("-A-")
  expect_equal(.fill_msa_gaps(ref_mapper, msa_ref), c("-1", "0", "0.1"))
})

test_that(".fill_msa_gaps errors on empty ref_mapper", {
  ref_mapper <- integer(0)
  msa_ref <- Biostrings::AAString("---")
  expect_error(.fill_msa_gaps(ref_mapper, msa_ref))
})

# .map_msa_to_ref --------------------------------------------------------------

test_that(".map_msa_to_ref works as expected", {
  mock_ref_mappers <- list(
    H3 = c("-1", as.character(0:25), "25.1", as.character(26:99), "99.1"),
    H4 = as.character(0:80)
  )
  locs <- list(list("10", "20"), list("25", "27"), list("27"), list("25/27"))
  families <- c("H4", "H3", "H3", "H3")

  expected <- list(c("10", "20"), c("24", "25.1"), "25.1", "24/25.1")
  expect_equal(.map_msa_to_ref(locs, families, mock_ref_mappers), expected)
})

test_that(".map_msa_to_ref handles NA inputs correctly", {
  mock_ref_mappers <- list(H3 = as.character(0:100))

  locs_na_loc <- list(list(NA_character_))
  families_na_loc <- "H3"
  expect_equal(.map_msa_to_ref(locs_na_loc, families_na_loc, mock_ref_mappers), list(list(NA_character_)))

  locs_na_fam <- list(list("25"))
  families_na_fam <- NA_character_
  expect_equal(.map_msa_to_ref(locs_na_fam, families_na_fam, mock_ref_mappers), list(list(NA_character_)))
})

test_that(".map_msa_to_ref handles empty lists", {
  expect_equal(.map_msa_to_ref(list(), character(), list()), list())
})

# .create_proforma -------------------------------------------------------------

test_that(".create_proforma handles single internal modifications correctly", {
  sequences <- "PEPTIDEK"
  locs <- "4"
  mods <- "Acetyl"
  charges <- 2L
  expected <- "PEPT[Acetyl]IDEK/2"
  expect_equal(.create_proforma(sequences, list(locs), list(mods), charges), expected)
})

test_that(".create_proforma handles multiple internal modifications correctly", {
  sequences <- "PEPTIDEK"
  locs <- c("4", "8")
  mods <- c("Acetyl", "Me3")
  charges <- 3L
  expected <- "PEPT[Acetyl]IDEK[Me3]/3"
  expect_equal(.create_proforma(sequences, list(locs), list(mods), charges), expected)
})

test_that(".create_proforma handles N-terminal modifications correctly", {
  sequences <- "PEPTIDEK"
  locs <- "N-term"
  mods <- "Propionyl"
  charges <- 2L
  expected <- "[Propionyl]-PEPTIDEK/2"
  expect_equal(.create_proforma(sequences, list(locs), list(mods), charges), expected)
})

test_that(".create_proforma handles C-terminal modifications correctly", {
  sequences <- "PEPTIDEK"
  locs <- "C-term"
  mods <- "Amidation"
  charges <- 1L
  expected <- "PEPTIDEK-[Amidation]/1"
  expect_equal(.create_proforma(sequences, list(locs), list(mods), charges), expected)
})

test_that(".create_proforma handles mixed (N-term, internal, C-term) modifications", {
  sequences <- "PEPTIDEK"
  locs <- c("N-term", "4", "C-term")
  mods <- c("Ac", "Me2", "Amide")
  charges <- 4L
  expected <- "[Ac]-PEPT[Me2]IDEK-[Amide]/4"
  expect_equal(.create_proforma(sequences, list(locs), list(mods), charges), expected)
})

test_that(".create_proforma handles cases with no modifications", {
  sequences <- "PEPTIDEK"
  locs <- NA
  mods <- NA
  charges <- 2L
  expected <- "PEPTIDEK/2"
  expect_equal(.create_proforma(sequences, list(locs), list(mods), charges), expected)
})

test_that(".create_proforma works with vectorized inputs", {
  sequences <- c("SEQONE", "SEQTWO")
  locs <- list(c("2", "C-term"), NA_character_)
  mods <- list(c("ModA", "ModB"), NA_character_)
  charges <- c(2L, 3L)
  expected <- c("SE[ModA]QONE-[ModB]/2", "SEQTWO/3")
  expect_equal(.create_proforma(sequences, locs, mods, charges), expected)
})

test_that(".create_proforma handles empty inputs", {
  expect_equal(.create_proforma(character(0), list(), list(), integer(0)), character(0))
})

test_that(".create_proforma handles modification names with special characters", {
  sequences <- "PEPTIDE"
  locs <- list("3")
  mods <- list("weird-mod(1) ^*!")
  charges <- 2L
  expected <- "PEP[weird-mod(1) ^*!]TIDE/2"
  expect_equal(.create_proforma(sequences, list(locs), list(mods), charges), expected)
})

# ------------------------------------------------------------------------------

test_that("processMods.SummarizedExperiment default works", {
  res_se <- processMods(mock_se, mock_msa, mod_format = "progenesis_sw")
  rd <- rowData(res_se)
  expect_true(all(c("mods_pep", "mods_var", "mods_msa", "mods_ref", "precursor") %in% names(rd)))
  expect_equal(rd$mods_pep[1], "Q|2|Ac;R|5|Me2")
  expect_equal(rd$mods_pep[2], "K|6|Me3;K|9|Unmod")
  expect_equal(rd$mods_var[1], "Q|5|Ac;R|8|Me2")
  expect_equal(rd$mods_msa[2], "K|14|Me3;K|9|Unmod")
  expect_equal(rd$mods_ref[1], "Q|5|Ac;R|8.1|Me2")
  expect_equal(rd$precursor[1], "K[Ac]QTA[Me2]RK/2")
})

test_that("processMods with unmods=NULL works", {
  res_se <- processMods(mock_se, mock_msa, mod_format = "progenesis_sw", unmods = NULL)
  rd <- rowData(res_se)
  expect_equal(rd$mods_pep[2], "K|6|Me3")
})

test_that("processMods with strip_mods=NULL works", {
  res_se <- processMods(mock_se, mock_msa, mod_format = "progenesis_sw", strip_mods = NULL)
  rd <- rowData(res_se)
  expect_equal(rd$mods_pep[2], "K|1|Propionyl;K|6|Me3;K|9|Unmod")
})

test_that("processMods with rename_mods works", {
  res_se <- processMods(
    mock_se,
    mock_msa,
    mod_format = "progenesis_sw",
    rename_mods = list("Ac" ~ "Acetyl", "Dimethyl" ~ "Me2")
  )
  rd <- rowData(res_se)
  expect_equal(rd$mods_pep[1], "Q|2|Acetyl;R|5|Me2")
  expect_equal(rd$mods_pep[6], "Q|2|Me2;K|1|Unmod;K|6|Unmod")
})

test_that("processMods without msa_ref works", {
  msa_no_ref <- mock_msa
  msa_no_ref$msa_ref <- NULL
  res_se <- processMods(mock_se, msa_no_ref, mod_format = "progenesis_sw")
  expect_false("mods_ref" %in% names(rowData(res_se)))
})

test_that("processMods handles co-extracts and missing mods", {
  res_se <- processMods(mock_se, mock_msa, mod_format = "progenesis_sw")
  rd <- rowData(res_se)
  expect_true(is.na(rd$mods_pep[4]))
  expect_true(is.na(rd$mods_var[4]))
  expect_true(is.na(rd$mods_msa[4]))
  expect_true(is.na(rd$mods_ref[4]))
  expect_equal(rd$precursor[4], "PEPTIDE/2")
  expect_equal(rd$mods_pep[7], "K|1|Unmod;K|6|Unmod")
})

test_that("processMods.QFeatures works", {
  res_qf <- processMods(mock_qf, mock_msa, "peptides", mod_format = "progenesis_sw")
  rd <- rowData(res_qf[["peptides"]])
  expect_s4_class(res_qf, "QFeatures")
  expect_true("mods_pep" %in% names(rd))
  expect_equal(rd$mods_pep[1], "Q|2|Ac;R|5|Me2")
})

test_that("processMods throws expected errors", {
  bad_se <- mock_se
  rowData(bad_se)$sequence <- NULL
  expect_error(processMods(bad_se, mock_msa))

  bad_msa <- mock_msa
  bad_msa$msa <- NULL
  expect_error(processMods(mock_se, bad_msa))
})

test_that("processMods handles multiple start indices", {
  res_se <- processMods(mock_se, mock_msa, mod_format = "progenesis_sw")
  rd <- rowData(res_se)
  expect_equal(rd$mods_var[7], "K|9/11|Unmod;K|14/16|Unmod")
})
