# helper functions -------------------------------------------------------------

## .parse_mods -----------------------------------------------------------------

test_that(".parse_mods works with 'progenesis' format", {
  mod_string <- c("[27] Me3 (K)", "[36] Ac (K)|[37] Me (K)")
  expected <- list(
    locs = list("27", c("36", "37")),
    mods = list("Me3", c("Ac", "Me"))
  )
  expect_equal(suppressMessages(.parse_mods(mod_string, "progenesis")), expected)
})

test_that(".parse_mods works with 'progenesis_sw' format", {
  mod_string <- c("[27] (K) Me3", "[36] (K) Ac|[37] (K) Me")
  expected <- list(
    locs = list("27", c("36", "37")),
    mods = list("Me3", c("Ac", "Me"))
  )
  expect_equal(suppressMessages(.parse_mods(mod_string, "progenesis_sw")), expected)
})

test_that(".parse_mods handles N-term and C-term locations", {
  expected <- list(
    locs = list(c("N-term", "C-term")),
    mods = list(c("Ac", "Me"))
  )

  expect_equal(suppressMessages(.parse_mods("[N-term] Ac (N-term)|[C-term] Me (C-term)", "progenesis")), expected)
  expect_equal(suppressMessages(.parse_mods("[N-term] (N-term) Ac|[C-term] (C-term) Me", "progenesis_sw")), expected)
})

test_that(".parse_mods handles empty or NA mod strings", {
  mod_string <- c(NA_character_, "")
  expected <- list(
    locs = list(NA_character_, NA_character_),
    mods = list(NA_character_, NA_character_)
  )
  expect_equal(suppressMessages(.parse_mods(mod_string, "progenesis")), expected)
  expect_equal(suppressMessages(.parse_mods(mod_string, "progenesis_sw")), expected)
})

test_that(".parse_mods handles malformed or partially malformed strings", {
  mod_string <- c("invalid string", "[27] Me3 (K)|invalid")
  expected <- list(
    locs = list(NA_character_, "27"),
    mods = list(NA_character_, "Me3")
  )
  expect_equal(suppressMessages(.parse_mods(mod_string, "progenesis")), expected)
})

test_that(".parse_mods handles empty input vector", {
  mod_string <- character(0)
  expected <- list(locs = list(), mods = list())
  expect_equal(suppressMessages(.parse_mods(mod_string, "progenesis")), expected)
})

test_that(".parse_mods prints correct message for unique mods", {
  mod_string <- c("[27] Me3 (K)|[36] Ac (K)", "[9] Ac (K)")
  expect_message(.parse_mods(mod_string, "progenesis"), "Unique modifications found: Ac; Me3")
})

test_that(".parse_mods correctly handles mods with special regex characters", {
  mod_string <- "[27] Some+Mod (K)"
  expected <- list(
    locs = list("27"),
    mods = list("Some+Mod")
  )
  expect_equal(suppressMessages(.parse_mods(mod_string, "progenesis")), expected)
})

## .rename_mods ----------------------------------------------------------------

test_that(".rename_mods returns original list when rename_map is NULL", {
  mods <- list(c("Ac", "Me2"))
  expect_no_message(.rename_mods(mods, NULL))
  expect_equal(.rename_mods(mods, NULL), mods)
})

test_that(".rename_mods renames modifications correctly", {
  mods <- list(c("Ac", "Me2", "Other"))
  rename_map <- list("Ac" ~ "Acetyl", "Me2" ~ "Dimethyl")

  expected <- list(c("Acetyl", "Dimethyl", "Other"))
  expect_message(
    renamed <- .rename_mods(mods, rename_map),
    "Renaming mods: \"Ac\" ~ \"Acetyl\"; \"Me2\" ~ \"Dimethyl\""
  )
  expect_equal(renamed, expected)
})

test_that(".rename_mods handles empty rename_map", {
  mods <- list(c("Ac", "Me2"))
  expect_error(.rename_mods(mods, list()))
})

test_that(".rename_mods works vectorized", {
  mods <- list(c("Ac"), c("Me2", "Ac"), character(0))
  rename_map <- list("Ac" ~ "Acetyl", "Me2" ~ "Dimethyl")
  expected <- list(c("Acetyl"), c("Dimethyl", "Acetyl"), character(0))
  expect_equal(suppressMessages(.rename_mods(mods, rename_map)), expected)
})

test_that(".rename_mods preserves NA values", {
  mods <- list(c("Ac", NA, "Me2"))
  rename_map <- list("Ac" ~ "Acetyl", "Me2" ~ "Dimethyl")
  expected <- list(c("Acetyl", NA, "Dimethyl"))
  expect_equal(suppressMessages(.rename_mods(mods, rename_map)), expected)
})

test_that(".rename_mods handles an empty list of mods", {
  mods <- list()
  rename_map <- list("Ac" ~ "Acetyl")
  expect_equal(suppressMessages(.rename_mods(mods, rename_map)), list())
})

## .process_modifications ------------------------------------------------------

test_that(".process_modifications returns original mods when no actions are specified", {
  locs <- list(c("4", "9"))
  mods <- list(c("Ac", "Me1"))
  sequences <- "PEPKTIDEK"

  expected <- list(
    loc = list(c("4", "9")),
    mod = list(c("Ac", "Me1")),
    aa = list(c("K", "K"))
  )
  expect_equal(
    .process_modifications(locs, mods, sequences, NULL, NULL, NULL),
    expected
  )
})

test_that(".process_modifications correctly strips specified modifications", {
  locs <- list(c("N-term", "4", "7", "8"), c("4"))
  mods <- list(c("Prop", "Ac", "Prop", "Prop"), c("Prop"))
  sequences <- c("KGGKGGKSK", "KGGK")
  strip_mods <- list(Prop = c("K", "N-term"))
  is_histone <- c(TRUE, TRUE)

  expect_message(
    result <- .process_modifications(locs, mods, sequences, strip_mods, NULL, is_histone),
    "Stripping mods: Prop on K, N-term"
  )
  expected <- list(
    loc = list(c("4", "8"), character(0)),
    mod = list(c("Ac", "Prop"), character(0)),
    aa = list(c("K", "S"), character(0))
  )
  expect_equal(result, expected)
})

test_that(".process_modifications correctly adds 'Unmod' to histone peptides", {
  locs <- list(c("4"), NA_character_)
  mods <- list(c("Ac"), NA_character_)
  sequences <- c("TGGKARTK", "TGGKARTK")
  unmods <- "K"
  is_histone <- c(TRUE, TRUE)

  expect_message(
    result <- .process_modifications(locs, mods, sequences, NULL, unmods, is_histone),
    "Adding Unmods to: K"
  )
  expected <- list(
    loc = list(c("4", "8"), c("4", "8")),
    mod = list(c("Ac", "Unmod"), c("Unmod", "Unmod")),
    aa = list(c("K", "K"), c("K", "K"))
  )
  expect_equal(result, expected)
})

test_that(".process_modifications does not add 'Unmod' to non-histone peptides", {
  locs <- list(NA_character_)
  mods <- list(NA_character_)
  sequences <- "TGGKARTK"
  unmods <- "K"
  is_histone <- FALSE

  expect_message(
    result <- .process_modifications(locs, mods, sequences, NULL, unmods, is_histone),
    "Adding Unmods to: K"
  )
  expected <- list(
    loc = list(NA_character_),
    mod = list(NA_character_),
    aa = list(NA_character_)
  )
  expect_equal(result, expected)
})

test_that(".process_modifications correctly adds 'Unmod' to N-term and C-term", {
  locs <- list(NA_character_)
  mods <- list(NA_character_)
  sequences <- "PEPTIDE"
  unmods <- c("N-term", "C-term")
  is_histone <- TRUE

  expect_message(
    result <- .process_modifications(locs, mods, sequences, NULL, unmods, is_histone),
    "Adding Unmods to: N-term, C-term"
  )
  expected <- list(
    loc = list(c("N-term", "C-term")),
    mod = list(c("Unmod", "Unmod")),
    aa = list(c("N-term", "C-term"))
  )
  expect_equal(result, expected)
})

test_that(".process_modifications strips mods before adding 'Unmod'", {
  locs <- list(c("4"))
  mods <- list(c("Prop"))
  sequences <- "TGGKARTK"
  strip_mods <- list(Prop = "K")
  unmods <- "K"
  is_histone <- TRUE

  expect_message(
    expect_message(
      result <- .process_modifications(locs, mods, sequences, strip_mods, unmods, is_histone),
      "Stripping mods: Prop on K"
    ),
    "Adding Unmods to: K"
  )
  expected <- list(
    loc = list(c("4", "8")),
    mod = list(c("Unmod", "Unmod")),
    aa = list(c("K", "K"))
  )
  expect_equal(result, expected)
})

test_that(".process_modifications handles empty inputs", {
  result <- .process_modifications(list(), list(), character(), NULL, NULL, logical())
  expected <- list(loc = list(), mod = list(), aa = list())
  expect_equal(result, expected)
})

## .aa_from_position -----------------------------------------------------------

test_that(".aa_from_position works with basic numeric locations", {
  locs <- list("3", "1")
  sequences <- list("ABCDE", "FGHIJ")
  expected <- list("C", "F")
  expect_equal(.aa_from_position(locs, sequences), expected)
})

test_that(".aa_from_position handles N-term and C-term", {
  locs <- list("N-term", "C-term")
  sequences <- list("ABCDE", "FGHIJ")
  expected <- list("N-term", "C-term")
  expect_equal(.aa_from_position(locs, sequences), expected)
})

test_that(".aa_from_position handles multiple locations per sequence", {
  locs <- list(c("1", "5", "7"))
  sequences <- list("PEPTIDE")
  expected <- list(c("P", "I", "E"))
  expect_equal(.aa_from_position(locs, sequences), expected)
})

test_that(".aa_from_position handles mixed location types", {
  locs <- list(c("N-term", "4", "C-term"))
  sequences <- list("SEQUENCE")
  expected <- list(c("N-term", "U", "C-term"))
  expect_equal(.aa_from_position(locs, sequences), expected)
})

test_that(".aa_from_position handles NA locations correctly", {
  locs <- list(NA_character_)
  sequences <- list("ANYSEQ")
  expected <- list(NA_character_)
  expect_equal(.aa_from_position(locs, sequences), expected)
})

test_that(".aa_from_position handles empty inputs", {
  expect_equal(.aa_from_position(list(), list()), list())
})

test_that(".aa_from_position handles empty sequence string", {
  locs <- list("1")
  sequences <- list("")
  expected <- list("")
  expect_equal(.aa_from_position(locs, sequences), expected)
})

## .add_unmods -----------------------------------------------------------------

test_that(".add_unmods adds Unmod to a single unmodified site", {
  res <- .add_unmods(
    locs = list(NA_character_),
    mods = list(NA_character_),
    sequences = "PEPTIDEK",
    unmods = "K",
    is_histone = TRUE
  )
  expect_equal(res$loc, list("8"))
  expect_equal(res$mod, list("Unmod"))
})

test_that(".add_unmods does not add Unmod to an already modified site", {
  res <- .add_unmods(
    locs = list("8"),
    mods = list("Ac"),
    sequences = "PEPTIDEK",
    unmods = "K",
    is_histone = TRUE
  )
  expect_equal(res$loc, list("8"))
  expect_equal(res$mod, list("Ac"))
})

test_that(".add_unmods adds Unmod to some sites but not modified sites", {
  res <- .add_unmods(
    locs = list("1"),
    mods = list("Ac"),
    sequences = "KGGK",
    unmods = "K",
    is_histone = TRUE
  )
  expect_equal(res$loc, list(c("1", "4")))
  expect_equal(res$mod, list(c("Ac", "Unmod")))
})

test_that(".add_unmods ignores non-histones (co-extracts)", {
  res <- .add_unmods(
    locs = list("1", NA_character_),
    mods = list("Ac", NA_character_),
    sequences = "KEPTIDEK",
    unmods = "K",
    is_histone = FALSE
  )
  expect_equal(res$loc, list("1", NA_character_))
  expect_equal(res$mod, list("Ac", NA_character_))
})

test_that(".add_unmods handles multiple sequences and mixed conditions", {
  res <- .add_unmods(
    locs = list("1", NA_character_, "8"),
    mods = list("Ac", NA_character_, "Me3"),
    sequences = c("KGGK", "SOMEPEPTIDE", "PEPTIDEK"),
    unmods = "K",
    is_histone = c(TRUE, FALSE, TRUE)
  )
  expect_equal(res$loc, list(c("1", "4"), NA_character_, "8"))
  expect_equal(res$mod, list(c("Ac", "Unmod"), NA_character_, "Me3"))
})

test_that(".add_unmods handles multiple unmod types", {
  res <- .add_unmods(
    locs = list(NA_character_),
    mods = list(NA_character_),
    sequences = "KGGRA",
    unmods = c("K", "R"),
    is_histone = TRUE
  )
  expect_equal(res$loc, list(c("1", "4")))
  expect_equal(res$mod, list(c("Unmod", "Unmod")))
})

test_that(".add_unmods handles terminal unmods", {
  res <- .add_unmods(
    locs = list(NA_character_, "N-term"),
    mods = list(NA_character_, "Ac"),
    sequences = c("KGGRA", "PEPTIDEK"),
    unmods = c("N-term", "K", "C-term"),
    is_histone = TRUE
  )
  expect_equal(res$loc, list(c("N-term", "1", "C-term"), c("N-term", "8", "C-term")))
  expect_equal(res$mod, list(c("Unmod", "Unmod", "Unmod"), c("Ac", "Unmod", "Unmod")))
})

test_that(".add_unmods correctly sorts with N-term and C-term mods", {
  res <- .add_unmods(
    locs = list(c("10", "N-term")),
    mods = list(c("Me3", "Ac")),
    sequences = "KPEPTIDERK",
    unmods = "K",
    is_histone = TRUE
  )
  expect_equal(res$loc, list(c("N-term", "1", "10")))
  expect_equal(res$mod, list(c("Ac", "Unmod", "Me3")))
})

test_that(".add_unmods handles empty loc/mod inputs", {
  res <- .add_unmods(
    locs = list(character(0)),
    mods = list(character(0)),
    sequences = "AK",
    unmods = "K",
    is_histone = TRUE
  )
  expect_equal(res$loc, list("2"))
  expect_equal(res$mod, list("Unmod"))
})

test_that(".add_unmods handles sequences with no potential unmod sites", {
  res <- .add_unmods(
    locs = list("1"),
    mods = list("Ac"),
    sequences = "PEPTIDE",
    unmods = "K",
    is_histone = TRUE
  )
  expect_equal(res$loc, list("1"))
  expect_equal(res$mod, list("Ac"))
})

test_that(".add_unmods handles case with no unmods to add", {
  res <- .add_unmods(
    locs = list(c("1", "4")),
    mods = list(c("Ac", "Me3")),
    sequences = "KGGK",
    unmods = "K",
    is_histone = TRUE
  )
  expect_equal(res$loc, list(c("1", "4")))
  expect_equal(res$mod, list(c("Ac", "Me3")))
})

## .create_proforma ------------------------------------------------------------

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

## .create_mod_string ----------------------------------------------------------

test_that(".create_mod_string handles a single modification", {
  aa <- "K"
  locs <- 27L
  mods <- "Me3"
  expect_equal(.create_mod_string(aa, locs, mods), "K|27|Me3")
})

test_that(".create_mod_string handles multiple modifications in one entry", {
  aa <- list(c("K", "K"))
  locs <- list(c(27L, 36L))
  mods <- list(c("Me3", "Ac"))
  expect_equal(.create_mod_string(aa, locs, mods), "K|27|Me3;K|36|Ac")
})

test_that(".create_mod_string handles multiple list entries (vectorizes)", {
  aa <- list(c("K", "K"), "R")
  locs <- list(c(27L, 36L), 17L)
  mods <- list(c("Me3", "Ac"), "Me2")
  expect_equal(.create_mod_string(aa, locs, mods), c("K|27|Me3;K|36|Ac", "R|17|Me2"))
})

test_that(".create_mod_string returns NA if there are no PTMs", {
  aa <- NA_character_
  locs <- NA_integer_
  mods <- NA_character_
  expect_equal(.create_mod_string(aa, locs, mods), NA_character_)
})

test_that(".create_mod_string handles special terminal locations", {
  aa <- list(c("N-term", "C-term"))
  locs <- list(c(1L, 8L))
  mods <- list(c("Propionyl", "Amidation"))
  expect_equal(.create_mod_string(aa, locs, mods), "N-term|1|Propionyl;C-term|8|Amidation")
})

test_that(".create_mod_string handles empty inputs", {
  expect_equal(.create_mod_string(list(), list(), list()), character(0))
})

test_that(".create_mod_string applies filter correctly", {
  aa <- list("K", "K")
  locs <- list(27L, 36L)
  mods <- list("Me3", "Ac")
  filter <- c(FALSE, TRUE)
  expect_equal(.create_mod_string(aa, locs, mods, filter = filter), c("K|27|Me3", NA_character_))
})

test_that(".create_mod_string combines filtering and NA mods", {
  aa <- list("K", "K", "R")
  locs <- list(27L, 36L, 17L)
  mods <- list("Me3", NA, "Me2")
  filter <- c(TRUE, FALSE, TRUE)
  expect_equal(.create_mod_string(aa, locs, mods, filter = filter), c(NA_character_, NA_character_, NA_character_))
})

## .map_pep_to_var -------------------------------------------------------------
# note the -1 offset due to initiator M getting index 0

test_that(".map_pep_to_var correctly maps a single numeric mod to one variant", {
  locs <- 5
  starts <- 27L
  seqs <- "PEPTIDEK"

  expected <- list(matrix(31L))
  expect_equal(.map_pep_to_var(locs, starts, seqs), expected)
})

test_that(".map_pep_to_var handles multiple numeric mods on one variant", {
  locs <- list(c(2, 8))
  starts <- 10L
  seqs <- "PEPTIDEKGG"

  expected <- list(matrix(c(11L, 17L), nrow = 1))
  expect_equal(.map_pep_to_var(locs, starts, seqs), expected)
})

test_that(".map_pep_to_var handles mapping to multiple variants", {
  locs <- 3
  starts <- list(c(27L, 10L))
  seqs <- "PEPTIDEK"

  expected <- list(matrix(c(29L, 12L), ncol = 1))
  expect_equal(.map_pep_to_var(locs, starts, seqs), expected)
})

test_that(".map_pep_to_var handles multiple mods on multiple variants", {
  locs <- list(c(3, 6))
  starts <- list(c(27L, 10L))
  seqs <- "PEPTIDEK"

  expected <- list(matrix(c(29L, 12L, 32L, 15L), nrow = 2))
  expect_equal(.map_pep_to_var(locs, starts, seqs), expected)
})

test_that(".map_pep_to_var returns NA when mods are NA", {
  locs <- NA_integer_
  starts <- 27L
  seqs <- "PEPTIDEK"

  expect_equal(.map_pep_to_var(locs, starts, seqs), list(NA_integer_))
})

test_that(".map_pep_to_var returns NA when start indices are NA (i.e., not a histone)", {
  locs <- 5
  starts <- NA_integer_
  seqs <- "PEPTIDEK"

  expect_equal(.map_pep_to_var(locs, starts, seqs), list(NA_integer_))
})

test_that(".map_pep_to_var works vectorized", {
  locs <- list(5, 5, NA_integer_, c(1, 3))
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
  locs <- 1
  starts <- 0L
  seqs <- "PEPTIDEK"

  expect_equal(
    .map_pep_to_var(locs, starts, seqs),
    list(matrix(0))
  )
})

## .mapper_from_msa ------------------------------------------------------------

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

## .map_var_to_msa -------------------------------------------------------------

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

## .fill_msa_gaps --------------------------------------------------------------

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

## .map_msa_to_ref -------------------------------------------------------------

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

# S4 Methods -------------------------------------------------------------------

## processMods.SummarizedExperiment --------------------------------------------

test_that("processMods works with default parameters", {
  mock_se <- new_mock_se()
  mock_msa <- new_mock_msa()

  expect_message(
    expect_message(
      expect_message(
        expect_message(
          expect_message(
            res <- processMods(mock_se, mock_msa),
            "Unique modifications found: Ac; Ph; Propionyl"
          ),
          "Stripping mods: Propionyl on K, N-term"
        ),
        "Adding Unmods to: K"
      ),
      "Modified sequences with ambiguous location post-alignment:\nFeat\tSequence"
    ),
    "4\tYQKSTELLIR"
  )
  # unchanged
  expect_s4_class(res, "SummarizedExperiment")
  expect_identical(colData(mock_se), colData(res))
  expect_identical(assays(mock_se, withDimnames = FALSE), assays(res, withDimnames = FALSE))
  # 5 new columns
  rd <- rowData(res, use.names = FALSE)
  expect_length(
    setdiff(
      names(rd),
      c(names(rowData(mock_se)), "mods_pep", "mods_var", "mods_msa", "mods_ref", "precursor")
    ),
    0
  )
  # of which the original columns are unchanged
  expect_identical(
    rd[c("feature_number", "sequence", "mods", "charge", "histone", "histone_family", "histone_group", "start_index")],
    rowData(mock_se, use.names = FALSE)
  )
  # new columns
  expect_equal(
    rd$mods_pep,
    c(
      "K|2|Ac",
      "K|1|Ac;K|10|Unmod;K|11|Unmod",
      "K|3|Unmod;S|4|Ph",
      "K|3|Unmod;S|4|Ph",
      "K|2|Unmod;K|5|Unmod;K|9|Unmod;K|13|Unmod",
      "K|8|Ac"
    )
  )
  expect_equal(
    rd$mods_var,
    c(
      "K|4|Ac",
      "K|27|Ac;K|36|Unmod;K|37|Unmod",
      "K|55/56|Unmod;S|56/57|Ph",
      "K|55/56|Unmod;S|56/57|Ph",
      "K|5|Unmod;K|8|Unmod;K|12|Unmod;K|16|Unmod",
      NA
    )
  )
  expect_equal(
    rd$mods_msa,
    c(
      "K|7|Ac",
      "K|34|Ac;K|43|Unmod;K|44|Unmod",
      "K|63|Unmod;S|64|Ph",
      "K|61/63|Unmod;S|62/64|Ph",
      "K|5|Unmod;K|8|Unmod;K|12|Unmod;K|16|Unmod",
      NA
    )
  )
  expect_equal(
    rd$mods_ref,
    c(
      "K|4|Ac",
      "K|27|Ac;K|36|Unmod;K|37|Unmod",
      "K|56|Unmod;S|57|Ph",
      "K|54/56|Unmod;S|55/57|Ph",
      "K|5|Unmod;K|8|Unmod;K|12|Unmod;K|16|Unmod",
      NA
    )
  )
  expect_equal(
    rd$precursor,
    c(
      "TK[Ac]QTAR/2",
      "K[Ac]SAPATGGVK[Unmod]K[Unmod]PHR/3",
      "YQK[Unmod]S[Ph]TELLIR/2",
      "YQK[Unmod]S[Ph]TELLIR/2",
      "GK[Unmod]GGK[Unmod]GLGK[Unmod]GGAK[Unmod]R/3",
      "PEPTIDEK[Ac]/2"
    )
  )
})

### argument checking ----------------------------------------------------------

test_that("processMods fails with missing rowData columns", {
  mock_msa <- new_mock_msa()
  mock_se <- new_mock_se()
  # drop mods column from rowData
  rowData(mock_se) <- rowData(mock_se)[-3]

  expect_error(processMods(mock_se, mock_msa))
})

test_that("processMods fails with invalid msa object", {
  mock_se <- new_mock_se()

  expect_error(processMods(mock_se, list(unaligned = new_mock_msa()$unaligned)))
  expect_error(processMods(mock_se, list(unaligned = "not_aassl", msa = "not_aassl")))
})

test_that("processMods fails with invalid argument types", {
  mock_se <- new_mock_se()
  mock_msa <- new_mock_msa()

  expect_error(processMods(mock_se, mock_msa, unmods = 123))
  expect_error(processMods(mock_se, mock_msa, strip_mods = list(c("K"))))
  expect_error(processMods(mock_se, mock_msa, rename_mods = list("Ac" = "Acetyl")))
})

### argument functionality -----------------------------------------------------

test_that("processMods correctly parses other formats", {
  mock_se <- new_mock_se(nrow = 1)
  mock_msa <- new_mock_msa()

  rowData(mock_se)$mods[[1]] <- "[2] (K) Ac"
  res <- suppressMessages(processMods(mock_se, mock_msa, mod_format = "progenesis_sw"))
  expect_equal(rowData(res)$mods_pep, "K|2|Ac")
})

test_that("processMods correctly adds unmods", {
  mock_se <- new_mock_se(nrow = 1)
  mock_msa <- new_mock_msa()

  res <- suppressMessages(processMods(mock_se, mock_msa, unmods = "R"))
  expect_equal(rowData(res)$mods_pep, "K|2|Ac;R|6|Unmod")
})

test_that("processMods correctly strips modifications", {
  mock_se <- new_mock_se(nrow = 1)
  mock_msa <- new_mock_msa()

  rowData(mock_se)$mods[[1]] <- "[N-term] Propionyl (N-term)|[2] Ac (K)|[6] Ac (R)"
  res <- suppressMessages(processMods(mock_se, mock_msa, strip_mods = list(Propionyl = "N-term", Ac = "K")))
  expect_equal(rowData(res)$mods_pep, "K|2|Unmod;R|6|Ac")
})

test_that("processMods correctly renames modifications", {
  mock_se <- new_mock_se()
  mock_msa <- new_mock_msa()

  res <- suppressMessages(processMods(mock_se, mock_msa, rename_mods = list("Ac" ~ "Acetyl")))
  expect_equal(
    rowData(res)$mods_pep,
    c(
      "K|2|Acetyl",
      "K|1|Acetyl;K|10|Unmod;K|11|Unmod",
      "K|3|Unmod;S|4|Ph",
      "K|3|Unmod;S|4|Ph",
      "K|2|Unmod;K|5|Unmod;K|9|Unmod;K|13|Unmod",
      "K|8|Acetyl"
    )
  )
})

test_that("processMods handles peptides with no modifications", {
  mock_se <- new_mock_se()[5, ]
  mock_msa <- new_mock_msa()

  res <- suppressMessages(processMods(mock_se, mock_msa, unmods = NULL))
  expect_true(
    is.na(rowData(res)$mods_pep) &&
      is.na(rowData(res)$mods_var) &&
      is.na(rowData(res)$mods_msa)
  )
  expect_equal(rowData(res)$precursor, "GKGGKGLGKGGAKR/3")
})

test_that("processMods handles non-histone peptides", {
  mock_se <- new_mock_se()[6, ]
  mock_msa <- new_mock_msa()

  res <- suppressMessages(processMods(mock_se, mock_msa))
  expect_equal(rowData(res)$mods_pep, "K|8|Ac")
  expect_true(is.na(rowData(res)$mods_var) && is.na(rowData(res)$mods_msa))
})

test_that("processMods does not add 'mods_ref' column when msa_ref is absent", {
  mock_se <- new_mock_se(nrow = 1)
  mock_msa <- new_mock_msa()[c("unaligned", "msa")]

  res <- suppressMessages(processMods(mock_se, mock_msa))
  expect_false("mods_ref" %in% colnames(rowData(res)))
})

test_that("strip_mods and unmods interact correctly", {
  mock_se <- new_mock_se()[3, ]
  mock_msa <- new_mock_msa()

  res <- suppressMessages(processMods(mock_se, mock_msa, strip_mods = list(Ph = "S"), unmods = "S"))
  expect_equal(rowData(res)$mods_pep, "S|4|Unmod")
})

## processMods.QFeatures -------------------------------------------------------

test_that("processMods.QFeatures works", {
  mock_msa <- new_mock_msa()
  mock_qf <- new_mock_qf()
  # repeat precursors are not allowed
  rowData(mock_qf[[1]])$mods[[4]] <- "[4] Su (S)"

  # warns that assayLinks are lost
  res_qf <- suppressWarnings(suppressMessages(processMods(mock_qf, mock_msa, "assay1")))
  # unchanged
  expect_s4_class(res_qf, "QFeatures")
  expect_equal(rowData(res_qf[[2]]), rowData(mock_qf[[2]]))
  # changed
  rd <- rowData(res_qf[["assay1"]], use.names = FALSE)
  expect_length(
    setdiff(
      names(rd),
      c(names(rowData(mock_qf[["assay1"]])), "mods_pep", "mods_var", "mods_msa", "mods_ref", "precursor")
    ),
    0
  )
  # of which the original columns are unchanged
  expect_identical(
    rd[c("feature_number", "sequence", "mods", "charge", "histone", "histone_family", "histone_group", "start_index")],
    rowData(mock_qf[["assay1"]], use.names = FALSE)
  )
  # new columns
  expect_equal(
    rd$mods_pep,
    c(
      "K|2|Ac",
      "K|1|Ac;K|10|Unmod;K|11|Unmod",
      "K|3|Unmod;S|4|Ph",
      "K|3|Unmod;S|4|Su",
      "K|2|Unmod;K|5|Unmod;K|9|Unmod;K|13|Unmod",
      "K|8|Ac"
    )
  )
  expect_equal(
    rd$mods_var,
    c(
      "K|4|Ac",
      "K|27|Ac;K|36|Unmod;K|37|Unmod",
      "K|55/56|Unmod;S|56/57|Ph",
      "K|55/56|Unmod;S|56/57|Su",
      "K|5|Unmod;K|8|Unmod;K|12|Unmod;K|16|Unmod",
      NA
    )
  )
  expect_equal(
    rd$mods_msa,
    c(
      "K|7|Ac",
      "K|34|Ac;K|43|Unmod;K|44|Unmod",
      "K|63|Unmod;S|64|Ph",
      "K|61/63|Unmod;S|62/64|Su",
      "K|5|Unmod;K|8|Unmod;K|12|Unmod;K|16|Unmod",
      NA
    )
  )
  expect_equal(
    rd$mods_ref,
    c(
      "K|4|Ac",
      "K|27|Ac;K|36|Unmod;K|37|Unmod",
      "K|56|Unmod;S|57|Ph",
      "K|54/56|Unmod;S|55/57|Su",
      "K|5|Unmod;K|8|Unmod;K|12|Unmod;K|16|Unmod",
      NA
    )
  )
  expect_equal(
    rd$precursor,
    c(
      "TK[Ac]QTAR/2",
      "K[Ac]SAPATGGVK[Unmod]K[Unmod]PHR/3",
      "YQK[Unmod]S[Ph]TELLIR/2",
      "YQK[Unmod]S[Su]TELLIR/2",
      "GK[Unmod]GGK[Unmod]GLGK[Unmod]GGAK[Unmod]R/3",
      "PEPTIDEK[Ac]/2"
    )
  )

  # res_qf_idx <- suppressMessages(processMods(mock_qf, mock_msa, 1))
  # expect_identical(res_qf, res_qf_idx)

  # res_qf_multi <- suppressMessages(processMods(mock_qf, mock_msa, 1:2))
  # expect_equal(rowData(res_qf_multi[[1]]), rowData(res_qf_multi[[2]]))
  # expect_equal(rowData(res_qf_multi[[1]]), rowData(res_qf[[1]]))
})
