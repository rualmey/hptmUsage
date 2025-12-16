## Import data -----------------------------------------------------------------
# Peptide ion export from the Progenesis QIP archive at https://www.ebi.ac.uk/pride/archive/projects/PXD009910
hptm_benchmark <- readProgenesis("./data-raw/hptm_benchmark.csv")
# Bovine histone samples (not negative controls) are not relevant, so discard
hptm_benchmark <- hptm_benchmark[, hptm_benchmark$group != "BH"]

# Manually edit colData (alternative to replaceColData)
hptm_benchmark$group <- hptm_benchmark$group |>
  forcats::fct_relevel("0min", "15min", "30min", "1h", "2h", "4h", "8h", "0min_neg", "8h_neg")
hptm_benchmark$include[hptm_benchmark$group == "BH"] <- FALSE
hptm_benchmark$time <- hptm_benchmark$group |>
  dplyr::case_match(
    c("0min", "0min_neg") ~ 0,
    "15min" ~ 0.25,
    "30min" ~ 0.5,
    "1h" ~ 1,
    "2h" ~ 2,
    "4h" ~ 4,
    c("8h", "8h_neg") ~ 8
  )
hptm_benchmark$treated <- hptm_benchmark$group |>
  dplyr::case_match(
    c("0min_neg", "8h_neg") ~ FALSE,
    c("0min", "15min", "30min", "1h", "2h", "4h", "8h") ~ TRUE
  )

## Process modifications -------------------------------------------------------
# need to retrieve bovine histones instead of human
histones <- histonesFromUniprot(query = c("organism_id:9913", "reviewed:true")) |>
  alignHistones()
hptm_benchmark <- matchHistones(hptm_benchmark, histones$unaligned, 1) |>
  processMods(
    histones,
    1,
    mod_format = "progenesis_sw",
    # remove propionylation, i.e., the fixed mod from derivatization
    strip_mods = list(`[+56.03]` = c("K", "N-term")),
    # examples of renaming two important hPTMs
    rename_mods = list("[+42.01]" ~ "Ac", "[+42.05]" ~ "Me3")
  )

## Global normalization --------------------------------------------------------
# VSN requires the non-log transformed values as it glog transforms the data
hptm_benchmark <- zeroIsNA(hptm_benchmark, i = "precursorRaw") |>
  normalize(
    i = "precursorRaw",
    method = "vsn",
    name = "precursor",
    verbose = FALSE,
    # while the below decreases the number of features used, it is more robust
    lts.quantile = .5
  )

## Filtering -------------------------------------------------------------------
# drop contaminant features
hptm_benchmark <- tagContaminants(hptm_benchmark, i = "precursor") |>
  filterFeatures(~ !contaminant, keep = TRUE)
# remove high missing features
# at least 3 samples should have the precursor quantified + at least 50% across
rowData(hptm_benchmark[["precursor"]]) <- rowData(hptm_benchmark[["precursor"]]) |>
  cbind(nNA(hptm_benchmark, i = "precursor")$nNArows[c("nNA", "pNA")])
hptm_benchmark <- hptm_benchmark |>
  filterFeatures(
    VariableFilter("nNA", ncols(hptm_benchmark)[["precursor"]] - 3, "<="),
    i = "precursor",
    keep = TRUE
  ) |>
  filterFeatures(~ pNA < .5, i = "precursor", keep = TRUE)

## Coextracts ------------------------------------------------------------------
# select co-extracts and aggregate to the protein level
hptm_benchmark <- hptm_benchmark |>
  addAssay(
    subset(hptm_benchmark[["precursor"]], subset = !histone),
    "precursorCoextr"
  ) |>
  addAssayLink("precursor", "precursorCoextr", "feature_number", "feature_number") |>
  aggregateFeatures(
    i = "precursorCoextr",
    fcol = "protein",
    name = "co-extracts",
    fun = MsCoreUtils::robustSummary
  )

## Usage calculation -----------------------------------------------------------
hptm_benchmark <- hptm_benchmark |>
  addAssay(
    subset(hptm_benchmark[["precursor"]], subset = histone),
    "precursorHistone"
  ) |>
  addAssayLink("precursor", "precursorHistone", "feature_number", "feature_number") |>
  # we will calculate hPTM usage against the entire nucleosome
  calculateUsage(
    i = "precursorHistone",
    name = "ptmVariantAgnostic"
  ) |>
  # against individual histone families (which we expect to be highly similar to the above)
  calculateUsage(
    i = "precursorHistone",
    name = "ptmFamilyCorrected",
    usage_level = "histone_family"
  ) |>
  # and against histone variant groups (resolved without any parsimony, see hptmUsage::matchHistones())
  calculateUsage(
    i = "precursorHistone",
    name = "ptmVariantCorrected",
    usage_level = "histone_group"
  ) |>
  # finally, calculate histone variant usage
  calculateUsage(
    i = "precursorHistone",
    name = "variant",
    target = "variant"
  )

## Data exploration ------------------------------------------------------------
plot(hptm_benchmark)
# manual visualization is possible and a differential usage analysis according to the msqrob2PTM framework
# both have been wrapped into the function hptmUsage::hptmAnalysis() for convenient analysis

# store the preprocessed dataset with this package for use in examples
usethis::use_data(hptm_benchmark, overwrite = TRUE)
