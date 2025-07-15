# Peptide ion export from the Progenesis QIP archive at https://www.ebi.ac.uk/pride/archive/projects/PXD028162
# Only a subset of naive (condition_A) and primed (condition_B) samples were kept for the example
ncbtoy <- readProgenesis("./data-raw/ncbtoy.csv")

usethis::use_data(ncbtoy, overwrite = TRUE)
