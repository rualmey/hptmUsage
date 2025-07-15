# Peptide ion export from the Progenesis QIP archive at https://www.ebi.ac.uk/pride/archive/projects/PXD009910
benchmark <- readProgenesis("./data-raw/benchmark.csv")

usethis::use_data(benchmark, overwrite = TRUE)
