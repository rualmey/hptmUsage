# Peptide ion export from the Progenesis QIP archive at https://www.ebi.ac.uk/pride/archive/projects/PXD009910
hPTM_benchmark <- readProgenesis("./data-raw/hPTM_benchmark.csv")

usethis::use_data(hPTM_benchmark, overwrite = TRUE)
