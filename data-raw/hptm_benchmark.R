# Peptide ion export from the Progenesis QIP archive at https://www.ebi.ac.uk/pride/archive/projects/PXD009910
hptm_benchmark <- readProgenesis("./data-raw/hptm_benchmark.csv")

colData(hptm_benchmark)

usethis::use_data(hptm_benchmark, overwrite = TRUE)
