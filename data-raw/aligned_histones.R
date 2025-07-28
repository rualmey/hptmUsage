# "2025-07-28 08:32:03 UTC"
Sys.time()

# Default call will retrieve all human, reviewed histones
# and add them to the corresponding HistoneDB 2.0 curated
# MSA profile.
aligned_histones <- alignHistones()

usethis::use_data(aligned_histones, overwrite = TRUE)
