# contaminant fasta files from "Frankenfield, A. M.; Ni, J.; Ahmed, M.; Hao, L.
# Protein Contaminants Matter: Building Universal Protein Contaminant Libraries
# for DDA and DIA Proteomics. J. Proteome Res. 2022, 21 (9), 2104–2113.
# https://doi.org/10.1021/acs.jproteome.2c00145.", retrieved from
# https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics

print(Sys.time())
# "2025-07-10 14:40:29 UTC"

# helper function to filter histones
remove_histones <- function(database) {
  mask <- stringr::str_detect(names(database), r"(^(sp|tr)\|Cont_[\w\d]+\|H(1|2A|2B|3|4)[\w\d]*_)")
  database[!mask]
}

# retrieve all contaminant databases and filter out histones
universal <- Biostrings::readAAStringSet(
  "https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/raw/refs/heads/main/Universal%20protein%20contaminant%20FASTA/0602_Universal%20Contaminants.fasta"
) |>
  remove_histones()
cell_culture <- Biostrings::readAAStringSet(
  "https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/raw/refs/heads/main/Sample-type%20specific%20contaminant%20FASTA/0602_Cell%20Culture%20Contaminants.fasta"
) |>
  remove_histones()
mouse_tissue <- Biostrings::readAAStringSet(
  "https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/raw/refs/heads/main/Sample-type%20specific%20contaminant%20FASTA/Aug2022_Mouse%20Tissue%20Contaminants.fasta"
) |>
  remove_histones()
rat_tissue <- Biostrings::readAAStringSet(
  "https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/raw/refs/heads/main/Sample-type%20specific%20contaminant%20FASTA/Aug2022_Rat%20Tissue%20Contaminants.fasta"
) |>
  remove_histones()
neuron_culture <- Biostrings::readAAStringSet(
  "https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/raw/refs/heads/main/Sample-type%20specific%20contaminant%20FASTA/Dec2022_Neuron%20Culture%20Contaminants.fasta"
) |>
  remove_histones()
stem_cell_culture <- Biostrings::readAAStringSet(
  "https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/raw/refs/heads/main/Sample-type%20specific%20contaminant%20FASTA/Dec2022_Stem%20Cell%20Culture%20Contaminants.fasta"
) |>
  remove_histones()

contaminants <- Biostrings::AAStringSetList(
  "universal" = universal,
  "cell_culture" = cell_culture,
  "mouse_tissue" = mouse_tissue,
  "rat_tissue" = rat_tissue,
  "neuron_culture" = neuron_culture,
  "stem_cell_culture" = stem_cell_culture
)

usethis::use_data(contaminants, internal = TRUE, overwrite = TRUE)
