#' Get Path to Data File(s)
#'
#' Convenience function for accessing data files (located in source at
#' `inst/extdata`) included with the hptmUsage package.
#'
#' @param file Package data file name. If `NULL` (default), all data files are
#'   listed.
#' @export
#' @examples
#' hptmUsageData()
#' hptmUsageData("H3.fasta")
hptmUsageData <- function(file = NULL) {
  if (is.null(file)) {
    dir(fs::path_package("extdata", package = "hptmUsage"))
  } else {
    fs::path_package("extdata", file, package = "hptmUsage")
  }
}

#' Example hPTM Data
#'
#' @description
#' A minimally preprocessed `QFeatures` of hPTM data for testing and
#' demonstration purposes. The samples consist of bottom-up histone extracts
#' from naive (condition_A) and primed (condition_B) stem cell cultures as
#' described in the reference below.
#'
#' @format
#' A QFeatures object with 1 assay named `precursorRaw`. This assay
#' contains 5663 features (peptide ions) and 10 samples from two conditions
#' ('condition_A' and 'condition_B').
#' \describe{
#'   \item{assay}{contains the raw peptide ion abundances}
#'   \item{rowData}{contains 21 variables, including (most relevant) "protein", "sequence", and "mods"}
#'   \item{colData}{contains 8 columns of sample metadata, including the "group" (biological factor), "ms_run"/"date_collected"/"prep_batch" (technical factors), and "include"/"outlier" (tags for analysis)}
#' }
#' @source <https://www.ebi.ac.uk/pride/archive/projects/PXD028162>
#' @references Zijlmans, D. W.; Talon, I.; Verhelst, S.; Bendall, A.;
#'   Van Nerum, K.; Javali, A.; Malcolm, A. A.; Van Knippenberg, S. S. F. A.;
#'   Biggins, L.; To, S. K.; Janiszewski, A.; Admiraal, D.; Knops, R.;
#'   Corthout, N.; Balaton, B. P.; Georgolopoulos, G.; Panda, A.; Bhanu, N. V.;
#'   Collier, A. J.; Fabian, C.; Allsop, R. N.; Chappell, J.; Pham, T. X. A.;
#'   Oberhuemer, M.; Ertekin, C.; Vanheer, L.; Athanasouli, P.; Lluis, F.;
#'   Deforce, D.; Jansen, J. H.; Garcia, B. A.; Vermeulen, M.; Rivron, N.;
#'   Dhaenens, M.; Marks, H.; Rugg-Gunn, P. J.; Pasque, V. Integrated
#'   Multi-Omics Reveal Polycomb Repressive Complex 2 Restricts Human
#'   Trophoblast Induction. Nat Cell Biol 2022, 24 (6), 858–871.
#'   <https://doi.org/10.1038/s41556-022-00932-w>.
#' @examples
#' ncbtoy |>
#'   show()
"ncbtoy"

#' hPTM Benchmark Data
#'
#' @description
#' A preprocessed `QFeatures` of hPTM data, functioning as a benchmark for
#' downstream analysis and statistics. The samples consist of commercial histone
#' extracts treated with a histone deacetylase (HDAC1) in a time-lapse design as
#' described in the reference below.
#'
#' @format
#' A QFeatures object containing 14 assays of 42 benchmark samples as described
#' in the reference.
#' \describe{
#'   \item{assay `precursorRaw`}{contains the raw (not log-transformed) peptide ion abundances of 677 features}
#'   \item{assay `precursor`}{contains the globally normalized peptide ion abundances of 578 filtered features}
#'   \item{assay `precursorCoextr`}{contains the globally normalized peptide ion abundances of 134 filtered features from co-extracts}
#'   \item{assay `co-extracts`}{contains the robustly summarized protein abundances of 45 co-extracts}
#'   \item{assay `precursorHistone`}{contains the globally normalized peptide ion abundances of 444 filtered features from histones}
#'   \item{assay `precursorHistoneNormHistone`}{contains the usage-normalized (usage defined against the nucleosome) peptide ion abundances of 444 histone features}
#'   \item{assay `precursorHistoneNormHistoneDeconv`}{contains the 1172 features deconvoluted from the assay `precursorHistoneNormHistone`}
#'   \item{assay `ptmVariantAgnostic`}{contains the robustly summarized hPTM usage of 347 histone groups}
#'   \item{assay `precursorHistoneNormHistone_family`}{contains the usage-normalized (usage defined against the corresponding histone family) peptide ion abundances of 444 histone features}
#'   \item{assay `precursorHistoneNormHistone_familyDeconv`}{contains the 1172 features deconvoluted from the assay `precursorHistoneNormHistone_family`}
#'   \item{assay `ptmFamilyCorrected`}{contains the robustly summarized hPTM usage of 347 histone groups}
#'   \item{assay `precursorHistoneNormHistone_group`}{contains the usage-normalized (usage defined against the corresponding histone group) peptide ion abundances of 444 histone features}
#'   \item{assay `precursorHistoneNormHistone_groupDeconv`}{contains the 1172 features deconvoluted from the assay `precursorHistoneNormHistone_group`}
#'   \item{assay `ptmVariantCorrected`}{contains the robustly summarized hPTM usage of 347 histone groups}
#'   \item{rowData}{contains up to 55 variables, possibly including:
#'     (I) variables exported from Progenesis QIP, e.g., "feature_number", "protein", "sequence", "mods", ...;
#'     (II) tags exported from Progenesis QIP, e.g., "histone_id_no_err_tol" to "no_protein_id_ra";
#'     (III) variables created during data processing using hptmUsage functions:
#'       "histone" = TRUE if the sequence matched any histone variant, else FALSE;
#'       "histone_family" = the histone family to which a sequence matched, NA is no match;
#'       "core_histone" = TRUE if the sequence matched to one of H2A/H2B/H3/H4, else FALSE;
#'       "histone_group" = all histone variants to which the sequence was matched;
#'       "start_index" = list of locations of the first amino acid within each matching variant;
#'       "end_index" = list of locations of the last amino acid within each matching variant;
#'       "ambiguous_match" = TRUE if the sequence has ambiguous histone family assignment or multiple possible positions per variant;
#'       "precursor" = precursor string in ProForma notation;
#'       "mods_pep" = clean mod string with location indexes on the peptide sequence;
#'       "mods_var" = clean mod string with location indexes on all matching histone variants;
#'       "mods_msa" = clean mod string with location indexes on all matching histone variants after MSA;
#'       "mods_ref" = clean mod string with "consensus" location index(es) on the histone family reference sequence after MSA;
#'       "contaminant" = TRUE if identified as a potential contaminant using hptmUsage::tagContaminants() else FALSE;
#'       "nNA" = the number of samples in which the precursor was missing;
#'       "pNA" = the fraction of samples in which the precursor was missing;
#'       ".n" = the number of precursors from which the abundance/usage was summarized;
#'       "hptm" = combination of histone variant + amino acid + location + PTM, e.g., H33_BOVIN#K|27|Me3, possibly in a "degenerated" hPTM group separated by ";".
#'   }
#'   \item{colData}{contains 7 columns of sample metadata, including "group", "time", and "treated" (experimental factors), and "include"/"outlier" (tags for analysis)}
#' }
#' @source <https://www.ebi.ac.uk/pride/archive/projects/PXD009910>
#' @references Clerck, L. D.; Willems, S.; Daled, S.; Puyvelde, B. V.;
#'   Verhelst, S.; Corveleyn, L.; Deforce, D.; Dhaenens, M. An Experimental
#'   Design to Extract More Information from MS-Based Histone Studies. Molecular
#'   Omics 2021, 17 (6), 929–938. <https://doi.org/10.1039/D1MO00201E>.
#' @examples
#' hptm_benchmark |>
#'   show()
"hptm_benchmark"

#' Aligned Human Histones
#'
#' @description
#' A default call to [alignHistones()] retrieved all human histone variants from
#' UniProt (reviewed entries only, "2025-07-28 08:32:03 UTC") and added them to
#' the predefined (HistoneDB 2.0 (1)) MSA profile.
#'
#' @format
#' A `list` of three `AAStringSetList` objects, each containing five elements
#' corresponding to the five histone families (H1, H2A, H2B, H3, H4):
#' \describe{
#'   \item{unaligned}{`AAStringSetList` containing the unaligned histone sequences}
#'   \item{msa}{`AAStringSetList` containing the aligned histone sequences after adding each sequence to the curated MSA profile of the corresponding family}
#'   \item{msa_ref}{`AAStringSetList` containing the aligned reference sequences for each family}
#' }
#' @references
#' (1) Draizen, E. J.; Shaytan, A. K.; Mariño-Ramírez, L.; Talbert, P. B.;
#' Landsman, D.; Panchenko, A. R. HistoneDB 2.0: A Histone Database with
#' Variants—an Integrated Resource to Explore Histones and Their Variants.
#' Database (Oxford) 2016, 2016, baw014. <https://doi.org/10.1093/database/baw014>.
#' @examples
#' # List all Histone H3 variants in the object
#' aligned_histones$unaligned$H3 |>
#'   names()
#' # Show the MSA of histone family H3
#' aligned_histones$msa$h3 |>
#'   print()
#' # Show the aligned reference sequences
#' aligned_histones$msa_ref |>
#'   purrr::walk(print)
"aligned_histones"
