skip_if_no_mafft <- function() {
  if (!nzchar(Sys.which("mafft"))) {
    testthat::skip("MAFFT not found")
  }
}

new_mock_msa <- function() {
  list(
    unaligned = Biostrings::AAStringSetList(
      H3 = Biostrings::AAStringSet(c(
        H31_HUMAN = "MARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEACEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA",
        H3C_HUMAN = "MARTKQTARKSTGGKAPRKQLATKAARKSTPSTCGVKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFNTDLRFQSAAVGALQEASEAYLVGLLEDTNLCAIHAKRVTIMPKDIQLARRIRGERA",
        H3C_MISALIGNED = "MARTKQTARKSTGGKAPRKQLATKAARKSTPSTCGVK-PHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFNTDLRFQSAAVGALQEASEAYLVGLLEDTNLCAIHAKRVTIMPKDIQLARRIRGERA"
      )),
      H4 = Biostrings::AAStringSet(c(
        H4_HUMAN = "MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG"
      ))
    ),
    msa = Biostrings::AAStringSetList(
      H3 = Biostrings::AAStringSet(c(
        H31_HUMAN = "---MARTKQTARKSTGGKAPRKQLATKAAR----KSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDF--KTDLRFQSSAVMALQEACEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA-----------",
        H3C_HUMAN = "---MARTKQTARKSTGGKAPRKQLATKAAR----KSTPSTCGVK-PHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDF--NTDLRFQSAAVGALQEASEAYLVGLLEDTNLCAIHAKRVTIMPKDIQLARRIRGERA-----------",
        H3C_MISALIGNED = "-MARTKQTARKSTGGKAPRKQLATKAAR----KSTPSTCGVK-PHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDF--NTDLRFQSAAVGALQEASEAYLVGLLEDTNLCAIHAKRVTIMPKDIQLARRIRGERA-----------"
      )),
      H4 = Biostrings::AAStringSet(c(
        H4_HUMAN = "MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG"
      ))
    ),
    msa_ref = Biostrings::AAStringSetList(
      H3 = Biostrings::AAStringSet(c(
        H3_REF = "---MARTKQTARKSTGGKAPRKQLATKAAR----KSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDF--KTDLRFQSSAVMALQEACEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA-----------"
      )),
      H4 = Biostrings::AAStringSet(c(
        H4_REF = "MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG"
      ))
    )
  )
}

new_mock_se <- function(nrow = 6) {
  SummarizedExperiment::SummarizedExperiment(
    assays = list(
      intensity = matrix(
        runif(nrow * 7),
        nrow = nrow,
        dimnames = list(c(as.character(1:nrow)), LETTERS[1:7])
      )
    ),
    #fmt: skip
    rowData = S4Vectors::DataFrame(
      feature_number = 1:6,
      sequence = c("TKQTAR", "KSAPATGGVKKPHR", "YQKSTELLIR", "YQKSTELLIR", "GKGGKGLGKGGAKR", "PEPTIDEK"),
      mods = c("[2] Ac (K)", "[1] Ac (K)|[10] Propionyl (K)", "[4] Ph (S)", "[4] Ph (S)", NA_character_, "[8] Ac (K)"),
      charge = c(2L, 3L, 2L, 2L, 3L, 2L),
      histone = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
      histone_family = c("H3", "H3", "H3", "H3", "H4", NA_character_),
      histone_group = c("H31_HUMAN", "H31_HUMAN", "H31_HUMAN/H3C_HUMAN", "H31_HUMAN/H3C_MISALIGNED", "H4_HUMAN", NA_character_),
      start_index = I(list(3L, 27L, c(54L, 53L), c(54L, 53L), 4L, NA_integer_)),
      row.names = as.character(1:6)
    )[1:nrow, ]
  )
}

new_mock_qf <- function(nrow = 6) {
  mock_se <- new_mock_se(nrow)
  QFeatures::QFeatures(
    list(assay1 = mock_se, assay2 = mock_se),
    colData = S4Vectors::DataFrame(Var1 = rnorm(7), Var2 = LETTERS[1:7], row.names = LETTERS[1:7])
  ) |>
    QFeatures::addAssayLinkOneToOne("assay1", "assay2")
}
