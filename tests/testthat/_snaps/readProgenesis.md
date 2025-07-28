# non-existent file throws error

    Code
      readProgenesis("foo")
    Condition
      Error in `readProgenesis()`:
      ! file.exists(file) is not TRUE

# readProgenesis catches invalid quant

    Code
      readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"), quant = NA)
    Condition
      Error in `readProgenesis()`:
      ! quant %in% c("Raw abundance", "Normalized abundance", "Intensity") &&  .... is not TRUE

---

    Code
      readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"), quant = -
      1)
    Condition
      Error in `readProgenesis()`:
      ! quant %in% c("Raw abundance", "Normalized abundance", "Intensity") &&  .... is not TRUE

---

    Code
      readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"), quant = Inf)
    Condition
      Error in `readProgenesis()`:
      ! quant %in% c("Raw abundance", "Normalized abundance", "Intensity") &&  .... is not TRUE

---

    Code
      readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"), quant = NULL)
    Condition
      Error in `readProgenesis()`:
      ! quant %in% c("Raw abundance", "Normalized abundance", "Intensity") &&  .... is not TRUE

---

    Code
      readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"), quant = "Spectral count")
    Condition
      Error in `readProgenesis()`:
      ! quant %in% c("Raw abundance", "Normalized abundance", "Intensity") &&  .... is not TRUE

---

    Code
      readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"), quant = c(
        "Raw abundance", "Intensity"))
    Condition
      Error in `quant %in% c("Raw abundance", "Normalized abundance", "Intensity") && length(
          quant) == 1`:
      ! 'length = 2' in coercion to 'logical(1)'

# readProgenesis catches invalid generate_metadata

    Code
      readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"),
      generate_metadata = NA)
    Condition
      Error in `readProgenesis()`:
      ! isFALSE(generate_metadata) || isTRUE(generate_metadata) || (is.character(generate_metadata) &&  .... is not TRUE

---

    Code
      readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"),
      generate_metadata = -1)
    Condition
      Error in `readProgenesis()`:
      ! isFALSE(generate_metadata) || isTRUE(generate_metadata) || (is.character(generate_metadata) &&  .... is not TRUE

---

    Code
      readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"),
      generate_metadata = Inf)
    Condition
      Error in `readProgenesis()`:
      ! isFALSE(generate_metadata) || isTRUE(generate_metadata) || (is.character(generate_metadata) &&  .... is not TRUE

---

    Code
      readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"),
      generate_metadata = NULL)
    Condition
      Error in `readProgenesis()`:
      ! isFALSE(generate_metadata) || isTRUE(generate_metadata) || (is.character(generate_metadata) &&  .... is not TRUE

---

    Code
      readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"),
      generate_metadata = c("./path_foo", "./path_bar"))
    Condition
      Error in `readProgenesis()`:
      ! isFALSE(generate_metadata) || isTRUE(generate_metadata) || (is.character(generate_metadata) &&  .... is not TRUE

# readProgenesis catches invalid headers

    Code
      readProgenesis(testthat::test_path("fixtures", "all_ion_export_no_quant.csv"))
    Condition
      Error in `readProgenesis()`:
      ! quant %in% header[1, ] is not TRUE

---

    Code
      readProgenesis(testthat::test_path("fixtures", "all_ion_export_no_quant.csv"),
      quant = "Intensity")
    Condition
      Error in `readProgenesis()`:
      ! quant %in% header[1, ] is not TRUE

---

    Code
      readProgenesis(testthat::test_path("fixtures", "all_ion_export_no_feat.csv"))
    Condition
      Error in `readProgenesis()`:
      ! all(c("#", "Charge", "Protein", "Sequence", "Variable modifications ([position] description)") %in%  .... is not TRUE

# messages and warnings work as expected

    Code
      invisible(readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv")))
    Message
      Some features had a note:
      * Feature 6: This feature has a note attached to it!
      * Feature 38342: This feature lost its ID, for example due to feature editing without redoing tags
      
    Condition
      Warning in `readProgenesis()`:
      Some features have no assigned sequence, please verify. These will be dropped: 38342

# metadata is generated as expected

    Code
      out <- readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"),
      generate_metadata = temp_path)
    Message
      Some features had a note:
      * Feature 6: This feature has a note attached to it!
      * Feature 38342: This feature lost its ID, for example due to feature editing without redoing tags
      
    Condition
      Warning in `readProgenesis()`:
      Some features have no assigned sequence, please verify. These will be dropped: 38342
    Message
      Metadata written at TEMP_FILEPATH
      

# overwrite works as expected

    Code
      out <- readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"),
      generate_metadata = temp_path, overwrite_metadata = FALSE)
    Message
      Some features had a note:
      * Feature 6: This feature has a note attached to it!
      * Feature 38342: This feature lost its ID, for example due to feature editing without redoing tags
      
    Condition
      Warning in `readProgenesis()`:
      Some features have no assigned sequence, please verify. These will be dropped: 38342
      Warning in `readProgenesis()`:
      File already exists at "TEMP_FILEPATH", see argument `overwrite_metadata`.

---

    Code
      out <- readProgenesis(testthat::test_path("fixtures", "all_ion_export.csv"),
      generate_metadata = temp_path, overwrite_metadata = TRUE)
    Message
      Some features had a note:
      * Feature 6: This feature has a note attached to it!
      * Feature 38342: This feature lost its ID, for example due to feature editing without redoing tags
      
    Condition
      Warning in `readProgenesis()`:
      Some features have no assigned sequence, please verify. These will be dropped: 38342
    Message
      Overwriting metadata... Metadata written at TEMP_FILEPATH
      

