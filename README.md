
# hptmUsage

<!-- badges: start -->

[![Project Status: Inactive – The project has reached a stable, usable
state but is no longer being actively developed; support/maintenance
will be provided as time
allows.](https://www.repostatus.org/badges/latest/inactive.svg)](https://www.repostatus.org/#inactive)
[![License](https://img.shields.io/badge/license-GPLv3-blue)](https://www.gnu.org/licenses/gpl-3.0.en.html)
<!-- badges: end -->

Preprocessing, visualization, and analysis (e.g., differential usage) of
histone post-translational modifications (hPTMs). This package builds on
the ‘msqrob2PTM’ workflow to enable robust and performant analysis of
hPTM data.

## Installation

Install the development version of hptmUsage from GitHub with:

``` r
# Install using "pak", alternatively use `devtools::install_github()` or `renv::install()`
# install.packages("pak")
pak::pak("rualmey/hptmUsage")
```

This package requires the following software to be installed (and to be
found on PATH):

- [MAFFT](https://mafft.cbrc.jp/alignment/software/) for multiple
  sequence alignment

- [Quarto](https://quarto.org/docs/get-started/) for reporting of
  results

## Usage

The easiest way to use this package is through the
`hptmUsage::generateReport()` function. A generic workflow is shown
below:

``` r
library(hptmUsage)

pe <- readProgenesis("./all-peptide-ion-export.csv", generate_metadata = "./metadata.csv")

# Now manually edit the metadata CSV and add it back to the dataset
pe |>
  replaceColData("./metadata.csv") |>
  # And perform the analysis + render the report
  generateReport(
    "./out/",
    # The following contrast will be positive (right side of the volcano) if the usage is higher in group B.
    contrasts = list(factor = c("A vs B" = "groupB - groupA")),
    # Note, if no contrasts are specified, the analysis is still done up to visualizing the design matrix
    # which can help in defining the contrasts
    generate_usageplots = "significant"
    # Note, generating usage plots takes quite some time so only do this once satisfied with the analysis
  )

# You should now be able to find your HTML report under ./out/date_hptmUsage.html
```

The parameters for `hptmUsage::generateReport()`, see
`help("generateReport")`, should allow a lot of freedom in how the
analysis is performed. Additionally, the fully processed dataset is
silently returned from `hptmUsage::generateReport()`, allowing for easy
post-processing.

If, however, more advanced changes to the workflow are required, you can
take a look at [this basic code example](/data-raw/hptm_benchmark.R)
that processes a [benchmark
dataset](https://doi.org/10.1039/D1MO00201E), as well as [the Quarto
report template](/inst/quarto/report_template.qmd) and [the package
reference manual](/REFERENCE_MANUAL.pdf).
