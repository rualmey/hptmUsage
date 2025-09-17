
# hptmUsage

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![License](https://img.shields.io/badge/license-GPLv3-blue)](https://www.gnu.org/licenses/gpl-3.0.en.html)
<!-- badges: end -->

Preprocessing, visualization, and analysis (e.g., differential usage) of
histone post-translational modifications (hPTMs). This package builds on
the ‘msqrob2PTM’ workflow to enable robust and performant analysis of
hPTM data.

Histone PTMs are parsed against consensus backbones, followed by a
bespoke msqrob2PTM workflow for differential hPTM usage analysis. Input
is, as of now, a peptide ion data export from Progenesis (all identified
features).

## Installation

Install the development version of hptmUsage from GitHub with:

``` r
# Install using "pak", alternatively use `devtools::install_github()` or `renv::install()`
# install.packages("pak")
pak::pak("rualmey/hptmUsage")
```

## Usage

[Example code](/data-raw/hptm_benchmark.R) on the processing of a
[benchmark dataset](https://doi.org/10.1039/D1MO00201E) is available in
this repo.

Below, some basic functionality is shown:

``` r
# hptmUsage automatically loads QFeatures for its infrastructure
library(hptmUsage)
#> Loading required package: QFeatures
#> Loading required package: MultiAssayExperiment
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: 'MatrixGenerics'
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> Loading required package: GenomeInfoDb
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians
#> 
#> Attaching package: 'QFeatures'
#> The following object is masked from 'package:base':
#> 
#>     sweep

# Reading a Progenesis QIP all ion export
hptmUsageData("all_ion_export.csv") |>
  readProgenesis()
#> Some features had a note:
#> * Feature 6: This feature has a note attached to it!
#> * Feature 38342: This feature lost its ID, for example due to feature editing without redoing tags
#> Warning in readProgenesis(hptmUsageData("all_ion_export.csv")): Some features
#> have no assigned sequence, please verify. These will be dropped: 38342
#> An instance of class QFeatures containing 1 set(s):
#>  [1] precursorRaw: SummarizedExperiment with 4 rows and 10 columns

# Retrieving histones from UniProt and aligning their sequences using MAFFT
# Note: this requires MAFFT to be on the system PATH
histonesFromUniprot() |>
 alignHistones()
#> $unaligned
#> AAStringSetList of length 5
#> [["H1"]] H10_HUMAN=MTENSTSAPAAKPKRAKASKKSTDHPKYSDMIVAAIQAEKNRAGSSRQSIQKYIKSHY...
#> [["H2A"]] H2AY_HUMAN=MSSRGGKKKSTKTSRSAKAGVIFPVGRMLRYIKKGHPKYRIGVGAPVYMAAVLEYL...
#> [["H2B"]] H2B1K_HUMAN=MPEPAKSAPAPKKGSKKAVTKAQKKDGKKRKRSRKESYSVYVYKVLKQVHPDTGI...
#> [["H3"]] CENPA_HUMAN=MGPRRRSRKPEAPRRRSPSPTPTPGPSRRGPSLGASSHQHSRRRQGWLKEIRKLQK...
#> [["H4"]] H4_HUMAN=MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVL...
#> 
#> $msa
#> AAStringSetList of length 5
#> [["H1"]] H10_HUMAN=-----------------------------------MTENST-SAPAA-----------...
#> [["H2A"]] H2A1A_HUMAN=-MSGR-----GK-QGGKARAKSKSRSSRAGLQFPVGRIHRLLRKGNYAE-RIGAG...
#> [["H2B"]] H2B1A_HUMAN=--------MPEVSSKGAT---ISKK-----G-FKKAVV--------KTQKK-EGK...
#> [["H3"]] H33_HUMAN=MARTKQTARKSTGGKAPRKQLATKAAR----KSAPSTGGVKKPHRYRPGTVALREIRR...
#> [["H4"]] H4_HUMAN=MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVL...
#> 
#> $msa_ref
#> AAStringSetList of length 5
#> [["H1"]] ref_H11_HUMAN=-----------------------------------MSETVP-PAPAASAAP---...
#> [["H2A"]] ref_H2A1B_HUMAN=-MSGR-----GK-QGGKARAKAKTRSSRAGLQFPVGRVHRLLRKGNYSE-R...
#> [["H2B"]] ref_H2B1J_HUMAN=--------MPE-PAKSAP---APKK-----G-SKKAVT--------KAQKK...
#> [["H3"]] ref_H31_HUMAN=MARTKQTARKSTGGKAPRKQLATKAAR----KSAPATGGVKKPHRYRPGTVALR...
#> [["H4"]] ref_H4_HUMAN=MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEET...

# Match peptide sequences to the retrieved histones
ncbtoy |>
  matchHistones(aligned_histones$unaligned, 1)
#> ⠙ 0/5 ETA: ? | Matching sequences
#> ⠹ 2/5 ETA:  4s | Matching sequences
#> An instance of class QFeatures containing 1 set(s):
#>  [1] precursorRaw: SummarizedExperiment with 5472 rows and 10 columns

# And more to be found in the example code above
```
