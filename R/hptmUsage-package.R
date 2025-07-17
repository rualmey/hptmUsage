#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dplyr across
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @importMethodsFrom SummarizedExperiment rowData colData
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom QFeatures QFeatures
## usethis namespace: end
setClassUnion("QFeatures_OR_SummarizedExperiment", c("QFeatures", "SummarizedExperiment"))
