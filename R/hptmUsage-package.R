#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importClassesFrom QFeatures QFeatures
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom dplyr across mutate
#' @importFrom methods as is
#' @importFrom rlang .data %||%
#' @importMethodsFrom SummarizedExperiment rowData rowData<- colData colData<-
## usethis namespace: end
setClassUnion("QFeatures_OR_SummarizedExperiment", c("QFeatures", "SummarizedExperiment"))
