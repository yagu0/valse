#' @include generateXY.R
#' @include EMGLLF.R
#' @include EMGrank.R
#' @include initSmallEM.R
#' @include computeGridLambda.R
#' @include constructionModelesLassoMLE.R
#' @include constructionModelesLassoRank.R
#' @include filterModels.R
#' @include selectVariables.R
#' @include main.R
#' @include plot.R
#'
#' @useDynLib valse
#'
#' @importFrom parallel makeCluster parLapply stopCluster clusterExport
#' @importFrom MASS ginv
NULL
