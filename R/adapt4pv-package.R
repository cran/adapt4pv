#' Simulated data for the adapt4pv package
#'
#' Simple simulated data, used to demonstrate the features of
#' functions from adapt4cv package.
#'
#' @name data_PV
#' @aliases X Y data_PV
#' @format
#' \describe{
#' \item{X}{large sparse and binary matrix with 117160 rows and 300 columns.
#' Drug matrix exposure: each row corresponds to an individual and
#' each column corresponds to a drug.}
#' \item{Y}{large spase and binary vector of length 117160. Indicator of
#' the presence/absence of an adverse event for ech individual.
#' Only the first 30 drugs (out of the 300)
#' are associated with the outcome.}
#' }
#' @keywords datasets
#' @import Matrix
#' @import glmnet
#' @import speedglm
#' @import xgboost
#' @import doParallel
#' @import foreach
#' @importFrom stats pnorm fisher.test quantile median deviance binomial predict p.adjust
#' @importFrom speedglm speedglm.wfit
#' @importFrom parallel makeCluster stopCluster
#'
#'
#'
#' @examples
#'
#' data(ExamplePvData)
#'
#'
NULL


#' Adaptive approaches for signal detection in PharmacoVigilance
#'
#' This package fits adaptive lasso approaches in high dimension for signal detection in
#' pharmacovigilance.
#' In addition to classical implementations found in the litterature, we implemented
#' two approaches particularly appropriated to variable selections framework, which
#' is the one that stands in pharmacovigilance.
#' We also supply in this package signal detection approaches based on lasso regression
#' and propensity score in high dimension.
#'
#'
#' @name adapt4pv-package
#' @docType package
#' @author Emeline Courtois \cr Maintainer:
#' Emeline Courtois <emeline.courtois@@inserm.fr>

NULL
