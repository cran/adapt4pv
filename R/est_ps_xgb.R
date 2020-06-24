#' propensity score estimation in high dimension using gradient tree boosting
#'
#' Estimate a propensity score to a given drug exposure (treatment)
#' with extreme gradient boosting.
#' Depends on \code{xgboost} package.
#' Internal function, not supposed to be used directly.
#'
#'
#' @param idx_expo Index of the column in \code{x} that corresponds to the
#' drug covariate for which we aim at estimating the PS.
#' @param x Input matrix, of dimension nobs x nvars. Each row is an
#' observation vector. Can be in sparse matrix format (inherit from class
#' \code{"sparseMatrix"} as in package \code{Matrix}).
#' @param parameters correspond to \code{params} in \code{xgb.train} function.
#' The complete list of parameters is available at
#' \url{http://xgboost.readthedocs.io/en/latest/parameter.html}.
#' Default is a list with \code{eta=0.1} (learning rate),
#' \code{max_depth = 6} (maximum length of a tree),
#' \code{objective = "binary:logistic"}
#' and  \code{nthread = 1} (number of threads for parallelization).
#' @param nrounds Maximum number of boosting iterations. Default is 200.
#' @param \dots Other arguments that can be passed to \code{xgb.train} function.
#'
#' @return An object with S3 class \code{"ps", "xgb"}.
#' \item{expo_name}{Character, name of the drug exposure for which the PS was
#' estimated. Correspond to \code{colnames(x)[idx_expo]}}.
#' \item{indicator_expo}{One-column Matrix object. Indicator of the drug
#' exposure for which the PS was estimated.
#' Defined by \code{x[, idx_expo]}.}.
#' \item{score_variables}{Character vector, names of covariates(s) used in
#' a at list one tree in the gradient tree boosting algorithm.
#' Obtained with \code{xgb.importance} function from \code{xgboost} package.}
#' \item{score}{One-column Matrix object, the estimated score.}
#' @examples
#'
#' set.seed(15)
#' drugs <- matrix(rbinom(100*20, 1, 0.2), nrow = 100, ncol = 20)
#' colnames(drugs) <- paste0("drugs",1:ncol(drugs))
#' ae <- rbinom(100, 1, 0.3)
#' psxgb2 <- est_ps_xgb(idx_expo = 2, x = drugs, nrounds = 100)
#' psxgb2$score_variables #selected variables to include in the PS model of drug_2
#'
#'
#' @author Emeline Courtois \cr Maintainer: Emeline Courtois
#' \email{emeline.courtois@@inserm.fr}
#' @export est_ps_xgb

est_ps_xgb <- function(idx_expo, x, parameters = list("eta" =0.1, max_depth = 6, objective = "binary:logistic", nthread = 1),
                       nrounds = 200, ...){


  # Formating data ----
  indic <- x[,idx_expo]
  data4xgb <- xgb.DMatrix(data = x[, -idx_expo], label = indic)


  # Gradient tree boosting ----
  bst <- xgb.train(params = parameters, data = data4xgb,
                   nrounds = nrounds, ...)

  # Variable present in a least one tree ----
  var <- xgb.importance(feature_names = NULL, model = bst)
  var <- var$Feature


  # PS estimation ----
  score <- predict(bst, x[,-idx_expo])
  score <- Matrix(score)

  indic <- Matrix(indic)

  # Results -----
  res <- list()
  res$expo_name <- colnames(x)[idx_expo]
  res$indicator_expo <- indic
  res$score_variables <- var
  res$score <- score

  class(res) = c("ps", "xgb")

  return(res)
}
