#' wrap function for \code{cv.glmnet}
#'
#' Fit a first cross-validation on lasso regression and return selected covariates.
#' Can deal with very large sparse data matrices.
#' Intended for binary reponse only (option \code{family = "binomial"} is forced).
#' Depends on the \code{cv.glmnet} function from the package \code{glmnet}.
#'
#'
#' @param x Input matrix, of dimension nobs x nvars. Each row is an observation
#' vector. Can be in sparse matrix format (inherit from class
#' \code{"sparseMatrix"} as in package \code{Matrix}).
#' @param y Binary response variable, numeric.
#' @param nfolds Number of folds - default is 5. Although \code{nfolds} can be
#' as large as the sample size (leave-one-out CV), it is not recommended for
#' large datasets. Smallest value allowable is \code{nfolds=3}.
#' @param foldid An optional vector of values between 1 and \code{nfolds}
#' identifying what fold each observation is in.
#' If supplied, \code{nfolds} can be missing.
#' @param betaPos Should the covariates selected by the procedure be positively
#' associated with the outcome ? Default is \code{TRUE}.
#' @param \dots Other arguments that can be passed to \code{cv.glmnet}
#' from package \code{glmnet} other than \code{nfolds}, \code{foldid},
#' and \code{family}.
#'
#' @return An object with S3 class \code{"log.lasso"}.
#' \item{beta}{Numeric vector of regression coefficients in the lasso.
#' In \code{lasso_cv} function, the regression coefficients are PENALIZED.
#' Length equal to nvars.}
#' \item{selected_variables}{Character vector, names of variable(s) selected with the
#' lasso-cv approach.
#' If \code{betaPos = TRUE}, this set is the covariates with a positive regression
#' coefficient in \code{beta}.
#' Else this set is the covariates with a non null regression coefficient in
#' \code{beta}.
#' Covariates are ordering according to magnitude of their regression
#' coefficients absolute value.}
#' @examples
#'
#' set.seed(15)
#' drugs <- matrix(rbinom(100*20, 1, 0.2), nrow = 100, ncol = 20)
#' colnames(drugs) <- paste0("drugs",1:ncol(drugs))
#' ae <- rbinom(100, 1, 0.3)
#' lcv <- lasso_cv(x = drugs, y = ae, nfolds = 3)
#'
#'
#' @author Emeline Courtois \cr Maintainer: Emeline Courtois
#' \email{emeline.courtois@@inserm.fr}
#' @export lasso_cv


lasso_cv <- function(x, y, nfolds = 5, foldid = NULL, betaPos = TRUE, ...){

  cv <- cv.glmnet(x = x, y = y, foldid = foldid, family = "binomial", ...)
  tmp <- match(cv$lambda.min, cv$glmnet.fit$lambda) # which lambda along the sequence is lambda min
  beta <- cv$glmnet.fit$beta[,tmp]

  if(betaPos){
    var <- names(sort(beta[which(beta>0)], decreasing = TRUE)) # ordre drecroissant
  }else{
    var <- names(sort(abs(beta[which(beta!=0)]), decreasing = TRUE)) # ordre decroissant en terme de valeur absolue
  }

  # Result -----
  res <- list()
  res$beta <- beta
  res$selected_variables <- var

  class(res) = "log.lasso"

  return(res)

}


