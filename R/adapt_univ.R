#' fit an adaptive lasso with adaptive weights derived from univariate coefficients
#'
#' Compute odd-ratios between each covariate of \code{x} and \code{y} then derived
#' adaptive weights to incorporate in an adaptive lasso.
#' BIC or cross-validation could either be used for the adaptive lasso for variable selection.
#' Can deal with very large sparse data matrices.
#' Intended for binary reponse only (option \code{family = "binomial"} is forced).
#' Depends on the \code{glmnet} and \code{relax.glmnet} function from the package
#' \code{glmnet}.
#'
#' The adaptive weight for a given covariate i is defined by
#' \deqn{w_i = 1/|\beta^univ_i|^\gamma} where
#'  \eqn{\beta^univ_i = log(OR_i)}, with
#'  \eqn{OR_i} is the odd-ratio associated to covariate \eqn{i}
#'  with the outcome.
#'
#' @param x Input matrix, of dimension nobs x nvars. Each row is an observation
#' vector. Can be in sparse matrix format (inherit from class
#' \code{"sparseMatrix"} as in package \code{Matrix}).
#' @param y Binary response variable, numeric.
#' @param gamma Tunning parameter to defined the penalty weights. See details below.
#' Default is set to 1.
#' @param criterion Character, indicates which criterion is used with the
#' adaptive lasso for variable selection. Could be either "bic" or "cv".
#' Default is "bic"
#' @param maxp Used only if \code{criterion} = "bic", ignored if \code{criterion} = "cv".
#' A limit on how many relaxed coefficients are allowed. Default is 50, in \code{glmnet}
#' option default is 'n-3', where 'n' is the sample size.
#' @param path Used only if \code{criterion} = "bic", ignored if \code{criterion} = "cv".
#' Since \code{glmnet} does not do stepsize optimization, the Newton
#' algorithm can get stuck and not converge, especially with relaxed fits.
#' With \code{path=TRUE}, each relaxed fit on a particular set of variables
#' is computed pathwise using the original sequence of lambda values
#' (with a zero attached to the end). Default is \code{path=TRUE}.
#' @param nfolds Used only if \code{criterion} = "cv", ignored if \code{criterion} = "bic".
#' Number of folds - default is 5. Although \code{nfolds} can be
#' as large as the sample size (leave-one-out CV), it is not recommended for
#' large datasets. Smallest value allowable is \code{nfolds=3}.
#' @param foldid Used only if \code{criterion} = "cv", ignored if \code{criterion} = "bic".
#' An optional vector of values between 1 and \code{nfolds}
#' identifying what fold each observation is in. If supplied, \code{nfolds} can
#' be missing.
#' @param betaPos Should the covariates selected by the procedure be
#' positively associated with the outcome ? Default is \code{TRUE}.
#' @param \dots Other arguments that can be passed to \code{glmnet}
#' from package \code{glmnet} other than \code{penalty.factor},
#' \code{family}, \code{maxp} and \code{path}.
#'
#' @return An object with S3 class \code{"adaptive"}.
#' \item{aws}{Numeric vector of penalty weights derived from odds-ratios.
#' Length equal to nvars.}
#' \item{criterion}{Character, same as input. Could be either "bic" or "cv".}
#' \item{beta}{Numeric vector of regression coefficients in the adaptive lasso.
#' If \code{criterion} = "cv" the regression coefficients are PENALIZED, if
#' \code{criterion} = "bic" the regression coefficients are UNPENALIZED.
#' Length equal to nvars. Could be NA if adaptive weights are all equal to infinity.}
#' \item{selected_variables}{Character vector, names of variable(s) selected
#' with this adaptive approach.
#' If \code{betaPos = TRUE}, this set is the covariates with a positive regression
#' coefficient in \code{beta}.
#' Else this set is the covariates with a non null regression coefficient in \code{beta}.
#' If \code{criterion} = "bic", covariates are ordering according to magnitude of their regression
#' coefficients absolute value in the adaptive lasso.
#' If \code{criterion} = "bic", covariates are ordering according to the p-values (two-sided if \code{betaPos = FALSE} ,
#' one-sided if \code{betaPos = TRUE}) in the classical multiple logistic regression
#' model that minimzes the BIC in the adaptive lasso.}
#' @examples
#'
#' set.seed(15)
#' drugs <- matrix(rbinom(100*20, 1, 0.2), nrow = 100, ncol = 20)
#' colnames(drugs) <- paste0("drugs",1:ncol(drugs))
#' ae <- rbinom(100, 1, 0.3)
#' au <- adapt_univ(x = drugs, y = ae, criterion ="cv", nfolds = 3)
#'
#' @author Emeline Courtois \cr Maintainer: Emeline Courtois
#' \email{emeline.courtois@@inserm.fr}
#' @export adapt_univ


adapt_univ <- function(x, y, gamma = 1, criterion = "bic", maxp = 50, path = TRUE,
                       nfolds = 5, foldid = NULL, betaPos = TRUE, ...){

  # Odds ratio ----
  tab22 <- contingence_all(data_cases = x[which(y==1),], data_controls = x[which(y==0),])
  or <- lapply(tab22, fisher.test)
  or <- sapply(or, function(x) unname(x$estimate))

  # Define adaptive weights ----
  # when OR = 0 => beta univ= INF => apw = 0 nope
  apws <- rep(Inf, length(or)) ; names(apws) <- names(or)
  ok <- which(or>0)
  apws[ok] <- 1/(abs(log(or[ok]))^gamma)

  # Adaptive lasso -----
  if(criterion == "bic"){
    # BIC
    ada.univ <- lasso_bic(x = x, y = y, maxp = maxp, path = path, betaPos = betaPos,
                          penalty.factor = apws, ...)

    beta2 <- ada.univ$beta
    select_var <- ada.univ$selected_variables

  }else if(criterion == "cv"){
    if(missing(foldid)){
      foldid <- sample(1:nfolds, size = nrow(x), replace = TRUE)
    }

    if(length(which(apws!=Inf))==0){  #no covariate is selected in the first step: do nothing IMPOSSIBLE BUT...
      beta2 <- NA
      select_var <- character(0)
    }else{
      top <- control_adaCV(dat = x, adaptive.weights = apws, id.folds = foldid)

      if(top>0){ # Not fine, new fold repartition until top = 0
        while (top>0) {
          foldid <- sample(1:nfolds, size = nrow(x), replace = TRUE)
          top <- control_adaCV(dat = x, adaptive.weights = apws, id.folds = foldid )
        }
      }

      ada.cv <- cv.glmnet(x = x, y = y, foldid = foldid, family = "binomial",
                          penalty.factor = apws, ...)
      tmp <- match(ada.cv$lambda.min, ada.cv$glmnet.fit$lambda)
      beta2 <- ada.cv$glmnet.fit$beta[,tmp]

      if(betaPos){ #positive
        select_var <- names(sort(beta2[which(beta2>0)], decreasing = TRUE))
      }else{ # all
        select_var <- names(sort(abs(beta2[which(beta2!=0)]), decreasing = TRUE))
      }

    }

  }else{
    #user get wrong
    stop("criterion parameter is not valid, it must be 'cv' or 'bic'.")
  }




  # Result -----
  res <- list()
  res$aws <- apws
  res$criterion <- criterion
  res$beta <- beta2
  res$selected_variables <- select_var

  class(res) = "adaptive"

  return(res)

}
