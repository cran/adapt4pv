#' propensity score estimation in high dimension with automated covariates selection using lasso-bic
#'
#' Estimate a propensity score to a given drug exposure by
#' (i) selecting among other drug covariates in \code{x} which ones to
#' include in the PS estimation model automatically using lasso-bic
#' approach,
#' (ii) estimating a score using a classical logistic regression
#' with the afore selected covariates.
#' Internal function, not supposed to be used directly.
#'
#' \code{betaPos} option of \code{lasso_bic} function is set to
#' \code{FALSE} and \code{maxp} is set to 20.
#' For optimal storage, the returned elements \code{indicator_expo} and
#' \code{score} are Matrix with ncol = 1.
#'
#' @param idx_expo Index of the column in \code{x} that corresponds to the
#' drug covariate for which we aim at estimating the PS.
#' @param x Input matrix, of dimension nobs x nvars. Each row is an observation
#' vector. Can be in sparse matrix format (inherit from class
#' \code{"sparseMatrix"} as in package \code{Matrix}).
#' @param penalty TEST OPTION penalty weights in the variable selection to
#' include in the PS.
#' @param \dots Other arguments that can be passed to \code{glmnet} from package
#' \code{glmnet} other than \code{penalty.factor},
#' \code{family}, \code{maxp} and \code{path}.
#'
#'
#' @return An object with S3 class \code{"ps", "bic"}.
#' \item{expo_name}{Character, name of the drug exposure for which the PS was
#' estimated. Correspond to \code{colnames(x)[idx_expo]}}.
#' \item{indicator_expo}{One-column Matrix object.
#' Indicator of the drug exposure for which the PS was estimated.
#' Defined by \code{x[, idx_expo]}.}.
#' \item{score_variables}{Character vector, names of covariates(s) selected
#' with the lasso-bic approach to include in the PS estimation model.
#' Could be empty.}
#' \item{score}{One-column Matrix object, the estimated score.}
#' @examples
#'
#' set.seed(15)
#' drugs <- matrix(rbinom(100*20, 1, 0.2), nrow = 100, ncol = 20)
#' colnames(drugs) <- paste0("drugs",1:ncol(drugs))
#' ae <- rbinom(100, 1, 0.3)
#' psb2 <- est_ps_bic(idx_expo = 2, x = drugs)
#' psb2$score_variables #selected variables to include in the PS model of drug_2
#'
#' @author Emeline Courtois \cr Maintainer: Emeline Courtois
#' \email{emeline.courtois@@inserm.fr}
#' @export est_ps_bic


est_ps_bic <- function(idx_expo, x, penalty = rep(1,nvars-1), ...){

  nvars <- ncol(x)

  # Drug exposure indicator ----
  indic <- x[,idx_expo]


  # Penalty specification ----
  # penalty : numeric vector of legth nvars - 1
  if(!missing(penalty)){
    le <- length(penalty)
    if(le!=(nvars-1) & le==nvars){ #pb
      warning("penalty option is a vector of size nvars, the idx_expo element will be discared")
      penalty <- penalty[-idx_expo]
    } else if(le!= (nvars-1) & le!=nvars){
      stop("penalty option not valid")
    }
  }





  # Covariate selection ----
  # direct use of lasso-bic option is not implement here because
  # of some bugs in relax.glmnet, here the data set tmp is set as a global
  # variable to fix it
  # PS  = P(T=1|X)
  tmp <-  x[,-idx_expo]
  ll <- glmnet(x = tmp, y = indic, family = "binomial", penalty.factor = penalty, ...)
  tmp <<- tmp
  rr <- relax.glmnet(fit = ll, x = tmp, y = indic, family = "binomial", path = TRUE, maxp = 20, penalty.factor = penalty, ... )
  rm(tmp)


  pseudo.bic <- as.data.frame(cbind(DF = rr$relaxed$df,
                                    Lambda = rr$relaxed$lambda,
                                    BIC = log(rr$relaxed$nobs) * rr$relaxed$df - (rr$relaxed$nulldev-deviance(rr$relaxed)) )  )


  beta <- rr$relaxed$beta[,which.min(pseudo.bic$BIC)]
  var <- names(which(beta!=0)) # not >0 : all the covariates associated with the exposure


  # Score estimation ----
  modele <- tryCatch(speedglm.wfit(y = as.vector(indic), X = cbind(1,x[,var]), family = binomial(), intercept = TRUE),error = function(e) NA)


  if(inherits(modele,"speedglm")){
    score <- predict_speedglm.wfit(speedglm = modele, newmatrix = x[, var])
    score <-  Matrix(score)
  }else{
    score <- NA
  }

  indic <- Matrix(indic)

  # Results -----
  res <- list()
  res$expo_name <- colnames(x)[idx_expo]
  res$indicator_expo <- indic
  res$score_variables <- var
  res$score <- score

  class(res) <- c("ps", "bic")

  return(res)

}
