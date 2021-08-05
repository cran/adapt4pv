#' Class Imbalanced Subsampling Lasso
#'
#' Implementation of CISL and the stability selection according to subsampling options.
#'
#' CISL is a variation of the stability method adapted to characteristics of pharmacovigilance databases.
#' Tunning \code{r = 4} and \code{replace = TRUE} are used to implement our CISL sampling.
#' For instance, \code{r = NULL} and \code{replace = FALSE} can be used to
#' implement the \eqn{n \over 2} sampling in Stability Selection.
#'
#'
#' @param x Input matrix, of dimension nobs x nvars. Each row is an
#' observation vector. Can be in sparse matrix format (inherit from class
#' \code{"sparseMatrix"} as in package \code{Matrix}).
#' @param y Binary response variable, numeric.
#' @param r Number of control in the CISL sampling. Default is 4.
#' See details below for other implementations.
#' @param nB Number of sub-samples. Default is 100.
#' @param dfmax Corresponds to the maximum size of the models visited with the
#' lasso (E in the paper). Default is 50.
#' @param nlambda Number of lambda values as is \code{glmnet} documentation.
#' Default is 250.
#' @param nMin Minimum number of events for a covariate to be considered.
#' Default is 0, all the covariates from \code{x} are considered.
#' @param replace Should sampling be with replacement? Default is TRUE.
#' @param betaPos  If \code{betaPos=TRUE}, variable selection is based on positive
#' regression coefficient.
#' Else, variable selection is based on non-zero regression coefficient.
#' Default is TRUE.
#' @param ncore The number of calcul units used for parallel computing.
#' This has to be set to 1 if the \code{parallel} package is not available.
#' Default is 1.
#' WARNING: parallel computing is not supported for windows machines!
#'
#' @return An object with S3 class \code{"cisl"}.
#' \item{prob}{Matrix of dimension nvars x \code{nB}.
#' Quantity compute by CISL for each covariate, for each subsample.}
#' \item{q05}{5 \eqn{\%} quantile of the CISL quantity for each covariates.
#' Numeric, length equal to nvars.}
#' \item{q10}{10 \eqn{\%} quantile of the CISL quantity for each covariates.
#' Numeric, length equal to nvars.}
#' \item{q15}{15 \eqn{\%} quantile of the CISL quantity for each covariates.
#' Numeric, length equal to nvars.}
#' \item{q20}{20 \eqn{\%} quantile of the CISL quantity for each covariates.
#' Numeric, length equal to nvars.}
#' @examples
#'
#' set.seed(15)
#' drugs <- matrix(rbinom(100*20, 1, 0.2), nrow = 100, ncol = 20)
#' colnames(drugs) <- paste0("drugs",1:ncol(drugs))
#' ae <- rbinom(100, 1, 0.3)
#' lcisl <- cisl(x = drugs, y = ae, nB = 50)
#'
#' @author Ismail Ahmed
#' @references Ahmed, I., Pariente, A., & Tubert-Bitter, P. (2018). "Class-imbalanced subsampling lasso algorithm for discovering adverse drug reactions".
#' \emph{Statistical Methods in Medical Research}. 27(3), 785–797, \doi{10.1177/0962280216643116}
#'
#' @export cisl

cisl <- function(x, y, r = 4, nB = 100, dfmax = 50, nlambda = 250, nMin = 0,
                 replace = TRUE, betaPos = TRUE, ncore = 1) {

  nObs <- nrow(x)
  nS <- round(nrow(x) / 2)
  y <- as.numeric(y)
  n <- t(x) %*% y
  nIdx <- which(n >= nMin)
  x <- x[, nIdx, drop = FALSE]

  if (length(nIdx) > 1){
    R <- 4 * length(nIdx)
    if (is.null(r)) {
      if (replace){
        idx <- as.data.frame(matrix(sample(1:nObs, nS * nB, replace = TRUE), ncol = nB))
      }else{
        idx <- matrix(ncol = nB, nrow = nS)
        for (b in 1:nB) idx[, b] <- sample(1:nObs, nS, replace = FALSE)
        idx <- as.data.frame(idx)
      }
    }else{
      idx <- sampleImb(y, nB = nB, r = r, nS0min = R, replace = replace)
    }

    if (ncore > 1){
      # # parallelisation UNIX
      # resStab  <- mclapply(idx, glmnetSub, x = x, y = y,
      #                      mc.cores = ncore, dfmax = dfmax,
      #                      nlambda = nlambda)

      # parallelization all OX
      i <- 0
      cl <- makeCluster(ncore)
      registerDoParallel(cl)
      resStab = foreach(i=1:nB, .packages = "glmnet", .export = "glmnetSub") %dopar% {
        glmnetSub(idx = idx[,i], x = x, y = y,
                  dfmax = dfmax, nlambda = nlambda)
      }
      stopCluster(cl)
      names(resStab) <- colnames(idx)
    }else{
      resStab  <- lapply(idx, glmnetSub, x = x, y = y,
                         dfmax = dfmax, nlambda = nlambda)
    }

    prob <- probStabMean(resStab, dfmax = dfmax, nCore = ncore, betaPos = betaPos)$prob
  }else{
    prob <- NA
  }
  res <- vector("list")
  res$prob <- prob
  res$q05 <- apply(prob, 1, quantile, 0.05)
  res$q10 <- apply(prob, 1, quantile, 0.10)
  res$q15 <- apply(prob, 1, quantile, 0.15)
  res$q20 <- apply(prob, 1, quantile, 0.20)

  class(res) = "cisl"

  return(res)
}

