#' propensity score estimation in high dimension with automated covariates selection using hdPS
#'
#' Estimate a propensity score to a given drug exposure by
#' (i) selecting among other drug covariates in \code{x} which ones to
#' include in the PS estimation model automatically using hdPS algorithm,
#' (ii) estimating a score using a classical logistic regression
#' with the afore selected covariates.
#' Internal function, not supposed to be used directly.
#'
#' Compared to the situation of the classic use of hdps
#' (i) there is only one dimension (the co-exposition matrix)
#' (ii) no need to expand covariates since they are already binary.
#' In other words, in our situation hdps consists in the "prioritize covariates"
#' step from the original algorithm, using Bross formula.
#' We consider the correction on the interpretation on this formula made
#' by Richard Wyss (drug epi).
#'
#' @param idx_expo Index of the column in \code{x} that corresponds to the
#' drug covariate for which we aim at estimating the PS.
#' @param x Input matrix, of dimension nobs x nvars. Each row is an observation
#' vector. Can be in sparse matrix format (inherit from class
#' \code{"sparseMatrix"} as in package \code{Matrix}).
#' @param y Binary response variable, numeric.
#' @param keep_total number of covariates to include in the PS estimation
#' model according to the hdps algorithm ordering. Default is 20.
#'
#' @return An object with S3 class \code{"ps", "hdps"}.
#' \item{expo_name}{Character, name of the drug exposure for which the PS was
#' estimated. Correspond to \code{colnames(x)[idx_expo]}}.
#' \item{indicator_expo}{One-column Matrix object. Indicator of the drug
#' exposure for which the PS was estimated.
#' Defined by \code{x[, idx_expo]}.}.
#' \item{score_variables}{Character vector, names of covariates(s)
#' selected with the hdPS algorithm to include in the PS estimation model.
#' Could be empty.}
#' \item{score}{One-column Matrix object, the estimated score.}
#' @examples
#'
#' set.seed(15)
#' drugs <- matrix(rbinom(100*20, 1, 0.2), nrow = 100, ncol = 20)
#' colnames(drugs) <- paste0("drugs",1:ncol(drugs))
#' ae <- rbinom(100, 1, 0.3)
#' pshdps2 <- est_ps_hdps(idx_expo = 2, x = drugs, y = ae, keep_total = 10)
#' pshdps2$score_variables #selected variables to include in the PS model of drug_2
#'
#' @author Emeline Courtois \cr Maintainer: Emeline Courtois
#' \email{emeline.courtois@@inserm.fr}
#'
#' @references Schneeweiss, S., Rassen, J. A., Glynn, R. J., Avorn, J., Mogun, H., Brookhart, M. A. (2009).
#' "High-dimensional propensity score adjustment in studies of treatment effects using health care claims data".
#' \emph{Epidemiology}. 20, 512–522, \doi{10.1097/EDE.0b013e3181a663cc}
#'
#' @export est_ps_hdps


est_ps_hdps <- function(idx_expo, x, y, keep_total = 20){

  # hdps algorithm in pv -----
  C <- x[,-idx_expo] #confusion
  E <- x[, idx_expo] #exposure
  D <- y #disease

  var <- hdps_pv(E = E, C = C, D = D, k = keep_total)

  indic <- Matrix(x[,idx_expo] ) #for output

  # Estimation PS -----
  ps <- tryCatch(speedglm.wfit(y = as.vector(indic), X = cbind(1,x[,var]), family = binomial(), intercept = TRUE),error = function(e) NA)
  if(inherits(ps, "speedglm")){
    ps <- predict_speedglm.wfit(speedglm = ps, newmatrix = x[, var])
    ps <- Matrix(ps)
  }else{
    ps <- NA
  }

  # Results -----
  res <- list()
  res$expo_name <- colnames(x)[idx_expo]
  res$indicator_expo <- indic
  res$score_variables <- var
  res$score <- ps

  class(res) = c("ps", "hdps")

  return(res)
}
