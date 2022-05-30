#' adjustment on propensity score for one drug exposure
#'
#' Implement the adjustment on propensity score for one drug exposure.
#' The binary outcome is regressed on the drug exposure of interest and
#' its estimated PS.
#' Internal function, not supposed to be used directly.
#'
#' The PS could be estimated in different ways: using lasso-bic approach,
#' the hdPS algorithm or gradient tree boosting using functions
#' \code{est_ps_bic}, \code{est_ps_hdps} and \code{est_ps_xgb}
#' respectivelly.
#'
#' @param ps_est An object of class \code{"ps", "*"} where \code{"*"} is
#' \code{"bic"}, \code{"hdps"} or \code{"xgb"} according on how the
#' score was estimated, respective outputs of internal functions
#' \code{est_ps_bic}, \code{est_ps_hdps}, \code{est_ps_xgb}.
#' It is a list with the following elements :
#' * score_type: character, name of the drug exposure for which the PS was
#' estimated.
#' * indicator_expo: indicator of the drugs exposure for which the
#' PS was estimated. One-column Matrix object.
#' * score_variables: Character vector, names of covariate(s) selected
#' to include in the PS estimation model. Could be empty.
#' *score: One-column Matrix object, the estimated score.
#' @param y Binary response variable, numeric.
#'
#' @return An object with S3 class \code{"ps","adjust" }
#' \item{expo_name}{Character, name of the drug exposure for which the PS
#' was estimated.}
#' \item{estimate}{Regression coefficient associated with the drug exposure
#' in adjustment on PS.}
#' \item{pval_1sided}{One sided p-value associated with the drug exposure
#' in adjustment on PS.}
#' \item{pval_2sided}{Two sided p-value associated with the drug exposure
#' in adjustment on PS.}
#' Could return NA if the adjustment on the PS did not converge.
#' @examples
#'
#' set.seed(15)
#' drugs <- matrix(rbinom(100*20, 1, 0.2), nrow = 100, ncol = 20)
#' colnames(drugs) <- paste0("drugs",1:ncol(drugs))
#' ae <- rbinom(100, 1, 0.3)
#' pshdps2 <- est_ps_hdps(idx_expo = 2, x = drugs, y = ae, keep_total = 10)
#' adjps2 <- ps_adjust_one(ps_est = pshdps2, y = ae)
#' adjps2$estimate #estimated strength of association between drug_2 and the outcome by PS adjustment
#'
#' @author Emeline Courtois \cr Maintainer: Emeline Courtois
#' \email{emeline.courtois@@inserm.fr}
#'
#' @export ps_adjust_one



ps_adjust_one <- function(ps_est, y){

  # Formating data -----
  mat <- cbind(ps_est$indicator_expo,ps_est$score)
  colnames(mat) <- c("expo", "score")

  # Adjustment ----
  res <- tryCatch(speedglm.wfit(y = as.vector(y), X = cbind(1,mat),
                                family = binomial(), intercept = TRUE),
                  error = function(e) NA)

  if(inherits(res, "speedglm")){
    pval.2sided <- as.numeric(as.character(summary(res)$coef["expo", "Pr(>|z|)"])) # 2sided pvalue
    pval.1sided <-  pnorm(summary(res)$coef["expo", "z value"], lower.tail = FALSE )

    final <- list()

    final$expo_name <- ps_est$expo_name
    final$estimate <- summary(res)$coef["expo", "Estimate"]
    final$pval_1sided <-  pval.1sided
    final$pval_2sided <- pval.2sided

    class(final) <- c("ps", "adjust")

    return(final)
  }else{
    return(NA)
  }


}
