#' weihting on propensity score for one drug exposure
#'
#' Implement the weighting on propensity score with Matching Weights (MW)
#' or the Inverse Probability of Treatment Weighting (IPTW) for
#' one drug exposure.
#' The binary outcome is regressed on the drug exposure of interest
#' through a classical weighted regression.
#' Internal function, not supposed to be used directly.
#'
#'
#' The MW are defined by
#' \deqn{mw_i = min(PS_i, 1-PS_i)/[(expo_i) * PS_i + (1-expo_i) * (1-PS_i) ] }
#' and weights from IPTW by
#'  \deqn{iptw_i = expo_i/PS_i + (1-expo_i)/(1-PS_i)  }
#' where \eqn{expo_i} is the drug exposure indicator.
#' The PS could be estimated in different ways: using lasso-bic approach,
#' the hdPS algorithm or gradient tree boosting using functions
#' \code{est_ps_bic}, \code{est_ps_hdps} and \code{est_ps_xgb}
#' respectivelly.
#'
#'
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
#' @param weights_type Character. Indicates which type of weighting
#' is implemented. Could be either "mw" or "iptw".
#' @param truncation Bouleen, should we do weight truncation?
#' Default is \code{FALSE}.
#' @param q If \code{truncation} is \code{TRUE}, quantile value for
#' weight truncation. Ignored if \code{truncation} is \code{FALSE}.
#' Default is 2.5 \eqn{\%}.
#'
#' @return An object with S3 class \code{"ps","*" },
#' where \code{"*"} is \code{"mw"} or \code{"iptw"}, same as the
#' input parameter \code{weights_type}
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
#' pondps2 <- ps_pond_one(ps_est = pshdps2, y = ae, weights_type = "iptw")
#' pondps2$estimate #estimated strength of association between drug_2 and the outcome by PS weighting
#'
#' @author Emeline Courtois \cr Maintainer: Emeline Courtois
#' \email{emeline.courtois@@inserm.fr}
#'
#' @export ps_pond_one


ps_pond_one <- function(ps_est, y, weights_type = c("mw", "iptw") ,
                        truncation = FALSE, q = 0.025){


  if(missing(weights_type)){
    stop("weighting type not specified")
  }

  if(missing(q) & truncation){#troncation quantile non precise
    warning("default quantile for weight truncation is 0.025")
  }

  # Weights calculation ------
  if(weights_type %in% c("mw", "MW")){
    numerateur <- cbind(ps_est$score, 1-ps_est$score)
    numerateur <- apply(numerateur, 1, min)
    denominateur <- (ps_est$indicator_expo * ps_est$score) + (1-ps_est$indicator_expo) * (1-ps_est$score)

    poids <- numerateur/denominateur
  }else if(weights_type %in% c("iptw", "IPTW")){
    poids <- ps_est$indicator_expo/ps_est$score + (1-ps_est$indicator_expo)/(1-ps_est$score)
  }

  if(truncation){
    val <- quantile(as.numeric(poids), probs = c(q, 1-q), na.rm = TRUE)

    poids2 <- as.numeric(poids)
    poids2[which(poids[,1]<val[1])] <- unname(val[1])
    poids2[which(poids[,1]>val[2])] <- unname(val[2])
    poids2 <- Matrix(poids2)
    poids <- poids2
  }


  # Weighted logistic regression -----
  colnames(ps_est$indicator_expo) <- "expo"

  res <- tryCatch(speedglm.wfit(y = as.vector(y), X = cbind(1,ps_est$indicator_expo),
                                family = binomial(), intercept = TRUE, weights = as.numeric(poids)),error = function(e) NA)


  if(inherits(res,"speedglm")){
    pval.2sided <- as.numeric(as.character(summary(res)$coef["expo", "Pr(>|z|)"])) #pvaleur bilaterale
    pval.1sided <-  pnorm(summary(res)$coef["expo", "z value"], lower.tail = FALSE )

    final <- list()

    final$expo_name <- ps_est$expo_name
    final$estimate <- summary(res)$coef["expo", "Estimate"]
    final$pval_1sided <-  pval.1sided
    final$pval_2sided <-  pval.2sided

    class(final) <- c("ps", weights_type)

    return(final)
  }else{
    return(NA)
  }

}
