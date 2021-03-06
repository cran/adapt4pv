#' adjustment on propensity score
#'
#' Implement the adjustment on propensity score for all the drug exposures
#' of the input drug matrix \code{x} which have more than a given
#' number of co-occurence with the outcome.
#' The binary outcome is regressed on a drug exposure and its
#' estimated PS, for each drug exposure considered after filtering.
#' With this approach, a p-value is obtained for each drug and a
#' variable selection is performed over the corrected for multiple
#' comparisons p-values.
#'
#'
#' The PS could be estimated in different ways: using lasso-bic approach,
#' the hdps algorithm or gradient tree boosting.
#' The scores are estimated using the default parameter values of
#' \code{est_ps_bic}, \code{est_ps_hdps} and \code{est_ps_xgb} functions
#' (see documentation for details).
#' We apply the same filter and the same multiple testing correction as in
#' the paper UPCOMING REFERENCE: first, PS are estimated only for drug covariates which have
#' more than \code{n_min} co-occurence with the outcome \code{y}.
#' Adjustment on the PS is performed for these covariates and
#' one sided or two-sided (depend on \code{betaPos} parameter)
#' p-values are obtained.
#' The p-values of the covariates not retained after filtering are set to 1.
#' All these p-values are then adjusted for multiple comparaison with the
#' Benjamini-Yekutieli correction.
#' COULD BE VERY LONG. Since this approach (i) estimate a score for several
#' drug covariates and (ii) perform an adjustment on these scores,
#' parallelization is highly recommanded.
#'
#'
#' @param x Input matrix, of dimension nobs x nvars. Each row is an observation
#' vector. Can be in sparse matrix format (inherit from class
#' \code{"sparseMatrix"} as in package \code{Matrix}).
#' @param y Binary response variable, numeric.
#' @param n_min Numeric, Minimal number of co-occurence between a drug
#' covariate and the outcome y to estimate its score. See details belows.
#' Default is 3.
#' @param betaPos Should the covariates selected by the procedure be
#' positively associated with the outcome ? Default is \code{TRUE}.
#' @param est_type Character, indicates which approach is used to estimate
#' the PS.
#' Could be either "bic", "hdps" or "xgb".
#' Default is "bic".
#' @param threshold Threshold for the p-values. Default is 0.05.
#' @param ncore The number of calcul units used for parallel computing.
#' Default is 1, no parallelization is implemented.
#'
#' @return An object with S3 class \code{"ps", "adjust", "*"}, where
#' \code{"*"} is \code{"bic"}, \code{"hdps"} or \code{"xgb"}according on how the
#' score were estimated.
#' \item{estimates}{Regression coefficients associated with
#' the drug covariates. Numeric, length equal to the number of selected
#' variables with this approach.
#' Some elements could be NA if
#' (i) the corresponding covariate was filtered out,
#' (ii) adjustment model did not converge. Trying to estimate the score in
#' a different way could help, but it's not insured.}
#' \item{corrected_pvals}{One sided p-values if \code{betaPos = TRUE},
#' two-sided p-values if \code{betaPos = FALSE} adjusted for multiple testing.
#' Numeric, length equal to nvars.}
#' \item{selected_variables}{Character vector, names of variable(s)
#' selected with the ps-adjust approach.
#' If \code{betaPos = TRUE}, this set is the covariates with a
#' corrected one-sided p-value lower than \code{threshold}.
#' Else this set is the covariates with a
#' corrected two-sided p-value lower than \code{threshold}.
#' Covariates are ordering according to their corrected p-value.}
#' @examples
#'
#' set.seed(15)
#' drugs <- matrix(rbinom(100*20, 1, 0.2), nrow = 100, ncol = 20)
#' colnames(drugs) <- paste0("drugs",1:ncol(drugs))
#' ae <- rbinom(100, 1, 0.3)
#' adjps <- ps_adjust(x = drugs, y = ae, n_min = 10)
#'
#' @author Emeline Courtois \cr Maintainer: Emeline Courtois
#' \email{emeline.courtois@@inserm.fr}
#'
#' @references Benjamini, Y., & Yekuteli, D. (2001). "The Control of the False Discovery Rate in Multiple Testing under Dependency".
#' \emph{The Annals of Statistics}. 29(4), 1165–1188, doi: \doi{10.1214/aos/1013699998}.
#'
#' @export ps_adjust



ps_adjust <- function(x, y, n_min = 3, betaPos = TRUE, est_type = "bic",
                      threshold = 0.05, ncore = 1){

  ## Step 1 : filter -----
  co_occurence <- t(t(y) %*% x)
  covar_ok_idx <- which(co_occurence[,1] >= n_min)
  covar_nope_idx <- which(co_occurence[,1] < n_min)


  ## Step 2 : multiple estimation of scores -----
  # scores are implemented only for the "ok" covariables,
  # all the variable of x are candidates covariables
  i <- 0
  if(est_type == "bic"){
    if(ncore == 1 ){ #no parallelisation
      step2 <- lapply(covar_ok_idx, est_ps_bic, x = x)
    }else if(ncore >1){
      # #parallelisation UNIX
      # step2 <- mclapply(covar_ok_idx, est_ps_bic, x = x,
      #                   mc.preschedule = FALSE, mc.cores = ncore)

      #parallelization all OX
      cl <- makeCluster(ncore)
      registerDoParallel(cl)
      step2 = foreach(i = covar_ok_idx, .packages = c("glmnet", "speedglm"),
                      .export = c("est_ps_bic", "predict_speedglm.wfit")) %dopar% {
        est_ps_bic(i, x= x)
      }
      stopCluster(cl)
      names(step2) <- names(covar_ok_idx)


    }

  }else if(est_type == "hdps"){
    if(ncore == 1 ){ #no parallelisation
      step2 <- lapply(covar_ok_idx, est_ps_hdps, x = x, y = y)
    }else if(ncore >1){
      # #parallelisation UNIX
      # step2 <- mclapply(covar_ok_idx, est_ps_hdps, x = x, y = y,
      #                   mc.preschedule = FALSE, mc.cores = ncore)

      #parallelization all OX
      cl <- makeCluster(ncore)
      registerDoParallel(cl)
      step2 = foreach(i = covar_ok_idx, .packages = c("speedglm"),
                      .export = c("est_ps_hdps", "hdps_pv", "predict_speedglm.wfit")) %dopar% {
        est_ps_hdps(i, x= x, y = y)
      }
      stopCluster(cl)
      names(step2) <- names(covar_ok_idx)

    }

  }else if(est_type == "xgb"){
    if(ncore == 1 ){ #no parallelisation
      step2 <- lapply(covar_ok_idx, est_ps_xgb, x = x)
    }else if(ncore >1){
      # #parallelisation UNIX
      # step2 <- mclapply(covar_ok_idx, est_ps_xgb, x = x,
      #                   mc.preschedule = FALSE, mc.cores = ncore)

      #parallelization all OX
      cl <- makeCluster(ncore)
      registerDoParallel(cl)
      step2 = foreach(i = covar_ok_idx, .packages = c("xgboost", "speedglm"),
                      .export = "est_ps_xgb") %dopar% {
        est_ps_xgb(i, x= x)
      }
      stopCluster(cl)
      names(step2) <- names(covar_ok_idx)
    }

  }else{
    stop("Invalide est_type argument")
  }


  # Step 3 : Adjustment -----
  # for the variables for which we estimated a score, adjustment
  if(ncore == 1 ){ #no parallelisation
    step3 <- lapply(step2, ps_adjust_one, y = y)
  }else if(ncore >1){
    # #parallelisation UNIX
    # step3 <- mclapply(step2, ps_adjust_one, y = y, mc.preschedule = FALSE,
    #                   mc.cores = ncore)

    #parallelization all OX
    cl <- makeCluster(ncore)
    registerDoParallel(cl)
    step3 = foreach(i = 1:length(step2), .packages = c("speedglm"),
                    .export = "ps_adjust_one") %dopar% {
      ps_adjust_one(ps_est = step2[[i]], y = y)
    }
    stopCluster(cl)
    names(step3) <- names(step2)

  }


  # Step 4 : extraction of p-values ----
  # Depend on betaPos
  if(betaPos){
    step4 <- unlist(lapply(step3, function(x) try_bis(x$pval_1sided)))
  }else{
    step4 <- unlist(lapply(step3, function(x) try_bis(x$pval_2sided)))
  }

  # Step 5 : correction ------
  #filtered covar ave their p-value set to 1, then BY correction
  step5 <- rep(1, ncol(x)) ; names(step5) <- colnames(x)
  step5[colnames(x)[covar_ok_idx]] <- step4[colnames(x)[covar_ok_idx]]
  if(sum(!is.na(step5)) > 20 ) step5 <- p.adjust(step5, method = "BY")


  #Step 6 : covariables selection + output ----
  variables <- names(sort(step5[which(step5<=threshold)], decreasing = FALSE))

  final <- list()
  final$estimate <- unlist(lapply(step3, function(x) try_bis(x$estimate)))
  final$estimate <- final$estimate[variables]
  final$corrected_pvals <- step5
  final$selected_variables <- variables

  class(final) <- c("ps", "adjust" , est_type)

  return(final)

}




