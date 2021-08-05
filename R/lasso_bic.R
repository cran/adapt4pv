#' fit a lasso regression and use standard BIC for variable selection
#'
#' Fit a lasso regression and use the Bayesian Information Criterion (BIC)
#' to select a subset of selected covariates.
#' Can deal with very large sparse data matrices.
#' Intended for binary reponse only (option \code{family = "binomial"} is forced).
#' Depends on the \code{glmnet} and \code{relax.glmnet} functions from the package \code{glmnet}.
#'
#' For each tested penalisation parameter \eqn{\lambda}, a standard version of the BIC
#' is implemented.
#' \deqn{BIC_\lambda = - 2 l_\lambda + df(\lambda) * ln (N)}
#' where \eqn{l_\lambda} is the log-likelihood of the non-penalized multiple logistic
#' regression model that includes the set of covariates with a non-zero coefficient
#' in the penalised regression coefficient vector associated to \eqn{\lambda},
#' and  \eqn{df(\lambda)} is the number of covariates with a non-zero coefficient
#' in the penalised regression coefficient vector associated to \eqn{\lambda},
#' The optimal set of covariates according to this approach is the one associated with
#' the classical multiple logistic regression model which minimizes the BIC.
#'
#'
#' @param x Input matrix, of dimension nobs x nvars. Each row is an observation
#' vector. Can be in sparse matrix format (inherit from class
#' \code{"sparseMatrix"} as in package \code{Matrix}).
#' @param y Binary response variable, numeric.
#' @param maxp A limit on how many relaxed coefficients are allowed.
#' Default is 50, in \code{glmnet} option default is 'n-3', where 'n' is the sample size.
#' @param path Since \code{glmnet} does not do stepsize optimization, the Newton
#' algorithm can get stuck and not converge, especially with relaxed fits. With \code{path=TRUE},
#' each relaxed fit on a particular set of variables is computed pathwise using the original sequence
#' of lambda values (with a zero attached to the end). Default is \code{path=TRUE}.
#' @param betaPos Should the covariates selected by the procedure be
#' positively associated with the outcome ? Default is \code{TRUE}.
#' @param \dots Other arguments that can be passed to \code{glmnet} from package
#' \code{glmnet} other than \code{family}, \code{maxp}
#' and \code{path}.
#'
#' @return An object with S3 class \code{"log.lasso"}.
#' \item{beta}{Numeric vector of regression coefficients in the lasso.
#' In \code{lasso_bic} function, the regression coefficients are UNPENALIZED.
#' Length equal to nvars.}
#' \item{selected_variables}{Character vector, names of variable(s) selected with the
#' lasso-bic approach.
#' If \code{betaPos = TRUE}, this set is the covariates with a positive regression
#' coefficient in \code{beta}.
#' Else this set is the covariates with a non null regression coefficient in \code{beta}.
#' Covariates are ordering according to the p-values (two-sided if \code{betaPos = FALSE} ,
#' one-sided if \code{betaPos = TRUE}) in the classical multiple logistic regression
#' model that minimzes the BIC.}
#' @examples
#'
#' set.seed(15)
#' drugs <- matrix(rbinom(100*20, 1, 0.2), nrow = 100, ncol = 20)
#' colnames(drugs) <- paste0("drugs",1:ncol(drugs))
#' ae <- rbinom(100, 1, 0.3)
#' lb <- lasso_bic(x = drugs, y = ae, maxp = 20)
#'
#'
#' @author Emeline Courtois \cr Maintainer: Emeline Courtois
#' \email{emeline.courtois@@inserm.fr}
#' @export lasso_bic

lasso_bic <- function(x, y, maxp = 50, path = TRUE, betaPos = TRUE, ...){

  # Lasso and relax lasso ------------
  # glmnet V3.0-2 bugg with Matrix object
  lasso <- glmnet(x = x, y = y, family = "binomial", ...)

  #x <<- x
  lasso.relax <- relax.glmnet(fit = lasso, x = x, y = y, path = path, maxp = maxp,
                              family = "binomial", ...)

  # BIC -------
  pseudo.bic <- as.data.frame(cbind(DF = lasso.relax$relaxed$df,
                                    Lambda = lasso.relax$relaxed$lambda,
                                    BIC = log(lasso.relax$relaxed$nobs) * lasso.relax$relaxed$df - (lasso.relax$relaxed$nulldev-deviance(lasso.relax$relaxed)) )  )

  #data frame with a quantity proportional to the classical BIC (constant : 2log likelihood of the null model)

  beta.np <- lasso.relax$relaxed$beta[,which.min(pseudo.bic$BIC)]
    #non penalyzed regression coefficient, same if df(lambda) remains unchange

  var1 <- which(beta.np!=0)
  # not selected variable : to refit a classical model betaPos doesn't matter


  # Classic logistic regression ----
  # In order to ordering covariates according to p values, re run a clasical model
  model.cla <- tryCatch(speedglm.wfit(y = as.vector(y),
                                      X = cbind(1,x[,var1]),
                                      family = binomial(), intercept = TRUE),
                        error = function(e) NA)


  # Variables and p-values ----
  coeffs.np <- summary(model.cla)$coefficients[-1,] # without intercept

  if(betaPos){ #positive, one sided p val
    coeffs.np <- coeffs.np[which(coeffs.np$Estimate>0), ] #coeff positifs
    pv <- pnorm(coeffs.np$`z value`, lower.tail = FALSE)
  }else{ # all
    pv <- as.numeric(as.character(coeffs.np$`Pr(>|z|)`))
  }

  select_var <- row.names(coeffs.np)[order(pv, decreasing = FALSE)]

  # Result -----
  res <- list()
  res$beta <- beta.np
  res$selected_variables <- select_var

  class(res) = "log.lasso"

  return(res)

}
