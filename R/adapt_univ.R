#' fit an adaptive lasso with adaptive weights derived from univariate coefficients
#'
#' Compute odd-ratios between each covariate of \code{x} and \code{y} then derived
#' adaptive weights to incorporate in an adaptive lasso.
#' BIC or cross-validation could either be used for the adaptive lasso for variable selection.
#' Two options for implementing cross-validation for the adaptive lasso are possible through the \code{type_cv} parameter (see bellow).
#' Can deal with very large sparse data matrices.
#' Intended for binary reponse only (option \code{family = "binomial"} is forced).
#' The cross-validation criterion used is deviance.
#' Depends on the \code{glmnet} and \code{relax.glmnet} function from the package
#' \code{glmnet}.
#'
#' The adaptive weight for a given covariate i is defined by
#' \deqn{w_i = 1/|\beta^{univ}_i|^\gamma} where
#'  \eqn{\beta^{univ}_i = log(OR_i)}, with
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
#' @param type_cv Used only if \code{criterion} = "cv", ignored if \code{criterion} = "bic".
#' Character, indicates which implementation of cross-validation is performed for the adaptive lasso:  a "naive" one,
#' where adaptive weights obtained on the full data are used, and a "proper" one, where adaptive weights are calculated for each training sets.
#' Could be either "naive" or "proper".
#' Default is "proper".
#' @param betaPos Should the covariates selected by the procedure be
#' positively associated with the outcome ? Default is \code{TRUE}.
#' @param \dots Other arguments that can be passed to \code{glmnet} from package
#' \code{glmnet} other than \code{family}, \code{maxp}, \code{standardize}, \code{intercept}
#'
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
                       nfolds = 5, foldid = NULL, type_cv = "proper", betaPos = TRUE, ...){

  # Odds ratio ----
  tab22 <- contingence_all(data_cases = x[which(y==1),], data_controls = x[which(y==0),])
  or <- lapply(tab22, fisher.test)
  or <- sapply(or, function(x) unname(x$estimate))

  # Define adaptive weights ----
  # when OR = 0 => beta univ= INF => apw = 0 nope
  apws <- rep(Inf, length(or)) ; names(apws) <- names(or)
  ok <- which(or>0)
  apws[ok] <- 1/(abs(log(or[ok]))^gamma)

  if(length(which(apws!=Inf))==0){  #no covariate is selected in the first step: do nothing IMPOSSIBLE BUT...
    beta2 <- NA
    select_var <- character(0)
  }
  else{ #Adaptive lasso is feasible
    if(criterion == "bic"){
      # BIC
      ada.univ <- lasso_bic(x = x, y = y, maxp = maxp, path = path, betaPos = betaPos,
                            penalty.factor = apws, ...)

      beta2 <- ada.univ$beta
      select_var <- ada.univ$selected_variables
    }
    else if(criterion == "cv"){
      if(missing(foldid) || is.null(foldid)){
        foldid <- sample(1:nfolds, size = nrow(x), replace = TRUE)
      }
      #Check folds
      top <- control_adaCV(dat = x, adaptive.weights = apws, id.folds = foldid)
      stop.cv <- 1

      if(top>0){ # Not fine, new fold repartition until top = 0
        while (top>0 & stop.cv<=1000) {
          foldid <- sample(1:nfolds, size = nrow(x), replace = TRUE)
          top <- control_adaCV(dat = x, adaptive.weights = apws, id.folds = foldid)
          stop.cv <- stop.cv+1
        }
      }

      if(stop.cv==1001){
        beta2 <- NA
        select_var <- character(0)
      }
      else{ #cross validation is feasible
        if(type_cv=="naive"){
          adaptive.cv.lasso.all <- cv.glmnet(x = x, y = y, family = "binomial", standardize = T,
                                             intercept = T,foldid = foldid, penalty.factor = apws, ...)
          beta2 <- as.numeric(predict.glmnet(object = adaptive.cv.lasso.all$glmnet.fit, s = adaptive.cv.lasso.all$lambda.min, type = "coefficients"))[-1]
          names(beta2) <- colnames(x)

        }
        else if(type_cv=="proper"){
          # CV proper adaptive lasso with the previous adaptive weights ALL DATA ----
          # use for lambda.sequence
          xall <- x[, which(apws!=Inf), drop = F] #bis
          xallstand <- standardize(xall)
          xall <- scale(xallstand$Xs, center=FALSE, scale=apws[apws!=Inf])

          adaptive.lasso.all <- glmnet(x = xall, y = y, family = "binomial", standardize = FALSE, intercept = T, ...)
          lambda.sequence <- adaptive.lasso.all$lambda

          devmat <- matrix(0, nrow = length(lambda.sequence), ncol = max(foldid))

          # Cross validation with adaptive weights determined on training sets
          for(cc in sort(unique(foldid))){
            Xtrain <- x[which(foldid!=cc), ]
            ytrain <- y[which(foldid!=cc)]
            Xtest <- x[which(foldid==cc), ]
            ytest <- y[which(foldid==cc)]

            tab22.inner <- contingence_all(data_cases = Xtrain[which(ytrain==1),], data_controls = Xtrain[which(ytrain==0),])
            or.inner <- lapply(tab22.inner, fisher.test)
            or.inner <- sapply(or.inner, function(x) unname(x$estimate))
            apws.inner <- 1/abs(log(or.inner[or.inner>0])^gamma)

            Xtrain <- Xtrain[ , names(apws.inner), drop=FALSE]
            Xtrainstand <- standardize(X = Xtrain)
            Xtrain <- scale(Xtrainstand$Xs, center=FALSE, scale=apws.inner)

            Xtest <- Xtest[ , names(apws.inner), drop=FALSE]
            Xteststand <- standardize(X = Xtest)
            Xtest <- scale(Xteststand$Xs, center=FALSE, scale=apws.inner)

            adaptive.lasso.inner <- glmnet(x = Xtrain, y = ytrain, family = "binomial", standardize = FALSE,
                                           intercept = T, ...)

            pred <- predict(object = adaptive.lasso.inner, newx = Xtest,
                            s = lambda.sequence,  type = "response")
            deviance <- cv.lognet.short(pred.matrix = pred, yy = ytest)
            devmat[,cc] <- apply(deviance, 2, mean, na.rm=TRUE)

          }

          cvm <- apply(devmat, 1, mean, na.rm=TRUE)
          lambda.min <- lambda.sequence[which.min(cvm)]

          beta2 <- as.numeric(predict.glmnet(object = adaptive.lasso.all, s = lambda.min,
                                             type = "coefficients"))[-1]
          names(beta2) <- names(apws[apws!=Inf])
          #manip to have length(beta) = ncol(x)
          beta2 <- beta2[colnames(x)]
          beta2[is.na(beta2)] <- 0
          names(beta2) <- colnames(x)
        }
        else{
          stop("type_cv parameter is not valid, it must be 'naive' or 'proper' (character).")
        }#end cv type

        if(betaPos){ #positive
          select_var <- names(sort(beta2[which(beta2>0)], decreasing = TRUE))
        }else{ # all
          select_var <- names(sort(abs(beta2[which(beta2!=0)]), decreasing = TRUE))
        }

      }#end stop.cv
    }
    else{
      stop("criterion parameter is not valid, it must be 'cv' or 'bic' (character).")
    }#end criterion


  }#end adaptive lasso feasible

  # Result -----
  res <- list()
  res$aws <- apws
  res$criterion <- criterion
  res$beta <- beta2
  res$selected_variables <- select_var

  class(res) = "adaptive"

  return(res)

}

