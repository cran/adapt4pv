#' fit an adaptive lasso with adaptive weights derived from lasso-cv
#'
#' Fit a first lasso regression with cross-validation to determine adaptive weights.
#' Run a cross-validation to determine an optimal lambda.
#' Two options for implementing cross-validation for the adaptive lasso are possible through the \code{type_cv} parameter (see bellow).
#' Can deal with very large sparse data matrices.
#' Intended for binary reponse only (option \code{family = "binomial"} is forced).
#' The cross-validation criterion used is deviance.
#' Depends on the \code{cv.glmnet} function from the package \code{glmnet}.
#'
#' The adaptive weight for a given covariate i is defined by
#' \deqn{w_i = 1/|\beta^{CV}_i|^\gamma} where
#' \eqn{\beta^{CV}_i} is the PENALIZED regression coefficient associated
#' to covariate \eqn{i} obtained with cross-validation.
#'
#' @param x Input matrix, of dimension nobs x nvars. Each row is an observation
#' vector. Can be in sparse matrix format (inherit from class
#' \code{"sparseMatrix"} as in package \code{Matrix}).
#' @param y Binary response variable, numeric.
#' @param gamma Tunning parameter to defined the penalty weights. See details below.
#' Default is set to 1.
#' @param nfolds Number of folds - default is 5. Although \code{nfolds} can be
#' as large as the sample size (leave-one-out CV), it is not recommended for
#' large datasets. Smallest value allowable is \code{nfolds=3}.
#' @param foldid An optional vector of values between 1 and \code{nfolds}
#' identifying what fold each observation is in. If supplied, \code{nfolds} can
#' be missing.
#' @param type_cv Character, indicates which implementation of cross-validation is performed for the adaptive lasso:  a "naive" one,
#' where adaptive weights obtained on the full data are used, and a "proper" one, where adaptive weights are calculated for each training sets.
#' Could be either "naive" or "proper".
#' Default is "proper".
#' @param betaPos Should the covariates selected by the procedure be
#' positively associated with the outcome ? Default is \code{TRUE}.
#' @param \dots Other arguments that can be passed to \code{glmnet}  from package \code{glmnet} other than \code{nfolds}, \code{foldid},
#' \code{penalty.factor}, \code{standardize}, \code{intercept}  and \code{family}.
#'
#' @return An object with S3 class \code{"adaptive"}.
#' \item{aws}{Numeric vector of penalty weights derived from cross-validation.
#' Length equal to nvars.}
#' \item{criterion}{Character, indicates which criterion is used with the
#' adaptive lasso for variable selection. For \code{adapt_cv} function, \code{criterion}
#' is "cv".}
#' \item{beta}{Numeric vector of regression coefficients in the adaptive lasso.
#' If \code{criterion} = "cv" the regression coefficients are PENALIZED, if
#' \code{criterion} = "bic" the regression coefficients are UNPENALIZED.
#' Length equal to nvars. Could be NA if adaptive weights are all equal to infinity.}
#' \item{selected_variables}{Character vector, names of variable(s) selected
#' with this adaptive approach.
#' If \code{betaPos = TRUE}, this set is the covariates with a positive regression
#' coefficient in \code{beta}.
#' Else this set is the covariates with a non null regression coefficient in \code{beta}.
#' Covariates are ordering according to magnitude of their regression
#' coefficients absolute value in the adaptive lasso.}
#' @examples
#'
#' set.seed(15)
#' drugs <- matrix(rbinom(100*20, 1, 0.2), nrow = 100, ncol = 20)
#' colnames(drugs) <- paste0("drugs",1:ncol(drugs))
#' ae <- rbinom(100, 1, 0.3)
#' acv <- adapt_cv(x = drugs, y = ae, nfolds = 5)
#'
#'
#' @author Emeline Courtois \cr Maintainer: Emeline Courtois
#' \email{emeline.courtois@@inserm.fr}
#' @export adapt_cv

adapt_cv <- function(x, y, gamma = 1, nfolds = 5, foldid = NULL, type_cv = "proper", betaPos = TRUE, ...){

  nfolds <- nfolds
  apws <- NULL

  # Check args fold  --------
  if(missing(foldid) || is.null(foldid)){
    foldid <- sample(1:nfolds, size = nrow(x), replace = TRUE)
  }

  # First step : lasso with cross-validation ALL DATA --------
  # use to determine weights (fold indifferent)
  cv <- cv.glmnet(x = x, y = y, nfolds = nfolds, family = "binomial", standardize = T,
                  intercept = T, ...)
  beta1 <- as.numeric(predict.glmnet(object = cv$glmnet.fit, s = cv$lambda.min, type = "coefficients"))[-1]
  names(beta1) <- colnames(x)
  apws <- 1/(abs(beta1)^gamma) #penalty weights

  if(length(which(beta1!=0))==0){  #no covariate is selected in the first step: do nothing
    beta2 <- NA
    select_var <- character(0)
    }
  else if(length(which(beta1!=0))==1){
      beta2 <- beta1
      select_var <- names(beta1[which(beta1>0)])
    }
  else if(length(which(beta1!=0))>1){
    ## Check fold for cross-validation -----
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
    }else{ #cross validation is feasible
      if(type_cv=="naive"){
        # CV naive ----
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

        #CV proper : cross-validation -----
        devmat <- matrix(0, nrow = length(lambda.sequence), ncol = max(foldid))

        for(cc in sort(unique(foldid))){
          Xtrain <- x[which(foldid!=cc), ]
          ytrain <- y[which(foldid!=cc)]
          Xtest <- x[which(foldid==cc), ]
          ytest <- y[which(foldid==cc)]

          id.cv.inner <- sample(1:5, size = nrow(Xtrain), replace = T)

          cv.inner <- cv.glmnet(x = Xtrain, y = ytrain,
                                family = "binomial", standardize = T, intercept = T,
                                foldid = id.cv.inner, ...)
          beta.cv.inner <- as.numeric(predict.glmnet(object = cv.inner$glmnet.fit,
                                                     s = cv.inner$lambda.min, type = "coefficients"))[-1]
          names(beta.cv.inner) <- colnames(x)
          apws.inner <- 1/abs(beta.cv.inner[ abs(beta.cv.inner)>0 ]^gamma)

          #Check weights
          if(length(apws.inner)==0){
            # no var selected, prediction with intercept
            intercept.inner <- as.numeric(predict.glmnet(object = cv.inner$glmnet.fit,
                                                         s = cv.inner$lambda.min, type = "coefficients"))[1]
            z <- Matrix(rep(intercept.inner, nrow(Xtest)), ncol = 1)
            row.names(z) <- row.names(Xtest)
            pred <- 1  / (1 + exp(-(z)) ) #proba du succes pour chaque individus de Dcc
            deviance <- cv.lognet.short(pred.matrix = pred, yy = ytest)
            devmat[,cc] <- apply(deviance, 2, mean, na.rm=TRUE)
          }
          else if(length(apws.inner)==1){
            #prediction avec seulement intercept et la variable selectionnee par cv inner
            intercept.inner <- as.numeric(predict.glmnet(object = cv.inner$glmnet.fit,
                                                         s = cv.inner$lambda.min, type = "coefficients"))[1]
            z <- intercept.inner + Xtest %*%beta.cv.inner
            pred <- 1  / (1 + exp(-(z)) )
            deviance <- cv.lognet.short(pred.matrix = pred, yy = ytest)
            devmat[,cc] <- apply(deviance, 2, mean, na.rm=TRUE)
          }
          else if (length(apws.inner)>1){
            Xtrain <- Xtrain[ , names(apws.inner), drop=FALSE]
            Xtrainstand <- standardize(X = Xtrain)
            Xtrain <- scale(Xtrainstand$Xs, center=FALSE, scale=apws.inner)

            Xtest <- Xtest[ , names(apws.inner), drop=FALSE]
            Xteststand <- standardize(X = Xtest)
            Xtest <- scale(Xteststand$Xs, center=FALSE, scale=apws.inner)

            #adaptive lasso sur donnees scale + sans pf + st =FALSE mais faite via Xtrain
            adaptive.lasso.inner <- glmnet(x = Xtrain, y = ytrain,
                                           family = "binomial", standardize = FALSE,
                                           intercept = T, ... )


            #prediction du Dcc, grille de lambda est la meme que sur le jeux de données entier
            pred <- predict(object = adaptive.lasso.inner, newx = Xtest,
                            type = "response", s = lambda.sequence)
            deviance <- cv.lognet.short(pred.matrix = pred, yy = ytest)
            devmat[,cc] <- apply(deviance, 2 , mean, na.rm=TRUE)
          }
        }#end for
        cvm <- apply(devmat, 1, mean, na.rm=TRUE)
        lambda.min <- lambda.sequence[which.min(cvm)]

        beta2 <- as.numeric(predict.glmnet(object = adaptive.lasso.all, s = lambda.min,
                                             type = "coefficients", ))[-1]
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

    }#end if stop cv

  }#end length(beta1)


  # Result -----
  res <- list()
  res$aws <- apws
  res$criterion <- "cv"
  res$beta <- beta2
  res$selected_variables <- select_var

  class(res) = "adaptive"

  return(res)
}









