#' fit a lasso regression and use standard permutation of the outcome for variable selection
#'
#' Performed K lasso logistic regression with K different permuted version of the outcome.
#' For earch of the lasso regression, the \eqn{\lambda_max}(i.e. the smaller
#' \eqn{\lambda} such as all penalized regression coefficients are shrunk to zero)
#' is obtained.
#' The median value of these K  \eqn{\lambda_max} is used to for variable selection
#' in the lasso regression with the non-permuted outcome.
#' Depends on the \code{glmnet} function from the package \code{glmnet}.
#'
#' The selected \eqn{\lambda} with this approach is defined as the closest
#' \eqn{\lambda} from the median value of the K \eqn{\lambda_max} obtained
#' with permutation of the outcome.
#'
#' @param x Input matrix, of dimension nobs x nvars. Each row is an observation
#' vector. Can be in sparse matrix format (inherit from class
#' \code{"sparseMatrix"} as in package \code{Matrix}).
#' @param y Binary response variable, numeric.
#' @param K Number of permutations of \code{y}. Default is 20.
#' @param keep Do some variables of \code{x} have to be permuted in the same way
#' as \code{y}? Default is NULL, means no.
#' If yes, must be a vector of covariates indices. TEST OPTION
#' @param betaPos Should the covariates selected by the procedure be positively
#' associated with the outcome ? Default is \code{TRUE}.
#' @param ncore The number of calcul units used for parallel computing.
#' Default is 1, no parallelization is implemented.
#' @param \dots Other arguments that can be passed to \code{glmnet}
#' from package \code{glmnet} other than \code{family}.
#'
#' @return An object with S3 class \code{"log.lasso"}.
#' \item{beta}{Numeric vector of regression coefficients in the lasso
#' In \code{lasso_perm} function, the regression coefficients are PENALIZED.
#' Length equal to nvars.}
#' \item{selected_variables}{Character vector, names of variable(s) selected with the
#' lasso-perm approach.
#' If \code{betaPos = TRUE}, this set is the covariates with a positive regression
#' coefficient in \code{beta}.
#' Else this set is the covariates with a non null regression coefficient in
#' \code{beta}.
#' Covariates are ordering according to magnitude of their regression
#' coefficients absolute value.}
#' @examples
#'
#' set.seed(15)
#' drugs <- matrix(rbinom(100*20, 1, 0.2), nrow = 100, ncol = 20)
#' colnames(drugs) <- paste0("drugs",1:ncol(drugs))
#' ae <- rbinom(100, 1, 0.3)
#' lp <- lasso_perm(x = drugs, y = ae, K = 10)
#'
#'
#' @author Emeline Courtois \cr Maintainer: Emeline Courtois
#' \email{emeline.courtois@@inserm.fr}
#'
#' @references Sabourin, J. A., Valdar, W., & Nobel, A. B. (2015). "A permutation approach for selecting the penalty parameter in penalized model selection".
#' \emph{Biometrics}. 71(4), 1185–1194, \doi{10.1111/biom.12359}
#'
#' @export lasso_perm


lasso_perm <- function(x, y, K = 20, keep = NULL, betaPos = TRUE, ncore = 1, ...){


  y <- as.numeric(y)

  if(K==0){
    stop("Error : K invalid")
  }

  ## Permutations ----
  perm.id <- list()

  for(ii in 1:K){
    perm.id[[ii]] <- sample(1:length(y), size = length(y), replace = FALSE) # random permutation
  }
  names(perm.id) <- paste("p",1:K, sep="")
  y.perm <- lapply(perm.id, function(a) y[a]) # K permuted outcomes

  # Keep is non null ------
  if(!is.null(keep)){
    # variables indexed by keep are permuted in the same way as y
    # new datasets.
    xbis <- lapply(perm.id, function(a){
      t <- x
      t[,keep] <- x[a,keep]
      return(t)
      })
    perm.names <- names(xbis)

    # lasso with permuted y and xbis
    if(ncore==1){
      #no parallelization
      lasso <- lapply(perm.names, function(x) glmnet(x = xbis[[x]], y = y.perm[[x]], ...))
    }else if(ncore>=1){
      # #parallelisation UNIX
      # lasso <- mclapply(perm.names, function(x) glmnet(x = xbis[[x]],
      #                                                  y = y.perm[[x]], ...),
      #                   mc.cores = ncore, mc.preschedule = FALSE)

      #parallelization all OX
      cl <- makeCluster(ncore)
      registerDoParallel(cl)
      lasso = foreach(i = perm.names, .packages = c("glmnet")) %dopar% {
        glmnet(x = xbis[[i]], y = y.perm[[i]],  ...)
        }
      stopCluster(cl)
      names(lasso) <- perm.names
    }

  # Keep is null (classical) -----
  }else if(is.null(keep)){
    if(ncore==1){
      #no parallelization
      lasso <- lapply(y.perm, glmnet, x= x,... )
    }else if(ncore>=1){
      # #parallelisation UNIX
      # lasso <- mclapply(y.perm, glmnet, x= x, ...,
      #                   mc.cores = ncore, mc.preschedule = FALSE)

      perm.names <- names(y.perm)
      #parallelization all OX
      i <- 0
      cl <- makeCluster(ncore)
      registerDoParallel(cl)
      lasso = foreach(i = perm.names, .packages = c("glmnet")) %dopar% {
        glmnet(x = x, y = y.perm[[i]], ...)
      }
      stopCluster(cl)
      names(lasso) <- perm.names
    }

  }

  # Lambda max in each of the K lasso regressions  -----
  lmax <- unlist(lapply(lasso, function(a) a$lambda[1]))

  # Lamba perm (median) + covariates ----
  lmed <- median(lmax)
  lasso <- glmnet(x = x, y = y, ...) #standard lasso

  closest <- which.min((lasso$lambda-lmed)^2)  #minimum of SCE
  beta <- lasso$beta[,closest]

  if(betaPos){
    var <- names(sort(beta[which(beta>0)], decreasing = TRUE)) # ordre drecroissant
  }else{
    var <- names(sort(abs(beta[which(beta!=0)]), decreasing = TRUE)) # ordre decroissant en terme de valeur absolue
  }


  # Result -----
  res <- list()
  res$beta <- beta
  res$selected_variables <- var

  class(res) = "log.lasso"

  return(res)

}

