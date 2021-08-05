
control_adaCV <- function(dat, id.folds, adaptive.weights){
  #Internal function for adapt_cv function
  #
  #Check if adapt-cv could be implemented on a given fold repartition. Return 0 if yes, return 1 if
  #on at least one fold, all the the non-penalyzed covariates are null (absent).
  #
  #dat : data matrix
  #id.folds : equivalent to foldid
  #adaptive.weights : adaptive weights to integrete to the adaptive lasso

  varCV <- names(which(adaptive.weights!=Inf)) # non excluded covariates

  top <- list()
  for(tt in 1:max(id.folds)){
    top[[tt]] <- sum(colSums(dat[which(id.folds==tt),varCV, drop = FALSE]))
    #on folds, colsums of non penalysed covariates
  }

  top <- unlist(top)
  # if at least one element of top is null: there is at least one fold where all non-penalyzed covariates are null on fold's individuals
  top <- sum(top == 0)
  return(top)
}


contingence_all <- function(data_cases, data_controls){
  #Internal function for adapt_univ function
  #
  #For a given dataset X and a outcome Y, after split X in two subdatasets, compute contigency table
  #for all the covariates in X with Y.
  #Expected for binary X and Y
  #Return a list
  #
  #data_cases : X when Y=1
  #data_controls : X when Y=0

  #test si on est bien sur les memes expos entre les cas et les temoins
  if(!all.equal(colnames(data_cases), colnames(data_controls))){
    stop("colnames differ from data_controls to data_case")
  }

  nb.expo <- ncol(data_cases)
  expo.name <- colnames(data_cases)

  cas_expose <- Matrix::colSums(data_cases)
  cas_non_expose <- nrow(data_cases) - cas_expose

  tem_expose <- Matrix::colSums(data_controls)
  tem_non_expose <- nrow(data_controls) - tem_expose

  res <- list()

  for(ii in 1:nb.expo){
    res[[ii]] <- matrix(data = c(cas_expose[ii], cas_non_expose[ii], tem_expose[ii], tem_non_expose[ii]), nrow = 2, ncol = 2)
    row.names(res[[ii]]) <- c("E+","E-")
    colnames(res[[ii]]) <- c("M+","M-")
  }

  names(res) <- expo.name

  return(res)
}



poidsCISL_v2 <- function(distrib){
  #Internal function for adapt_cisl function
  #
  #Compute apdaptive penalty weights derived from the quantity computed in CISL
  #
  #distrib : one row of the matrice return by cisl function. Correspond to the
  # distribution of the CISL quantity for one covariate. Length equal to nB.


  res <- length(which(distrib==0))/length(distrib)

  if(res==0){ # quantity always positive : apw no null
    res <- 1/length(distrib)
  }else if(res == 1){  # quantity always null : apw infty
    res <- Inf
  }
  return(res)

}


hdps_pv <- function(E, C, D, k){
  #internal function for ps_hdps
  #
  #Compute the step 4 (prioritize covariates) of hdps in the pharmacovigilance context
  #
  #E : exposure indicator
  #C : confonders, Matrix
  #D : outcome indicator
  #k : number of covariates to keep in the PS estimation model

  Tot_C <- colSums(C)
  Tot_E <- sum(E)
  Tot_D <- sum(D)
  Tot <- nrow(C)

  #Pc1 and Pc0
  a <- E %*% C
  PC1 <- t(a/Tot)
  PC0 <- t((Tot_C - a)/(Tot - Tot_E))

  #RRCD
  e <- D %*% C
  RRCD <- t( (e/Tot_C) / ( (Tot_D - e)/(Tot - Tot_C) )  )

  #Bias M - correction
  BiasM <- (PC1 *(RRCD -1) +1) / (PC0 * (RRCD -1)+1)
  BiasM <- abs(log(BiasM))

  #Ordering
  res <- row.names(BiasM)[order(BiasM, decreasing = TRUE)][1:k]

  return(res)

}



predict_speedglm.wfit <- function(speedglm, newmatrix){
  #Internal function for summary score estimation
  #
  #Make a prediction for a speedglm object were newmatrix is not a data.frame but a sparse matrix.
  #Return a numeric vector of lenght nrows(newmatrix), correspond to prediction of class "response"
  #
  #speedglm : an object of class "speedglm"
  #newmatrix a matrix (could be sparse) with new data or the original data

  newmatrix<-cbind(1,newmatrix)
  eta <- newmatrix %*% speedglm$coefficients
  eta <- speedglm$family$linkinv(as.vector(eta))

  return(eta)
}


try_bis <- function(x, classe = "numeric"){
  #Internal function to retrieve p-value element in ps analysis
  #
  #Use to extract the pval_1sided or the pval_2sided element
  #from an object of class adjust, mw or iptw
  #
  #
  #x : step 3 object in the ps analysis
  #classe : class of the desired element

  res <- try(x)
  if(class(res)!=classe) res <- NA
  return(res)
}



standardize <- function(X){
  #Internal function for standardization
  #
  #Given a data matrix X, perform standardization (different from scale, use the biase estimation of the standard deviation)
  #
  #
  #X: data matrix, could be sparse
  #Return a list xith elements
  #Xs : X matrixstandardized
  #X_Mean and X_sd : mean and sd of X's columns

  X_centered <- apply(X, 2, function(x) x - mean(x))
  Xs <- apply(X_centered, 2, function(x) x/sqrt(sum(x^2)/nrow(X)))

  X_mean <- colMeans(X)
  X_sd <- apply(X_centered, 2, function(x) sqrt(sum(x^2) / nrow(X)))

  return(list(Xs = Xs, X_mean = X_mean, X_sd = X_sd))
}


cv.lognet.short <- function(pred.matrix, yy){
  #Internal function for prediction error for CV
  #
  #Given a the binary outcome yy and the prediction pred.matrix (one column corresponding to a lambda value) return the deviance
  #Code from cv.lognet from glmnet package

  prob_min = 1e-05
  prob_max = 1 - prob_min
  nc = dim(yy)
  if (is.null(nc)) {
    yy = as.factor(yy)
    ntab = table(yy)
    nc = as.integer(length(ntab))
    yy = diag(nc)[as.numeric(yy), ]
  }
  N = nrow(yy)
  nlambda=ncol(pred.matrix)
  ywt = apply(yy, 1, sum)
  yy = yy/ywt
  N = nrow(yy) - apply(is.na(pred.matrix), 2, sum)
  # mse = (yy[, 1] - (1 - pred.matrix))^2 + (yy[, 2] - pred.matrix)^2
  # mae = abs(yy[, 1] - (1 - pred.matrix)) + abs(yy[, 2] - pred.matrix)
  deviance = {pred.matrix = pmin(pmax(pred.matrix, prob_min), prob_max)
  lp = yy[, 1] * log(1 - pred.matrix) + yy[, 2] * log(pred.matrix)
  ly = log(yy)
  ly[yy == 0] = 0
  ly = drop((yy * ly) %*% c(1, 1))
  2 * (ly - lp)
  }
  # class = yy[, 1] * (pred.matrix > 0.5) + yy[, 2] * (pred.matrix <= 0.5)

  return(deviance)

}

