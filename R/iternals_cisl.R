
sampleImb <- function(y, nB, r, nS0min = NULL, replace = TRUE){
  #Internal function for cisl function
  #
  #Compute a CISL sub sampling for r!=NULL.
  #
  #y : same as cisl
  #nB : same as cisl
  #r : same as cisl
  #nS0min : use to determine number of controls to sample. Equivalent to 4P in the paper
  #replace : same as cisl

  idx1 <- which(y == 1)
  idx0 <- which(y == 0)
  n1 <- length(idx1)

  if (replace) vec1 <- sample(idx1, size = nB * n1, replace = TRUE)
  else vec1 <- rep(idx1, times = nB)
  mat1 <- matrix(vec1, ncol = nB)

  n0 <- r * n1
  if (!is.null(nS0min)) n0 <- max(n0, nS0min)
  n0 <- min(n0, length(idx0))

  if (replace) {
    vec0 <- sample(idx0, nB * n0, replace = TRUE)
    mat0 <- matrix(vec0, ncol = nB)
  } else {
    mat0 <- matrix(ncol = nB, nrow = n0)
    for (b in 1:nB) mat0[, b] <- sample(idx0, n0, replace = FALSE)
  }
  as.data.frame(rbind(mat1, mat0))
}


glmnetSub <- function(idx, x, y, ...){
  #Internal function for cisl function
  #
  #Compute a CISL sub sampling for r!=NULL.
  #
  #idx : indexes of subsampled rows
  #x : x in cisl
  #y : y in cisl
  #... other arguments that can be passed to glmnet  from package glmnet.
  #Here it is dfmax and nlambda

  y <- y[idx]
  x <- x[idx,]
  res <- glmnet(x, y, family = "binomial", standardize = TRUE, ...)
}




probInc <-function(resGlmnet, dfmax, betaPos = TRUE) {
  #Internal function for probStabMean function
  #
  #Compute the CISL quantity
  #
  #resGlmnet : an object of class glmnet
  #dfmax : same as cisl
  #betaPos : same as cisl


  idxDf <- !duplicated(resGlmnet$df)
  if (!missing(dfmax)) idxDf <- idxDf & (resGlmnet$df <= dfmax)
  if (betaPos){
    betaCoef <- resGlmnet$beta[, idxDf, drop = FALSE] > 0
  }else{
    betaCoef <- resGlmnet$beta[, idxDf, drop = FALSE] != 0
  }

  prob <- apply(betaCoef, 1, mean)
  return(prob)
}



probStabMean <- function(obj, dfmax = 50, nCore = 1, betaPos = TRUE){
  #Internal function for cisl function
  #
  #Compute the CISL quantity
  #
  #obj : list of length nB of glmnet objects
  #dfmax : same as cisl
  #nCore : ncore option in cisl
  #betaPos : same as cisl

  prob <- 0
  nB <- length(obj)

  if (nCore > 1){
    # #parallelisation UNIX
    # prob <- mclapply(obj, probInc, mc.cores = nCore, dfmax = dfmax,
    #                  betaPos = betaPos)

    # parallelization all OX
    i <- 0
    cl <- makeCluster(nCore)
    registerDoParallel(cl)
    prob = foreach(i=1:nB, .packages = "glmnet",
                       .export = "probInc") %dopar% {
      probInc(resGlmnet = obj[[i]], dfmax = dfmax,
              betaPos = betaPos)
    }
    stopCluster(cl)
    names(prob) <- names(obj)

  }else{ # no parallelization
    prob <- lapply(obj, probInc, dfmax, betaPos)
  }
  prob <- as.data.frame(prob)
  probMean <- apply(prob, 1, mean)

  res <- vector("list")
  res$prob <- prob
  res$probMean <- probMean
  return(res)
}
