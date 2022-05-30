#' Summary statistics for main adapt4pv package functions
#'
#' Return the Sensitivity and the False Discovery Rate of an approach
#' implemeted by the main functions of adapt4pv package.
#'
#' @param object An object of class \code{"log.lasso"},
#' \code{"cisl"}, \code{"adaptive"} and \code{"*", "ps","**" } where
#' \code{"*"} is either \code{"adjust"}, \code{"iptw"} or \code{"mw"}
#' and \code{"**"} is either \code{"bic"}, \code{"hdps"} or \code{"xgb"}.
#' @param true_pos Character vector, names of the true positives
#' controls
#' @param q Quantile value for variable selection with
#' an object of class \code{"cisl"}.
#' Possible values are 5, 10, 15, 20. Default is 10
#'
#' @return A data frame wich details for the signal detection method
#' implemented in \code{object}: its number of generated signals, its
#' sensitivity and its false discovery rate.
#' @examples
#'
#' set.seed(15)
#' drugs <- matrix(rbinom(100*20, 1, 0.2), nrow = 100, ncol = 20)
#' colnames(drugs) <- paste0("drugs",1:ncol(drugs))
#' ae <- rbinom(100, 1, 0.3)
#' lcv <- lasso_cv(x = drugs, y = ae, nfolds = 3)
#' summary_stat(object = lcv, true_pos = colnames(drugs)[1:10])
#' # the data are not simulated in such a way that there are true positives
#'
#'
#' @author Emeline Courtois \cr Maintainer: Emeline Courtois
#' \email{emeline.courtois@@inserm.fr}
#'
#' @export summary_stat

summary_stat <- function(object, true_pos, q = 10){

  ## Test class ----
  temp <- class(object)[1]
  if(length(class(object))==1 && inherits(object, "try-error")){
    var <- character(0)
  }else if(length(class(object))==1 && inherits(object,c("log.lasso", "adaptive"))){
    var <- object$selected_variables
  }else if(length(class(object))==3 && temp == "ps"){
    var <- object$selected_variables
  }else if(length(class(object))==1 && inherits(object, "cisl")){
    if(!(q %in% c(5,10,15,20))) stop("invalid q value")
    tmp <- paste0("q",q)
    if(tmp=="q5") tmp <- "q05"
    var <- sort(object[[tmp]][which(object[[tmp]]>0)], decreasing = TRUE)
    var <- names(var)
  }else{
    stop("object does not have a valid class")
  }

  # Test ----
  if(missing(true_pos)) stop("true_pos is not specified")

  # Summary statistics ----
  a <- length(var)
  b <- sum(var %in% true_pos) / length(true_pos) #Sensitivity

  if(a==0){ #pas de signaux detectes, donc FDR = 0
    c <- 0
  }else if(a>0){
    c <-  sum(!(var %in% true_pos)) / a #FDR
  }

  res <- c(a,b,c)
  names(res) <- c("nb.sig", "Sensitivity", "FDR")
  res <- round(res, 2)

  return(res)

}


