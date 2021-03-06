% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_ps_bic.R
\name{est_ps_bic}
\alias{est_ps_bic}
\title{propensity score estimation in high dimension with automated covariates selection using lasso-bic}
\usage{
est_ps_bic(idx_expo, x, penalty = rep(1, nvars - 1), ...)
}
\arguments{
\item{idx_expo}{Index of the column in \code{x} that corresponds to the
drug covariate for which we aim at estimating the PS.}

\item{x}{Input matrix, of dimension nobs x nvars. Each row is an observation
vector. Can be in sparse matrix format (inherit from class
\code{"sparseMatrix"} as in package \code{Matrix}).}

\item{penalty}{TEST OPTION penalty weights in the variable selection to
include in the PS.}

\item{\dots}{Other arguments that can be passed to \code{glmnet} from package
\code{glmnet} other than \code{penalty.factor},
\code{family}, \code{maxp} and \code{path}.}
}
\value{
An object with S3 class \code{"ps", "bic"}.
\item{expo_name}{Character, name of the drug exposure for which the PS was
estimated. Correspond to \code{colnames(x)[idx_expo]}}.
\item{indicator_expo}{One-column Matrix object.
Indicator of the drug exposure for which the PS was estimated.
Defined by \code{x[, idx_expo]}.}.
\item{score_variables}{Character vector, names of covariates(s) selected
with the lasso-bic approach to include in the PS estimation model.
Could be empty.}
\item{score}{One-column Matrix object, the estimated score.}
}
\description{
Estimate a propensity score to a given drug exposure by
(i) selecting among other drug covariates in \code{x} which ones to
include in the PS estimation model automatically using lasso-bic
approach,
(ii) estimating a score using a classical logistic regression
with the afore selected covariates.
Internal function, not supposed to be used directly.
}
\details{
\code{betaPos} option of \code{lasso_bic} function is set to
\code{FALSE} and \code{maxp} is set to 20.
For optimal storage, the returned elements \code{indicator_expo} and
\code{score} are Matrix with ncol = 1.
}
\examples{

set.seed(15)
drugs <- matrix(rbinom(100*20, 1, 0.2), nrow = 100, ncol = 20)
colnames(drugs) <- paste0("drugs",1:ncol(drugs))
ae <- rbinom(100, 1, 0.3)
psb2 <- est_ps_bic(idx_expo = 2, x = drugs)
psb2$score_variables #selected variables to include in the PS model of drug_2

}
\author{
Emeline Courtois \cr Maintainer: Emeline Courtois
\email{emeline.courtois@inserm.fr}
}
