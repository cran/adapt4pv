% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adapt_bic.R
\name{adapt_bic}
\alias{adapt_bic}
\title{fit an adaptive lasso with adaptive weights derived from lasso-bic}
\usage{
adapt_bic(x, y, gamma = 1, maxp = 50, path = TRUE, betaPos = TRUE, ...)
}
\arguments{
\item{x}{Input matrix, of dimension nobs x nvars. Each row is an observation
vector. Can be in sparse matrix format (inherit from class
\code{"sparseMatrix"} as in package \code{Matrix}).}

\item{y}{Binary response variable, numeric.}

\item{gamma}{Tunning parameter to defined the penalty weights. See details below.
Default is set to 1.}

\item{maxp}{A limit on how many relaxed coefficients are allowed.
Default is 50, in \code{glmnet} option default is 'n-3', where 'n' is the sample size.}

\item{path}{Since \code{glmnet} does not do stepsize optimization, the Newton
algorithm can get stuck and not converge, especially with relaxed fits.
With \code{path=TRUE}, each relaxed fit on a particular set of variables
is computed pathwise using the original sequence of lambda values
(with a zero attached to the end). Default is \code{path=TRUE}.}

\item{betaPos}{Should the covariates selected by the procedure be
positively associated with the outcome ? Default is \code{TRUE}.}

\item{\dots}{Other arguments that can be passed to \code{glmnet} from package
\code{glmnet} other than \code{penalty.factor},
\code{family}, \code{maxp} and \code{path}.}
}
\value{
An object with S3 class \code{"adaptive"}.
\item{aws}{Numeric vector of penalty weights derived from lasso-bic.
Length equal to nvars.}
\item{criterion}{Character, indicates which criterion is used with the
adaptive lasso for variable selection. For \code{adapt_bic} function, \code{criterion}
is "bic".}
\item{beta}{Numeric vector of regression coefficients in the adaptive lasso.
If \code{criterion} = "cv" the regression coefficients are PENALIZED, if
\code{criterion} = "bic" the regression coefficients are UNPENALIZED.
Length equal to nvars. Could be NA if adaptive weights are all equal to infinity.}
\item{selected_variables}{Character vector, names of variable(s) selected
with this adaptive approach.
If \code{betaPos = TRUE}, this set is the covariates with a positive regression
coefficient in \code{beta}.
Else this set is the covariates with a non null regression coefficient in \code{beta}.
Covariates are ordering according to the p-values (two-sided if \code{betaPos = FALSE} ,
one-sided if \code{betaPos = TRUE}) in the classical multiple logistic regression
model that minimzes the BIC in the adaptive lasso.}
}
\description{
Fit a first lasso regression and use Bayesian Information Criterion to determine `
adaptive weights (see \code{lasso_bic} function for more details),
then run an adaptive lasso with this penalty weighting.
BIC is used for the adaptive lasso for variable selection.
Can deal with very large sparse data matrices.
Intended for binary reponse only (option \code{family = "binomial"} is forced).
Depends on the \code{glmnet} and \code{relax.glmnet} function from the package
\code{glmnet}.
}
\details{
The adaptive weight for a given covariate i is defined by
\deqn{w_i = 1/|\beta^{BIC}_i|^\gamma} where
 \eqn{\beta^{BIC}_i} is the NON PENALIZED regression coefficient
 associated to covariate \eqn{i} obtained with lasso-bic.
}
\examples{

set.seed(15)
drugs <- matrix(rbinom(100*20, 1, 0.2), nrow = 100, ncol = 20)
colnames(drugs) <- paste0("drugs",1:ncol(drugs))
ae <- rbinom(100, 1, 0.3)
ab <- adapt_bic(x = drugs, y = ae, maxp = 50)


}
\author{
Emeline Courtois \cr Maintainer: Emeline Courtois
\email{emeline.courtois@inserm.fr}
}
