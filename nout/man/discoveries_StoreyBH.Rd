% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/discoveries_StoreyBH.R
\name{discoveries_StoreyBH}
\alias{discoveries_StoreyBH}
\title{discoveries_StoreyBH}
\usage{
discoveries_StoreyBH(S_X, S_Y, alpha = 0.1, lambda = 0.5)
}
\arguments{
\item{S_X}{: score vector for the calibration set}

\item{S_Y}{: score vector for the test set}

\item{alpha}{: significance level. Default level is set equal to 0.1}

\item{lambda}{: parameter involved in the computation of Storey estimator. Default value is set equal to 0.5}
}
\value{
A vector corresponding to indices of test observations which are rejected
by Benjamini-Hochberg procedure applied to conformal \emph{p}-values and using Storey
estimator for the proportion of true null hypoteses.
The selection set is trivial, i.e., we are interested in testing all the observations in the test set by default.
Storey estimator is computed as
\deqn{\hat\pi_0 = \frac{1+\sum_{i=1}^n \mathbb{1}\{p_i>\lambda\}}{m(1-\lambda)}}
where \eqn{\lambda\in(0,1)} is a tuning parameter, \eqn{m} and \eqn{n} are the size of the calibration set and of the test set
respectively and \eqn{p_i} is the \emph{p}-value related to the \emph{i}th observation in the test set.
}
\description{
Applies Benjamini-Hochberg procedure to conformal \emph{p}-values using
Storey estimator for the proportion of true null hypotheses.
}
\examples{
set.seed(321)
Sxy = sample(x=1:1000, size=100)
Sx = sample(Sxy, size=70)
Sy = setdiff(Sxy, Sx)
discoveries_StoreyBH(Sx, Sy)


}
