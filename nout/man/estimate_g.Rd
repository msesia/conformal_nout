% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_g.R
\name{estimate_g}
\alias{estimate_g}
\title{estimate_g}
\usage{
estimate_g(
  scores_reference,
  scores_pooled,
  method = "betamix",
  monotone = FALSE
)
}
\arguments{
\item{scores_reference}{: calibration score vector of inliers.}

\item{scores_pooled}{: pooled score vector of calibration and test scores.}

\item{method}{: character value indicating the method to be used to estimate the PDF.
It can be either "betamix" or "mixmodel".}

\item{monotone}{: logical value. If \code{TRUE} the estimated PDF via the specified method
is made monotone. The direction (increasing or decreasing) is automatically chosen from the data.}
}
\value{
A list of three elements which are the estimated PDF and CDF and the monotonicity
used in the estimation process.
}
\description{
Given a two-component mixture model where the first distribution is a Standard Uniform
and the second one is a Beta distribution, it estimates the parameter of the Beta distribution
and returns the PDF and the CDF.
}
