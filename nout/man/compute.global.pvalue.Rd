% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_global_functions.R
\name{compute.global.pvalue}
\alias{compute.global.pvalue}
\title{compute.global.pvalue}
\usage{
compute.global.pvalue(
  T.obs,
  m,
  n,
  local_test = "wmw",
  k = NULL,
  stats_G_vector,
  n_perm = 0,
  B = 100,
  seed = 321
)
}
\arguments{
\item{T.obs}{: observed value of the chosen test statistic.}

\item{m}{: calibration sample size.}

\item{n}{: test sample size.}

\item{local_test}{: local test to be used in the closed testing procedure. Default value is Wilcoxon sum-rank test.}

\item{k}{: order of the higher order Wilcoxon sum-rank  test statistic.}

\item{stats_G_vector}{: vector of elementary test statistics to perform the test in Shiraishi (1985).}

\item{n_perm}{: if \eqn{min(m,n)\leq n_perm} the \emph{p}-value for the global null will be computed via permutation. Default value is 0.}

\item{B}{: number of permutations.}

\item{seed}{: seed to ensure reproducible results.}
}
\value{
A number, the \emph{p}-value for the global null according to the chosen test statistic.
}
\description{
It computes the \emph{p}-value for the global null according to the chosen
test statistic. The \emph{p}-value is computed via permutation if either the calibration
sample size or the test sample size is smaller than \code{n_perm}.
Otherwise, it is computed using the asymptotic distribution.
}
