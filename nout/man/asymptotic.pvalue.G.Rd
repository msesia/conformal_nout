% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_g_old.R
\name{asymptotic.pvalue.G}
\alias{asymptotic.pvalue.G}
\title{asymptotic.pvalue.G}
\usage{
asymptotic.pvalue.G(m, n, stats_G_vector, T.obs)
}
\arguments{
\item{m}{: calibration size.}

\item{n}{: test size.}

\item{stats_G_vector}{: vector of elementary test statistics to perform the test in Shiraishi (1985).}

\item{T.obs}{: observed value of the test statistic.}
}
\value{
A numeric value that is the asymptotic \emph{p}-value corresponding to the global null.
}
\description{
asymptotic.pvalue.G
}
