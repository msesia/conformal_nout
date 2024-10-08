% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_g_old.R
\name{asymptotic.critical.G}
\alias{asymptotic.critical.G}
\title{asymptotic.critical.G}
\usage{
asymptotic.critical.G(m, n, stats_G_vector, alpha = 0.1)
}
\arguments{
\item{m}{: calibration size.}

\item{n}{: test size.}

\item{stats_G_vector}{: vector of elementary test statistics to perform the test in Shiraishi (1985).}

\item{alpha}{: significance level.}
}
\value{
A numeric value that is the asymptotic critical value for the Shiraishi test statistic at level \eqn{\alpha}.
}
\description{
asymptotic.critical.G
}
