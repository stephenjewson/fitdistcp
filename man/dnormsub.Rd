% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/30b_norm_libs.R
\name{dnormsub}
\alias{dnormsub}
\title{Densities from MLE and RHP}
\usage{
dnormsub(x, y, aderivs = TRUE)
}
\arguments{
\item{x}{a vector of training data values}

\item{y}{a vector of values at which to calculate the density and distribution functions}

\item{aderivs}{logical for whether to use analytic derivatives (instead of numerical)}
}
\value{
A vector of parameter estimates, two pdf vectors, two cdf vectors
}
\description{
Densities from MLE and RHP
}
