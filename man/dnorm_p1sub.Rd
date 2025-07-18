% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/060b_norm_p1_libs.R
\name{dnorm_p1sub}
\alias{dnorm_p1sub}
\title{Densities from MLE and RHP}
\usage{
dnorm_p1sub(x, t, y, t0)
}
\arguments{
\item{x}{a vector of training data values}

\item{t}{a vector or matrix of predictors}

\item{y}{a vector of values at which to calculate the density and distribution functions}

\item{t0}{a single value of the predictor (specify either \code{t0} or \code{n0} but not both)}
}
\value{
A vector of parameter estimates, two pdf vectors, two cdf vectors
}
\description{
Densities from MLE and RHP
}
