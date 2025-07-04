% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/53b_gev_k3_libs.R
\name{dgev_k3sub}
\alias{dgev_k3sub}
\title{Densities from MLE and RHP}
\usage{
dgev_k3sub(x, y, d1 = 0.01, fd2 = 0.01, kshape, aderivs = TRUE)
}
\arguments{
\item{x}{a vector of training data values}

\item{y}{a vector of values at which to calculate the density and distribution functions}

\item{d1}{the delta used in the numerical derivatives with respect to the parameter}

\item{fd2}{the fractional delta used in the numerical derivatives with respect to the parameter}

\item{kshape}{the known shape parameter}

\item{aderivs}{logical for whether to use analytic derivatives (instead of numerical)}
}
\value{
A vector of parameter estimates, two pdf vectors, two cdf vectors
}
\description{
Densities from MLE and RHP
}
