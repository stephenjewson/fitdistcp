% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/074c_gev_p12k3_derivs.R
\name{gev_p12k3_lddda}
\alias{gev_p12k3_lddda}
\title{The third derivative of the normalized log-likelihood}
\usage{
gev_p12k3_lddda(x, t1, t2, v1, v2, v3, v4, kshape)
}
\arguments{
\item{x}{a vector of training data values}

\item{t1}{a vector of predictors for the mean}

\item{t2}{a vector of predictors for the sd}

\item{v1}{first parameter}

\item{v2}{second parameter}

\item{v3}{third parameter}

\item{v4}{fourth parameter}

\item{kshape}{the known shape parameter}
}
\value{
3d array
}
\description{
The third derivative of the normalized log-likelihood
}
