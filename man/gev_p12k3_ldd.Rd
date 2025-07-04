% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/151b_gev_p12_libs.R
\name{gev_p12k3_ldd}
\alias{gev_p12k3_ldd}
\title{Second derivative matrix of the normalized log-likelihood, with fixed shape parameter}
\usage{
gev_p12k3_ldd(x, t1, t2, v1, d1, v2, d2, v3, d3, v4, d4, v5)
}
\arguments{
\item{x}{a vector of training data values}

\item{t1}{a vector of predictors for the mean}

\item{t2}{a vector of predictors for the sd}

\item{v1}{first parameter}

\item{d1}{the delta used in the numerical derivatives with respect to the parameter}

\item{v2}{second parameter}

\item{d2}{the delta used in the numerical derivatives with respect to the parameter}

\item{v3}{third parameter}

\item{d3}{the delta used in the numerical derivatives with respect to the parameter}

\item{v4}{fourth parameter}

\item{d4}{the delta used in the numerical derivatives with respect to the parameter}

\item{v5}{fifth parameter}
}
\value{
Square scalar matrix
}
\description{
Second derivative matrix of the normalized log-likelihood, with fixed shape parameter
}
