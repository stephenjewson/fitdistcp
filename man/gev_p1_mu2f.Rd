% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/150b_gev_p1_libs.R
\name{gev_p1_mu2f}
\alias{gev_p1_mu2f}
\title{GEVD-with-p1: DMGS equation 3.3 mu2 term}
\usage{
gev_p1_mu2f(alpha, t0, v1, d1, v2, d2, v3, fd3, v4, d4)
}
\arguments{
\item{alpha}{a vector of values of alpha (one minus probability)}

\item{t0}{a single value of the predictor (specify either \code{t0} or \code{n0} but not both)}

\item{v1}{first parameter}

\item{d1}{the delta used in the numerical derivatives with respect to the parameter}

\item{v2}{second parameter}

\item{d2}{the delta used in the numerical derivatives with respect to the parameter}

\item{v3}{third parameter}

\item{fd3}{the fractional delta used in the numerical derivatives with respect to the parameter}

\item{v4}{fourth parameter}

\item{d4}{the delta used in the numerical derivatives with respect to the parameter}
}
\value{
3d array
}
\description{
GEVD-with-p1: DMGS equation 3.3 mu2 term
}
