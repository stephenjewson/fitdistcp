% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/032b_gnorm_k3_libs.R
\name{gnorm_k3_ldd}
\alias{gnorm_k3_ldd}
\title{Second derivative matrix of the normalized log-likelihood}
\usage{
gnorm_k3_ldd(x, v1, d1, v2, fd2, kbeta)
}
\arguments{
\item{x}{a vector of training data values}

\item{v1}{first parameter}

\item{d1}{the delta used in the numerical derivatives with respect to the parameter}

\item{v2}{second parameter}

\item{fd2}{the fractional delta used in the numerical derivatives with respect to the parameter}

\item{kbeta}{the known beta parameter}
}
\value{
Square scalar matrix
}
\description{
Second derivative matrix of the normalized log-likelihood
}
