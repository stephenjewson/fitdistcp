% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/063b_lst_p1k3_libs.R
\name{lst_p1k3_lddd}
\alias{lst_p1k3_lddd}
\title{Third derivative tensor of the normalized log-likelihood}
\usage{
lst_p1k3_lddd(x, t, v1, d1, v2, d2, v3, fd3, kdf)
}
\arguments{
\item{x}{a vector of training data values}

\item{t}{a vector or matrix of predictors}

\item{v1}{first parameter}

\item{d1}{the delta used in the numerical derivatives with respect to the parameter}

\item{v2}{second parameter}

\item{d2}{the delta used in the numerical derivatives with respect to the parameter}

\item{v3}{third parameter}

\item{fd3}{the fractional delta used in the numerical derivatives with respect to the parameter}

\item{kdf}{the known degrees of freedom parameter}
}
\value{
Cubic scalar array
}
\description{
Third derivative tensor of the normalized log-likelihood
}
