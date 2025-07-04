% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/56b_pareto_p1k2_libs.R
\name{pareto_p1k2_lmnp}
\alias{pareto_p1k2_lmnp}
\title{One component of the second derivative of the normalized log-likelihood}
\usage{
pareto_p1k2_lmnp(x, t, v1, d1, v2, d2, kscale, mm, nn, rr)
}
\arguments{
\item{x}{a vector of training data values}

\item{t}{a vector or matrix of predictors}

\item{v1}{first parameter}

\item{d1}{the delta used in the numerical derivatives with respect to the parameter}

\item{v2}{second parameter}

\item{d2}{the delta used in the numerical derivatives with respect to the parameter}

\item{kscale}{the known scale parameter}

\item{mm}{an index for which derivative to calculate}

\item{nn}{an index for which derivative to calculate}

\item{rr}{an index for which derivative to calculate}
}
\value{
Scalar value
}
\description{
One component of the second derivative of the normalized log-likelihood
}
