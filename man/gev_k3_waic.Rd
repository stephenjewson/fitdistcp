% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/53b_gev_k3_libs.R
\name{gev_k3_waic}
\alias{gev_k3_waic}
\title{Waic}
\usage{
gev_k3_waic(
  waicscores,
  x,
  v1hat,
  d1,
  v2hat,
  fd2,
  kshape,
  lddi,
  lddd,
  lambdad,
  aderivs
)
}
\arguments{
\item{waicscores}{logical that indicates whether to return estimates for the waic1 and waic2 scores (longer runtime)}

\item{x}{a vector of training data values}

\item{v1hat}{first parameter}

\item{d1}{the delta used in the numerical derivatives with respect to the parameter}

\item{v2hat}{second parameter}

\item{fd2}{the fractional delta used in the numerical derivatives with respect to the parameter}

\item{kshape}{the known shape parameter}

\item{lddi}{inverse observed information matrix}

\item{lddd}{third derivative of log-likelihood}

\item{lambdad}{derivative of the log prior}

\item{aderivs}{logical for whether to use analytic derivatives (instead of numerical)}
}
\value{
Two numeric values.
}
\description{
Waic
}
