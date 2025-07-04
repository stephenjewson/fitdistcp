% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/73b_weibull_p2_libs.R
\name{weibull_p2_waic}
\alias{weibull_p2_waic}
\title{Waic}
\usage{
weibull_p2_waic(
  waicscores,
  x,
  t,
  v1hat,
  fd1,
  v2hat,
  d2,
  v3hat,
  d3,
  lddi,
  lddd,
  lambdad,
  aderivs = TRUE
)
}
\arguments{
\item{waicscores}{logical that indicates whether to return estimates for the waic1 and waic2 scores (longer runtime)}

\item{x}{a vector of training data values}

\item{t}{a vector or matrix of predictors}

\item{v1hat}{first parameter}

\item{fd1}{the fractional delta used in the numerical derivatives with respect to the parameter}

\item{v2hat}{second parameter}

\item{d2}{the delta used in the numerical derivatives with respect to the parameter}

\item{v3hat}{third parameter}

\item{d3}{the delta used in the numerical derivatives with respect to the parameter}

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
