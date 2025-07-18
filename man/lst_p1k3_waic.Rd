% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/063b_lst_p1k3_libs.R
\name{lst_p1k3_waic}
\alias{lst_p1k3_waic}
\title{Waic}
\usage{
lst_p1k3_waic(
  waicscores,
  x,
  t,
  v1hat,
  d1,
  v2hat,
  d2,
  v3hat,
  fd3,
  kdf,
  lddi,
  lddd,
  lambdad,
  aderivs
)
}
\arguments{
\item{waicscores}{logical that indicates whether to return estimates for the waic1 and waic2 scores (longer runtime)}

\item{x}{a vector of training data values}

\item{t}{a vector or matrix of predictors}

\item{v1hat}{first parameter}

\item{d1}{the delta used in the numerical derivatives with respect to the parameter}

\item{v2hat}{second parameter}

\item{d2}{the delta used in the numerical derivatives with respect to the parameter}

\item{v3hat}{third parameter}

\item{fd3}{the fractional delta used in the numerical derivatives with respect to the parameter}

\item{kdf}{the known degrees of freedom parameter}

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
