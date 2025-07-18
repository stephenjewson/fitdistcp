% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/151b_gev_p12_libs.R
\name{dgev_p12sub}
\alias{dgev_p12sub}
\title{Densities for 5 predictions}
\usage{
dgev_p12sub(
  x,
  t1,
  t2,
  y,
  t01,
  t02,
  ics,
  minxi,
  maxxi,
  debug,
  extramodels = FALSE
)
}
\arguments{
\item{x}{a vector of training data values}

\item{t1}{a vector of predictors for the mean}

\item{t2}{a vector of predictors for the sd}

\item{y}{a vector of values at which to calculate the density and distribution functions}

\item{t01}{a single value of the predictor (specify either \code{t01} or \code{n01} but not both)}

\item{t02}{a single value of the predictor (specify either \code{t02} or \code{n02} but not both)}

\item{ics}{initial conditions for the maximum likelihood search}

\item{minxi}{minimum value of shape parameter xi}

\item{maxxi}{maximum value of shape parameter xi}

\item{debug}{debug flag}

\item{extramodels}{logical that indicates whether to add three additional prediction models}
}
\value{
A vector of parameter estimates, two pdf vectors, two cdf vectors
}
\description{
Densities for 5 predictions
}
