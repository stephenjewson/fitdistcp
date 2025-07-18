% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/152b_gev_p123_libs.R
\name{gev_p123_means}
\alias{gev_p123_means}
\title{Analytical expressions for Predictive Means
RHP mean based on the expectation of DMGS equation 2.1}
\usage{
gev_p123_means(means, t01, t02, t03, ml_params, nx)
}
\arguments{
\item{means}{logical that indicates whether to return analytical estimates for the distribution means (longer runtime)}

\item{t01}{a single value of the predictor (specify either \code{t01} or \code{n01} but not both)}

\item{t02}{a single value of the predictor (specify either \code{t02} or \code{n02} but not both)}

\item{t03}{a single value of the predictor (specify either \code{t03} or \code{n03} but not both)}

\item{ml_params}{parameters}

\item{nx}{length of training data}
}
\value{
Two scalars
}
\description{
Analytical expressions for Predictive Means
RHP mean based on the expectation of DMGS equation 2.1
}
