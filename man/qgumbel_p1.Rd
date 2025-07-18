% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/070b_gumbel_p1_libs.R
\name{qgumbel_p1}
\alias{qgumbel_p1}
\title{Gumbel-with-p1 quantile function}
\usage{
qgumbel_p1(p, t0, ymn, slope, sigma)
}
\arguments{
\item{p}{a vector of probabilities at which to generate predictive quantiles}

\item{t0}{a single value of the predictor (specify either \code{t0} or \code{n0} but not both)}

\item{ymn}{the location parameter of the function of the predictor}

\item{slope}{the slope of the function of the predictor}

\item{sigma}{the sigma parameter of the distribution}
}
\value{
Vector
}
\description{
Gumbel-with-p1 quantile function
}
