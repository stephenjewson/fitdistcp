% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/073b_weibull_p2_libs.R
\name{dweibull_p2}
\alias{dweibull_p2}
\title{Weibull-with-p1 density function}
\usage{
dweibull_p2(x, t0, shape, ymn, slope, log = FALSE)
}
\arguments{
\item{x}{a vector of training data values}

\item{t0}{a single value of the predictor (specify either \code{t0} or \code{n0} but not both)}

\item{shape}{the shape parameter of the distribution}

\item{ymn}{the location parameter of the function of the predictor}

\item{slope}{the slope of the function of the predictor}

\item{log}{logical for the density evaluation}
}
\value{
Vector
}
\description{
Weibull-with-p1 density function
}
