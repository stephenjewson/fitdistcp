% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/151b_gev_p12_libs.R
\name{dgev_p12}
\alias{dgev_p12}
\title{GEVD-with-p1: Density function}
\usage{
dgev_p12(x, t1, t2, ymn, slope, sigma1, sigma2, xi, log = FALSE)
}
\arguments{
\item{x}{a vector of training data values}

\item{t1}{a vector of predictors for the mean}

\item{t2}{a vector of predictors for the sd}

\item{ymn}{the location parameter of the function of the predictor}

\item{slope}{the slope of the function of the predictor}

\item{sigma1}{first coefficient for the sigma parameter of the distribution}

\item{sigma2}{second coefficient for the sigma parameter of the distribution}

\item{xi}{the shape parameter of the distribution}

\item{log}{logical for the density evaluation}
}
\value{
Vector
}
\description{
GEVD-with-p1: Density function
}
