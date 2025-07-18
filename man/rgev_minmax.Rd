% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/110b_gev_libs.R
\name{rgev_minmax}
\alias{rgev_minmax}
\title{rgev but with maxlik xi guaranteed within bounds}
\usage{
rgev_minmax(nx, mu, sigma, xi, minxi = -1, maxxi = 1)
}
\arguments{
\item{nx}{length of training data}

\item{mu}{the location parameter of the distribution}

\item{sigma}{the sigma parameter of the distribution}

\item{xi}{the shape parameter of the distribution}

\item{minxi}{minimum value of shape parameter xi}

\item{maxxi}{maximum value of shape parameter xi}
}
\value{
Vector
}
\description{
rgev but with maxlik xi guaranteed within bounds
}
