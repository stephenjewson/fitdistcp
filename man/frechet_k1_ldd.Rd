% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/51b_frechet_k1_libs.R
\name{frechet_k1_ldd}
\alias{frechet_k1_ldd}
\title{Second derivative matrix of the normalized log-likelihood}
\usage{
frechet_k1_ldd(x, v1, fd1, v2, fd2, kloc)
}
\arguments{
\item{x}{a vector of training data values}

\item{v1}{first parameter}

\item{fd1}{the fractional delta used in the numerical derivatives with respect to the parameter}

\item{v2}{second parameter}

\item{fd2}{the fractional delta used in the numerical derivatives with respect to the parameter}

\item{kloc}{the known location parameter}
}
\value{
Square scalar matrix
}
\description{
Second derivative matrix of the normalized log-likelihood
}
