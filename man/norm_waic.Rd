% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/30b_norm_libs.R
\name{norm_waic}
\alias{norm_waic}
\title{Waic}
\usage{
norm_waic(waicscores, x, v1hat, d1, v2hat, fd2, aderivs)
}
\arguments{
\item{waicscores}{logical that indicates whether to return estimates for the waic1 and waic2 scores (longer runtime)}

\item{x}{a vector of training data values}

\item{v1hat}{first parameter}

\item{d1}{the delta used in the numerical derivatives with respect to the parameter}

\item{v2hat}{second parameter}

\item{fd2}{the fractional delta used in the numerical derivatives with respect to the parameter}

\item{aderivs}{logical for whether to use analytic derivatives (instead of numerical)}
}
\value{
Two numeric values.
}
\description{
Waic
}
