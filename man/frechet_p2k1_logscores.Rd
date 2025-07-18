% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/071b_frechet_p2k1_libs.R
\name{frechet_p2k1_logscores}
\alias{frechet_p2k1_logscores}
\title{Log scores for MLE and RHP predictions calculated using leave-one-out}
\usage{
frechet_p2k1_logscores(logscores, x, t, kloc)
}
\arguments{
\item{logscores}{logical that indicates whether to return leave-one-out estimates estimates of the log-score (much longer runtime)}

\item{x}{a vector of training data values}

\item{t}{a vector or matrix of predictors}

\item{kloc}{the known location parameter}
}
\value{
Two scalars
}
\description{
Log scores for MLE and RHP predictions calculated using leave-one-out
}
