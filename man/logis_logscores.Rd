% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/040b_logis_libs.R
\name{logis_logscores}
\alias{logis_logscores}
\title{Log scores for MLE and RHP predictions calculated using leave-one-out}
\usage{
logis_logscores(logscores, x)
}
\arguments{
\item{logscores}{logical that indicates whether to return leave-one-out estimates estimates of the log-score (much longer runtime)}

\item{x}{a vector of training data values}
}
\value{
Two scalars
}
\description{
Log scores for MLE and RHP predictions calculated using leave-one-out
}
