% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/003b_reltest2_libs.R
\name{reltest2_simulate}
\alias{reltest2_simulate}
\title{Random training data from one model}
\usage{
reltest2_simulate(model = "gev", nx = 50, tt, params)
}
\arguments{
\item{model}{which distribution to test. Possibles values are
"gev",
"gpd_k1",
"gev_pred1".}

\item{nx}{the length of the training data.}

\item{tt}{the predictor}

\item{params}{values for the parameters for the specified distribution}
}
\value{
Vector
}
\description{
Random training data from one model
}
