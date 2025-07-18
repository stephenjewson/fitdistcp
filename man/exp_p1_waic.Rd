% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/055b_exp_p1_libs.R
\name{exp_p1_waic}
\alias{exp_p1_waic}
\title{Waic}
\usage{
exp_p1_waic(waicscores, x, t, v1hat, v2hat, lddi, lddd, lambdad)
}
\arguments{
\item{waicscores}{logical that indicates whether to return estimates for the waic1 and waic2 scores (longer runtime)}

\item{x}{a vector of training data values}

\item{t}{a vector or matrix of predictors}

\item{v1hat}{first parameter}

\item{v2hat}{second parameter}

\item{lddi}{inverse observed information matrix}

\item{lddd}{third derivative of log-likelihood}

\item{lambdad}{derivative of the log prior}
}
\value{
Two numeric values.
}
\description{
Waic
}
