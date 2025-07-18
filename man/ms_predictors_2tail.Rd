% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/004d_ms_predictors_2tail.R
\name{ms_predictors_2tail}
\alias{ms_predictors_2tail}
\title{Model Selection Among 6 Distributions with predictors from the \code{fitdistcp} Package}
\usage{
ms_predictors_2tail(x, t)
}
\arguments{
\item{x}{data vector}

\item{t}{predictor vector}
}
\value{
Plots QQ plots to the screen, for each of the 6 models,
and returns a data frame containing
\itemize{
\item AIC scores, AIC weights
\item WAIC1 scores, WAIC1 weights
\item WAIC2 scores, WAIC2 weights
\item logscores and logscore weights
}
}
\description{
Applies model selection using AIC, WAIC1, WAIC2 and leave-one-out logscore
to the input data \eqn{x,t},
for 6 two tail models with predictors in the \code{fitdistcp} packages
(although for the GEV, the logscore is NA for mathematical reasons).

The code is straightforward, and the point is to illustrate what is
possible using the model selection outputs from the \code{fitdistcp} routines.

GEVD is temperamental in that
it doesn't work if the shape parameter is extreme.
}
\details{
The 11 models are:
\code{norm_p1},
\code{gumbel_p1},
\code{logis_p1},
\code{lst_k3_p1},
\code{cauchy_p1} and
\code{gev_p1}.
}
\examples{
\donttest{ # because it's too slow for CRAN
set.seed(2)
nx=100
predictor=c(1:nx)/nx
x=rnorm(nx,mean=predictor,sd=1)
print(ms_predictors_2tail(x,predictor))
}



}
\author{
Stephen Jewson \email{stephen.jewson@gmail.com}
}
