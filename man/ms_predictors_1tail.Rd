% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/004c_ms_predictors_1tail.R
\name{ms_predictors_1tail}
\alias{ms_predictors_1tail}
\title{Model Selection Among 5 Distributions with predictors from the \code{fitdistcp} Package}
\usage{
ms_predictors_1tail(x, t)
}
\arguments{
\item{x}{data vector}

\item{t}{predictor vector}
}
\value{
Plots QQ plots to the screen, for each of the 5 models,
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
for 5 one tailed models with predictors in the \code{fitdistcp} package.

The code is straightforward, and the point is to illustrate what is
possible using the model selection outputs from the \code{fitdistcp} routines.

The input data may be automatically shifted so that the minimum value is positive.

For the Pareto, the data is so that the minimum value is slightly greater than 1.
}
\details{
The 5 models are:
\code{exp_p1},
\code{pareto_p1k2},
\code{lnorm_p1},
\code{frechet_p2k1},
\code{weibull_p2}.
}
\examples{
\donttest{ # because it's too slow for CRAN
set.seed(3)
nx=100
predictor=c(1:nx)/nx
x=rlnorm(nx,meanlog=predictor,sdlog=0.1)
print(ms_predictors_1tail(x,predictor))
}



}
\author{
Stephen Jewson \email{stephen.jewson@gmail.com}
}
