% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04b_modelselection_predictors.R
\name{modelselection_predictors}
\alias{modelselection_predictors}
\title{Model Selection Among 6 Distributions with predictors from the \code{fitdistcp} Package}
\usage{
modelselection_predictors(x, t)
}
\arguments{
\item{x}{data vector}

\item{t}{predictor vector}
}
\value{
Plots QQ plots to the screen, for each of the 11 models,
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
for 11 models with predictors in the \code{fitdistcp} packages
(although for the GEV and GPD, the logscore is NA).

The code is straightforward, and the point is to illustrate what is
possible using the outputs from the \code{fitdistcp} routines.

The input data is automatically shifted so that the minimum value is positive,
to avoid errors from the models that require positive random variable.

For the Pareto, the data is further shifted so that the minimum value is slightly greater than 1.

GEVD is temperamental in that
it doesn't work if the shape parameter is extreme.
}
\details{
The 11 models are:
\code{exp_p1}
\code{pareto_p1k2}
\code{norm_p1},
\code{lnorm_p1},
\code{gumbel_p1},
\code{frechet_p2k1},
\code{weibull_p2},
\code{logis_p1},
\code{lst_k3_p1},
\code{cauchy_p1} and
\code{gev_p1}.
}
\examples{
set.seed(2)
nx=100
predictor=c(1:nx)/nx
x=rlnorm(nx,meanlog=predictor,sdlog=0.1)
print(modelselection_predictors(x,predictor))



}
\author{
Stephen Jewson \email{stephen.jewson@gmail.com}
}
