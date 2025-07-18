% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/300_ppm_temporary.R
\name{qgev_p1_ppm}
\alias{qgev_p1_ppm}
\title{Temporary dummy for one of the ppm models}
\usage{
qgev_p1_ppm(x, t, n0, p)
}
\arguments{
\item{x}{a vector of training data values}

\item{t}{a vector of predictors, such that \code{length(t)=length(x)}}

\item{n0}{an index for the predictor
(specify either \code{t0} or \code{n0} but not both)}

\item{p}{a vector of probabilities at which to generate predictive quantiles}
}
\value{
\code{q****} returns a list containing at least the following:

\itemize{
\item \code{ml_params:} maximum likelihood estimates for the parameters.
\item \code{ml_value:} the value of the log-likelihood at the maximum.
\item \code{standard_errors:} estimates of the standard errors on the parameters,
from the inverse observed information matrix.
\item \code{ml_quantiles:} quantiles calculated using maximum likelihood.
\item \code{cp_quantiles:} predictive quantiles calculated using a calibrating prior.
\item \code{maic:} the AIC score for the maximum likelihood model, times -1/2.
\item \code{cp_method:} a comment about the method used to generate the \code{cp} prediction.
}

For models with predictors, \code{q****} additionally returns:
\itemize{
\item \code{predictedparameter:} the estimated value for parameter,
as a function of the predictor.
\item \code{adjustedx:} the detrended values of \code{x}
}

\code{r****} returns a list containing the following:

\itemize{
\item \code{ml_params:} maximum likelihood estimates for the parameters.
\item \code{ml_deviates:} random deviates calculated using maximum likelihood.
\item \code{cp_deviates:} predictive random deviates calculated using a calibrating prior.
\item \code{cp_method:} a comment about the method used to generate the \code{cp} prediction.
}

\code{d****} returns a list containing the following:

\itemize{
\item \code{ml_params:} maximum likelihood estimates for the parameters.
\item \code{ml_pdf:} density function from maximum likelihood.
\item \code{cp_pdf:} predictive density function calculated using a calibrating prior
(not available in EVT routines, for mathematical reasons, unless using RUST).
\item \code{cp_method:} a comment about the method used to generate the \code{cp} prediction.
}

\code{p***} returns a list containing the following:

\itemize{
\item \code{ml_params:} maximum likelihood estimates for the parameters.
\item \code{ml_cdf:} distribution function from maximum likelihood.
\item \code{cp_cdf:} predictive distribution function  calculated using a calibrating prior
(not available in EVT routines, for mathematical reasons, unless using RUST).
\item \code{cp_method:} a comment about the method used to generate the \code{cp} prediction.
}

\code{t***} returns a list containing the following:

\itemize{
\item \code{theta_samples:} random samples from the parameter posterior.
}
}
\description{
Temporary dummy for one of the ppm models
}
