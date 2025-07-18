% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/151a_gev_p12.R
\name{gev_p12_cp}
\alias{gev_p12_cp}
\alias{qgev_p12_cp}
\alias{rgev_p12_cp}
\alias{dgev_p12_cp}
\alias{pgev_p12_cp}
\alias{tgev_p12_cp}
\title{Generalized Extreme Value Distribution with Two Predictors, Predictions based on a Calibrating Prior}
\usage{
qgev_p12_cp(
  x,
  t1,
  t2,
  t01 = NA,
  t02 = NA,
  n01 = NA,
  n02 = NA,
  p = seq(0.1, 0.9, 0.1),
  ics = c(0, 0, 0, 0, 0),
  fdalpha = 0.01,
  minxi = -1,
  maxxi = 1,
  means = FALSE,
  waicscores = FALSE,
  extramodels = FALSE,
  pdf = FALSE,
  dmgs = TRUE,
  rust = FALSE,
  nrust = 1e+05,
  predictordata = TRUE,
  centering = TRUE,
  debug = FALSE
)

rgev_p12_cp(
  n,
  x,
  t1,
  t2,
  t01 = NA,
  t02 = NA,
  n01 = NA,
  n02 = NA,
  ics = c(0, 0, 0, 0, 0),
  minxi = -1,
  maxxi = 1,
  extramodels = FALSE,
  rust = FALSE,
  mlcp = TRUE,
  centering = TRUE,
  debug = FALSE
)

dgev_p12_cp(
  x,
  t1,
  t2,
  t01 = NA,
  t02 = NA,
  n01 = NA,
  n02 = NA,
  y = x,
  ics = c(0, 0, 0, 0, 0),
  minxi = -1,
  maxxi = 1,
  extramodels = FALSE,
  rust = FALSE,
  nrust = 10,
  centering = TRUE,
  debug = FALSE
)

pgev_p12_cp(
  x,
  t1,
  t2,
  t01 = NA,
  t02 = NA,
  n01 = NA,
  n02 = NA,
  y = x,
  ics = c(0, 0, 0, 0, 0),
  minxi = -1,
  maxxi = 1,
  extramodels = FALSE,
  rust = FALSE,
  nrust = 1000,
  centering = TRUE,
  debug = FALSE
)

tgev_p12_cp(
  n,
  x,
  t1,
  t2,
  ics = c(0, 0, 0, 0, 0),
  extramodels = FALSE,
  debug = FALSE
)
}
\arguments{
\item{x}{a vector of training data values}

\item{t1}{a vector of predictors for the mean,
such that \code{length(t1)=length(x)}}

\item{t2}{a vector of predictors for the sd,
such that \code{length(t2)=length(x)}}

\item{t01}{a single value of the predictor
(specify either \code{t01} or \code{n01} but not both)}

\item{t02}{a single value of the predictor
(specify either \code{t02} or \code{n02} but not both)}

\item{n01}{an index for the predictor
(specify either \code{t01} or \code{n01} but not both)}

\item{n02}{an index for the predictor
(specify either \code{t02} or \code{n02} but not both)}

\item{p}{a vector of probabilities at which to generate predictive quantiles}

\item{ics}{initial conditions for the maximum likelihood search}

\item{fdalpha}{if \code{pdf=TRUE}, the fractional delta used for numerical derivatives
with respect to probability, for calculating the pdf as a function of quantiles}

\item{minxi}{the minimum allowed value of the shape parameter (decrease with caution)}

\item{maxxi}{the maximum allowed value of the shape parameter (increase with caution)}

\item{means}{logical that indicates whether to run additional
calculations and return analytical estimates for the distribution means (longer runtime)}

\item{waicscores}{logical that indicates whether to run additional
calculations and return estimates for the WAIC1 and WAIC2 scores (longer runtime)}

\item{extramodels}{logical that indicates whether to run additional
calculations and add three additional prediction models (longer runtime)}

\item{pdf}{logical that indicates whether to run additional
calculations and return density functions evaluated at quantiles specified by
the input probabilities (longer runtime)}

\item{dmgs}{logical that indicates whether DMGS calculations
should be run or not (longer run time)}

\item{rust}{logical that indicates whether RUST-based posterior
sampling calculations should be run or not (longer run time)}

\item{nrust}{the number of posterior samples used in the RUST calculations}

\item{predictordata}{logical that indicates whether predictordata should be calculated}

\item{centering}{logical that indicates whether the predictor should be centered}

\item{debug}{logical for turning on debug messages}

\item{n}{the number of random samples required}

\item{mlcp}{logical that indicates whether maxlik and parameter
uncertainty calculations should be performed (turn off to speed up RUST)}

\item{y}{a vector of values at which to calculate the density and
distribution functions}
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
The \code{fitdistcp} package contains functions that generate predictive distributions
for various statistical models,
with and without parameter uncertainty.
Parameter uncertainty is included by using Bayesian prediction with a type of objective
prior known as a calibrating prior.
Calibrating priors are chosen to give predictions that give good reliability
(i.e., are well calibrated), for any
underlying true parameter values.

There are five functions for each model,
each of which uses training data \code{x}.
For model \code{****} the five functions are as follows:

\itemize{
\item \code{q****_cp} returns predictive quantiles at the specified probabilities \code{p},
and various other diagnostics.
\item \code{r****_cp} returns \code{n} random deviates from the predictive distribution.
\item \code{d****_cp} returns the predictive density function at the specified values \code{y}
\item \code{p****_cp} returns the predictive distribution function at the specified values \code{y}
\item \code{t****_cp} returns \code{n} random deviates from the posterior distribution
of the model parameters.
}
The \code{q, r, d, p} routines return two sets of results,
one based on maximum likelihood, and the other based on a calibrating prior.
The prior used depends on the model, and
is given under Details below.

Where possible, the Bayesian prediction integral is solved analytically.
Otherwise, DMGS asymptotic expansions are used.
Optionally, a third set of results is returned that integrates the prediction
integral by sampling the parameter posterior
distribution using the RUST rejection sampling algorithm.
}
\section{Details of the Model}{

The GEV distribution with two predictors has distribution function
\deqn{F(x;a_1,b_1,a_2,b_2,\xi)=\exp{(-t(x;\mu(a_1,b_1),\sigma(a_2,b_2),\xi))}}
where
\deqn{t(x;\mu(a_1,b_1),\sigma(a_2,b_2),\xi) =
   \begin{cases}
     {\left[1+\xi\left(\frac{x-\mu(a_1,b_1)}{\sigma(a_2,b_2)}\right)\right]}^{-1/\xi} & \text{if $\xi \ne 0$}\\
     \exp{\left(-\frac{x-\mu(a_1,b_1)}{\sigma(a_2,b_2)}\right)} & \text{if $\xi=0$}
   \end{cases}}
where
\eqn{x} is the random variable,
\eqn{\mu=a_1+b_1t_1} is the location parameter,
modelled as a function of parameters \eqn{a_1,b_1} and predictor \eqn{t_1},
\eqn{\sigma=e^{a_2+b_2t_2}} is the scale parameter,
modelled as a function of parameters \eqn{a_2,b_2} and predictor \eqn{t_2},
and \eqn{\xi} is the shape parameter.

The calibrating prior we use is given by
\deqn{\pi(a_1,b_1,a_2,b_2,\xi) \propto 1}
as given in Jewson et al. (2025).

The code will stop with an error if the
input data gives a maximum likelihood
value for the shape parameter that lies outside the range \code{(minxi,maxxi)},
since outside this range there may be numerical problems.
Such values seldom occur
in real observed data for maxima.
}

\section{Optional Return Values}{


\code{q****} optionally returns the following:

If \code{rust=TRUE}:
\itemize{
\item \code{ru_quantiles:} predictive quantiles calculated using a calibrating prior,
using posterior sampling with the RUST algorithm, based on inverting
an empirical CDF based on \code{nrust} samples.
}
If \code{waicscores=TRUE}:
\itemize{
\item \code{waic1:} the WAIC1 score for the calibrating prior model.
\item \code{waic2:} the WAIC2 score for the calibrating prior model.
}
If \code{logscores=TRUE}:
\itemize{
\item \code{ml_oos_logscore:} the leave-one-out logscore for the
maximum likelihood prediction
(not available in EVT routines, for mathematical reasons)
\item \code{cp_oos_logscore:} the leave-one-out logscore for the
parameter uncertainty model
available in EVT routines, for mathematical reasons)
}
If \code{means=TRUE}:
\itemize{
\item \code{ml_mean:} analytic estimate of the mean of the MLE predictive distribution,
where possible
\item \code{cp_mean:} analytic estimate of the mean of the calibrating prior
predictive distribution, where mathematically possible.
Can be compared with the mean estimated from random deviates.
}

\code{r****} optionally returns the following:

If \code{rust=TRUE}:
\itemize{
\item \code{ru_deviates:} \code{nrust} predictive random deviatives calculated using a
calibrating prior,
using posterior sampling with RUST.
}

\code{d****} optionally returns the following:

If \code{rust=TRUE}:
\itemize{
\item \code{ru_pdf:} predictive density calculated using a calibrating prior,
using posterior sampling with RUST,
averaging over \code{nrust} density functions.
}

\code{p****} optionally returns the following:

If \code{rust=TRUE}:
\itemize{
\item \code{ru_cdf:} predictive probability calculated using a calibrating prior,
using posterior sampling with RUST,
averaging over \code{nrust} distribution functions.
}
Selecting these additional outputs increases runtime.
They are optional so that runtime for the basic outputs is minimised.
This facilitates repeated experiments that evaluate reliability
over many thousands of repeats.
}

\section{Optional Return Values (EVT models only)}{


\code{q****} optionally returns the following, for EVT models only:

\itemize{
\item \code{cp_pdf:} the density function at quantiles corresponding to input
probabilities \code{p}.
We provide this for EVD models,
because direct estimation of the density function using the DMGS density
equation is not possible.
}
}

\section{Details (non-homogeneous models)}{

This model is not homogeneous, i.e. it
does not have a transitive transformation group,
and so there is no right Haar prior and no method
for generating exactly reliable predictions.
The \code{cp} outputs are generated using a prior that has been shown
in tests to give reasonable reliability.
See Jewson et al. (2024) for discussion of the prior and test results.
For non-homogeneous models, reliability is generally poor for small sample sizes (<20),
but is still much better than maximum likelihood.
For small sample sizes, it is advisable to check the level of reliability
using the routine \code{reltest}.
}

\section{Details (DMGS integration)}{

For this model, the Bayesian prediction equation cannot be solved analytically,
and is approximated using the DMGS asymptotic
expansions given by Datta et al. (2000).
This approximation seems to work well for medium and large sample sizes,
but may not work well for small sample sizes (e.g., <20 data points).
For small sample sizes, it may be preferable to use posterior
sampling by setting \code{rust=TRUE} and looking at the \code{ru} outputs.
The performance for any sample size, in terms of reliability,
can be tested using \code{reltest}.
}

\section{Details (RUST)}{

The Bayesian prediction equation can also be integrated using
ratio-of-uniforms-sampling-with-transformation (RUST), using the option
\code{rust=TRUE}.
\code{fitdistcp} then calls Paul Northrop's \code{rust} package (Northrop, 2023).
The RUST calculations are slower than the DMGS calculations.

For small sample sizes (e.g., n<20), and the very extreme tail, the DMGS approximation is somewhat poor
(although always better than maximum likelihood) and
it may be better to use RUST.
For medium sample sizes (30+), DMGS is reasonably accurate, except for the very far tail.

It is advisable to check the RUST results for convergence versus the number
of RUST samples.

It may be interesting to compare the DMGS and RUST results.
}

\examples{
# example 1
x=fitdistcp::d151gev_p12_example_data_v1_x
tt=fitdistcp::d151gev_p12_example_data_v1_t
t1=tt[,1]
t2=tt[,2]
p=c(1:9)/10
n01=10
n02=10
q=qgev_p12_cp(x=x,t1=t1,t2=t2,n01=n01,n02=n02,t01=NA,t02=NA,p=p,rust=TRUE,nrust=1000)
xmin=min(q$ml_quantiles,q$cp_quantiles);
xmax=max(q$ml_quantiles,q$cp_quantiles);
plot(q$ml_quantiles,p,xlab="quantile estimates",xlim=c(xmin,xmax),
	sub="(from qgev_p12_cp)",
	main="GEVD w/ p12: quantile estimates");
points(q$cp_quantiles,p,col="red",lwd=2)
points(q$ru_quantiles,p,col="blue",lwd=2)
cat(" ml_params=",q$ml_params,"\n")
}
\references{
If you use this package, we would be grateful if you would cite the following reference,
which gives the various calibrating priors, and tests them for reliability:

\itemize{
\item Jewson S., Sweeting T. and Jewson L. (2024): Reducing Reliability Bias in
Assessments of Extreme Weather Risk using Calibrating Priors;
ASCMO Advances in Statistical Climatology, Meteorology and Oceanography),
\url{https://ascmo.copernicus.org/articles/11/1/2025/}.
}
}
\seealso{
An introduction to \code{fitdistcp}, with more examples,
is given \href{http://www.fitdistcp.info/index.html}{on this webpage}.

The \code{fitdistcp} package currently includes the following models (in alphabetical order):
\itemize{
\item Cauchy (\code{cauchy}),
\item Cauchy with linear predictor on the mean (\code{cauchy_p1}),
\item Exponential (\code{exp}),
\item Exponential with log-linear predictor on the scale (\code{exp_p1}),
\item Frechet with known location parameter (\code{frechet_k1}),
\item Frechet with log-linear predictor on the scale and known location parameter
(\code{frechet_p2k1}),
\item Gamma (\code{gamma}),
\item Generalized normal (\code{gnorm}),
\item GEV (\code{gev}),
\item GEV with linear predictor on the location (\code{gev_p1}),
\item GEV with linear predictor on the location and log-linear prediction on the scale
(\code{gev_p12}),
\item GEV with linear predictor on the location, log-linear prediction on the scale,
and linear predictor on the shape (\code{gev_p123}),
\item GEV with linear predictor on the location and known shape (\code{gev_p1k3}),
\item GEV with known shape (\code{gev_k3}),
\item GPD with known location (\code{gpd_k1}),
\item Gumbel (\code{gumbel}),
\item Gumbel with linear predictor on the mean(\code{gumbel_p1}),
\item Half-normal (\code{halfnorm}),
\item Inverse gamma (\code{invgamma}),
\item Inverse Gaussian (\code{invgauss}),
\item t distribution with unknown location and scale and known DoF (\code{lst_k3}),
\item t distribution with unknown location and scale, linear predictor on the location,
and known DoF (\code{lst_p1k3}),
\item Logistic (\code{logis}),
\item Logistic with linear predictor on the location (\code{logis_p1}),
\item Log-normal (\code{lnorm}),
\item Log-normal with linear predictor on the location (\code{lnorm_p1}),
\item Normal (\code{norm}),
\item Normal with linear predictor on the mean (\code{norm_p1}),
\item Pareto with known scale (\code{pareto_k2}),
\item Pareto with log-linear predictor on the shape and known scale (\code{pareto_p1k2}),
\item Uniform (\code{unif}),
\item Weibull (\code{weibull}),
\item Weibull with linear predictor on the scale (\code{weibull_p2}),
}

The level of predictive probability matching achieved
by the maximum likelihood and calibrating prior quantiles, for any model,
sample size and true parameter values, can be demonstrated using the
routine \code{reltest}.

Model selection among models can be demonstrated using the routines
\code{ms_flat_1tail},
\code{ms_flat_2tail},
\code{ms_predictors_1tail},
and \code{ms_predictors_2tail},
}
\author{
Stephen Jewson \email{stephen.jewson@gmail.com}
}
