#' A blank function I use for setting up the man page information
#'
# #' 1) Generic description ##########################################################
# #'
#' @description
#' The \code{fitdistcp} package contains functions that generate predictive distributions
#' for various statistical models,
#' with and without parameter uncertainty.
#' Parameter uncertainty is included by using Bayesian prediction with a type of objective
#' prior known as a calibrating prior.
#'
#' There are five functions for each model,
#' each of which uses training data \code{x}.
#' For model \code{****} the five functions are as follows:
#'
#'\itemize{
#'\item \code{q****_cp} returns predictive quantiles at the specified probabilities \code{p},
#'and various other diagnostics.
#'\item \code{r****_cp} returns \code{n} random deviates from the predictive distribution.
#'\item \code{d****_cp} returns the predictive density function at the specified values \code{y}
#'\item \code{p****_cp} returns the predictive distribution function at the specified values \code{y}
#'\item \code{t****_cp} returns \code{n} random deviates from the posterior distribution
#'of the model parameters.
#'}
#' The \code{q, r, d, p} routines return two sets of results,
#' one based on maximum likelihood, and the other based on a calibrating prior.
#' The prior used depends on the model, and
#' is given under Details below.
#'
#' Where possible, the Bayesian prediction integral is solved analytically.
#' Otherwise, the DMGS asymptotic expansions are used.
#' Optionally, a third set of results is returned that integrates the prediction
#' integral by sampling the parameter posterior
#' distribution using the RUST rejection sampling algorithm.
#'
# #' 2) All possible parameters ########################################################
#' @param x 								a vector of training data values
#' @param t 								a vector of predictors, such that \code{length(t)=length(x)}
#' @param t1 								a vector of predictors for the mean,
#' such that \code{length(t1)=length(x)}
#' @param t2 								a vector of predictors for the sd,
#' such that \code{length(t2)=length(x)}
#' @param t3 								a vector of predictors for the shape,
#' such that \code{length(t3)=length(x)}
#' @param t0 								a single value of the predictor
#' (specify either \code{t0} or \code{n0} but not both)
#' @param t01 							a single value of the predictor
#' (specify either \code{t01} or \code{n01} but not both)
#' @param t02 							a single value of the predictor
#' (specify either \code{t02} or \code{n02} but not both)
#' @param t03 							a single value of the predictor
#' (specify either \code{t03} or \code{n03} but not both)
#' @param t10 							a single value of the predictor for the mean
#' (specify either \code{t10} or \code{n10} but not both)
#' @param t20 							a single value of the predictor for the sd
#' (specify either \code{t20} or \code{n20} but not both)
#' @param n0 								an index for the predictor
#' (specify either \code{t0} or \code{n0} but not both)
#' @param n01 							an index for the predictor
#' (specify either \code{t01} or \code{n01} but not both)
#' @param n02 							an index for the predictor
#' (specify either \code{t02} or \code{n02} but not both)
#' @param n03 							an index for the predictor
#' (specify either \code{t03} or \code{n03} but not both)
#' @param n10 							an index for the predictor for the mean
#' (specify either \code{t10} or \code{n10} but not both)
#' @param n20 							an index for the predictor for the sd
#' (specify either \code{t20} or \code{n20} but not both)
#' @param p 								a vector of probabilities at which to generate predictive quantiles
#' @param n									the number of random samples required
#' @param y									a vector of values at which to calculate the density and
#' distribution functions
#' @param ics								initial conditions for the maximum likelihood search
#' @param kloc							the known location parameter
#' @param kscale						the known scale parameter
#' @param kshape						the known shape parameter
#' @param kdf								the known degrees of freedom parameter
#' @param kbeta							the known beta parameter
#' @param d1								if \code{aderivs=FALSE}, the delta used for numerical derivatives
#' with respect to the first parameter
#' @param fd1								if \code{aderivs=FALSE}, the fractional delta used for numerical derivatives
#' with respect to the first parameter
#' @param d2								if \code{aderivs=FALSE}, the delta used for numerical derivatives
#' with respect to the second parameter
#' @param fd2								if \code{aderivs=FALSE}, the fractional delta used for numerical derivatives
#' with respect to the second parameter
#' @param d3								if \code{aderivs=FALSE}, the delta used for numerical derivatives
#' with respect to the third parameter
#'
#' @param fd3								if \code{aderivs=FALSE}, the fractional delta used for numerical derivatives
#' with respect to the third parameter
#' @param d4								if \code{aderivs=FALSE}, the delta used for numerical derivatives
#' with respect to the fourth parameter
#' @param fd4								if \code{aderivs=FALSE}, the fractional delta used for numerical derivatives
#' with respect to the fourth parameter
#' @param d5								if \code{aderivs=FALSE}, the delta used for numerical derivatives
#' with respect to the fifth parameter
#' @param fd5								if \code{aderivs=FALSE}, the fractional delta used for numerical derivatives
#' with respect to the fourth parameter
#' @param d6								if \code{aderivs=FALSE}, the delta used for numerical derivatives
#' with respect to the sixth parameter
#' @param fd6								if \code{aderivs=FALSE}, the fractional delta used for numerical derivatives
#' with respect to the fourth parameter
#' @param fdalpha						if \code{pdf=TRUE}, the fractional delta used for numerical derivatives
#' with respect to probability, for calculating the pdf as a function of quantiles
#' @param minxi							the minimum allowed value of the shape parameter (decrease with caution)
#' @param maxxi							the maximum allowed value of the shape parameter (increase with caution)
#' @param dlogpi						gradient of the log prior
#' @param means							a logical that indicates whether to run additional
#' calculations and return analytical estimates for the distribution means (longer runtime)
#' @param waicscores				a logical that indicates whether to run additional
#' calculations and return estimates for the WAIC1 and WAIC2 scores (longer runtime)
#' @param logscores					a logical that indicates whether to run additional
#' calculations and return leave-one-out estimates of the log-score
#' (much longer runtime, non-EVT models only)
#' @param extramodels				a logical that indicates whether to run additional
#' calculations and add three additional prediction models (longer runtime)
#' @param pdf								a logical that indicates whether to run additional
#' calculations and return density functions evaluated at quantiles specified by
#' the input probabilities (longer runtime)
#' @param customprior				a custom value for the slope of the log prior at the
#' maxlik estimate
#' @param dmgs							a logical that indicates whether DMGS calculations
#' should be run or not (longer run time)
#' @param mlcp							a logical that indicates whether maxlik and parameter
#' uncertainty calculations should be performed (turn off to speed up RUST)
#' @param predictordata		 	a logical that indicates whether predictordata should be calculated
#' @param centering					a logical that indicates whether the predictor should be centered
#' @param nonnegslopesonly	a logical that indicates whether to disallow non-negative slopes
#' @param rnonnegslopesonly	a logical that indicates whether to disallow non-negative slopes
#' @param	prior							a logical indicating which prior to use
#' @param debug							a logical for turning on debug messages
#' @param rust							a logical that indicates whether RUST-based posterior
#' sampling calculations should be run or not (longer run time)
#' @param nrust							the number of posterior samples used in the RUST calculations
#' @param pwm								a logical for whether to include PWM results (longer runtime)
#' @param unbiasedv					a logical for whether to include unbiased variance results in norm
#' @param aderivs						(for code testing only) a logical for whether to use
#' analytic derivatives (instead of numerical). By default almost all models now
#' use analytical derivatives.
#'
# #' 3a) Default Returns for all cases #################################################
#'
#' @section Default Return Values:
#'
#' \code{q****} returns a list containing the following:
#'
#' \itemize{
#' \item \code{ml_params:} maximum likelihood estimates for the parameters.
#' \item \code{ml_value:} the value of the log-likelihood at the maximum.
#' \item \code{standard_errors:} estimates of the standard errors on the parameters,
#' from the inverse observed information matrix.
#' \item \code{ml_quantiles:} quantiles calculated using maximum likelihood.
#' \item \code{cp_quantiles:} predictive quantiles calculated using a calibrating prior.
#' \item \code{maic:} the AIC score for the maximum likelihood model, times -1/2.
#' \item \code{cp_method:} a comment about the method used to generate the \code{cp} prediction.
#' }
#'
#' For models with predictors, \code{q****} additionally returns:
#' \itemize{
#' \item \code{predictedparameter:} the estimated value for parameter,
#' as a function of the predictor.
#' \item \code{adjustedx:} the detrended values of \code{x}
#' }
#'
#' \code{r****} returns a list containing the following:
#'
#' \itemize{
#' \item \code{ml_params:} maximum likelihood estimates for the parameters.
#' \item \code{ml_deviates:} random deviates calculated using maximum likelihood.
#' \item \code{cp_deviates:} predictive random deviates calculated using a calibrating prior.
#' \item \code{cp_method:} a comment about the method used to generate the \code{cp} prediction.
#' }
#'
#' \code{d****} returns a list containing the following:
#'
#' \itemize{
#' \item \code{ml_params:} maximum likelihood estimates for the parameters.
#' \item \code{ml_pdf:} density function from maximum likelihood.
#' \item \code{cp_pdf:} predictive density function calculated using a calibrating prior
#' (not available in EVT routines, for mathematical reasons, but available using RUST).
#' \item \code{cp_method:} a comment about the method used to generate the \code{cp} prediction.
#' }
#'
#' \code{p***} returns a list containing the following:
#'
#' \itemize{
#' \item \code{ml_params:} maximum likelihood estimates for the parameters.
#' \item \code{ml_cdf:} distribution function from maximum likelihood.
#' \item \code{cp_cdf:} predictive distribution function  calculated using a calibrating prior
#' (not available in EVT routines, for mathematical reasons, but available using RUST).
#' \item \code{cp_method:} a comment about the method used to generate the \code{cp} prediction.
#' }
#'
#' \code{t***} returns a list containing the following:
#'
#' \itemize{
#' \item \code{theta_samples:} random samples from the parameter posterior.
#' }
#'
# #' 3b) Optional returns for all cases #################################################
#'
#' @section Optional Return Values:
#'
#' \code{q****} optionally returns the following:
#'
#' If \code{rust=TRUE}:
#' \itemize{
#' \item \code{ru_quantiles:} predictive quantiles calculated using a calibrating prior,
#' using posterior sampling with the RUST algorithm, based on inverting
#' an empirical CDF based on \code{nrust} samples.
#' }
#' If \code{waicscores=TRUE}:
#' \itemize{
#' \item \code{waic1:} the WAIC1 score for the calibrating prior model.
#' \item \code{waic2:} the WAIC2 score for the calibrating prior model.
#' }
#' If \code{logscores=TRUE}:
#' \itemize{
#' \item \code{ml_oos_logscore:} the leave-one-out logscore for the
#' maximum likelihood prediction
#' (not available in EVT routines, for mathematical reasons)
#' \item \code{cp_oos_logscore:} the leave-one-out logscore for the
#' parameter uncertainty model
#'  available in EVT routines, for mathematical reasons)
#' }
#' If \code{means=TRUE}:
#' \itemize{
#' \item \code{ml_mean:} analytic estimate of the mean of the MLE predictive distribution,
#' where possible
#' \item \code{cp_mean:} analytic estimate of the mean of the calibrating prior
#' predictive distribution, where mathematically possible.
#' Can be compared with the mean estimated from random deviates.
#' }
#'
#' \code{r****} optionally returns the following:
#'
#' If \code{rust=TRUE}:
#' \itemize{
#' \item \code{ru_deviates:} \code{nrust} predictive random deviatives calculated using a
#' calibrating prior,
#' using posterior sampling with RUST.
#' }
#'
#' \code{d****} optionally returns the following:
#'
#' If \code{rust=TRUE}:
#' \itemize{
#' \item \code{ru_pdf:} predictive density calculated using a calibrating prior,
#' using posterior sampling with RUST,
#' averaging over \code{nrust} density functions.
#' }
#'
#' \code{p****} optionally returns the following:
#'
#' If \code{rust=TRUE}:
#' \itemize{
#' \item \code{ru_cdf:} predictive probability calculated using a calibrating prior,
#' using posterior sampling with RUST,
#' averaging over \code{nrust} distribution functions.
#' }
#'Selecting these additional outputs increases runtime.
#'They are optional so that runtime for the basic outputs is minimised.
#'This facilitates repeated experiments that evaluate reliability
#'over many thousands of repeats.
#'
# #' 3c) Returns EVD models #################################################
#'
#' @section Optional Return Values (EVT models only):
#'
#' \code{q****} optionally returns the following, for EVT models only:
#'
#' \itemize{
#' \item \code{cp_pdf:} the density function at quantiles corresponding to input
#' probabilities \code{p}.
#' We provide this for EVD models,
#' because direct estimation of the density function using the DMGS density
#' equation is not possible.
#' }
# #' 3d) Returns for some EVT models #################################################
#'
#' @section Optional Return Values (some EVT models only):
#'
#' {q****} optionally returns the following, for some EVT models only:
#'
#' If \code{extramodels=TRUE}:
#' \itemize{
#' \item \code{flat_quantiles:} predictive quantiles
#' calculated from Bayesian integration with a flat prior.
#' \item \code{rh_ml_quantiles:} predictive quantiles
#' calculated from Bayesian integration with the calibrating prior, and the maximmum
#' likelihood estimate for the shape parameter.
#' \item \code{jp_quantiles:} predictive quantiles
#' calculated from Bayesian integration with Jeffreys' prior.
#' }
#'
#' \code{r****} additionally returns the following, for some EVT models only:
#'
#' If \code{extramodels=TRUE}:
#' \itemize{
#' \item \code{flat_deviates:} predictive random deviates calculated using a
#' Bayesian analysis with a flat prior.
#' \item \code{rh_ml_deviates:} predictive random deviates calculated using a
#' Bayesian analysis with the RHP-MLE prior.
#' \item \code{jp_deviates:} predictive random deviates calculated using a
#' Bayesian analysis with the JP.
#' }
#'
#' \code{d****} additionally returns the following, for some EVT models only:
#'
#' If \code{extramodels=TRUE}:
#' \itemize{
#' \item \code{flat_pdf:} predictive density function from a Bayesian analysis
#' with the flat prior.
#' \item \code{rh_ml_pdf:} predictive density function from a Bayesian analysis
#' with the RHP-MLE prior.
#' \item \code{jp_pdf:} predictive density function from a Bayesian analysis
#' with the JP.
#' }
#'
#' \code{p****} additionally returns the following, for some EVT models only:
#'
#' If \code{extramodels=TRUE}:
#' \itemize{
#' \item \code{flat_cdf:} predictive distribution function from a Bayesian analysis
#' with the flat prior.
#' \item \code{rh_ml_cdf:} predictive distribution function from a Bayesian analysis
#' with the RHP-MLE prior.
#' \item \code{jp_cdf:} predictive distribution function from a Bayesian analysis
#' with the JP.
#' }

#'These additional predictive distributions are included for comparison
#'with the calibrating prior model.
#'They generally give less good reliability than the calibrating prior.
#'
# #' 4a) Details for homogeneous models ##################################################
#'
#' @section Details (homogeneous models):
#' This model is a homogeneous model, and the \code{cp} results are based on the
#' right Haar prior.
#' For homogeneous models
#' (models with sharply transitive transformations),
#' a Bayesian prediction based on the right Haar prior
#' gives exact predictive probability matching,
#' as shown by Severini et al. (2002).
#' Exact predictive probability matching then implies exact reliability,
#' even when the true parameters are unknown.
#' This means that probabilities in the prediction will correspond to frequencies
#' of future outcomes in repeated trials (if the model is correct).
#'
#' Maximum likelihood prediction does not give reliable predictions,
#' even when the model is correct, because it does not account for
#' parameter uncertainty.
#' In particular, maximum likelihood predictions typically underestimate the tail in repeated trials.
#'
#' The reliability of the maximum likelihood and the calibrating prior
#' predictive quantiles produced by
#' the \code{q****_cp} routines in \code{fitdistcp} can be quantified
#' using repeated simulations with the routine \code{reltest}.
#'
# #' 4b) Details for non-homogeneous models ##################################################

#' @section Details (non-homogeneous models):
#' This model is not homogeneous, i.e. it
#' does not have a transitive transformation group,
#' and so there is no right Haar prior and no method
#' for generating exactly reliable predictions.
#' The \code{cp} outputs are generated using a prior that has been shown
#' in tests to give reasonable reliability.
#' See Jewson et al. (2024) for discussion of the prior and test results.
#' For non-homogeneous models, reliability is generally poor for small sample sizes (<20),
#' but is still much better than maximum likelihood.
#' For small sample sizes, it is advisable to check the level of reliability
#' using the routine \code{reltest}.
#'
# #' 5a) Details for analytic models ###############################################
#'
#' @section Details (analytic integration):
#' For this model, the Bayesian prediction equation is integrated analytically.
#'
# #' 5b) Details for DMGS models ###############################################
#'
#' @section Details (DMGS integration):
#' For this model, the Bayesian prediction equation cannot be solved analytically,
#' and is approximated using the DMGS asymptotic
#' expansions given by Datta et al. (2000).
#' This approximation seems to work well for medium and large sample sizes,
#' but may not work well for small sample sizes (e.g., <20 data points).
#' For small sample sizes, it may be preferable to use posterior
#' sampling by setting \code{rust=TRUE} and looking at the \code{ru} outputs.
#' The performance for any sample size, in terms of reliability,
#' can be tested using \code{reltest}.
#'
# #' 5c) Details for RUST #################################################
#' @section Details (RUST):
#' The Bayesian prediction equation can also be integrated using
#' ratio-of-uniforms-sampling-with-transformation (RUST), using the option
#' \code{rust=TRUE}.
#' \code{fitdistcp} then calls Paul Northrop's \code{rust} package (Northrop, 2023).
#' The RUST calculations are slower than the DMGS calculations.
#'
#' For small sample sizes (e.g., n<20), and the very extreme tail, the DMGS approximation is somewhat poor
#' (although always better than maximum likelihood) and
#' it may be better to use RUST.
#' For medium sample sizes (30+), DMGS is reasonably accurate, except for the very far tail.
#'
#' It is advisable to check the RUST results for convergence versus the number
#' of RUST samples. It may also be interesting to compare the DMGS and RUST results.
#
# #' 6) author, seealso and references ##############################################
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @references
#'
#' If you use this package, we would be grateful if you would cite the following reference,
#' which gives the various calibrating priors, and tests them for reliability:
#'
#'\itemize{
#' \item Jewson S., Sweeting T. and Jewson L. (2024): Reducing Reliability Bias in
#' Assessments of Extreme Weather Risk using Calibrating Priors;
#' ASCMO (Advances in Statistical Climatology, Meteorology and Oceanography)
#' }
#'
#' The proof that using the right Haar prior gives exact probability
#' matching for homogeneous models is given in:
#'\itemize{
#' \item Severini, T., Mukerjee, R. and Ghosh, M. (2002); Biometrika
#' }
#'
#' The DMGS asymptotic expansions come from:
#'\itemize{
#' \item Datta, G., Mukerjee, R., Ghosh, M. and Sweeting, T. (2000); Annals of Statistics
#' }
#'
#' The RUST library is described here:
#'\itemize{
#' \item Northrop, P. (2023); Ratio-of-uniforms simulation with transformation, R package.
#' }

#'
#' @seealso
#'
#' An introduction to \code{fitdistcp}
#' is given [on this webpage](http://www.fitdistcp.info/index.html).
#'
#' The \code{fitdistcp} package currently includes the following models (in alphabetical order):
#' \itemize{
#' \item Cauchy (\code{cauchy}),
#' \item Cauchy with linear predictor on the mean (\code{cauchy_p1}),
#' \item Exponential (\code{exp}),
#' \item Exponential with log-linear predictor on the scale (\code{exp_p1}),
#' \item Frechet with known location parameter (\code{frechet_k1}),
#' \item Frechet with log-linear predictor on the scale and known location parameter
#' (\code{frechet_p2k1}),
#' \item Gamma (\code{gamma}),
#' \item Generalized normal (\code{gnorm}),
#' \item GEV (\code{gev}),
#' \item GEV with linear predictor on the location (\code{gev_p1}),
#' \item GEV with linear predictor on the location and log-linear prediction on the scale
#' (\code{gev_p12}),
#' \item GEV with linear predictor on the location, log-linear prediction on the scale,
#' and linear predictor on the shape (\code{gev_p123}),
#' \item GEV with linear predictor on the location and known shape (\code{gev_p1k3}),
#' \item GEV with known shape (\code{gev_k3}),
#' \item GPD with known location (\code{gpd_k1}),
#' \item Gumbel (\code{gumbel}),
#' \item Gumbel with linear predictor on the mean(\code{gumbel_p1}),
#' \item Half-normal (\code{halfnorm}),
#' \item Inverse gamma (\code{invgamma}),
#' \item Inverse Gaussian (\code{invgauss}),
#' \item t distribution with unknown location and scale and known DoF (\code{lst_k3}),
#' \item t distribution with unknown location and scale, linear predictor on the location,
#' and known DoF (\code{lst_p1k3}),
#' \item Logistic (\code{logis}),
#' \item Logistic with linear predictor on the location (\code{logis_p1}),
#' \item Log-normal (\code{lnorm}),
#' \item Log-normal with linear predictor on the location (\code{lnorm_p1}),
#' \item Normal (\code{norm}),
#' \item Normal with linear predictor on the mean (\code{norm_p1}),
#' \item Pareto with known scale (\code{pareto_k2}),
#' \item Pareto with log-linear predictor on the shape and known scale (\code{pareto_p1k2}),
#' \item Uniform (\code{unif}),
#' \item Weibull (\code{weibull}),
#' \item Weibull with linear predictor on the scale (\code{weibull_p2}),
#' }
#'
#' The level of predictive probability matching achieved
#' by the maximum likelihood and calibrating prior quantiles, for any model,
#' sample size and true parameter values, can be demonstrated using the
#' routine \code{reltest}.
#'
#' Model selection among models can be demonstrated using the routines \code{modelselection_flat}
#' and \code{modelselection_predictors}.
#'
#' @name man
#' @export
man=function(x,t,t1,t2,t3,t0,t01,t02,t03,t10,t20,n0,n01,n02,n03,n10,n20,p,n,y,ics,
	kloc,kscale,kshape,kdf,
	kbeta,d1,fd1,d2,fd2,d3,fd3,d4,fd4,d5,fd5,d6,fd6,fdalpha,minxi,maxxi,dlogpi,
	means,waicscores,logscores,extramodels,pdf,customprior,dmgs,mlcp,
	predictordata,centering,nonnegslopesonly,rnonnegslopesonly,
	prior,debug,rust,nrust,pwm,unbiasedv,aderivs){}
