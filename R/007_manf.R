#' Blank function I use for setting up the man page information for the functions
#' @param vv								parameters
#' @param ml_params					parameters
#' @param	nx								length of training data
#' @param	nxx								length of training data
#' @param x 								a vector of training data values
#' @param xx 								a vector of training data values
#' @param t 								a vector or matrix of predictors
#' @param nt 								the number of columns in \code{t}
#' @param ta								a vector of predictors for the mean (first column)
#' @param tb								a vector of predictors for the mean (second column)
#' @param tc								a vector of predictors for the mean (third column)
#' @param t1 								a vector of predictors for the mean
#' @param t2 								a vector of predictors for the sd
#' @param t3 								a vector of predictors for the shape
#' @param tt 								a vector of predictors
#' @param tt1 							a vector of predictors for the mean
#' @param tt2 							a vector of predictors for the sd
#' @param tt3 							a vector of predictors for the shape
#' @param tt2d 							a matrix of predictors (nx by 2)
#' @param tt3d 							a matrix of predictors (nx by 3)
#' @param t0 								a single value of the predictor (specify either \code{t0} or \code{n0} but not both)
#' @param t0a								a single value of the predictor, for the first column of the predictor	(specify either \code{t0a} or \code{n0a} but not both)
#' @param t0b								a single value of the predictor, for the second column of the predictor (specify either \code{t0b} or \code{n0b} but not both)
#' @param t0c								a single value of the predictor, for the third column of the predictor	(specify either \code{t0c} or \code{n0c} but not both)
#' @param t01 							a single value of the predictor (specify either \code{t01} or \code{n01} but not both)
#' @param t02 							a single value of the predictor (specify either \code{t02} or \code{n02} but not both)
#' @param t03 							a single value of the predictor (specify either \code{t03} or \code{n03} but not both)
#' @param t10 							a single value of the predictor for the mean (specify either \code{t10} or \code{n10} but not both)
#' @param t20 							a single value of the predictor for the sd (specify either \code{t20} or \code{n20} but not both)
#' @param t30 							a single value of the predictor for the shape (specify either \code{t30} or \code{n30} but not both)
#' @param n0 								an index for the predictor (specify either \code{t0} or \code{n0} but not both)
#' @param n10 							an index for the predictor for the mean (specify either \code{t10} or \code{n10} but not both)
#' @param n20 							an index for the predictor for the sd (specify either \code{t10} or \code{n10} but not both)
#' @param p 								a vector of probabilities at which to generate predictive quantiles
#' @param n									number of random samples required
#' @param y									a vector of values at which to calculate the density and distribution functions
#' @param ics								initial conditions for the maximum likelihood search
#' @param tresid						predictor residuals
#' @param tresid0						predictor residual at the point being predicted
#' @param muhat0						muhat at the point being predicted
#' @param kscale						the known scale parameter
#' @param kloc 							the known location parameter
#' @param kshape						the known shape parameter
#' @param kdf								the known degrees of freedom parameter
#' @param	kbeta							the known beta parameter
#' @param vhat							vector of all parameters
#' @param v1								first parameter
#' @param v1hat							first parameter
#' @param v1h								first parameter
#' @param d1								the delta used in the numerical derivatives with respect to the parameter
#' @param fd1								the fractional delta used in the numerical derivatives with respect to the parameter
#' @param v2								second parameter
#' @param v2hat							second parameter
#' @param v2h								second parameter
#' @param d2								the delta used in the numerical derivatives with respect to the parameter
#' @param fd2								the fractional delta used in the numerical derivatives with respect to the parameter
#' @param v3								third parameter
#' @param v3hat							third parameter
#' @param v3h								third parameter
#' @param d3								the delta used in the numerical derivatives with respect to the parameter
#' @param fd3								the fractional delta used in the numerical derivatives with respect to the parameter
#' @param v4								fourth parameter
#' @param v4hat							fourth parameter
#' @param v4h								fourth parameter
#' @param d4								the delta used in the numerical derivatives with respect to the parameter
#' @param fd4								the fractional delta used in the numerical derivatives with respect to the parameter
#' @param v5								fifth parameter
#' @param v5hat							fifth parameter
#' @param v5h								fifth parameter
#' @param d5								the delta used in the numerical derivatives with respect to the parameter
#' @param v6								sixth parameter
#' @param v6hat							sixth parameter
#' @param v6h								sixth parameter
#' @param d6								the delta used in the numerical derivatives with respect to the parameter
#' @param minxi							minimum value of shape parameter xi
#' @param maxxi							maximum value of shape parameter xi
#' @param ximin							minimum value of shape parameter xi
#' @param ximax							maximum value of shape parameter xi
#' @param alpha							a vector of values of alpha (one minus probability)
#' @param fdalpha						the fractional delta used in the numerical derivatives with respect to probability, for calculating the pdf as a function of quantiles
#' @param means							logical that indicates whether to return analytical estimates for the distribution means (longer runtime)
#' @param waicscores				logical that indicates whether to return estimates for the waic1 and waic2 scores (longer runtime)
#' @param logscores					logical that indicates whether to return leave-one-out estimates estimates of the log-score (much longer runtime)
#' @param extramodels				logical that indicates whether to add three additional prediction models
#' @param pdf								logical that indicates whether to return density functions evaluated at quantiles specified by input probabilities
#' @param predictordata			logical that indicates whether to calculate and return predictordata
#' @param nonnegslopesonly	logical that indicates whether to disallow non-negative slopes
#' @param rnonnegslopesonly	logical that indicates whether to disallow non-negative slopes
#' @param aderivs						logical for whether to use analytic derivatives (instead of numerical)
#' @param	ymn								the location parameter of the function of the predictor
#' @param slope							the slope of the function of the predictor
#' @param mu								the location parameter of the distribution
#' @param sigma							the sigma parameter of the distribution
#' @param sigma1						first coefficient for the sigma parameter of the distribution
#' @param sigma2						second coefficient for the sigma parameter of the distribution
#' @param scale							the scale parameter of the distribution
#' @param shape							the shape parameter of the distribution
#' @param xi								the shape parameter of the distribution
#' @param xi1								first coefficient for the shape parameter of the distribution
#' @param xi2								second coefficient for the shape parameter of the distribution
#' @param lambda						the lambda parameter of the distribution
#' @param log								logical for the density evaluation
#' @param mm								an index for which derivative to calculate
#' @param nn								an index for which derivative to calculate
#' @param rr								an index for which derivative to calculate
#' @param lddi							inverse observed information matrix
#' @param lddi_k2						inverse observed information matrix, fixed shape parameter
#' @param lddi_k3						inverse observed information matrix, fixed shape parameter
#' @param lddi_k4						inverse observed information matrix, fixed shape parameter
#' @param lddd							third derivative of log-likelihood
#' @param lddd_k2						third derivative of log-likelihood, fixed shape parameter
#' @param lddd_k3						third derivative of log-likelihood, fixed shape parameter
#' @param lddd_k4						third derivative of log-likelihood, fixed shape parameter
#' @param lambdad						derivative of the log prior
#' @param lambdad_cp				derivative of the log prior
#' @param lambdad_rhp 			derivative of the log RHP prior
#' @param lambdad_flat			derivative of the log flat prior
#' @param lambdad_rh_mle		derivative of the log CRHP-MLE prior
#' @param lambdad_rh_flat		derivative of the log CRHP-FLAT prior
#' @param lambdad_jp				derivative of the log JP prior
#' @param lambdad_custom		custom value of the derivative of the log prior
#' @param dim								number of parameters
#' @param customprior				a custom value for the slope of the log prior at the maxlik estimate
#' @param	prior							logical indicating which prior to use
#' @param	params						model parameters for calculating logf
#' @param yy								vector of samples
#' @param pp								vector of probabilities
#' @param	dlogpi						gradient of the log prior
#' @param	debug							debug flag
#' @param	centering					indicates whether the routine should center the data or not
#' @return No return value
#' @name manf
#' @export
manf=function(dim,vv,ml_params,nx,nxx,x,xx,
	t,nt,ta,tb,tc,t1,t2,t3,tt,tt1,tt2,tt3,tt2d,tt3d,
	t0,t0a,t0b,t0c,
	t01,t02,t03,
	t10,t20,t30,
	n0,n10,n20,p,n,y,ics,
	tresid,tresid0,muhat0,
	vhat,
	v1,v1hat,v1h,d1,fd1,v2,v2hat,v2h,d2,fd2,v3,v3hat,v3h,d3,fd3,
	v4,v4hat,v4h,d4,fd4,v5,v5hat,v5h,d5,v6,v6hat,v6h,d6,minxi,maxxi,ximin,ximax,fdalpha,
	kscale,kloc,kshape,kdf,kbeta,alpha,
	ymn,slope,mu,sigma,sigma1,sigma2,scale,shape,xi,xi1,xi2,lambda,log,
	mm,nn,rr,lddi,lddi_k2,lddi_k3,lddi_k4,
	lddd,lddd_k2,lddd_k3,lddd_k4,
	lambdad,lambdad_cp,lambdad_rhp,lambdad_flat,lambdad_rh_mle,
	lambdad_rh_flat,lambdad_jp,lambdad_custom,
	means,waicscores,logscores,extramodels,pdf,predictordata,
	nonnegslopesonly,rnonnegslopesonly,customprior,prior,params,
	yy,pp,dlogpi,debug,centering,aderivs){}
