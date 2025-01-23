#' Normal Distribution with a Predictor, Predictions Based on a Calibrating Prior
#'
#' @inherit man description author references seealso
#' @inheritParams man
#'
#' @inheritSection man Default Return Values
#' @inheritSection man Optional Return Values
# #' @inheritSection man Optional Return Values (EVD models only)
# #' @inheritSection man Optional Return Values (non-RHP models only)
#' @inheritSection man Details (homogeneous models)
# #' @inheritSection man Details (non-homogeneous models)
#' @inheritSection man Details (analytic integration)
# #' @inheritSection man Details (DMGS integration)
#' @inheritSection man Details (RUST)
#'
#' @section Details of the Model:
#' The normal distribution with a predictor has probability density function
#' \deqn{f(x;a,b,\sigma)=\frac{1}{\sqrt{2\pi}\sigma}e^{-(x-\mu(a,b))^2/(2\sigma^2)}}
#' where
#' \eqn{x} is the random variable,
#' \eqn{\mu=a+bt} is the location parameter, modelled as a function
#' of parameters \eqn{a,b} and predictor \eqn{t},
#' and \eqn{\sigma>0} is the scale parameter.
#'
#' The calibrating prior is given by the right Haar prior, which is
#' \deqn{\pi(a,b,\sigma) \propto \frac{1}{\sigma}}
#' as given in Jewson et al. (2024).
#'
#' @example man/examples/example_60_norm_p1.R
#'
#' @name norm_p1_cp
NULL
#' @rdname norm_p1_cp
#' @inheritParams man
#' @export
#'
qnorm_p1_cp=function(x,t,t0=NA,n0=NA,p=seq(0.1,0.9,0.1),d1=0.01,d2=0.01,fd3=0.01,
	means=FALSE,waicscores=FALSE,logscores=FALSE,rust=FALSE,nrust=100000,
	centering=TRUE,
	debug=FALSE,aderivs=TRUE){
#
# 1 intro
#
	debug=TRUE
	debug=FALSE
	stopifnot(	is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1)
	alpha=1-p
	nx=length(x)
#
# 2 centering
#
	if(centering){
		meant=mean(t)
		t=t-meant
		t0=t0-meant
	}
	ta=t-mean(t)
	t0=maket0(t0,n0,t)
	ta0=maketa0(t0,n0,t)
#
# 3 ml estimates
#
  ml_params=norm_p1_mlparams(x,t)
  v1hat=ml_params[1]
  v2hat=ml_params[2]
  v3hat=ml_params[3]
	muhat0=makemuhat0(t0,n0,t,ml_params)
  if(debug)cat("  v1hat,v2hat,v3hat=",v1hat,v2hat,v3hat,"\n")
#
# 4 predictordata
#
	prd=norm_p1_predictordata(x,t,t0,ml_params)
	predictedparameter=prd$predictedparameter
	adjustedx=prd$adjustedx
#
# 5 aic
#
 	ml_value=norm_p1_loglik(ml_params,x,t)
	maic=make_maic(ml_value,nparams=3)

# 6 ml quantiles at muhat0
#
	ml_quantiles=qnorm((1-alpha),mean=muhat0,sd=v3hat)

	ldd="only relevant for DMGS models, not analytic models"
	lddi="only relevant for DMGS models, not analytic models"
	expinfmat=matrix(0,3,3)
	expinfmat[1,1]=nx/(v3hat*v3hat)
	expinfmat[2,2]=sum(ta*ta)/(v3hat*v3hat)
	expinfmat[3,3]=2*nx/(v3hat*v3hat)
	expinfmat[1,2]=0
	expinfmat[2,1]=0
	expinfmat[1,3]=0
	expinfmat[3,1]=0
	expinfmati=solve(expinfmat)
	standard_errors=matrix(0,3)
	standard_errors[1]=sqrt(expinfmati[1,1])
	standard_errors[2]=sqrt(expinfmati[2,2])
	standard_errors[3]=sqrt(expinfmati[3,3])
#
# normal for comparison
#
#	mu=v1hat
#	fact=sqrt(1+(2/(nx-1)))
#	sg=v2hat*fact
#	temp=qt((1-alpha),df=nx-1)
#	rh_quantiles=mu+temp*sg
#
# 7 rhp quantiles at t0 and muhat0
#
	rh_quantiles=qnorm_p1_formula(alpha,ta,ta0,nx,muhat0,v3hat)
#
# 8 means (might as well always calculate)
#
	ml_mean=muhat0
	rh_mean=ml_mean
#
# 9 waicscores
#
	waic=norm_p1_waic(waicscores,x,t,v1hat,d1,v2hat,d2,v3hat,fd3,aderivs=aderivs)
	waic1=waic$waic1
	waic2=waic$waic2
#
# 10 logscores
#
	logscores=norm_p1_logscores(logscores,x,t,aderivs)
	ml_oos_logscore=logscores$ml_oos_logscore
	rh_oos_logscore=logscores$rh_oos_logscore
#
# 11 rust
#
	ru_quantiles="rust not selected"
	if(rust){
		rustsim=rnorm_p1_cp(nrust,x=x,t=t,t0=t0,rust=TRUE,mlcp=FALSE)
		ru_quantiles=makeq(rustsim$ru_deviates,p)
	}
#
# 12 decentering
#
	if(centering){
		ml_params[1]=ml_params[1]-ml_params[2]*meant
	}

# return
	list(	ml_params=ml_params,
				ml_value=ml_value,
				predictedparameter=predictedparameter,
				adjustedx=adjustedx,
#				ldd=ldd,
#				lddi=lddi,
#				expinfmat=expinfmat,
#				expinfmati=expinfmati,
				standard_errors=standard_errors,
				ml_quantiles=ml_quantiles,
				cp_quantiles=rh_quantiles,
				ru_quantiles=ru_quantiles,
				maic=maic,
				waic1=waic1,
				waic2=waic2,
				ml_oos_logscore=ml_oos_logscore,
				cp_oos_logscore=rh_oos_logscore,
				ml_mean=ml_mean,
				cp_mean=rh_mean,
				cp_method=analytic_cpmethod())

}
#' @rdname norm_p1_cp
#' @inheritParams man
#' @export
rnorm_p1_cp=function(n,x,t,t0=NA,n0=NA,rust=FALSE,mlcp=TRUE,debug=FALSE,aderivs=TRUE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),is.finite(t),!is.na(t))
	stopifnot(is.finite(x),!is.na(x),is.finite(t),!is.na(t))

	t0=maket0(t0,n0,t)
#
# 2 centering
#
	meant=mean(t)
	t=t-meant
	t0=t0-meant

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qnorm_p1_cp(x,t,t0=t0,n0=NA,p=runif(n),aderivs=aderivs)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tnorm_p1_cp(n,x,t)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			mu=th[i,1]+t0*th[i,2]
			ru_deviates[i]=rnorm(1,mean=mu,sd=th[i,3])
		}
	}
#
# decentering
#
	if(mlcp)ml_params[1]=ml_params[1]-ml_params[2]*meant

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
			 cp_method=analytic_cpmethod())

	return(op)

}
#' @rdname norm_p1_cp
#' @inheritParams man
#' @export
dnorm_p1_cp=function(x,t,t0=NA,n0=NA,y=x,rust=FALSE,nrust=1000,centering=TRUE,
	debug=FALSE,aderivs=TRUE){
# why is rust an option, given that we have exact solutions?

	t0=maket0(t0,n0,t)

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),is.finite(t),!is.na(t))

#
# centering
#
	if(centering){
		meant=mean(t)
		t=t-meant
		t0=t0-meant
	}

	dd=dnorm_p1sub(x=x,t=t,y=y,t0=t0,aderivs=aderivs)
	ml_params=dd$ml_params

	ru_pdf="rust not selected"
	if(rust){
		th=tnorm_p1_cp(nrust,x,t)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+t0*th[ir,2]
			ru_pdf=ru_pdf+dnorm(y,mean=mu,sd=th[ir,3])
		}
		ru_pdf=ru_pdf/nrust
	}

#
# decentering
#
	if(centering)ml_params[1]=ml_params[1]-ml_params[2]*meant

	op=list(	ml_params=ml_params,
					ml_pdf=dd$ml_pdf,
					cp_pdf=dd$rh_pdf,
					ml_cdf=dd$ml_cdf,
					cp_cdf=dd$rh_cdf,
					ru_pdf=ru_pdf,
					cp_method=analytic_cpmethod())
	return(op)
}
#' @rdname norm_p1_cp
#' @inheritParams man
#' @export
pnorm_p1_cp=function(x,t,t0=NA,n0=NA,y=x,rust=FALSE,nrust=1000,centering=TRUE,
	debug=FALSE,aderivs=TRUE){

	t0=maket0(t0,n0,t)

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),is.finite(t),!is.na(t))

#
# centering
#
	if(centering){
		meant=mean(t)
		t=t-meant
		t0=t0-meant
	}

	dd=dnorm_p1sub(x=x,t=t,y=y,t0=t0,aderivs=aderivs)
	ml_params=dd$ml_params

	ru_cdf="rust not selected"
	if(rust){
		th=tnorm_p1_cp(nrust,x,t)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+t0*th[ir,2]
			ru_cdf=ru_cdf+pnorm(y,mean=mu,sd=th[ir,3])
		}
		ru_cdf=ru_cdf/nrust
	}

#
# decentering
#
	if(centering)ml_params[1]=ml_params[1]-ml_params[2]*meant

	op=list(	ml_params=ml_params,
					ml_cdf=dd$ml_cdf,
					cp_cdf=dd$rh_cdf,
					ru_cdf=ru_cdf,
					cp_method=analytic_cpmethod())
	return(op)
}
#' @rdname norm_p1_cp
#' @inheritParams man
#' @export
tnorm_p1_cp=function(n,x,t,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),is.finite(t),!is.na(t))
	stopifnot(is.finite(x),!is.na(x),is.finite(t),!is.na(t))

#
# centering
#
	meant=mean(t)
	t=t-meant

	th=ru(norm_p1_logf,x=x,t=t,n=n,d=3,init=c(0,0,1))
	theta_samples=th$sim_vals

#
# decentering
#
	theta_samples[,1]=theta_samples[,1]-theta_samples[,2]*meant

	list(theta_samples=theta_samples)

}


