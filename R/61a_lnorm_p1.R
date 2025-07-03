#' Log-normal Distribution with a Predictor, Predictions Based on a Calibrating Prior
#'
#' @inherit man description author references seealso return
#' @inheritParams man
#'
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
#' The log normal distribution with a predictor has probability density function
#' \deqn{f(x;a,b,\sigma)=\frac{1}{\sqrt{2\pi}x\sigma}e^{-(\log(x)-\mu(a,b))^2/(2\sigma^2)}}
#' where
#' \eqn{x} is the random variable,
#' \eqn{\mu=a+bt} is the location parameter of the log of the random variable,
#' modelled as a function
#' of parameters \eqn{a,b} and predictor \eqn{t},
#' and \eqn{\sigma>0} is the scale parameter of the log of the random variable.
#'
#' The calibrating prior is given by the right Haar prior, which is
#' \deqn{\pi(a,b,\sigma) \propto \frac{1}{\sigma}}
#' as given in Jewson et al. (2025).
#'
#' @example man/examples/example_61_lnorm_p1.R
#'
#' @name lnorm_p1_cp
NULL
#' @rdname lnorm_p1_cp
#' @inheritParams man
#' @export
#'
qlnorm_p1_cp=function(x,t,t0=NA,n0=NA,p=seq(0.1,0.9,0.1),d1=0.01,d2=0.01,fd3=0.01,
	means=FALSE,waicscores=FALSE,logscores=FALSE,rust=FALSE,nrust=100000,
	centering=TRUE,
	debug=FALSE,aderivs=TRUE){
#
# 1 intro
#
	stopifnot(	is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1,
		!x<0)
	alpha=1-p
	z=log(x)
	nx=length(x)
#
# 2 centering
#
  if(centering){
	  meant=mean(t)
    t=t-meant
    t0=t0-meant
  }
	ta=t-mean(t) #whether or not we center, the equations use this ta
	t0=maket0(t0,n0,t)
	ta0=maketa0(t0,n0,t)
#
# 3 ml estimates
#
  ml_params=norm_p1_mlparams(z,t)
  v1hat=ml_params[1]
  v2hat=ml_params[2]
  v3hat=ml_params[3]
	muhatz0=makemuhat0(t0,n0,t,ml_params)
#	muhatz=ml_params[1]+ml_params[2]*t
#	muhatx=exp(muhatz)
#	residuals=x-muhatx
  if(debug)message("  v1hat,v2hat,v3hat=",v1hat,v2hat,v3hat)
#
# 4 predictordata
#
	prd=lnorm_p1_predictordata(x,t,t0,ml_params)
	predictedparameter=prd$predictedparameter
	adjustedx=prd$adjustedx
#
# 5 aic
#
 	ml_value=lnorm_p1_loglik(ml_params,x,t)
	maic=make_maic(ml_value,nparams=3)

# 6 ml quantiles at muhatz0
#
	ml_quantiles=qlnorm((1-alpha),meanlog=muhatz0,sdlog=v3hat) #basic, no parameter uncertainty

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
# 7 rhp quantiles at n0
#
	rh_quantiles=exp(qnorm_p1_formula(alpha,ta,ta0,nx,muhatz0,v3hat))
#
# 8 means (might as well always calculate)
#
	ml_mean=exp(muhatz0+0.5*(ml_params[3])**2)
	rh_mean="no analytic expression"
#
# 9 waicscores
#
# x is log-normal, y is normal
	waic=lnorm_p1_waic(waicscores,x,t,v1hat,d1,v2hat,d2,v3hat,fd3,aderivs=aderivs)
	waic1=waic$waic1
	waic2=waic$waic2
#
# 10 logscores
#
# x is log-normal, y is normal
	logscores=lnorm_p1_logscores(logscores,x,t)
	ml_oos_logscore=logscores$ml_oos_logscore
	rh_oos_logscore=logscores$rh_oos_logscore
#
# 11 rust
#
	ru_quantiles="rust not selected"
	if(rust){
		rustsim=rlnorm_p1_cp(nrust,x,t=t,t0=t0,rust=TRUE,mlcp=FALSE)
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
#' @rdname lnorm_p1_cp
#' @inheritParams man
#' @export
rlnorm_p1_cp=function(n,x,t,t0=NA,n0=NA,rust=FALSE,mlcp=TRUE,debug=FALSE,aderivs=TRUE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),
#		is.finite(t),!is.na(t),!x<0)
	stopifnot(is.finite(x),!is.na(x),
		is.finite(t),!is.na(t),!x<0)

	t0=maket0(t0,n0,t)

#
# centering
#
  meant=mean(t)
  t=t-meant
  t0=t0-meant

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qlnorm_p1_cp(x,t,t0=t0,n0=NA,p=runif(n),aderivs=aderivs)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tlnorm_p1_cp(n,x,t)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			mu=th[i,1]+t0*th[i,2]
			ru_deviates[i]=rlnorm(1,meanlog=mu,sdlog=th[i,3])
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
#' @rdname lnorm_p1_cp
#' @inheritParams man
#' @export
dlnorm_p1_cp=function(x,t,t0=NA,n0=NA,y=x,rust=FALSE,nrust=1000,centering=TRUE,
	debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
		is.finite(t),!is.na(t),!x<0,!y<0)

	t0=maket0(t0,n0,t)

#
# centering
#
  if(centering){
    meant=mean(t)
    t=t-meant
    t0=t0-meant
  }

	dd=dlnorm_p1sub(x=x,t=t,y=y,t0=t0,debug=debug,aderivs=aderivs)
	ru_pdf="rust not selected"
	ml_params=dd$ml_params
	if(rust){
		th=tlnorm_p1_cp(nrust,x,t)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+t0*th[ir,2]
			ru_pdf=ru_pdf+dlnorm(y,meanlog=mu,sdlog=th[ir,3])
		}
		ru_pdf=ru_pdf/nrust
	}
#
# decentering
#
	if(debug)message("before decentering:",ml_params)
	if(centering)ml_params[1]=ml_params[1]-ml_params[2]*meant
	if(debug)message("after decentering:",ml_params)

	op=list(	ml_params=ml_params,
				ml_pdf=dd$ml_pdf,
				cp_pdf=dd$rh_pdf,
				ru_pdf=ru_pdf,
				cp_method=analytic_cpmethod())
	return(op)
}
#' @rdname lnorm_p1_cp
#' @inheritParams man
#' @export
plnorm_p1_cp=function(x,t,t0=NA,n0=NA,y=x,rust=FALSE,nrust=1000,centering=TRUE,
	debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
		is.finite(t),!is.na(t),!x<0,!y<0)

	t0=maket0(t0,n0,t)

#
# centering
#
  if(centering){
    meant=mean(t)
    t=t-meant
    t0=t0-meant
  }
	dd=dlnorm_p1sub(x=x,t=t,y=y,t0=t0,debug=debug,aderivs=aderivs)
	ru_cdf="rust not selected"
	ml_params=dd$ml_params
	if(rust){
		th=tlnorm_p1_cp(nrust,x,t)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+t0*th[ir,2]
			ru_cdf=ru_cdf+plnorm(y,meanlog=mu,sdlog=th[ir,3])
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
#' @rdname lnorm_p1_cp
#' @inheritParams man
#' @export
tlnorm_p1_cp=function(n,x,t,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),
#		is.finite(t),!is.na(t),!x<0)
	stopifnot(is.finite(x),!is.na(x),
		is.finite(t),!is.na(t),!x<0)

#
# centering
#
  meant=mean(t)
  t=t-meant

	th=ru(lnorm_p1_logf,x=x,t=t,n=n,d=3,init=c(0,0,1))
  theta_samples=th$sim_vals
#
# decentering
#
  theta_samples[,1]=theta_samples[,1]-theta_samples[,2]*meant

	list(theta_samples=theta_samples)

}
