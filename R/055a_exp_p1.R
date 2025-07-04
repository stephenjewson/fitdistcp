#' Exponential Distribution with a Predictor, Predictions Based on a Calibrating Prior
#'
#' @inherit man description author references seealso return
#' @inheritParams man
#'
#' @inheritSection man Optional Return Values
# #' @inheritSection man Optional Return Values (EVD models only)
# #' @inheritSection man Optional Return Values (non-RHP models only)
#' @inheritSection man Details (homogeneous models)
# #' @inheritSection man Details (non-homogeneous models)
# #' @inheritSection man Details (analytic integration)
#' @inheritSection man Details (DMGS integration)
#' @inheritSection man Details (RUST)
#'
#' @section Details of the Model:
#' The exponential distribution with a predictor has exceedance distribution function
#' \deqn{S(x;a,b)=\exp(-x \lambda (a,b))}
#' where
#' \eqn{x \ge 0} is the random variable
#' and
#' \eqn{\lambda(a,b)=e^{-a-bt}} is the rate parameter, modelled as a function
#' of the parameters \eqn{a,b} and a predictor \eqn{t}.
#'
#' The calibrating prior is given by the right Haar prior, which is
#' \deqn{\pi(a,b) \propto 1}.
#' as given in Jewson et al. (2025).
#'
#' @example man/examples/example_055_exp_p1.R
#'
#' @name exp_p1_cp
NULL
#' @rdname exp_p1_cp
#' @inheritParams man
#' @export
#'
qexp_p1_cp=function(x,t,t0=NA,n0=NA,p=seq(0.1,0.9,0.1),d1=0.01,d2=0.01,
	means=FALSE,waicscores=FALSE,logscores=FALSE,
	dmgs=TRUE,rust=FALSE,nrust=100000,predictordata=TRUE,centering=TRUE,
	debug=FALSE,aderivs=TRUE){
#
# 1 intro
#
	debug=TRUE
	debug=FALSE
	if(debug)message("inside qexp_p1")
#	stopifnot(	is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1,length(t)==length(x),!x<0)
	stopifnot(	is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1,!x<0)
	alpha=1-p
	nx=length(x)
	nalpha=length(alpha)
	t0=maket0(t0,n0,t)
#
# 2 centering
#
  if(centering){
	  meant=mean(t)
    t=t-meant
    t0=t0-meant
  }
#
# 3 ml param estimate
#
	if(debug)message("calc ml param estimate")
	lm=lm(x~t)
	v1start=lm$coefficients[1]
	v2start=lm$coefficients[2]
#	xhat=v1start+v2start*t
	v1start=0
	v2start=0
	opt=optim(c(v1start,v2start),exp_p1_loglik,x=x,t=t,control=list(fnscale=-1))
	v1hat=opt$par[1]
	v2hat=opt$par[2]
	ml_params=c(v1hat,v2hat)
#
# 4 predictordata
#
	prd=exp_p1_predictordata(predictordata,x,t,t0,ml_params)
	predictedparameter=prd$predictedparameter
	adjustedx=prd$adjustedx
	if(debug)message("  ml_params=",ml_params)
#
# 5 aic
#
	if(debug)message("calc aic")
 	ml_value=opt$val
	maic=make_maic(ml_value,nparams=2)

#
# 6 mle quantiles
#
	if(debug)message("calc mle quantiles")
	ml_quantiles=qexp_p1((1-alpha),t0,ymn=v1hat,slope=v2hat)
#
# dmgs
#
	standard_errors="dmgs not selected"
	rh_quantiles="dmgs not selected"
	ru_quantiles="dmgs not selected"
	waic1="dmgs not selected"
	waic2="dmgs not selected"
	ml_oos_logscore="dmgs not selected"
	rh_oos_logscore="dmgs not selected"
	cp_oos_logscore="dmgs not selected"
	ml_mean="dmgs not selected"
	rh_mean="dmgs not selected"
	cp_mean="dmgs not selected"
	cp_method="dmgs not selected"
	if(dmgs){
#
# 7 lddi
#
		if(debug)message("calc ldd")
		if(aderivs)	ldd=exp_p1_ldda(x,t,v1hat,v2hat)
		if(!aderivs)ldd=exp_p1_ldd(x,t,v1hat,d1,v2hat,d2)
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
#
# 8 lddd
#
		if(debug)message("calculate lddd")
		if(aderivs) lddd=exp_p1_lddda(x,t,v1hat,v2hat)
		if(!aderivs)lddd=exp_p1_lddd(x,t,v1hat,d1,v2hat,d2)
#
# 9 mu1
#
		if(debug)message("calculate mu1")
		if(aderivs) mu1=exp_p1_mu1fa(alpha,t0,v1hat,v2hat)
		if(!aderivs)mu1=exp_p1_mu1f(alpha,t0,v1hat,d1,v2hat,d2)
#
# 10 mu2
#
		if(debug)message("calculate mu2")
		if(aderivs) mu2=exp_p1_mu2fa(alpha,t0,v1hat,v2hat)
		if(!aderivs)mu2=exp_p1_mu2f(alpha,t0,v1hat,d1,v2hat,d2)
#
# 11 rhp
#
		if(debug)message("  rhp")
		lambdad_rhp=c(0,0) #this is rhp
#
# 12 rhp quantiles
#
		if(debug)message("  rhp quantiles")
		fhat=dexp_p1(ml_quantiles,t0,ymn=v1hat,slope=v2hat,log=FALSE)
		dq=dmgs(lddi,lddd,mu1,lambdad_rhp,mu2,dim=2)
		rh_quantiles=ml_quantiles+dq/(nx*fhat)
#
# 13 means
#
		means=exp_p1_means(means,t0,ml_params,lddi,lddd,lambdad_rhp,nx,dim=2)
		ml_mean=means$ml_mean
		rh_mean=means$rh_mean
#
# 14 waicscores
#
		waic=exp_p1_waic(waicscores,x,t,v1hat,d1,v2hat,d2,lddi,lddd,lambdad_rhp,aderivs)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 15 logscores
#
		logscores=exp_p1_logscores(logscores,x,t,d1,d2,aderivs)
		ml_oos_logscore=logscores$ml_oos_logscore
		rh_oos_logscore=logscores$rh_oos_logscore
#
# 16 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rexp_p1_cp(nrust,x,t=t,t0=t0,rust=TRUE,mlcp=FALSE)
			ru_quantiles=makeq(rustsim$ru_deviates,p)
		}
	} #end of if(dmgs)

#
# 17 decentering
#
  if(centering){
    ml_params[1]=ml_params[1]-ml_params[2]*meant
    if(predictordata)predictedparameter=predictedparameter-ml_params[2]*meant
  }

# return
	list(	ml_params=ml_params,
				ml_value=ml_value,
				predictedparameter=predictedparameter,
				adjustedx=adjustedx,
#				ldd=ldd,
#				lddi=lddi,
#				lddd=lddd,
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
				cp_method=rhp_dmgs_cpmethod())

}
#' @rdname exp_p1_cp
#' @inheritParams man
#' @export
rexp_p1_cp=function(n,x,t,t0=NA,n0=NA,d1=0.01,d2=0.01,rust=FALSE,mlcp=TRUE,
	debug=FALSE,aderivs=TRUE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),
#		is.finite(t),!is.na(t),length(t)==length(x),!x<0)
#	stopifnot(is.finite(x),!is.na(x),
#		is.finite(t),!is.na(t),length(t)==length(x),!x<0)
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
		q=qexp_p1_cp(x,t,t0=t0,n0=NA,p=runif(n),d1=d1,d2=d2,aderivs=aderivs)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=texp_p1_cp(n,x,t)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			mu=exp(th[i,1]+t0*th[i,2])
			ru_deviates[i]=rexp(1,rate=1/mu)
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
				cp_method=rhp_dmgs_cpmethod())

	return(op)

}
#' @rdname  exp_p1_cp
#' @inheritParams 	man
#' @param d1				the fractional delta used in the numerical derivatives with respect to the location parameter
#' @param d2				the fractional delta used in the numerical derivatives with respect to the slope parameter
#' @export
dexp_p1_cp=function(x,t,t0=NA,n0=NA,y=x,d1=0.01,d2=0.01,rust=FALSE,nrust=1000,
	centering=TRUE,debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),is.finite(t),!is.na(t),!x<0,!y<0)

	t0=maket0(t0,n0,t)

#
# centering
#
  if(centering){
    meant=mean(t)
    t=t-meant
    t0=t0-meant
  }

	dd=dexp_p1sub(x=x,t=t,y=y,t0=t0,d1,d2,aderivs=aderivs)
	ru_pdf="rust not selected"
	ml_params=dd$ml_params
	if(rust){
		th=texp_p1_cp(nrust,x,t)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=exp(th[ir,1]+t0*th[ir,2])
			ru_pdf=ru_pdf+dexp(y,rate=1/mu)
		}
		ru_pdf=ru_pdf/nrust
	}
#
# decentering
#
 if(centering){
   ml_params[1]=ml_params[1]-ml_params[2]*meant
 }

	op=list(	ml_params=ml_params,
				ml_pdf=dd$ml_pdf,
				cp_pdf=dd$rh_pdf,
				ru_pdf=ru_pdf,
				cp_method=rhp_dmgs_cpmethod())
	return(op)
}
#' @rdname  exp_p1_cp
#' @inheritParams 	man
#' @export
pexp_p1_cp=function(x,t,t0=NA,n0=NA,y=x,d1=0.01,d2=0.01,rust=FALSE,nrust=1000,
	centering=TRUE,debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),is.finite(t),!is.na(t),!x<0,!y<0)

	t0=maket0(t0,n0,t)

#
# centering
#
  if(centering){
    meant=mean(t)
    t=t-meant
    t0=t0-meant
  }

	dd=dexp_p1sub(x=x,t=t,y=y,t0=t0,d1,d2,aderivs=aderivs)
	ru_cdf="rust not selected"
	ml_params=dd$ml_params
	if(rust){
		th=texp_p1_cp(nrust,x,t)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=exp(th[ir,1]+t0*th[ir,2])
			ru_cdf=ru_cdf+pexp(y,rate=1/mu)
		}
		ru_cdf=ru_cdf/nrust
	}
#
# decentering
#
 if(centering){
   ml_params[1]=ml_params[1]-ml_params[2]*meant
 }
	op=list(	ml_params=ml_params,
				ml_cdf=dd$ml_cdf,
				cp_cdf=dd$rh_cdf,
				ru_cdf=ru_cdf,
				cp_method=rhp_dmgs_cpmethod())
	return(op)
}
#' @rdname exp_p1_cp
#' @inheritParams man
#' @export
texp_p1_cp=function(n,x,t,d1=0.01,d2=0.01,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),
#		is.finite(t),!is.na(t),length(t)==length(x),!x<0)
#	stopifnot(is.finite(x),!is.na(x),
#		is.finite(t),!is.na(t),length(t)==length(x),!x<0)
	stopifnot(is.finite(x),!is.na(x),
		is.finite(t),!is.na(t),!x<0)
#
# centering
#
  meant=mean(t)
  t=t-meant

	th=ru(exp_p1_logf,x=x,t=t,n=n,d=2,init=c(0,0))
  theta_samples=th$sim_vals

#
# decentering
#
  theta_samples[,1]=theta_samples[,1]-theta_samples[,2]*meant

	list(theta_samples=theta_samples)

}
