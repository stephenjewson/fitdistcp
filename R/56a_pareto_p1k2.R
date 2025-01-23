# extraDistr is very confusing:
# a = shape parameter (they call scale), that I'm varying here, so which is parameters 1 and 2
# b = scale parameter (they call location), that I'm calling kscale, which is parameter 3
# so it's p1k3 (shape has a predictor, scale is known)
#' Pareto Distribution with a Predictor, Predictions Based on a Calibrating Prior
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
# #' @inheritSection man Details (analytic integration)
#' @inheritSection man Details (DMGS integration)
#' @inheritSection man Details (RUST)
#'
#' @section Details of the Model:
#' The Pareto distribution with a predictor has various forms.
#' The form we are using has exceedance distribution function
#' \deqn{S(x;a,b)={\left(\frac{\sigma}{x}\right)^{\alpha(a,b)}}}
#' where
#' \eqn{x \ge \sigma} is the random variable,
#' \eqn{\alpha=\exp(-a-bt)} is the shape parameter, modelled
#' as a function of parameters \eqn{a,b}, and \eqn{\sigma}
#' is the scale parameter.
#' We consider the scale parameter \eqn{\sigma} to be known
#' (hence the \code{k2} in the name).
#'
#' The calibrating prior is given by the right Haar prior, which is
#' \deqn{\pi(a,b) \propto 1}
#' as given in Jewson et al. (2024).
#' Note that others authors have referred to the shape and scale parameters
#' as the scale and location parameters, respectively.
#'
#' @example man/examples/example_56_pareto_p1k2.R
#'
#' @name pareto_p1k2_cp
NULL
#' @rdname pareto_p1k2_cp
#' @inheritParams man
#' @export
#'
qpareto_p1k2_cp=function(x,t,t0=NA,n0=NA,p=seq(0.1,0.9,0.1),d1=0.01,d2=0.01,kscale=1,
	means=FALSE,waicscores=FALSE,logscores=FALSE,dmgs=TRUE,rust=FALSE,nrust=100000,
	predictordata=TRUE,centering=TRUE,
	debug=FALSE,aderivs=TRUE){
#
# 1 intro
#
	if(debug)cat("inside qpareto_p1k2\n")
	stopifnot(	is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1,!x<kscale)
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
	if(debug)cat("calc ml param estimate\n")
	lm=lm(x~t)
	v1start=lm$coefficients[1]
	v2start=lm$coefficients[2]
	v1start=0
	v2start=1
	xhat=v1start+v2start*t
	opt=optim(c(v1start,v2start),pareto_p1k2_loglik,x=x,t=t,
		kscale=kscale,control=list(fnscale=-1))
	v1hat=opt$par[1]
	v2hat=opt$par[2]
	ml_params=c(v1hat,v2hat)
	muhat=ml_params[1]+ml_params[2]*t
#	muhat0=makemuhat0(t0,n0,t,ml_params)
	residuals=x-muhat
	if(debug)cat("  inside q: ml_params=",ml_params,"\n")
#
# 4 predictordata
#
	prd=pareto_p1k2_predictordata(predictordata,x,t,t0,ml_params,kscale)
	predictedparameter=prd$predictedparameter
	adjustedx=prd$adjustedx
#
# 5 aic
#
	if(debug)cat("calc aic\n")
 	ml_value=opt$val
	maic=make_maic(ml_value,nparams=2)

#
# 6 mle quantiles
#
	if(debug)cat("calc mle quantiles\n")
	ml_quantiles=qpareto_p1k2((1-alpha),t0,ymn=v1hat,slope=v2hat,kscale=kscale)
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
		if(debug)cat("calc ldd\n")
		if(aderivs) ldd=pareto_p1k2_ldda(x,t,v1hat,v2hat,kscale)
		if(!aderivs)ldd=pareto_p1k2_ldd(x,t,v1hat,d1,v2hat,d2,kscale)
		if(debug)cat("ldd=",ldd,"\n")
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
#
# 8 lddd
#
		if(debug)cat("calculate lddd\n")
		if(aderivs) lddd=pareto_p1k2_lddda(x,t,v1hat,v2hat,kscale)
		if(!aderivs)lddd=pareto_p1k2_lddd(x,t,v1hat,d1,v2hat,d2,kscale)
#
# 9 mu1
#
		if(debug)cat("calculate mu1\n")
		if(aderivs) mu1=pareto_p1k2_mu1fa(alpha,t0,v1hat,v2hat,kscale)
		if(!aderivs)mu1=pareto_p1k2_mu1f(alpha,t0,v1hat,d1,v2hat,d2,kscale)
#
# 10 mu2
#
		if(debug)cat("calculate mu2\n")
		if(aderivs) mu2=pareto_p1k2_mu2fa(alpha,t0,v1hat,v2hat,kscale)
		if(!aderivs)mu2=pareto_p1k2_mu2f(alpha,t0,v1hat,d1,v2hat,d2,kscale)
#
# 11 rhp
#
		if(debug)cat("  rhp\n")
		lambdad_rhp=c(0,0) #this is rhp
#
# 12 rhp quantiles
#
		if(debug)cat("  rhp quantiles\n")
		fhat=dpareto_p1k2(ml_quantiles,t0,ymn=v1hat,slope=v2hat,kscale=kscale,log=FALSE)
		dq=dmgs(lddi,lddd,mu1,lambdad_rhp,mu2,dim=2)
		rh_quantiles=ml_quantiles+dq/(nx*fhat)
#
# 13 means
#
		means=pareto_p1k2_means(means,t0,ml_params,lddi,lddd,lambdad_rhp,nx,dim=2,kscale)
		ml_mean=means$ml_mean
		rh_mean=means$rh_mean
#
# 14 waicscores
#
		waic=pareto_p1k2_waic(waicscores,x,t,v1hat,d1,v2hat,d2,kscale,lddi,lddd,lambdad_rhp)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 15 logscores
#
		if(debug)cat("calc logscores\n")
		logscores=pareto_p1k2_logscores(logscores,x,t,d1,d2,kscale,aderivs,debug)
		ml_oos_logscore=logscores$ml_oos_logscore
		rh_oos_logscore=logscores$rh_oos_logscore
#
# 16 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rpareto_p1k2_cp(nrust,x,t=t,t0=t0,kscale=kscale,rust=TRUE,mlcp=FALSE)
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
	if(debug)cat("return\n")
	list(	ml_params=ml_params,
				ml_value=ml_value,
				predictedparameter=predictedparameter,
				adjustedx=adjustedx,
				ldd=ldd,
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
				cp_method=rhp_dmgs_cpmethod())

}
#' @rdname pareto_p1k2_cp
#' @inheritParams man
#' @export
rpareto_p1k2_cp=function(n,x,t,t0=NA,n0=NA,d1=0.01,d2=0.01,
	kscale=1,rust=FALSE,mlcp=TRUE,centering=TRUE,
	debug=FALSE,aderivs=TRUE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),
#		is.finite(t),!is.na(t),!x<kscale)
	stopifnot(is.finite(x),!is.na(x),
		is.finite(t),!is.na(t),!x<kscale)

	t0=maket0(t0,n0,t)

#
# centering
#
	if(centering){
   meant=mean(t)
	 t=t-meant
	 t0=t0-meant
	}

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qpareto_p1k2_cp(x,t,t0=t0,n0=NA,p=runif(n),d1,d2,kscale=kscale,
			centering=centering,aderivs=aderivs)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
		if(debug)cat("  inside r: ml_params=",ml_params,"\n")
	}

	if(rust){
		th=tpareto_p1k2_cp(n,x,t,kscale)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			mu=th[i,1]+t0*th[i,2]
			ru_deviates[i]=rpareto(1,a=1/exp(mu),b=kscale)
		}
	}

#
# decentering
#
  if(mlcp&centering)ml_params[1]=ml_params[1]-ml_params[2]*meant

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
				cp_method=rhp_dmgs_cpmethod())

	return(op)

}
#' @rdname pareto_p1k2_cp
#' @inheritParams man
#' @export
dpareto_p1k2_cp=function(x,t,t0=NA,n0=NA,y=x,d1=0.01,d2=0.01,
	kscale=1,rust=FALSE,nrust=1000,centering=TRUE,
	debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),
		is.finite(y),!is.na(y),is.finite(t),!is.na(t),!x<kscale,!y<kscale)

	t0=maket0(t0,n0,t)
#
# centering
#
  if(centering){
    meant=mean(t)
    t=t-meant
    t0=t0-meant
  }

	dd=dpareto_p1k2sub(x=x,t=t,y=y,t0=t0,d1,d2,kscale=kscale,aderivs=aderivs)
	ru_pdf="rust not selected"
	ml_params=dd$ml_params
	if(rust){
		th=tpareto_p1k2_cp(nrust,x,t,kscale=kscale)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=exp(th[ir,1]+t0*th[ir,2])
			ru_pdf=ru_pdf+dpareto(y,a=1/mu,b=kscale)
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
					ru_pdf=ru_pdf,
					cp_method=rhp_dmgs_cpmethod())
	return(op)
}
#' @rdname pareto_p1k2_cp
#' @inheritParams man
#' @export
ppareto_p1k2_cp=function(x,t,t0=NA,n0=NA,y=x,d1=0.01,d2=0.01,kscale=1,
				rust=FALSE,nrust=1000,centering=TRUE,
	debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
		is.finite(t),!is.na(t),!x<kscale,!y<kscale)

	t0=maket0(t0,n0,t)
#
# centering
#
  if(centering){
    meant=mean(t)
    t=t-meant
    t0=t0-meant
  }

	dd=dpareto_p1k2sub(x=x,t=t,y=y,t0=t0,d1,d2,kscale=kscale,aderivs=aderivs)
	ru_cdf="rust not selected"
	ml_params=dd$ml_params
	if(rust){
		th=tpareto_p1k2_cp(nrust,x,t,kscale=kscale)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=exp(th[ir,1]+t0*th[ir,2])
			ru_cdf=ru_cdf+ppareto(y,a=1/mu,b=kscale)
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
					cp_method=rhp_dmgs_cpmethod())
	return(op)
}
#' @rdname pareto_p1k2_cp
#' @inheritParams man
#' @export
tpareto_p1k2_cp=function(n,x,t,d1=0.01,d2=0.01,kscale=1,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),
#		is.finite(t),!is.na(t),!x<kscale)
	stopifnot(is.finite(x),!is.na(x),
		is.finite(t),!is.na(t),!x<kscale)

#
# centering
#
  meant=mean(t)
  t=t-meant

	th=ru(pareto_p1k2_logf,x=x,t=t,kscale=kscale,n=n,d=2,init=c(0,0))
  theta_samples=th$sim_vals
#
# decentering
#
  theta_samples[,1]=theta_samples[,1]-theta_samples[,2]*meant

	list(theta_samples=theta_samples)

}
