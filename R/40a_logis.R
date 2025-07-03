#' Logistic Distribution Predictions Based on a Calibrating Prior
#'
#' @inherit man description author references seealso return
#' @inheritParams man
#'
#' @inheritSection man Optional Return Values
# #' @inheritSection man Optional Return Values (EVD models only)
# #' @inheritSection man Optional Return Values (non-RHP models only)
# #' @inheritSection man Details (homogeneous models)
#' @inheritSection man Details (non-homogeneous models)
# #' @inheritSection man Details (analytic integration)
#' @inheritSection man Details (DMGS integration)
#' @inheritSection man Details (RUST)
#'
#' @section Details of the Model:
#' The logistic distribution has distribution function
#' \deqn{f(x;\mu,\sigma)=\frac{1}{1+e^{-(x-\mu)/\sigma}}}
#' where
#' \eqn{x} is the random variable
#' and
#' \eqn{\mu,\sigma>0} are the parameters.
#'
#' The calibrating prior is given by the right Haar prior, which is
#' \deqn{\pi(\sigma) \propto \frac{1}{\sigma}}
#' as given in Jewson et al. (2025).
#'
#' @example man/examples/example_40_logis.R
#'
#' @name logis_cp
NULL
#' @rdname logis_cp
#' @inheritParams man
#' @export
#'
qlogis_cp=function(x,p=seq(0.1,0.9,0.1),d1=0.01,fd2=0.01,
	means=FALSE,waicscores=FALSE,logscores=FALSE,dmgs=TRUE,rust=FALSE,nrust=100000,
	debug=FALSE,aderivs=TRUE){
#
# 1 intro
#
	debug=TRUE
	debug=FALSE
	stopifnot(is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1)
	alpha=1-p
	nx=length(x)
	nalpha=length(alpha)
#
# 2 ml param estimate
#
	opt=optim(c(mean(x),sd(x)),logis_loglik,x=x,control=list(fnscale=-1))
	v1hat=opt$par[1]
	v2hat=opt$par[2]
	ml_params=c(v1hat,v2hat)
	if(debug)message("  v1hat,v2hat=",v1hat,v2hat,"")
#
# 3 aic
#
	ml_value=opt$val
	maic=make_maic(ml_value,nparams=2)

#
# 4 ml quantiles (vectorized over alpha)
#
	ml_quantiles=qlogis((1-alpha),location=v1hat,scale=v2hat)
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
# 5 lddi
#
		if(aderivs)	ldd=logis_ldda(x,v1hat,v2hat)
		if(!aderivs)ldd=logis_ldd(x,v1hat,d1,v2hat,fd2)
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
#
# 6 lddd
#
		if(debug)message("  calculate lddd")
		if(aderivs)	lddd=logis_lddda(x,v1hat,v2hat)
		if(!aderivs)lddd=logis_lddd(x,v1hat,d1,v2hat,fd2)
#
# 7 mu1
#
		if(debug)message("calculate mu1")
		if(aderivs) mu1=logis_mu1fa(alpha,v1hat,v2hat)
		if(!aderivs)mu1=logis_mu1f(alpha,v1hat,d1,v2hat,fd2)
#
# 8 mu2
#
		if(debug)message("calculate mu2")
		if(aderivs) mu2=logis_mu2fa(alpha,v1hat,v2hat)
		if(!aderivs)mu2=logis_mu2f(alpha,v1hat,d1,v2hat,fd2)
#
# 9 q rhp
#
		if(debug)message("  rhp")
		lambdad_rhp=c(0,-1/v2hat)
#
# 10 derive the bayesian dq based on v2hat
#
		if(debug)message("  fhat, dq and rhp quantiles")
		fhat=dlogis(ml_quantiles,location=v1hat,scale=v2hat)
		dq=dmgs(lddi,lddd,mu1,lambdad_rhp,mu2,dim=2)
		rh_quantiles=ml_quantiles+dq/(nx*fhat)
#
# 11 means (might as well always calculate)
#
		ml_mean=ml_params[1]
		rh_mean=ml_params[1]
#
# 12 waicscores
#
		waic=logis_waic(waicscores,x,v1hat,d1,v2hat,fd2,lddi,lddd,lambdad_rhp,aderivs)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 13 logscores
#
		logscores=logis_logscores(logscores,x,d1,fd2,aderivs)
		ml_oos_logscore=logscores$ml_oos_logscore
		rh_oos_logscore=logscores$rh_oos_logscore
#
# 14 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rlogis_cp(nrust,x,rust=TRUE,mlcp=FALSE)
			ru_quantiles=makeq(rustsim$ru_deviates,p)
		}
	}	# end of if(dmgs)

# return
	list(	ml_params=ml_params,
				ml_value=ml_value,
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
				cp_method=rhp_dmgs_cpmethod())

}
#' @rdname logis_cp
#' @inheritParams man
#' @export
rlogis_cp=function(n,x,d1=0.01,fd2=0.01,rust=FALSE,mlcp=TRUE,
	debug=FALSE,aderivs=TRUE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x))
	stopifnot(is.finite(x),!is.na(x))

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qlogis_cp(x,runif(n),d1=d1,fd2=fd2,aderivs=aderivs)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tlogis_cp(n,x)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			ru_deviates[i]=rlogis(1,location=th[i,1],scale=th[i,2])
		}
	}
	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
				cp_method=rhp_dmgs_cpmethod())
	return(op)

}
#' @rdname logis_cp
#' @inheritParams man
#' @export
dlogis_cp=function(x,y=x,d1=0.01,fd2=0.01,rust=FALSE,nrust=1000,
	debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y))

	dd=dlogis2sub(x=x,y=y,d1,fd2,aderivs=aderivs)
	ru_pdf="rust not selected"
	if(rust){
		th=tlogis_cp(nrust,x)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_pdf=ru_pdf+dlogis(y,location=th[ir,1],scale=th[ir,2])
		}
		ru_pdf=ru_pdf/nrust
	}
	op=list(	ml_params=dd$ml_params,
					ml_pdf=dd$ml_pdf,
					cp_pdf=dd$rh_pdf,
					ru_pdf=ru_pdf,
					cp_method=rhp_dmgs_cpmethod())
	return(op)
}
#' @rdname logis_cp
#' @inheritParams man
#' @export
plogis_cp=function(x,y=x,d1=0.01,fd2=0.01,rust=FALSE,nrust=1000,
	debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y))

	dd=dlogis2sub(x=x,y=y,d1,fd2,aderivs=aderivs)
	ru_cdf="rust not selected"
	if(rust){
		th=tlogis_cp(nrust,x)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_cdf=ru_cdf+plogis(y,location=th[ir,1],scale=th[ir,2])
		}
		ru_cdf=ru_cdf/nrust
	}
	op=list(	ml_params=dd$ml_params,
					ml_cdf=dd$ml_cdf,
					cp_cdf=dd$rh_cdf,
					ru_cdf=ru_cdf,
					cp_method=rhp_dmgs_cpmethod())
	return(op)
}
#' @rdname logis_cp
#' @inheritParams man
#' @export
tlogis_cp=function(n,x,d1=0.01,fd2=0.01,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x))
	stopifnot(is.finite(x),!is.na(x))

	th=ru(logis_logf,x=x,n=n,d=2,init=c(mean(x),sd(x)))

	list(theta_samples=th$sim_vals)
}
