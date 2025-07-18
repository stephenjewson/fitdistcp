#' Weibull Distribution Predictions Based on a Calibrating Prior
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
#' The Weibull distribution has exceedance distribution function
#' \deqn{S(x;k,\sigma)=\exp\left(-\left(\frac{x}{\sigma}\right)^{k}\right)}
#' where
#' \eqn{x \ge 0} is the random variable and
#' \eqn{k>0,\sigma>0} are the parameters.
#'
#' The calibrating prior is given by the right Haar prior, which is
#' \deqn{\pi(k,\sigma) \propto \frac{1}{k \sigma}}
#' as given in Jewson et al. (2025).
#'
#' @example man/examples/example_052_weibull.R
#'
#' @name weibull_cp
NULL
#' @rdname weibull_cp
#' @inheritParams man
#' @export
#'
qweibull_cp=function(x,p=seq(0.1,0.9,0.1),
	means=FALSE,waicscores=FALSE,logscores=FALSE,dmgs=TRUE,rust=FALSE,nrust=100000,
	debug=FALSE){
#
# 1 intro
#
	debug=TRUE
	debug=FALSE
	stopifnot(is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1,!x<0)
	alpha=1-p
	nx=length(x)
	nalpha=length(alpha)
#
# 2 ml param estimate
#
	opt=optim(c(1,1),weibull_loglik,x=x,control=list(fnscale=-1))
	v1hat=opt$par[1]
	v2hat=opt$par[2]
	ml_params=c(v1hat,v2hat)
	if(debug)message("  v1hat,v2hat=",v1hat,v2hat)
#
# 3 aic
#
	ml_value=opt$val
	maic=make_maic(ml_value,nparams=2)
#
# 4 ml quantiles (vectorized over alpha)
#
	ml_quantiles=qweibull((1-alpha),shape=v1hat,scale=v2hat)
#
# 5 dmgs
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
# 6 lddi
#
		ldd=weibull_ldda(x,v1hat,v2hat)
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
#
# 7 lddd
#
		if(debug)message("  calculate lddd")
		lddd=weibull_lddda(x,v1hat,v2hat)
#
# 7 mu1
#
		if(debug)message("calculate mu1")
		mu1=weibull_mu1fa(alpha,v1hat,v2hat)

# 8 mu2
#
		if(debug)message("calculate mu2")
		mu2=weibull_mu2fa(alpha,v1hat,v2hat)

# 9 q rhp
#
		if(debug)message("  rhp")
		lambdad_rhp=c(-1/v1hat,-1/v2hat)
#
# 10 derive the bayesian dq based on v2hat
#
		if(debug)message("  fhat, dq and rhp quantiles")
		fhat=dweibull(ml_quantiles,shape=v1hat,scale=v2hat)
		dq=dmgs(lddi,lddd,mu1,lambdad_rhp,mu2,dim=2)
		rh_quantiles=ml_quantiles+dq/(nx*fhat)
#
# 11 means
#
		means=weibull_means(means,ml_params,lddi,lddd,lambdad_rhp,nx,dim=2)
		ml_mean=means$ml_mean
		rh_mean=means$rh_mean
#
# 12 waicscores
#
		waic=weibull_waic(waicscores,x,v1hat,v2hat,lddi,lddd,
			lambdad_rhp)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 13 logscores
#
		logscores=weibull_logscores(logscores,x)
		ml_oos_logscore=logscores$ml_oos_logscore
		rh_oos_logscore=logscores$rh_oos_logscore
#
# 14 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rweibull_cp(nrust,x,rust=TRUE,mlcp=FALSE)
			ru_quantiles=makeq(rustsim$ru_deviates,p)
		}
	} #end of if dmgs

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
#' @rdname weibull_cp
#' @inheritParams man
#' @export
rweibull_cp=function(n,x,rust=FALSE,mlcp=TRUE,
	debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),!x<0)
	stopifnot(is.finite(x),!is.na(x),!x<0)

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qweibull_cp(x,runif(n))
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tweibull_cp(n,x)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			ru_deviates[i]=rweibull(1,shape=th[i,1],scale=th[i,2])
		}
	}

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
				cp_method=rhp_dmgs_cpmethod())
	return(op)

}
#' @rdname weibull_cp
#' @inheritParams man
#' @export
dweibull_cp=function(x,y=x,rust=FALSE,nrust=1000,
	debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),!x<0,!y<0)

	dd=dweibullsub(x=x,y=y)
	ru_pdf="rust not selected"
	if(rust){
		th=tweibull_cp(nrust,x)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_pdf=ru_pdf+dweibull(y,shape=th[ir,1],scale=th[ir,2])
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
#' @rdname weibull_cp
#' @inheritParams man
#' @export
pweibull_cp=function(x,y=x,rust=FALSE,nrust=1000,
	debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),!x<0,!y<0)

	dd=dweibullsub(x=x,y=y)
	ru_cdf="rust not selected"
	if(rust){
		th=tweibull_cp(nrust,x)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_cdf=ru_cdf+pweibull(y,shape=th[ir,1],scale=th[ir,2])
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
#' @rdname weibull_cp
#' @inheritParams man
#' @export
tweibull_cp=function(n,x,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),!x<0)
	stopifnot(is.finite(x),!is.na(x),!x<0)

	t=ru(weibull_logf,x=x,n=n,d=2,init=c(1,1))

	list(theta_samples=t$sim_vals)

}

