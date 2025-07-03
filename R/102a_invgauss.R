#' Inverse Gauss Distribution, Predictions Based on a Calibrating Prior
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
#' The Inverse Gaussian distribution has probability density function
#' \deqn{f(x;\mu,\phi)=\left(\frac{1}{2\pi\phi x^3}\right)^{1/2}
#' \exp\left(-\frac{(x-\mu)^2}{2\mu^2\phi x}\right)}
#' where
#' \eqn{x \ge 0} is the random variable and
#' \eqn{\mu>0,\phi>0} are the parameters.
#'
#' The calibrating prior we use by default is
#' \deqn{\pi(\alpha,\sigma) \propto \frac{1}{\phi}}
#' The prior
#' \deqn{\pi(\alpha,\sigma) \propto \frac{1}{\mu\phi}}
#' is also available as an option with \code{prior="type 2"}.
#'
#' @example man/examples/example_102_invgauss.R
#'
#' @name invgauss_cp
NULL
#' @rdname invgauss_cp
#' @inheritParams man
#' @export
#'
qinvgauss_cp=function(x,p=seq(0.1,0.9,0.1),fd1=0.01,fd2=0.01,
	means=FALSE,waicscores=FALSE,logscores=FALSE,dmgs=TRUE,
	rust=FALSE,nrust=100000,prior="type 1",debug=FALSE,aderivs=TRUE){
#
# 1 intro
#
	stopifnot(is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1,!x<0)
	alpha=1-p
	nx=length(x)
	nalpha=length(alpha)
#
# 2 ml param estimate
#
	opt=optim(c(1,1),invgauss_loglik,x=x,control=list(fnscale=-1))
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
	ml_quantiles=qinvgauss((1-alpha),mean=v1hat,shape=v2hat)
#
# 5 dmgs
#
	standard_errors="dmgs not selected"
	cp_quantiles="dmgs not selected"
	ru_quantiles="dmgs not selected"
	waic1="dmgs not selected"
	waic2="dmgs not selected"
	ml_oos_logscore="dmgs not selected"
	cp_oos_logscore="dmgs not selected"
	cp_oos_logscore="dmgs not selected"
	ml_mean="dmgs not selected"
	cp_mean="dmgs not selected"
	cp_mean="dmgs not selected"
	cp_method="dmgs not selected"
	if(dmgs){
#
# 6 lddi
#
		if(aderivs) ldd=invgauss_ldda(x,v1hat,v2hat)
		if(!aderivs)ldd=invgauss_ldd(x,v1hat,fd1,v2hat,fd2)
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
#
# 7 lddd
#
		if(debug)message("  calculate lddd")
		if(aderivs) lddd=invgauss_lddda(x,v1hat,v2hat)
		if(!aderivs)lddd=invgauss_lddd(x,v1hat,fd1,v2hat,fd2)
#
# 8 mu1
#
		if(debug)message("calculate mu1")
		mu1=invgauss_mu1f(alpha,v1hat,fd1,v2hat,fd2)
#
# 9 mu2
#
		if(debug)message("calculate mu2")
		mu2=invgauss_mu2f(alpha,v1hat,fd1,v2hat,fd2)
#
# 10 q cp
# the tests I did for the actuary paper show that prior1 does better
# for sample sizes of 20, but prior2 does better for sample sizes of 80
# so I let the user choose both, pending a better system
#
		if(debug)message("  cp")
		if(prior=="type 1"){
			lambdad_cp=c(-1/v1hat,-1/v2hat)
		} else if (prior=="type 2"){
			lambdad_cp=c(0,-1/v2hat)
		} else {
			message("invalid prior choice.")
			stop()
		}
#
# 11 derive the bayesian dq based on v2hat
#
		if(debug)message("  fhat, dq and cp quantiles")
		fhat=dinvgauss(ml_quantiles,mean=v1hat,shape=v2hat)
		dq=dmgs(lddi,lddd,mu1,lambdad_cp,mu2,dim=2)
		cp_quantiles=ml_quantiles+dq/(nx*fhat)
#
# 12 means
#
		means=invgauss_means(means,ml_params,lddi,lddd,lambdad_cp,nx,dim=2)
		ml_mean=means$ml_mean
		cp_mean=means$cp_mean
#
# 13 waicscores
#
		waic=invgauss_waic(waicscores,x,v1hat,fd1,v2hat,fd2,lddi,lddd,
			lambdad_cp,aderivs)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 14 logscores
#
		logscores=invgauss_logscores(logscores,x,prior,fd1,fd2,aderivs)
		ml_oos_logscore=logscores$ml_oos_logscore
		cp_oos_logscore=logscores$cp_oos_logscore
#
# 15 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rinvgauss_cp(nrust,x,rust=TRUE,prior=prior,mlcp=FALSE)
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
				cp_quantiles=cp_quantiles,
				ru_quantiles=ru_quantiles,
				maic=maic,
				waic1=waic1,
				waic2=waic2,
				ml_oos_logscore=ml_oos_logscore,
				cp_oos_logscore=cp_oos_logscore,
				ml_mean=ml_mean,
				cp_mean=cp_mean,
				cp_method=adhoc_dmgs_cpmethod())

}
#' @rdname invgauss_cp
#' @inheritParams man
#' @export
rinvgauss_cp=function(n,x,fd1=0.01,fd2=0.01,rust=FALSE,prior="type 1",mlcp=TRUE,
	debug=FALSE,aderivs=TRUE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),!x<0)
	stopifnot(is.finite(x),!is.na(x),!x<0)

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qinvgauss_cp(x,runif(n),fd1=fd1,fd2=fd2,prior=prior,aderivs=aderivs)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tinvgauss_cp(n,x,prior=prior)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			ru_deviates[i]=rinvgauss(1,mean=th[i,1],shape=th[i,2])
		}
	}

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
				cp_method=adhoc_dmgs_cpmethod())
	return(op)

}
#' @rdname invgauss_cp
#' @inheritParams man
#' @export
dinvgauss_cp=function(x,y=x,fd1=0.01,fd2=0.01,rust=FALSE,nrust=1000,
	prior="type 1",debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),!x<0,!y<0)

	if(debug)message("before dd")
	dd=dinvgausssub(x=x,y=y,prior,fd1,fd2,aderivs=aderivs)
	ru_pdf="rust not selected"
	if(debug)message("before rust")
	if(rust){
		if(debug)message("before th=")
		th=tinvgauss_cp(nrust,x)$theta_samples
		ru_pdf=numeric(length(y))
		if(debug)message("before ru_pdf=")
		for (ir in 1:nrust){
			ru_pdf=ru_pdf+dinvgauss(y,mean=th[ir,1],shape=th[ir,2])
		}
		ru_pdf=ru_pdf/nrust
	}
	op=list(	ml_params=dd$ml_params,
					ml_pdf=dd$ml_pdf,
					cp_pdf=dd$cp_pdf,
					ru_pdf=ru_pdf,
					cp_method=adhoc_dmgs_cpmethod())
	return(op)
}
#' @rdname invgauss_cp
#' @inheritParams man
#' @export
pinvgauss_cp=function(x,y=x,fd1=0.01,fd2=0.01,rust=FALSE,nrust=1000,
	prior="type 1",debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),!x<0,!y<0)

	dd=dinvgausssub(x=x,y=y,prior,fd1,fd2,aderivs=aderivs)
	ru_cdf="rust not selected"
	if(rust){
		th=tinvgauss_cp(nrust,x)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_cdf=ru_cdf+pinvgauss(y,mean=th[ir,1],shape=th[ir,2])
		}
		ru_cdf=ru_cdf/nrust
	}
	op=list(	ml_params=dd$ml_params,
					ml_cdf=dd$ml_cdf,
					cp_cdf=dd$cp_cdf,
					ru_cdf=ru_cdf,
					cp_method=adhoc_dmgs_cpmethod())
	return(op)
}
#' @rdname invgauss_cp
#' @inheritParams man
#' @export
tinvgauss_cp=function(n,x,fd1=0.01,fd2=0.01,prior="type 1",debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),!x<0)
	stopifnot(is.finite(x),!is.na(x),!x<0)

	th=ru(invgauss_logf,x=x,prior=prior,n=n,d=2,init=c(1,1))

	list(theta_samples=th$sim_vals)

}

