#' Gamma Distribution Predictions Based on a Calibrating Prior
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
#' The Gamma distribution has probability density function
#' \deqn{f(x;\alpha,\sigma)=\frac{1}{\sigma^\alpha\Gamma(\alpha)}x^{\alpha-1}e^{-x/\sigma}}
#' where
#' \eqn{x \ge 0} is the random variable and
#' \eqn{\alpha>0,\sigma>0} are the parameters.
#'
#' The calibrating prior we use is
#' \deqn{\pi(\alpha,\sigma) \propto \frac{1}{\alpha \sigma}}
#'
#' @example man/examples/example_100_gamma.R
#'
#' @name gamma_cp
NULL
#' @rdname gamma_cp
#' @inheritParams man
#' @export
#'
qgamma_cp=function(x,p=seq(0.1,0.9,0.1),fd1=0.01,fd2=0.01,
	means=FALSE,waicscores=FALSE,logscores=FALSE,dmgs=TRUE,
	rust=FALSE,nrust=100000,prior="type 1",debug=FALSE,aderivs=TRUE){
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
	opt=optim(c(1,1),gamma_loglik,x=x,control=list(fnscale=-1))
	v1hat=opt$par[1]
	v2hat=opt$par[2]
	ml_params=c(v1hat,v2hat)
	if(debug)cat("  v1hat,v2hat=",v1hat,v2hat,"\n")
#
# 3 aic
#
	ml_value=opt$val
	maic=make_maic(ml_value,nparams=2)
#
# 4 ml quantiles (vectorized over alpha)
#
	ml_quantiles=qgamma((1-alpha),shape=v1hat,scale=v2hat)
#
# test of gg code (for future implementation of mpd theory, as a test of the mpd code)
#
#	gamma_gg(v1hat,fd1,v2hat,fd2)
#
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
		if(aderivs) ldd=gamma_ldda(x,v1hat,v2hat)
		if(!aderivs)ldd=gamma_ldd(x,v1hat,fd1,v2hat,fd2)
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
#
# 7 lddd
#
		if(debug)cat("  calculate lddd\n")
		if(aderivs) lddd=gamma_lddda(x,v1hat,v2hat)
		if(!aderivs)lddd=gamma_lddd(x,v1hat,fd1,v2hat,fd2)
#
# 8 mu1
#
		if(debug)cat("calculate mu1\n")
		mu1=gamma_mu1f(alpha,v1hat,fd1,v2hat,fd2)
#
# 9 mu2
#
		if(debug)cat("calculate mu2\n")
		mu2=gamma_mu2f(alpha,v1hat,fd1,v2hat,fd2)
#
# 10 q cp
# (I compared these two priors, and the double prior worked better)
# (so I made the better one type 1)
# (see the actuary paper)
		if(debug)cat("  cp\n")
		if(prior=="type 1"){
			lambdad_cp=c(-1/v1hat,-1/v2hat)
		} else if (prior=="type 2"){
			lambdad_cp=c(0,-1/v2hat)
		} else {
			cat("invalid prior choice.\n")
			stop()
		}
##		lambdad_cp=c(0,-1/v2hat) 				#this worked ok
##		lambdad_cp=c(-1/v1hat,-1/v2hat)	#but this worked better
#
# 11 derive the bayesian dq based on v2hat
#
		if(debug)cat("  fhat, dq and cp quantiles\n")
		fhat=dgamma(ml_quantiles,shape=v1hat,scale=v2hat)
		dq=dmgs(lddi,lddd,mu1,lambdad_cp,mu2,dim=2)
		cp_quantiles=ml_quantiles+dq/(nx*fhat)
#
# 12 means
#
		means=gamma_means(means,ml_params,lddi,lddd,lambdad_cp,nx,dim=2)
		ml_mean=means$ml_mean
		cp_mean=means$cp_mean
#
# 13 waicscores
#
		waic=gamma_waic(waicscores,x,v1hat,fd1,v2hat,fd2,lddi,lddd,
			lambdad_cp,aderivs)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 14 logscores
#
		logscores=gamma_logscores(logscores,x,fd1,fd2,aderivs)
		ml_oos_logscore=logscores$ml_oos_logscore
		cp_oos_logscore=logscores$cp_oos_logscore
#
# 15 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rgamma_cp(nrust,x,rust=TRUE,mlcp=FALSE)
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
#' @rdname gamma_cp
#' @inheritParams man
#' @export
rgamma_cp=function(n,x,fd1=0.01,fd2=0.01,rust=FALSE,mlcp=TRUE,
	debug=FALSE,aderivs=TRUE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),!x<0)
	stopifnot(is.finite(x),!is.na(x),!x<0)

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qgamma_cp(x,runif(n),fd1=fd1,fd2=fd2,aderivs=aderivs)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tgamma_cp(n,x)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			ru_deviates[i]=rgamma(1,shape=th[i,1],scale=th[i,2])
		}
	}

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
				cp_method=adhoc_dmgs_cpmethod())
	return(op)

}
#' @rdname gamma_cp
#' @inheritParams man
#' @export
dgamma_cp=function(x,y=x,fd1=0.01,fd2=0.01,rust=FALSE,nrust=1000,
	debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),!x<0,!y<0)

	dd=dgammasub(x=x,y=y,fd1,fd2,aderivs=aderivs)
	ru_pdf="rust not selected"
	if(rust){
		th=tgamma_cp(nrust,x)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_pdf=ru_pdf+dgamma(y,shape=th[ir,1],scale=th[ir,2])
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
#' @rdname gamma_cp
#' @inheritParams man
#' @export
pgamma_cp=function(x,y=x,fd1=0.01,fd2=0.01,rust=FALSE,nrust=1000,
	debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),!x<0,!y<0)

	dd=dgammasub(x=x,y=y,fd1,fd2,aderivs=aderivs)
	ru_cdf="rust not selected"
	if(rust){
		th=tgamma_cp(nrust,x)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_cdf=ru_cdf+pgamma(y,shape=th[ir,1],scale=th[ir,2])
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
#' @rdname gamma_cp
#' @inheritParams man
#' @export
tgamma_cp=function(n,x,fd1=0.01,fd2=0.01,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),!x<0)
	stopifnot(is.finite(x),!is.na(x),!x<0)

	t=ru(gamma_logf,x=x,n=n,d=2,init=c(1,1))

	list(theta_samples=t$sim_vals)

}

