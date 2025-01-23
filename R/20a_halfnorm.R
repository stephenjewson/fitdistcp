#' Half-Normal Distribution Predictions Based on a Calibrating Prior
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
#' The half-normal distribution has probability density function
#' \deqn{f(x;\theta)=\frac{2\theta}{\pi}e^{-\theta^2 x^2/\pi}}
#' where
#' \eqn{x \ge 0} is the random variable
#' and
#' \eqn{\theta >0} is the parameter.
#'
#' The calibrating prior is given by the right Haar prior, which is
#' \deqn{\pi(\theta) \propto \frac{1}{\theta}}
#' as given in Jewson et al. (2024).
#' Some other authors may parametrize the half-normal differently.
#'
#' @example man/examples/example_20_halfnorm.R
#'
#' @name halfnorm_cp
NULL
#' @rdname halfnorm_cp
#' @inheritParams man
#' @export
#'
qhalfnorm_cp=function(x,p=seq(0.1,0.9,0.1),fd1=0.01,
	means=FALSE,waicscores=FALSE,logscores=FALSE,dmgs=TRUE,rust=FALSE,nrust=100000,debug=FALSE,
	aderivs=TRUE){
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
	if(debug)cat("2 calc ml param estimate")
	v1start=sqrt((sum(x*x))/nx)
	opt=optim(c(v1start),halfnorm_loglik,x=x,method="Brent",
		lower=.Machine$double.eps,upper=999999999999,control=list(fnscale=-1))
	v1hat=opt$par[1]
	ml_params=c(v1hat)
	if(debug)cat("  v1start,v1hat=",v1start,v1hat,"//")
#
# 3 aic
#
	ml_value=opt$val
	maic=make_maic(ml_value,nparams=1)
#
# 4 ml quantiles (vectorized over alpha)
#
	ml_quantiles=qhalfnorm((1-alpha),theta=v1hat)
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
		if(debug)cat("  calculate ldd,lddi\n")
		if(aderivs)	ldd=halfnorm_ldda(x,v1hat)
		if(!aderivs)ldd=halfnorm_ldd(x,v1hat,fd1)

		lddi=solve(ldd)
		expinfmat=halfnorm_gg(v1hat,fd1)
		expinfmati=solve(expinfmat)
		standard_errors=sqrt(expinfmati/nx)
#
# 6 lddd
#
		if(debug)cat("  calculate lddd\n")
		if(aderivs)	lddd=halfnorm_lddda(x,v1hat)
		if(!aderivs)lddd=halfnorm_lddd(x,v1hat,fd1)
#
# 7 mu1
#
		if(debug)cat("  calculate mu1\n")
		mu1=halfnorm_mu1f(alpha,v1hat,fd1)
#
# 8 mu2
#
		if(debug)cat("  calculate mu2\n")
		mu2=halfnorm_mu2f(alpha,v1hat,fd1)
#
# 9 rhp
#
		lambdad_rhp=matrix(-1/v1hat,1)
#
# 10 fhat, dq and quantiles
#
		if(debug)cat("  fhat, dq and quantiles\n")
		fhat=dhalfnorm(ml_quantiles,theta=v1hat)
		dq=dmgs(lddi,lddd,mu1,lambdad_rhp,mu2,dim=1)
		rh_quantiles=ml_quantiles+dq/(nx*fhat)
#
# 11 means
#
		means=halfnorm_means(means,ml_params,lddi,lddd,lambdad_rhp,nx,dim=1)
		ml_mean=means$ml_mean
		rh_mean=means$rh_mean
#
# 12 waicscores
#
		waic=halfnorm_waic(waicscores,x,v1hat,fd1,lddi,lddd,lambdad_rhp,aderivs)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 13 logscores
#
		logscores=halfnorm_logscores(logscores,x,fd1,aderivs)
		ml_oos_logscore=logscores$ml_oos_logscore
		rh_oos_logscore=logscores$rh_oos_logscore
#
# 14 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rhalfnorm_cp(nrust,x,rust=TRUE,mlcp=FALSE)
			ru_quantiles=makeq(rustsim$ru_deviates,p)
		}
	} #end of if(dmgs)

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
#' @rdname halfnorm_cp
#' @inheritParams man
#' @export
rhalfnorm_cp=function(n,x,fd1=0.01,rust=FALSE,mlcp=TRUE,debug=FALSE,aderivs=TRUE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),!x<0)
	stopifnot(is.finite(x),!is.na(x),!x<0)

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qhalfnorm_cp(x,runif(n),aderivs=aderivs)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=thalfnorm_cp(n,x)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			ru_deviates[i]=rhalfnorm(1,theta=th[i])
		}
	}

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
			cp_method=rhp_dmgs_cpmethod())

	return(op)

}
#' @rdname halfnorm_cp
#' @inheritParams man
#' @export
dhalfnorm_cp=function(x,y=x,fd1=0.01,rust=FALSE,nrust=1000,debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),!x<0,!y<0)

	dd=dhalfnormsub(x=x,y=y,fd1=fd1,aderivs=aderivs)
	ru_pdf="rust not selected"
	if(rust){
		th=thalfnorm_cp(nrust,x)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_pdf=ru_pdf+dhalfnorm(y,theta=th[ir])
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
#' @rdname halfnorm_cp
#' @inheritParams man
#' @export
phalfnorm_cp=function(x,y=x,fd1=0.01,rust=FALSE,nrust=1000,debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),!x<0,!y<0)

	dd=dhalfnormsub(x=x,y=y,fd1=fd1,aderivs)
	ru_cdf="rust not selected"
	if(rust){
		th=thalfnorm_cp(nrust,x)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_cdf=ru_cdf+phalfnorm(y,theta=th[ir])
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
#' @rdname halfnorm_cp
#' @inheritParams man
#' @export
thalfnorm_cp=function(n,x,fd1=0.01,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),!x<0)
	stopifnot(is.finite(x),!is.na(x),!x<0)

	t=ru(halfnorm_logf,x=x,n=n,d=1,init=1)

	list(theta_samples=t$sim_vals)

}
