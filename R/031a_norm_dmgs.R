#' Normal Distribution Predictions Based on a Calibrating Prior, using DMGS (for testing only)
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
#' The normal distribution has probability density function
#' \deqn{f(x;\mu,\sigma)=\frac{1}{\sqrt{2\pi}\sigma}e^{-(x-\mu)^2/(2\sigma^2)}}
#' where
#' \eqn{x} is the random variable
#' and
#' \eqn{\mu,\sigma>0} are the parameters.
#'
#' The calibrating prior is given by the right Haar prior, which is
#' \deqn{\pi(\sigma) \propto \frac{1}{\sigma}}
#' as given in Jewson et al. (2025).
#'
#' @example man/examples/example_031_norm_dmgs.R
#'
#' @name norm_dmgs_cp
NULL
#' @rdname norm_dmgs_cp
#' @inheritParams man
#' @export
#'
qnorm_dmgs_cp=function(x,p=seq(0.1,0.9,0.1),
	means=FALSE,waicscores=FALSE,logscores=FALSE,dmgs=TRUE,debug=FALSE){
#
# 1 intro
#
	stopifnot(is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1)
	alpha=1-p
	nx=length(x)
	nalpha=length(alpha)
#
# 2 ml param estimate
#
	if(debug)message("2 calc ml param estimate")
	v1start=mean(x)
	v2start=sd(x)
	opt=optim(c(v1start,v2start),norm_dmgs_loglik,x=x,control=list(fnscale=-1))
	v1hat=opt$par[1]
	v2hat=opt$par[2]
	ml_params=c(v1hat,v2hat)
	if(debug)message("  v1hat,v2hat=",v1hat,v2hat,"//")
#
# 3 aic
#
	ml_value=opt$val
	maic=make_maic(ml_value,nparams=2)

#
# 4 ml quantiles (vectorized over alpha)
#
	ml_quantiles=qnorm((1-alpha),mean=v1hat,sd=v2hat)
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
		if(debug)message("  calculate ldd,lddi\n")
		ldd=norm_ldda(x,v1hat,v2hat)
		if(debug)message("ldd=",ldd,"\n")
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
#
# 6 lddd
#
		if(debug)message("  calculate lddd\n")
		lddd=norm_lddda(x,v1hat,v2hat)
#
# 7 mu1
#
		if(debug)message("  calculate mu1\n")
		mu1=norm_mu1fa(alpha,v1hat,v2hat)
#
# 8 mu2
#
		if(debug)message("  calculate mu2\n")
		mu2=norm_mu2fa(alpha,v1hat,v2hat)
#
# 9 rhp
#
		lambdad_rhp=c(0,-1/v2hat)
#
# 10 fhat, dq and quantiles
#
		if(debug)message("  fhat, dq and quantiles\n")
		fhat=dnorm(ml_quantiles,mean=v1hat,sd=v2hat)
		dq=dmgs(lddi,lddd,mu1,lambdad_rhp,mu2,dim=2)
		rh_quantiles=ml_quantiles+dq/(nx*fhat)
#
# 11 means
#
		means=norm_dmgs_means(means,ml_params,lddi,lddd,lambdad_rhp,nx,dim=2)
		ml_mean=means$ml_mean
		rh_mean=means$rh_mean
#
# 12 waicscores
#
		waic=norm_dmgs_waic(waicscores,x,v1hat,v2hat,lddi,lddd,
			lambdad_rhp)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 13 logscores
#
		logscores=norm_dmgs_logscores(logscores,x)
		ml_oos_logscore=logscores$ml_oos_logscore
		rh_oos_logscore=logscores$rh_oos_logscore
	} # end of if(dmgs)

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
				maic=maic,
				waic1=waic1,
				waic2=waic2,
				ml_oos_logscore=ml_oos_logscore,
				cp_oos_logscore=rh_oos_logscore,
				ml_mean=ml_mean,
				cp_mean=rh_mean,
				cp_method=rhp_dmgs_cpmethod())

}
#' @rdname norm_dmgs_cp
#' @inheritParams man
#' @export
rnorm_dmgs_cp=function(n,x,mlcp=TRUE,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x))
	stopifnot(is.finite(x),!is.na(x))

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"

	if(mlcp){
		q=qnorm_dmgs_cp(x,runif(n))
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
				cp_method=rhp_dmgs_cpmethod())

}
#' @rdname norm_dmgs_cp
#' @inheritParams man
#' @export
dnorm_dmgs_cp=function(x,y=x,debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y))

	dd=dnorm_dmgssub(x=x,y=y)

	list(	ml_params=dd$ml_params,
				ml_pdf=dd$ml_pdf,
				cp_pdf=dd$rh_pdf,
				cp_pdf1=dd$rh_pdf1,
				cp_pdf2=dd$rh_pdf2,
				cp_method=rhp_dmgs_cpmethod())
}
#' @rdname norm_dmgs_cp
#' @inheritParams man
#' @export
pnorm_dmgs_cp=function(x,y=x,debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y))

	dd=dnorm_dmgssub(x=x,y=y)

	list(	ml_params=dd$ml_params,
				ml_cdf=dd$ml_cdf,
				cp_cdf=dd$rh_cdf,
				cp_cdf1=dd$rh_cdf1,
				cp_cdf2=dd$rh_cdf2,
				cp_method=rhp_dmgs_cpmethod())
}
