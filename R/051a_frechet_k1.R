#' Frechet Distribution Predictions Based on a Calibrating Prior
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
#' The Frechet distribution has distribution function
#' \deqn{F(x;\sigma,\lambda)=\exp\left(-\left(\frac{x-\mu}{\sigma}\right)^{-\lambda}\right)}
#' where
#' \eqn{x>\mu} is the random variable,
#' \eqn{\sigma>0,\lambda>0} are the parameters
#' and we consider \eqn{\mu} to be known (hence the \code{k1} in the name).
#'
#' The calibrating prior is given by the right Haar prior, which is
#' \deqn{\pi(\sigma,\lambda) \propto \frac{1}{\sigma \lambda}}
#' as given in Jewson et al. (2025).
#'
#' @example man/examples/example_051_frechet_k1.R
#'
#' @name frechet_k1_cp
NULL
#' @rdname frechet_k1_cp
#' @inheritParams man
#' @export
#'
qfrechet_k1_cp=function(x,p=seq(0.1,0.9,0.1),kloc=0,
	means=FALSE,waicscores=FALSE,logscores=FALSE,dmgs=TRUE,rust=FALSE,nrust=100000,
	debug=FALSE){
#
# 1 intro
#
	debug=TRUE
	debug=FALSE
	stopifnot(is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1,!x<kloc)
	alpha=1-p
	nx=length(x)
	nalpha=length(alpha)
#
# 2 calc ml param estimate
#
	opt=optim(c(1,1),frechet_loglik,x=x,kloc=kloc,control=list(fnscale=-1))
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
	ml_quantiles=qfrechet((1-alpha),mu=kloc,sigma=v1hat,lambda=v2hat)
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
		ldd=frechet_k1_ldda(x,v1hat,v2hat,kloc)
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
#
# 6 lddd
#
	  if(debug)message("  calculate lddd")
		lddd=frechet_k1_lddda(x,v1hat,v2hat,kloc)
#
# 7 mu1
#
		if(debug)message("calculate mu1")
		mu1=frechet_k1_mu1fa(alpha,v1hat,v2hat,kloc)
#
# 8 mu2
#
		if(debug)message("calculate mu2")
		mu2=frechet_k1_mu2fa(alpha,v1hat,v2hat,kloc)
#
# 9 rhp
#
		lambdad_rhp=c(-1/v1hat,-1/v2hat)
#
# 10 bayesian dq
#
		if(debug)message("  fhat, dq and rhp quantiles")
		fhat=dfrechet(ml_quantiles,mu=kloc,sigma=v1hat,lambda=v2hat)
		dq=dmgs(lddi,lddd,mu1,lambdad_rhp,mu2,dim=2)
		rh_quantiles=ml_quantiles+dq/(nx*fhat)
#
# 11 means
#
		means=frechet_means(means,ml_params,lddi,lddd,lambdad_rhp,nx,dim=2,kloc)
		ml_mean=means$ml_mean
		rh_mean=Inf
#
# 12 waicscores
#
		waic=frechet_k1_waic(waicscores,x,v1hat,v2hat,kloc,lddi,lddd,
			lambdad_rhp)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 13 logscores
#
		logscores=frechet_logscores(logscores,x,kloc)
		ml_oos_logscore=logscores$ml_oos_logscore
		rh_oos_logscore=logscores$rh_oos_logscore
#
# 14 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rfrechet_k1_cp(nrust,x,rust=TRUE,mlcp=FALSE)
			ru_quantiles=makeq(rustsim$ru_deviates,p)
		}
	}  #end of if(dmgs)

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
#' @rdname frechet_k1_cp
#' @inheritParams man
#' @export
rfrechet_k1_cp=function(n,x,kloc=0,rust=FALSE,mlcp=TRUE,
	debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),!x<kloc)
	stopifnot(is.finite(x),!is.na(x),!x<kloc)

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qfrechet_k1_cp(x,runif(n),kloc=kloc)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tfrechet_k1_cp(n,x)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			ru_deviates[i]=rfrechet(1,mu=kloc,sigma=th[i,1],lambda=th[i,2])
		}
	}

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
				cp_method=rhp_dmgs_cpmethod())

	return(op)

}
#' @rdname frechet_k1_cp
#' @inheritParams man
#' @export
dfrechet_k1_cp=function(x,y=x,kloc=0,rust=FALSE,nrust=1000,
	debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),!x<kloc,!y<kloc)

	dd=dfrechetsub(x=x,y=y,kloc=0)
	ru_pdf="rust not selected"
	if(rust){
		th=tfrechet_k1_cp(nrust,x,kloc)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_pdf=ru_pdf+dfrechet(y,mu=kloc,sigma=th[ir,1],lambda=th[ir,2])
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
#' @rdname frechet_k1_cp
#' @inheritParams man
#' @export
pfrechet_k1_cp=function(x,y=x,kloc=0,rust=FALSE,nrust=1000,
	debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),!x<kloc,!y<kloc)

	dd=dfrechetsub(x=x,y=y,kloc=0)
	ru_cdf="rust not selected"
	if(rust){
		th=tfrechet_k1_cp(nrust,x,kloc)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_cdf=ru_cdf+pfrechet(y,mu=kloc,sigma=th[ir,1],lambda=th[ir,2])
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
#' @rdname frechet_k1_cp
#' @inheritParams man
#' @export
tfrechet_k1_cp=function(n,x,kloc=0,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),!x<kloc)
	stopifnot(is.finite(x),!is.na(x),!x<kloc)

	t=ru(frechet_k1_logf,x=x,kloc=kloc,n=n,d=2,init=c(1,1))

	list(theta_samples=t$sim_vals)

}

