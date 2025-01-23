#' Generalized Normal Distribution Predictions Based on a Calibrating Prior
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
#' The generalized normal distribution has probability density function
#' \deqn{f(x;\mu,\alpha)=\frac{\beta}{2\alpha\Gamma(1/\beta)}e^{-(|x-\mu|/\alpha)^\beta}}
#' where
#' \eqn{x} is the random variable,
#' \eqn{\mu,\alpha>0} are the parameters and we consider
#' \eqn{\beta} to be known
#' (hence the \code{k3} in the name).
#'
#' The calibrating prior is given by the right Haar prior, which is
#' \deqn{\pi(\alpha) \propto \frac{1}{\alpha}}
#' as given in Jewson et al. (2024).
#'
#' @example man/examples/example_32_gnorm_k3.R
#'
#' @name gnorm_k3_cp
NULL
#' @rdname gnorm_k3_cp
#' @inheritParams man
#' @export
#'
qgnorm_k3_cp=function(x,p=seq(0.1,0.9,0.1),kbeta=4,d1=0.01,fd2=0.01,
	means=FALSE,waicscores=FALSE,logscores=FALSE,dmgs=TRUE,rust=FALSE,nrust=100000,
	debug=FALSE,aderivs=TRUE){
#
# 1 intro
#
	debug=TRUE
	debug=FALSE
	stopifnot(is.finite(x),
						!is.na(x),
						is.finite(p),
						!is.na(p),
						p>0,
						p<1,
						!is.na(kbeta))
	alpha=1-p
	nx=length(x)
	nalpha=length(alpha)
#
# 2 ml param estimate
#
	if(debug)cat("2 calc ml param estimate")
	v1start=mean(x)
	v2start=sd(x)
	opt=optim(c(v1start,v2start),gnorm_k3_loglik,x=x,kbeta=kbeta,control=list(fnscale=-1))
	v1hat=opt$par[1]
	v2hat=opt$par[2]
	ml_params=c(v1hat,v2hat)
	if(debug)cat("  v1hat,v2hat=",v1hat,v2hat,"//")
#
# 3 aic
#
	ml_value=opt$val
	maic=make_maic(ml_value,nparams=2)

#
# 4 ml quantiles (vectorized over alpha)
#
	ml_quantiles=qgnorm((1-alpha),mu=v1hat,alpha=v2hat,beta=kbeta)
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
		if(kbeta<=1){
			cat("dmgs cannot support gnorm shape parameter <=1,
				because the density is not differentiable\n")
			stop()
		}
#
# 5 lddi
#
		if(debug)cat("  calculate ldd,lddi\n")
		if(aderivs)	ldd=gnorm_k3_ldda(x,v1hat,v2hat,kbeta)
		if(!aderivs)ldd=gnorm_k3_ldd(x,v1hat,d1,v2hat,fd2,kbeta)
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
#
# 6 lddd
#
		if(debug)cat("  calculate lddd\n")
		if(aderivs)	lddd=gnorm_k3_lddda(x,v1hat,v2hat,kbeta)
		if(!aderivs)lddd=gnorm_k3_lddd(x,v1hat,d1,v2hat,fd2,kbeta)
#
# 7 mu1
#
		if(debug)cat("  calculate mu1\n")
		mu1=gnorm_k3_mu1f(alpha,v1hat,d1,v2hat,fd2,kbeta)
#
# 8 mu2
#
		if(debug)cat("  calculate mu2\n")
		mu2=gnorm_k3_mu2f(alpha,v1hat,d1,v2hat,fd2,kbeta)
#
# 9 rhp
#
		lambdad_rhp=c(0,-1/v2hat)
#
# 10 fhat, dq and quantiles
#
		if(debug)cat("  fhat, dq and quantiles\n")
		fhat=dgnorm(ml_quantiles,mu=v1hat,alpha=v2hat,beta=kbeta)
		dq=dmgs(lddi,lddd,mu1,lambdad_rhp,mu2,dim=2)
		rh_quantiles=ml_quantiles+dq/(nx*fhat)
#
# 11 means
#
		if(means){
			ml_mean=ml_params[1]
			rh_mean=ml_params[1]
		}else{
			ml_mean="means not selected"
			rh_mean="means not selected"
		}
#
# 12 waicscores
#
		waic=gnorm_waic(waicscores,x,v1hat,d1,v2hat,fd2,kbeta,lddi,lddd,lambdad_rhp,aderivs)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 13 logscores
#
		logscores=gnorm_k3_logscores(logscores,x,d1,fd2,kbeta,aderivs)
		ml_oos_logscore=logscores$ml_oos_logscore
		rh_oos_logscore=logscores$rh_oos_logscore
#
# 14 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rgnorm_k3_cp(nrust,x,rust=TRUE,mlcp=FALSE)
			ru_quantiles=makeq(rustsim$ru_deviates,p)
		}
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
#' @rdname gnorm_k3_cp
#' @inheritParams man
#' @export
rgnorm_k3_cp=function(n,x,d1=0.01,fd2=0.01,kbeta=4,rust=FALSE,mlcp=TRUE,
	debug=FALSE,aderivs=TRUE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x))
	stopifnot(is.finite(x),!is.na(x))

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qgnorm_k3_cp(x,runif(n),d1=d1,fd2=fd2,kbeta=kbeta,aderivs)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tgnorm_k3_cp(n,x,kbeta=kbeta)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			ru_deviates[i]=rgnorm(1,mu=th[i,1],alpha=th[i,2],beta=kbeta)
		}
	}

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
				 ru_deviates=ru_deviates,
		 		 cp_method=rhp_dmgs_cpmethod())
	return(op)


}
#' @rdname gnorm_k3_cp
#' @inheritParams man
#' @export
dgnorm_k3_cp=function(x,y=x,d1=0.01,fd2=0.01,kbeta=4,rust=FALSE,nrust=1000,
	debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y))

	dd=dgnorm_k3sub(x=x,y=y,d1,fd2,kbeta,aderivs)
	ru_pdf="rust not selected"
	if(rust){
		th=tgnorm_k3_cp(nrust,x,kbeta)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_pdf=ru_pdf+dgnorm(y,mu=th[ir,1],alpha=th[ir,2],beta=kbeta)
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
#' @rdname gnorm_k3_cp
#' @inheritParams man
#' @export
pgnorm_k3_cp=function(x,y=x,d1=0.01,fd2=0.01,kbeta=4,rust=FALSE,nrust=1000,
	debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y))

	dd=dgnorm_k3sub(x=x,y=y,d1,fd2,kbeta,aderivs)
	ru_cdf="rust not selected"
	if(rust){
		th=tgnorm_k3_cp(nrust,x,kbeta)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_cdf=ru_cdf+pgnorm(y,mu=th[ir,1],alpha=th[ir,2],beta=kbeta)
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
#' @rdname gnorm_k3_cp
#' @inheritParams man
#' @export
tgnorm_k3_cp=function(n,x,d1=0.01,fd2=0.01,kbeta=4,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x))
	stopifnot(is.finite(x),!is.na(x))

	t=ru(gnorm_k3_logf,x=x,kbeta=kbeta,n=n,d=2,init=c(mean(x),sd(x)))

	list(theta_samples=t$sim_vals)

}
