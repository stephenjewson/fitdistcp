#' t Distribution Predictions Based on a Calibrating Prior
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
#' The t distribution
#' (also known as the location-scale t distribution, hence the name \code{lst}),
#' has probability density function
#' \deqn{f(x;\mu,\sigma)
#' =\frac{\Gamma((\nu+1)/2)}{\sqrt{\pi\nu}\sigma\Gamma(\nu/2)}
#' \left(1+\frac{(x-\mu)^2}{\sigma^2\nu}\right)^{(\nu+1)/2}}
#' where
#' \eqn{x} is the random variable,
#' \eqn{\mu,\sigma>0} are the parameters,
#' and we consider the degrees of freedom \eqn{\nu} to be known
#' (hence the \code{k3} in the name).
#'
#' The calibrating prior is given by the right Haar prior, which is
#' \deqn{\pi(\sigma) \propto \frac{1}{\sigma}}
#' as given in Jewson et al. (2025).
#'
#' @example man/examples/example_41_lst_k3.R
#'
#' @name lst_k3_cp
NULL
#' @rdname lst_k3_cp
#' @inheritParams man
#' @export
#'
qlst_k3_cp=function(x,p=seq(0.1,0.9,0.1),kdf=5,d1=0.01,fd2=0.01,
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
						!is.na(kdf))
	alpha=1-p
	nx=length(x)
	nalpha=length(alpha)
#
# 2 ml param estimate
#
	if(debug)message("2 calc ml param estimate")
	v1start=mean(x)
	v2start=sd(x)
	opt=optim(c(v1start,v2start),lst_k3_loglik,x=x,kdf=kdf,control=list(fnscale=-1))
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
	ml_quantiles=qlst((1-alpha),mu=v1hat,sigma=v2hat,df=kdf)
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
		if(debug)message("  calculate ldd,lddi")
		if(aderivs)	ldd=lst_k3_ldda(x,v1hat,v2hat,kdf)
		if(!aderivs)ldd=lst_k3_ldd(x,v1hat,d1,v2hat,fd2,kdf)
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
#
# 6 lddd
#
		if(debug)message("  calculate lddd")
		if(aderivs)	lddd=lst_k3_lddda(x,v1hat,v2hat,kdf)
		if(!aderivs)lddd=lst_k3_lddd(x,v1hat,d1,v2hat,fd2,kdf)
#
# 7 mu1
#
		if(debug)message("  calculate mu1")
		mu1=lst_k3_mu1f(alpha,v1hat,d1,v2hat,fd2,kdf)
#
# 8 mu2
#
		if(debug)message("  calculate mu2")
		mu2=lst_k3_mu2f(alpha,v1hat,d1,v2hat,fd2,kdf)
#
# 9 rhp
#
		lambdad_rhp=c(0,-1/v2hat)
#
# 10 fhat, dq and quantiles
#
		if(debug)message("  fhat, dq and quantiles")
		fhat=dlst(ml_quantiles,mu=v1hat,sigma=v2hat,df=kdf)
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
		waic=lst_k3_waic(waicscores,x,v1hat,d1,v2hat,fd2,kdf,lddi,lddd,lambdad_rhp,aderivs)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 13 logscores
#
		logscores=lst_k3_logscores(logscores,x,d1,fd2,kdf,aderivs)
		ml_oos_logscore=logscores$ml_oos_logscore
		rh_oos_logscore=logscores$rh_oos_logscore
#
# 14 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rlst_k3_cp(nrust,x,rust=TRUE,mlcp=FALSE)
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
#' @rdname lst_k3_cp
#' @inheritParams man
#' @export
rlst_k3_cp=function(n,x,d1=0.01,fd2=0.01,kdf=5,rust=FALSE,mlcp=TRUE,
	debug=FALSE,aderivs=TRUE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x))
	stopifnot(is.finite(x),!is.na(x))

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qlst_k3_cp(x,runif(n),d1=d1,fd2=fd2,kdf=kdf,aderivs=aderivs)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tlst_k3_cp(n,x,kdf=kdf)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			ru_deviates[i]=rlst(1,mu=th[i,1],sigma=th[i,2],df=kdf)
		}
	}

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
				 ru_deviates=ru_deviates,
		 		 cp_method=rhp_dmgs_cpmethod())
	return(op)


}
#' @rdname lst_k3_cp
#' @inheritParams man
#' @export
dlst_k3_cp=function(x,y=x,d1=0.01,fd2=0.01,kdf=5,rust=FALSE,nrust=1000,
	debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y))

	dd=dlst_k3sub(x=x,y=y,d1,fd2,kdf,aderivs=aderivs)
	ru_pdf="rust not selected"
	if(rust){
		th=tlst_k3_cp(nrust,x,kdf)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_pdf=ru_pdf+dlst(y,mu=th[ir,1],sigma=th[ir,2],df=kdf)
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
#' @rdname lst_k3_cp
#' @inheritParams man
#' @export
plst_k3_cp=function(x,y=x,d1=0.01,fd2=0.01,kdf=5,rust=FALSE,nrust=1000,
	debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y))

	dd=dlst_k3sub(x=x,y=y,d1,fd2,kdf,aderivs=aderivs)
	ru_cdf="rust not selected"
	if(rust){
		th=tlst_k3_cp(nrust,x,kdf)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_cdf=ru_cdf+plst(y,mu=th[ir,1],sigma=th[ir,2],df=kdf)
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
#' @rdname lst_k3_cp
#' @inheritParams man
#' @export
tlst_k3_cp=function(n,x,d1=0.01,fd2=0.01,kdf=5,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x))
	stopifnot(is.finite(x),!is.na(x))

	t=ru(lst_k3_logf,x=x,kdf=kdf,n=n,d=2,init=c(mean(x),sd(x)))

	list(theta_samples=t$sim_vals)

}
