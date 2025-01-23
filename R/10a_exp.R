#' Exponential Distribution Predictions Based on a Calibrating Prior
#'
#' @inherit man description author references seealso
#' @inheritParams man
#'
#' @inheritSection man Default Return Values
#' @inheritSection man Optional Return Values
# #' @inheritSection man Optional Return Values (EVT models only)
# #' @inheritSection man Optional Return Values (non-homogeneous models only)
#' @inheritSection man Details (homogeneous models)
# #' @inheritSection man Details (non-homogeneous models)
#' @inheritSection man Details (analytic integration)
# #' @inheritSection man Details (DMGS integration)
#' @inheritSection man Details (RUST)
#'
#' @section Details of the Model:
#' The exponential distribution has exceedance distribution function
#' \deqn{S(x;\lambda)=\exp(-\lambda x)}
#' where
#' \eqn{x \ge 0} is the random variable
#' and
#' \eqn{\lambda >0} is the rate parameter.
#'
#' The calibrating prior is given by the right Haar prior, which is
#' \deqn{\pi(\lambda) \propto \frac{1}{\lambda}}
#' as given in Jewson et al. (2024).
#'
#' @example man/examples/example_10_exp.R
#'
#' @name exp_cp
NULL
#' @rdname exp_cp
#' @inheritParams man
#' @export
qexp_cp=function(x,p=seq(0.1,0.9,0.1),fd1=0.01,
	means=FALSE,waicscores=FALSE,logscores=FALSE,rust=FALSE,nrust=100000,
	debug=FALSE,aderivs=TRUE){
#
# 1 intro
#
	stopifnot(is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1,!x<0)
	alpha=1-p
	nx=length(x)
	sumx=sum(x)
#
# 2 ml param estimate
#
	v1hat=nx/sumx
	ml_params=v1hat
#
# 3 aic
#
	ml_value=sum(dexp(x,rate=v1hat,log=TRUE))
	maic=make_maic(ml_value,nparams=1)
#
# 4 ml quantiles (vectorized over alpha)
#
	ml_quantiles=qexp((1-alpha),rate=v1hat)
	ldd="only relevant for DMGS models, not analytic models"
	lddi="only relevant for DMGS models, not analytic models"
	expinfmat=nx/(v1hat*v1hat)
	expinfmati=1/expinfmat
	standard_errors=sqrt(expinfmati)
#	if(extras){
#		expinfmat=nx/(v1hat*v1hat)
#		expinfmati=1/expinfmat
#		standard_errors=sqrt(expinfmati)
#	}else{
#		expinfmat="extras not selected"
#		expinfmati="extras not selected"
#		standard_errors="extras not selected"
#	}

#
# 5 rh_quantiles (vectorized over alpha)
#
	rh_quantiles=sumx*((alpha^(-1/nx))-1)
#
# 6 means (might as well always calculate)
#
	ml_mean=1/v1hat
	rh_mean=sumx/(nx-1)
#
# 7 waicscores
#
	waic=exp_waic(waicscores,x,v1hat,fd1,aderivs)
	waic1=waic$waic1
	waic2=waic$waic2
#
# 8 logscores
#
	logscores=exp_logscores(logscores,x)
	ml_oos_logscore=logscores$ml_oos_logscore
	rh_oos_logscore=logscores$rh_oos_logscore
#
# 9 rust
#
	ru_quantiles="rust not selected"
	if(rust){
		rustsim=rexp_cp(nrust,x,rust=TRUE,mlcp=FALSE)
		ru_quantiles=makeq(rustsim$ru_deviates,p)
	}

# return
	list(	ml_params=ml_params,
				ml_value=ml_value,
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
				cp_method=analytic_cpmethod())

}
#' @rdname exp_cp
#' @inheritParams man
#' @export
rexp_cp=function(n,x,rust=FALSE,mlcp=TRUE,debug=FALSE,aderivs=TRUE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),!x<0)
	stopifnot(is.finite(x),!is.na(x),!x<0)

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qexp_cp(x,runif(n),aderivs=aderivs)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=texp_cp(n,x)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			ru_deviates[i]=rexp(1,rate=th[i])
		}
	}

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
	 		 cp_method=analytic_cpmethod())

	return(op)
}
#' @rdname exp_cp
#' @inheritParams man
#' @export
dexp_cp=function(x,y=x,rust=FALSE,nrust=1000,debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),!x<0,!y<0)

	dd=dexpsub(x=x,y=y,aderivs=aderivs)
	ru_pdf="rust not selected"
	if(rust){
		th=texp_cp(nrust,x)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_pdf=ru_pdf+dexp(y,rate=th[ir])
		}
		ru_pdf=ru_pdf/nrust
	}
	op=list(	ml_params=dd$ml_params,
				ml_pdf=dd$ml_pdf,
				cp_pdf=dd$rh_pdf,
				ru_pdf=ru_pdf,
				cp_method=analytic_cpmethod())
	return(op)

}
#' @rdname exp_cp
#' @inheritParams man
#' @export
pexp_cp=function(x,y=x,rust=FALSE,nrust=1000,debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),!x<0,!y<0)

	dd=dexpsub(x=x,y=y,aderivs=aderivs)
	ru_cdf="rust not selected"
	if(rust){
		th=texp_cp(nrust,x)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_cdf=ru_cdf+pexp(y,rate=th[ir])
		}
		ru_cdf=ru_cdf/nrust
	}
	op=list(	ml_params=dd$ml_params,
				ml_cdf=dd$ml_cdf,
				cp_cdf=dd$rh_cdf,
				ru_cdf=ru_cdf,
				cp_method=analytic_cpmethod())
	return(op)
}
#' @rdname exp_cp
#' @inheritParams man
#' @export
texp_cp=function(n,x,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),!x<0)
	stopifnot(is.finite(x),!is.na(x),!x<0)

	t=ru(exp_logf,x=x,n=n,d=1,init=mean(x))

	list(theta_samples=t$sim_vals)

}
