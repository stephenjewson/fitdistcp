# extraDistr is very confusing:
# a = shape parameter (they call scale), that I'm varying here
# b = scale parameter (they call location), that I'm calling kscale
#
#' Pareto Distribution Predictions Based on a Calibrating Prior
#'
#' @inherit man description author references seealso return
#' @inheritParams man
#'
#' @inheritSection man Optional Return Values
# #' @inheritSection man Optional Return Values (EVD models only)
# #' @inheritSection man Optional Return Values (non-RHP models only)
#' @inheritSection man Details (homogeneous models)
# #' @inheritSection man Details (non-homogeneous models)
#' @inheritSection man Details (analytic integration)
# #' @inheritSection man Details (DMGS integration)
#' @inheritSection man Details (RUST)
#'
#' @section Details of the Model:
#' The Pareto distribution has various forms.
#' The form we are using has exceedance distribution function
#' \deqn{S(x;\alpha)={\left(\frac{\sigma}{x}\right)^\alpha}}
#' where
#' \eqn{x \ge \sigma} is the random variable
#' and
#' \eqn{\alpha>0, \sigma>0} are the shape and scale parameters.
#' We consider the scale parameter \eqn{\sigma} to be known
#' (hence the \code{k2} in the name).
#'
#' The calibrating prior is given by the right Haar prior, which is
#' \deqn{\pi(\alpha) \propto \frac{1}{\alpha}}
#' as given in Jewson et al. (2025).
#' Some others authors may refer to the shape and scale parameters
#' as the scale and location parameters, respectively.
#'
#' @example man/examples/example_11_pareto_k2.R
#'
#' @name pareto_k2_cp
NULL
#' @rdname pareto_k2_cp
#' @inheritParams man
#' @export
#'
qpareto_k2_cp=function(x,p=seq(0.1,0.9,0.1),kscale=1,fd1=0.01,
	means=FALSE,waicscores=FALSE,logscores=FALSE,rust=FALSE,nrust=100000,
	debug=FALSE,aderivs=TRUE){
#
# 1 intro
#
	stopifnot(is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1,!x<kscale)
	alpha=1-p
	nx=length(x)
	lnx=log(x/kscale)
	slnx=sum(lnx)
	nalpha=length(alpha)
#
# 2 ml param estimate
#
	v1hat=pareto_k2_ml_params(x,kscale)
	ml_params=v1hat
#
# 3 aic
#
	ml_value=sum(log(v1hat)+v1hat*log(kscale)-(v1hat+1)*log(x))
	maic=make_maic(ml_value,nparams=2)
#
# 4 ml quantiles (vectorized over alpha)
#
	logalpha=log(alpha)
	ml_quantiles=exp(log(kscale)-logalpha/v1hat)
	ldd="only relevant for DMGS models, not analytic models"
	lddi="only relevant for DMGS models, not analytic models"
	expinfmat=nx/(v1hat*v1hat)
	expinfmati=1/expinfmat
	standard_errors=sqrt(expinfmati)
#
# 5 rhp quantiles (Vectorized over alpha)
#
	alphan=alpha**(1/nx)
	logq=(slnx/alphan)-slnx
	rh_quantiles=exp(log(kscale)+logq)
#
# 5.5 calc mpd quantiles (vectorized over alpha)
# -left here for now, but not returned
	logalphasq=logalpha*logalpha
	eps=0.5*alpha*logalphasq/nx
	ftop=v1hat
	fb1=(v1hat+1)/v1hat
	fbot=exp(-fb1*logalpha)
	ff=ftop/fbot
	delta=eps/ff
	mpd_quantiles=ml_quantiles+delta
#
# 6 means (might as well always calculate)
#
	if(v1hat>1){
		ml_mean=v1hat*kscale/(v1hat-1)
	} else{
		ml_mean=Inf
	}
	rh_mean=Inf
#
# 7 waicscores
#
	waic=pareto_k2_waic(waicscores,x,v1hat,fd1,kscale,aderivs)
	waic1=waic$waic1
	waic2=waic$waic2
#
# 8 logscores
#
	logscores=pareto_k2_logscores(logscores,x,kscale)
	ml_oos_logscore=logscores$ml_oos_logscore
	rh_oos_logscore=logscores$rh_oos_logscore
#
# 9 rust
#
	ru_quantiles="rust not selected"
	if(rust){
		rustsim=rpareto_k2_cp(nrust,x,kscale,rust=TRUE,mlcp=FALSE)
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
#' @rdname pareto_k2_cp
#' @inheritParams man
#' @export
rpareto_k2_cp=function(n,x,kscale=1,rust=FALSE,mlcp=TRUE,debug=FALSE,aderivs=TRUE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),!x<kscale)
	stopifnot(is.finite(x),!is.na(x),!x<kscale)

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qpareto_k2_cp(x,runif(n),kscale=kscale,aderivs=aderivs)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tpareto_k2_cp(n,x)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			ru_deviates[i]=rpareto(1,a=th[i],b=kscale)
		}
	}
	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
	 		 cp_method=analytic_cpmethod())

	return(op)

}
#' @rdname pareto_k2_cp
#' @inheritParams man
#' @export
dpareto_k2_cp=function(x,y=x,kscale=1,rust=FALSE,nrust=1000,debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),!x<kscale,!y<kscale)

	dd=dpareto_k2_sub(x=x,y=y,kscale=kscale,aderivs=aderivs)
	ru_pdf="rust not selected"
	if(rust){
		th=tpareto_k2_cp(nrust,x,kscale)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_pdf=ru_pdf+dpareto(y,a=th[ir],b=kscale)
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
#' @rdname pareto_k2_cp
#' @inheritParams man
#' @export
ppareto_k2_cp=function(x,y=x,kscale=1,rust=FALSE,nrust=1000,debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),x>=kscale)
#	stopifnot(is.finite(y),!is.na(y),y>=kscale)

	dd=dpareto_k2_sub(x=x,y=y,kscale=kscale,aderivs=aderivs)
	ru_cdf="rust not selected"
	if(rust){
		th=tpareto_k2_cp(nrust,x,kscale)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_cdf=ru_cdf+ppareto(y,a=th[ir],b=kscale)
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
#' @rdname pareto_k2_cp
#' @inheritParams man
#' @export
tpareto_k2_cp=function(n,x,kscale=1,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),!x<kscale)
	stopifnot(is.finite(x),!is.na(x),!x<kscale)

	t=ru(pareto_k2_logf,x=x,kscale=kscale,n=n,d=1,init=1)

	list(theta_samples=t$sim_vals)

}
