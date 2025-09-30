#' Normal Distribution Predictions Based on a Calibrating Prior
#'
#' @inherit man description author references seealso return
#' @inheritParams man
#'
# #' @inheritSection man Optional Return Values (EVD models only)
# #' @inheritSection man Optional Return Values (non-RHP models only)
#' @inheritSection man Details (homogeneous models)
# #' @inheritSection man Details (non-homogeneous models)
#' @inheritSection man Details (analytic integration)
# #' @inheritSection man Details (DMGS integration)
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
#' @example man/examples/example_030_norm.R
#'
#' @name norm_cp
NULL
#' @rdname norm_cp
#' @inheritParams man
#' @export
#'
qnorm_cp=function(x,p=seq(0.1,0.9,0.1),
	means=FALSE,waicscores=FALSE,logscores=FALSE,rust=FALSE,nrust=100000,
	unbiasedv=FALSE,debug=FALSE){
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
	ml_params=norm_ml_params(x)
	v1hat=ml_params[1]
	v2hat=ml_params[2]
	uv_params=rep("unbiasedv not selected",2)
	if(unbiasedv)uv_params=norm_unbiasedv_params(x)
#
# 3 aic
#
	ml_value=sum(dnorm(x,mean=v1hat,sd=v2hat,log=TRUE))
	maic=make_maic(ml_value,nparams=2)

#
# 4 ml quantiles (vectorized over alpha)
#
	ml_quantiles=qnorm((1-alpha),mean=v1hat,sd=v2hat)
	uv_quantiles="unbiasedv not selected"
	if(unbiasedv)uv_quantiles=qnorm((1-alpha),mean=uv_params[1],sd=uv_params[2])
#
# 5 rhp quantiles (vectorized over alpha)
#
	mu=v1hat

# first, convert sigma from maxlik to unbiased
	sgu=v2hat*sqrt(nx/(nx-1))
# then, convert sigma to predictive sigma
	sg=sgu*sqrt((nx+1)/nx)

	temp=qt((1-alpha),df=nx-1)
	rh_quantiles=mu+temp*sg

	ldd="only relevant for DMGS models, not analytic models"
	lddi="only relevant for DMGS models, not analytic models"
	expinfmat=matrix(0,2,2)
	expinfmat[1,1]=nx/(v2hat*v2hat)
	expinfmat[2,2]=2*nx/(v2hat*v2hat)
	expinfmati=solve(expinfmat)
	standard_errors=matrix(0,2)
	standard_errors[1]=sqrt(expinfmati[1,1])
	standard_errors[2]=sqrt(expinfmati[2,2])

#
# test of gg code (for future implementation of mpd theory, as a test of the mpd code)
#
#	norm_gg(nx,v1hat,v2hat)
#

# 6 means (might as well always calculate)
#
	ml_mean=v1hat
	rh_mean=v1hat
#
# 7 waicscores
#
	waic=norm_waic(waicscores,x,v1hat,v2hat)
	waic1=waic$waic1
	waic2=waic$waic2
#
# 8 logscores
#
	logscores=norm_logscores(logscores,x)
	ml_oos_logscore=logscores$ml_oos_logscore
	rh_oos_logscore=logscores$rh_oos_logscore
#
# 9 rust
#
	ru_quantiles="rust not selected"
	if(rust){
		rustsim=rnorm_cp(nrust,x,method="rust",rust=TRUE,mlcp=FALSE)
		ru_quantiles=makeq(rustsim$ru_deviates,p)
	}

	list(	ml_params=ml_params,
				ml_value=ml_value,
				uv_params=uv_params,
#				ldd=ldd,
#				lddi=lddi,
#				expinfmat=expinfmat,
#				expinfmati=expinfmati,
				standard_errors=standard_errors,
				ml_quantiles=ml_quantiles,
				cp_quantiles=rh_quantiles,
				ru_quantiles=ru_quantiles,
				uv_quantiles=uv_quantiles,
				maic=maic,
				waic1=waic1,
				waic2=waic2,
				ml_oos_logscore=ml_oos_logscore,
				cp_oos_logscore=rh_oos_logscore,
				ml_mean=ml_mean,
				cp_mean=rh_mean,
				cp_method=analytic_cpmethod())

}
#' @rdname norm_cp
#' @inheritParams man
#' @export
rnorm_cp=function(n,x,method="rust",rust=FALSE,mlcp=TRUE,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x))
	stopifnot(is.finite(x),!is.na(x))

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qnorm_cp(x,runif(n))
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tnorm_cp(method=method,n,x)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			ru_deviates[i]=rnorm(1,mean=th[i,1],sd=th[i,2])
		}
	}

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
				cp_method=analytic_cpmethod())
	return(op)

}
#' @rdname norm_cp
#' @inheritParams man
#' @export
dnorm_cp=function(x,y=x,
	rust=FALSE,nrust=1000,
	boot=FALSE,nboot=1000,
	debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y))

	dd=dnormsub(x=x,y=y)
	ru_pdf="rust not selected"
	bs_pdf="boot not selected"

	if(rust){
		th=tnorm_cp("rust",nrust,x)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_pdf=ru_pdf+dnorm(y,mean=th[ir,1],sd=th[ir,2])
		}
		ru_pdf=ru_pdf/nrust
	}

	if(boot){
		th=tnorm_cp("boot",nboot,x)$theta_samples
		bs_pdf=numeric(length(y))
		for (ir in 1:nboot){
			bs_pdf=bs_pdf+dnorm(y,mean=th[ir,1],sd=th[ir,2])
		}
		bs_pdf=bs_pdf/nboot
	}

	op=list(	ml_params=dd$ml_params,
					ml_pdf=dd$ml_pdf,
					cp_pdf=dd$rh_pdf,
					ru_pdf=ru_pdf,
					bs_pdf=bs_pdf,
					cp_method=analytic_cpmethod())
	return(op)
}
#' @rdname norm_cp
#' @inheritParams man
#' @export
pnorm_cp=function(x,y=x,
	rust=FALSE,nrust=1000,
	boot=FALSE,nboot=1000,
	debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y))

	dd=dnormsub(x=x,y=y)
	ru_cdf="rust not selected"
	bs_cdf="rust not selected"

	if(rust){
		th=tnorm_cp("rust",nrust,x)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_cdf=ru_cdf+pnorm(y,mean=th[ir,1],sd=th[ir,2])
		}
		ru_cdf=ru_cdf/nrust
	}

	if(boot){
		th=tnorm_cp("boot",nboot,x)$theta_samples
		bs_cdf=numeric(length(y))
		for (ir in 1:nboot){
			bs_cdf=bs_cdf+pnorm(y,mean=th[ir,1],sd=th[ir,2])
		}
		bs_cdf=bs_cdf/nboot
	}

		op=list(	ml_params=dd$ml_params,
					ml_cdf=dd$ml_cdf,
					cp_cdf=dd$rh_cdf,
					ru_cdf=ru_cdf,
					bs_cdf=bs_cdf,
					cp_method=analytic_cpmethod())
	return(op)
}
#' @rdname norm_cp
#' @inheritParams man
#' @export
tnorm_cp=function(method,n,x,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x))
	stopifnot(is.finite(x),!is.na(x))

	if(method=="rust"){
		th=ru(norm_logf,x=x,n=n,d=2,init=c(mean(x),sd(x)))
	} else if (method=="boot"){
		th=norm_boot(x=x,n=n)
	} else{
		message("tnorm method not valid so stopping.\n")
		stop()
	}

	list(theta_samples=th$sim_vals)

}


