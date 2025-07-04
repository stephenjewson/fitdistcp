#' Log-normal Distribution Predictions Based on a Calibrating Prior
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
#' The log normal distribution has probability density function
#' \deqn{f(x;\mu,\sigma)=\frac{1}{\sqrt{2\pi}\sigma x}e^{-(\log(x)-\mu)^2/(2\sigma^2)}}
#' where
#' \eqn{x} is the random variable
#' and
#' \eqn{\mu,\sigma>0} are the parameters.
#'
#' The calibrating prior is given by the right Haar prior, which is
#' \deqn{\pi(\sigma) \propto \frac{1}{\sigma}}
#' as given in Jewson et al. (2025).
#'
#' @example man/examples/example_035_lnorm.R
#'
#' @name lnorm_cp
NULL
#' @rdname lnorm_cp
#' @inheritParams man
#' @export
#'
qlnorm_cp=function(x,p=seq(0.1,0.9,0.1),d1=0.01,fd2=0.01,
	means=FALSE,waicscores=FALSE,logscores=FALSE,rust=FALSE,nrust=100000,
	debug=FALSE,aderivs=TRUE){
#
# 1 intro
#
	stopifnot(is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1,!x<0)
	alpha=1-p
	y=log(x)
	nx=length(x)
	nalpha=length(alpha)
#
# 2 ml param estimate
#
	ml_params=norm_ml_params(y) #note that it uses y, and the normal routine
	v1hat=ml_params[1]
	v2hat=ml_params[2]
#
# 3 aic
#
	ml_value=sum(dlnorm(x,meanlog=v1hat,sdlog=v2hat,log=TRUE))
	maic=make_maic(ml_value,nparams=2)
#
# 4 ml quantiles (vectorized over alpha)
#
	ml_quantiles=qlnorm((1-alpha),meanlog=v1hat,sdlog=v2hat)

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
# 5 rhp quantiles
#
	mu=mean(y)
# calculate the unbiased variance
	s1=sqrt(var(y))
	temp=qt((1-alpha),df=nx-1)
# convert the unbiased to predictive
	rh_quantiles=exp(mu+temp*s1*sqrt((1+1/nx)))
#
# 6 means (might as well always calculate)
#
	ml_mean=exp(v1hat+0.5*v2hat*v2hat)
	rh_mean="no analytic expression"
#
# 7 waicscores
#
	waic=lnorm_waic(waicscores,x,v1hat,d1,v2hat,fd2,aderivs)
	waic1=waic$waic1
	waic2=waic$waic2
#
# 8 logscores
#
	logscores=lnorm_logscores(logscores,x)
	ml_oos_logscore=logscores$ml_oos_logscore
	rh_oos_logscore=logscores$rh_oos_logscore
#
# 9 rust
#
	ru_quantiles="rust not selected"
	if(rust){
		rustsim=rlnorm_cp(nrust,x,rust=TRUE,mlcp=FALSE)
		ru_quantiles=makeq(rustsim$ru_deviates,p)
	}


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
				cp_method=analytic_cpmethod())

}
#' @rdname lnorm_cp
#' @inheritParams man
#' @export
rlnorm_cp=function(n,x,rust=FALSE,mlcp=TRUE,debug=FALSE,aderivs=TRUE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),!x<0)
	stopifnot(is.finite(x),!is.na(x),!x<0)


	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qlnorm_cp(x,runif(n),aderivs=aderivs)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tlnorm_cp(n,x)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			ru_deviates[i]=rlnorm(1,meanlog=th[i,1],sdlog=th[i,2])
		}
	}
	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
				cp_method=analytic_cpmethod())

	return(op)


}
#' @rdname lnorm_cp
#' @inheritParams man
#' @export
dlnorm_cp=function(x,y=x,rust=FALSE,nrust=1000,debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),!x<0,!y<0)

	dd=dlnormsub(x=x,y=y,aderivs=aderivs)
	ru_pdf="rust not selected"
	if(rust){
		th=tlnorm_cp(nrust,x)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_pdf=ru_pdf+dlnorm(y,meanlog=th[ir,1],sdlog=th[ir,2])
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
#' @rdname lnorm_cp
#' @inheritParams man
#' @export
plnorm_cp=function(x,y=x,rust=FALSE,nrust=1000,debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),!x<0,!y<0)

	dd=dlnormsub(x=x,y=y,aderivs=aderivs)
	ru_cdf="rust not selected"
	if(rust){
		th=tlnorm_cp(nrust,x)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_cdf=ru_cdf+plnorm(y,meanlog=th[ir,1],sdlog=th[ir,2])
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
#' @rdname lnorm_cp
#' @inheritParams man
#' @export
tlnorm_cp=function(n,x,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),!x<0)
	stopifnot(is.finite(x),!is.na(x),!x<0)

	t=ru(lnorm_logf,x=x,n=n,d=2,init=c(mean(exp(x)),sd(exp(x))))

	list(theta_samples=t$sim_vals)
}

