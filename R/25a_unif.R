#' Uniform Distribution Predictions Based on a Calibrating Prior
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
#' @inheritSection man Details (analytic integration)
# #' @inheritSection man Details (DMGS integration)
#' @inheritSection man Details (RUST)
#'
#' @section Details of the Model:
#' The uniform distribution has probability density function
#' \deqn{f(x;min,max)=\frac{1}{max-min}}
#' and zero otherwise,
#' where
#' \eqn{min \le x \le max } is the random variable
#' and
#' \eqn{min, max} are the parameters.
#'
#' The calibrating prior is given by the right Haar prior, which is
#' \deqn{\pi(\lambda) \propto \frac{1}{max-min}}
#' as given in Jewson et al. (2024).
#'
#' @example man/examples/example_25_unif.R
#'
#' @name unif_cp
NULL
#' @rdname unif_cp
#' @inheritParams man
#' @export
#'
qunif_cp=function(x,p=seq(0.1,0.9,0.1),
	means=FALSE,debug=FALSE,aderivs=TRUE){
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
	v1hat=min(x)
	v2hat=max(x)
	ml_params=c(v1hat,v2hat)
#
# 3 aic
#
	ml_value=sum(dunif(x,min=v1hat,max=v2hat,log=TRUE))
	maic=make_maic(ml_value,nparams=2)
#
# 4 all quantile calculations
#
	qq=qunif_formula(x,1-alpha)
	ml_quantiles=qq$ml_quantiles
	rh_quantiles=qq$rh_quantiles

	ldd="only relevant for DMGS models, not analytic models"
	lddi="only relevant for DMGS models, not analytic models"
#
# 6 means (might as well always calculate)
#
	ml_mean=qq$mean
	rh_mean=qq$mean

	list(	ml_params=ml_params,
				ml_value=ml_value,
				ml_quantiles=ml_quantiles,
				cp_quantiles=rh_quantiles,
				maic=maic,
				ml_mean=ml_mean,
				cp_mean=rh_mean,
				cp_method=analytic_cpmethod())

}
#' @rdname unif_cp
#' @inheritParams man
#' @export
runif_cp=function(n,x,mlcp=TRUE,debug=FALSE,aderivs=TRUE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x))
	stopifnot(is.finite(x),!is.na(x))

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"

	if(mlcp){
		q=qunif_cp(x,runif(n),aderivs=aderivs)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
				cp_method=analytic_cpmethod())
	return(op)

}
#' @rdname unif_cp
#' @inheritParams man
#' @export
dunif_cp=function(x,y=x,debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y))

	dd=dunif_formula(x,y)
#	cat("dd$ml_pdf=",dd$ml_pdf,"\n")

	op=list(	ml_params=dd$ml_params,
					ml_pdf=dd$ml_pdf,
					cp_pdf=dd$rh_pdf,
					cp_method=analytic_cpmethod())
	return(op)
}
#' @rdname unif_cp
#' @inheritParams man
#' @export
punif_cp=function(x,y=x,debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y))

	dd=punif_formula(x,y)

	op=list(	ml_params=dd$ml_params,
					ml_cdf=dd$ml_cdf,
					cp_cdf=dd$rh_cdf,
					cp_method=analytic_cpmethod())
	return(op)
}


