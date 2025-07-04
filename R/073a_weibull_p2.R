#' weibull Distribution with a Predictor on the Scale Parameter, Predictions Based on a Calibrating Prior
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
#' The Weibull distribution with predictor on the scale parameter
#' has exceedance distribution function
#' \deqn{S(x;k,a,b)=\exp\left(-\left(\frac{x}{\sigma(a,b)}\right)^{k}\right)}
#' where
#' \eqn{x \ge 0} is the random variable,
#' \eqn{k>0} is the shape parameter and
#' \eqn{\sigma=e^{a+bt}} is the scale parameter, modelled as a function of
#' parameters \eqn{a,b} and predictor \eqn{t}.
#'
#' The calibrating prior is given by the right Haar prior, which is
#' \deqn{\pi(k,\sigma) \propto \frac{1}{k}}
#' as given in Jewson et al. (2025).
#'
#' @example man/examples/example_073_weibull_p2.R
#'
#' @name weibull_p2_cp
NULL
#' @rdname weibull_p2_cp
#' @inheritParams man
#' @export
#'
qweibull_p2_cp=function(x,t,t0=NA,n0=NA,p=seq(0.1,0.9,0.1),fd1=0.01,d2=0.01,d3=0.01,
	means=FALSE,waicscores=FALSE,logscores=FALSE,
	dmgs=TRUE,rust=FALSE,nrust=100000,predictordata=TRUE,centering=TRUE,
	debug=FALSE,aderivs=TRUE){
#
# 1 intro
#
	if(debug)message("inside qweibull_p2")
	stopifnot(	is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1,
		length(t)==length(x),!x<0)
	alpha=1-p
	nx=length(x)
	nalpha=length(alpha)
	t0=maket0(t0,n0,t)
#
# 2 centering
#
  if(centering){
	  meant=mean(t)
    t=t-meant
    t0=t0-meant
  }
#
# 3 ml param estimate
#
	if(debug)message("calc ml param estimate")
	opt=optim(c(1,0,0),weibull_p2_loglik,x=x,t=t,
		control=list(fnscale=-1))
	v1hat=opt$par[1]
	v2hat=opt$par[2]
	v3hat=opt$par[3]
	ml_params=c(v1hat,v2hat,v3hat)
	if(debug)message("  ml_params=",ml_params)
#
# 4 predictordata
#
	prd=weibull_p2_predictordata(predictordata,x,t,t0,ml_params)
	predictedparameter=prd$predictedparameter
	adjustedx=prd$adjustedx
#
# 5 aic
#
	if(debug)message("calc aic")
 	ml_value=opt$val
	maic=make_maic(ml_value,nparams=3)
#
# 6 mle quantiles
#
	if(debug)message("calc mle quantiles")
	ml_quantiles=qweibull_p2((1-alpha),t0,shape=v1hat,ymn=v2hat,slope=v3hat)
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
# 7 lddi
#
		if(debug)message("calc ldd")
		if(aderivs) ldd=weibull_p2_ldda(x,t,v1hat,v2hat,v3hat)
		if(!aderivs)ldd=weibull_p2_ldd(x,t,v1hat,fd1,v2hat,d2,v3hat,d3)
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
#
# 8 lddd
#
		if(debug)message("calculate lddd")
		if(aderivs) lddd=weibull_p2_lddda(x,t,v1hat,v2hat,v3hat)
		if(!aderivs)lddd=weibull_p2_lddd(x,t,v1hat,fd1,v2hat,d2,v3hat,d3)
#
# 9 mu1
#
		if(debug)message("calculate mu1")
		if(aderivs) mu1=weibull_p2_mu1fa(alpha,t0,v1hat,v2hat,v3hat)
		if(!aderivs)mu1=weibull_p2_mu1f(alpha,t0,v1hat,fd1,v2hat,d2,v3hat,d3)
#
# 10 mu2
#
		if(debug)message("calculate mu2")
		if(aderivs) mu2=weibull_p2_mu2fa(alpha,t0,v1hat,v2hat,v3hat)
		if(!aderivs)mu2=weibull_p2_mu2f(alpha,t0,v1hat,fd1,v2hat,d2,v3hat,d3)
#
# 11 rhp
#
		if(debug)message("  rhp")
		lambdad_rhp=c(-1/v1hat,0,0) #this is rhp
#
# 12 rhp quantiles
#
		if(debug)message("  rhp quantiles")
		fhat=dweibull_p2(ml_quantiles,t0,shape=v1hat,ymn=v2hat,slope=v3hat,log=FALSE)
		dq=dmgs(lddi,lddd,mu1,lambdad_rhp,mu2,dim=3)
		rh_quantiles=ml_quantiles+dq/(nx*fhat)
#
# 13 means
#
		means=weibull_p2_means(means,t0,ml_params,lddi,lddd,lambdad_rhp,nx,dim=3)
		ml_mean=means$ml_mean
		rh_mean=means$rh_mean
#
# 14 waicscores
#
		waic=weibull_p2_waic(waicscores,x,t,v1hat,fd1,v2hat,d2,v3hat,d3,
			lddi,lddd,lambdad_rhp,aderivs)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 15 logscores
#
		logscores=weibull_p2_logscores(logscores,x,t,fd1,d2,d3,aderivs)
		ml_oos_logscore=logscores$ml_oos_logscore
		rh_oos_logscore=logscores$rh_oos_logscore
#
# 16 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rweibull_p2_cp(n=nrust,x=x,t=t,t0=t0,rust=TRUE,mlcp=FALSE)
			ru_quantiles=makeq(rustsim$ru_deviates,p)
		}
	} #end of if(dmgs)
#
# 17 decentering
#
  if(centering){
    ml_params[2]=ml_params[2]-ml_params[3]*meant
    if(predictordata)predictedparameter=predictedparameter-ml_params[3]*meant
  }

# return
	list(	ml_params=ml_params,
				ml_value=ml_value,
				predictedparameter=predictedparameter,
				adjustedx=adjustedx,
				ldd=ldd,
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
#' @rdname weibull_p2_cp
#' @inheritParams man
#' @export
rweibull_p2_cp=function(n,x,t,t0=NA,n0=NA,fd1=0.01,d2=0.01,d3=0.01,rust=FALSE,
	mlcp=TRUE,debug=FALSE,aderivs=TRUE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),
#		is.finite(t),!is.na(t),!x<0)
	stopifnot(is.finite(x),!is.na(x),
		is.finite(t),!is.na(t),!x<0)

	t0=maket0(t0,n0,t)

#
# centering
#
  meant=mean(t)
  t=t-meant
  t0=t0-meant

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qweibull_p2_cp(x=x,t=t,t0=t0,n0=NA,p=runif(n),fd1,d2,d3,aderivs=aderivs)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tweibull_p2_cp(n=n,x=x,t=t)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			sc=exp(th[i,2]+t0*th[i,3])
			ru_deviates[i]=rweibull(1,shape=th[i,1],scale=sc)
		}
	}

#
# decentering
#
  if(mlcp)ml_params[2]=ml_params[2]-ml_params[3]*meant

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
				cp_method=rhp_dmgs_cpmethod())

	return(op)

}
#' @rdname weibull_p2_cp
#' @inheritParams 	man
#' @export
dweibull_p2_cp=function(x,t,t0=NA,n0=NA,y=x,fd1=0.01,d2=0.01,d3=0.01,
	rust=FALSE,nrust=1000,centering=TRUE,debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
		is.finite(t),!is.na(t),!x<0,!y<0)

	t0=maket0(t0,n0,t)
#
# centering
#
  if(centering){
    meant=mean(t)
    t=t-meant
    t0=t0-meant
  }

	dd=dweibull_p2sub(x=x,t=t,y=y,t0=t0,fd1,d2,d3,aderivs=aderivs)
	ru_pdf="rust not selected"
	ml_params=dd$ml_params
	if(rust){
		th=tweibull_p2_cp(nrust,x,t)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			sc=exp(th[ir,2]+t0*th[ir,3])
			ru_pdf=ru_pdf+dweibull(y,shape=th[ir,1],scale=sc)
		}
		ru_pdf=ru_pdf/nrust
	}

#
# decentering
#
 if(centering)ml_params[2]=ml_params[2]-ml_params[3]*meant

	op=list(	ml_params=ml_params,
					ml_pdf=dd$ml_pdf,
					cp_pdf=dd$rh_pdf,
					ru_pdf=ru_pdf,
					cp_method=rhp_dmgs_cpmethod())
	return(op)

}
#' @rdname weibull_p2_cp
#' @inheritParams 	man
#' @export
pweibull_p2_cp=function(x,t,t0=NA,n0=NA,y=x,fd1=0.01,d2=0.01,d3=0.01,
	rust=FALSE,nrust=1000,centering=TRUE,debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
		is.finite(t),!is.na(t),!x<0,!y<0)

	t0=maket0(t0,n0,t)
#
# centering
#
  if(centering){
    meant=mean(t)
    t=t-meant
    t0=t0-meant
  }

	dd=dweibull_p2sub(x=x,t=t,y=y,t0=t0,fd1,d2,d3,aderivs=aderivs)
	ru_cdf="rust not selected"
	ml_params=dd$ml_params
	if(rust){
		th=tweibull_p2_cp(nrust,x,t)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			sc=exp(th[ir,2]+t0*th[ir,3])
			ru_cdf=ru_cdf+pweibull(y,shape=th[ir,1],scale=sc)
		}
		ru_cdf=ru_cdf/nrust
	}
#
# decentering
#
 if(centering)ml_params[2]=ml_params[2]-ml_params[3]*meant

	op=list(	ml_params=ml_params,
					ml_cdf=dd$ml_cdf,
					cp_cdf=dd$rh_cdf,
					ru_cdf=ru_cdf,
					cp_method=rhp_dmgs_cpmethod())
	return(op)
}
#' @rdname weibull_p2_cp
#' @inheritParams man
#' @export
tweibull_p2_cp=function(n,x,t,fd1=0.01,d2=0.01,d3=0.01,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),
#		is.finite(t),!is.na(t),!x<0)
	stopifnot(is.finite(x),!is.na(x),
		is.finite(t),!is.na(t),!x<0)

#
# centering
#
  meant=mean(t)
  t=t-meant

	th=ru(weibull_p2_logf,x=x,t=t,n=n,d=3,init=c(1,0,0))
  theta_samples=th$sim_vals

#
# decentering
#
  theta_samples[,1]=theta_samples[,1]-theta_samples[,2]*meant

	list(theta_samples=theta_samples)

}
