#' Frechet Distribution with Predictor, Predictions Based on a Calibrating Prior
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
#' The Frechet distribution with predictor has distribution function
#' \deqn{F(x;a,b,\lambda)=\exp\left(-\left(\frac{x-\mu}{\sigma(a,b)}\right)^{-\lambda}\right)}
#' where
#' \eqn{x>\mu} is the random variable,
#' \eqn{\sigma=e^{a+bt}} is the scale parameter, modelled as a function of
#' parameters \eqn{a,b} and predictor \eqn{t},
#' and \eqn{\lambda>0} is the shape parameter.
#' We consider \eqn{\mu} to be known (hence the \code{k1} in the name).
#'
#' The calibrating prior is given by the right Haar prior, which is
#' \deqn{\pi(a,b) \propto 1}
#' as given in Jewson et al. (2025).
#'
#' @example man/examples/example_071_frechet_p2k1.R
#'
#' @name frechet_p2k1_cp
NULL
#' @rdname frechet_p2k1_cp
#' @inheritParams man
#' @export
#'
qfrechet_p2k1_cp=function(x,t,t0=NA,n0=NA,p=seq(0.1,0.9,0.1),
	means=FALSE,waicscores=FALSE,logscores=FALSE,kloc=0,
	dmgs=TRUE,rust=FALSE,nrust=100000,predictordata=TRUE,centering=TRUE,
	debug=FALSE){
#
# 1 intro
#
	if(debug)message("inside qfrechet_p2k1")
	stopifnot(	is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1,
		length(t)==length(x),!x<kloc)
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
	lm=lm(x~t)
	v1start=lm$coefficients[1]
	v2start=lm$coefficients[2]
	xhat=v1start+v2start*t
	v3start=sd(lm$residuals)
	v1start=kloc
	v2start=0
	v3start=2
	opt=optim(c(v1start,v2start,v3start),frechet_p2k1_loglik,x=x,t=t,
		kloc=kloc,control=list(fnscale=-1))
	v1hat=opt$par[1]
	v2hat=opt$par[2]
	v3hat=opt$par[3]
	ml_params=c(v1hat,v2hat,v3hat)
	muhat=ml_params[1]+ml_params[2]*t
#	muhat0=makemuhat0(t0,n0,t,ml_params)
	residuals=x-muhat
	if(debug)message("  ml_params=",ml_params)
#
# 4 predictordata
#
	prd=frechet_p2k1_predictordata(predictordata,x,t,t0,ml_params,kloc)
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
	ml_quantiles=qfrechet_p2k1((1-alpha),t0,ymn=v1hat,slope=v2hat,lambda=v3hat,kloc=kloc)
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
		ldd=frechet_p2k1_ldda(x,t,v1hat,v2hat,v3hat,kloc)
		if(debug)message("ldd=",ldd)
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
#
# 8 lddd
#
		if(debug)message("calculate lddd")
		lddd=frechet_p2k1_lddda(x,t,v1hat,v2hat,v3hat,kloc)
#
# 9 mu1
#
		if(debug)message("calculate mu1")
		mu1=frechet_p2k1_mu1fa(alpha,t0,v1hat,v2hat,v3hat,kloc)
#
# 10 mu2
#
		if(debug)message("calculate mu2")
		mu2=frechet_p2k1_mu2fa(alpha,t0,v1hat,v2hat,v3hat,kloc)
#
# 11 rhp
#
		if(debug)message("  rhp")
		lambdad_rhp=c(0,0,-1/v3hat) #this is rhp
#
# 12 rhp quantiles
#
		if(debug)message("  rhp quantiles")
		fhat=dfrechet_p2k1(ml_quantiles,t0,ymn=v1hat,slope=v2hat,lambda=v3hat,log=FALSE,kloc=kloc)
		dq=dmgs(lddi,lddd,mu1,lambdad_rhp,mu2,dim=3)
		rh_quantiles=ml_quantiles+dq/(nx*fhat)
#
# 13 means
#
		means=frechet_p2k1_means(means,t0,ml_params,lddi,lddd,lambdad_rhp,nx,dim=3,kloc=kloc)
		ml_mean=means$ml_mean
		rh_mean=means$rh_mean
#
# 14 waicscores
#
		waic=frechet_p2k1_waic(waicscores,x,t,v1hat,v2hat,v3hat,kloc=kloc,
			lddi,lddd,lambdad_rhp)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 15 logscores
#
		logscores=frechet_p2k1_logscores(logscores,x,t,kloc=kloc)
		ml_oos_logscore=logscores$ml_oos_logscore
		rh_oos_logscore=logscores$rh_oos_logscore
#
# 16 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rfrechet_p2k1_cp(nrust,x,t=t,t0=t0,kloc=kloc,rust=TRUE,mlcp=FALSE)
			ru_quantiles=makeq(rustsim$ru_deviates,p)
		}
	} #end of if(dmgs)
#
# 17 decentering
#
  if(centering){
    ml_params[1]=ml_params[1]-ml_params[2]*meant
    if(predictordata)predictedparameter=predictedparameter-ml_params[2]*meant
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
#' @rdname frechet_p2k1_cp
#' @inheritParams man
#' @export
rfrechet_p2k1_cp=function(n,x,t,t0=NA,n0=NA,
	kloc=0,rust=FALSE,mlcp=TRUE,centering=TRUE,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),
#		is.finite(t),!is.na(t),!x<kloc)
	stopifnot(is.finite(x),!is.na(x),
		is.finite(t),!is.na(t),!x<kloc)

	t0=maket0(t0,n0,t)

#
# centering
#
  if(centering){
  	meant=mean(t)
  	t=t-meant
  	t0=t0-meant
  }

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qfrechet_p2k1_cp(x,t,t0=t0,n0=NA,p=runif(n),kloc=kloc,
			centering=centering)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tfrechet_p2k1_cp(n,x,t)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			sg=exp(th[i,1]+t0*th[i,2])
			ru_deviates[i]=rfrechet(1,mu=kloc,sigma=sg,lambda=th[i,3])
		}
	}

#
# decentering
#
  if(mlcp&centering)ml_params[1]=ml_params[1]-ml_params[2]*meant

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
				cp_method=rhp_dmgs_cpmethod())

	return(op)

}
#' @rdname frechet_p2k1_cp
#' @inheritParams 	man
#' @export
dfrechet_p2k1_cp=function(x,t,t0=NA,n0=NA,y=x,
	kloc=0,rust=FALSE,nrust=1000,centering=TRUE,debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
		is.finite(t),!is.na(t),!x<kloc,!y<kloc)

	t0=maket0(t0,n0,t)
#
# centering
#
  if(centering){
    meant=mean(t)
    t=t-meant
    t0=t0-meant
  }

	dd=dfrechet_p2k1sub(x=x,t=t,y=y,t0=t0,kloc=kloc)
	ru_pdf="rust not selected"
	ml_params=dd$ml_params
	if(rust){
		th=tfrechet_p2k1_cp(nrust,x,t)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			sg=exp(th[ir,1]+t0*th[ir,2])
			ru_pdf=ru_pdf+dfrechet(y,mu=kloc,sigma=sg,lambda=th[ir,3])
		}
		ru_pdf=ru_pdf/nrust
	}
#
# decentering
#
	if(centering)ml_params[1]=ml_params[1]-ml_params[2]*meant

	op=list(	ml_params=ml_params,
				ml_pdf=dd$ml_pdf,
				cp_pdf=dd$rh_pdf,
				ru_pdf=ru_pdf,
				cp_method=rhp_dmgs_cpmethod())
	return(op)
}
#' @rdname frechet_p2k1_cp
#' @inheritParams 	man
#' @export
pfrechet_p2k1_cp=function(x,t,t0=NA,n0=NA,y=x,
	kloc=0,rust=FALSE,nrust=1000,centering=TRUE,debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
		is.finite(t),!is.na(t),!x<kloc,!y<kloc)

	t0=maket0(t0,n0,t)
#
# centering
#
  if(centering){
    meant=mean(t)
    t=t-meant
    t0=t0-meant
  }

	dd=dfrechet_p2k1sub(x=x,t=t,y=y,t0=t0,kloc=kloc)
	ru_cdf="rust not selected"
	ml_params=dd$ml_params
	if(rust){
		th=tfrechet_p2k1_cp(nrust,x,t)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			sg=exp(th[ir,1]+t0*th[ir,2])
			ru_cdf=ru_cdf+pfrechet(y,mu=kloc,sigma=sg,lambda=th[ir,3])
		}
		ru_cdf=ru_cdf/nrust
	}
#
# decentering
#
  if(centering)ml_params[1]=ml_params[1]-ml_params[2]*meant

	op=list(	ml_params=ml_params,
				ml_cdf=dd$ml_cdf,
				cp_cdf=dd$rh_cdf,
				ru_cdf=ru_cdf,
				cp_method=rhp_dmgs_cpmethod())
	return(op)
}
#' @rdname frechet_p2k1_cp
#' @inheritParams man
#' @export
tfrechet_p2k1_cp=function(n,x,t,kloc=0,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),
#		is.finite(t),!is.na(t),!x<kloc)
	stopifnot(is.finite(x),!is.na(x),
		is.finite(t),!is.na(t),!x<kloc)

#
# centering
#
  meant=mean(t)
  t=t-meant

	th=ru(frechet_p2k1_logf,x=x,t=t,kloc=kloc,n=n,d=3,init=c(0,0,1))
  theta_samples=th$sim_vals

#
# decentering
#
  theta_samples[,1]=theta_samples[,1]-theta_samples[,2]*meant

	list(theta_samples=theta_samples)

}
