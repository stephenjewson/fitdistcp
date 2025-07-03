#' t Distribution with a Predictor, Predictions Based on a Calibrating Prior
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
#' The t distribution with a predictor
#' (also known as the location-scale t distribution with a predictor,
#' hence the name \code{lst}),
#' has probability density function
#' \deqn{f(x;a,b,\sigma)
#' =\frac{\Gamma((\nu+1)/2)}{\sqrt{\pi\nu}\sigma\Gamma(\nu/2)}
#' \left(1+\frac{(x-\mu(a,b))^2}{\sigma^2\nu}\right)^{(\nu+1)/2}}
#' where
#' \eqn{x} is the random variable,
#' \eqn{\mu=a+bt} is the location parameter,
#' and \eqn{\sigma>0} is the scale parameter.
#' We consider the degrees of freedom \eqn{\nu} to be known
#' (hence the \code{k3} in the name).
#'
#' The calibrating prior is given by the right Haar prior, which is
#' \deqn{\pi(\sigma) \propto \frac{1}{\sigma}}
#' as given in Jewson et al. (2025).
#'
#' @example man/examples/example_63_lst_p1k3.R
#'
#' @name lst_p1k3_cp
NULL
#' @rdname lst_p1k3_cp
#' @inheritParams man
#' @export
#'
qlst_p1k3_cp=function(x,t,t0=NA,n0=NA,p=seq(0.1,0.9,0.1),d1=0.01,d2=0.01,fd3=0.01,kdf=10,
	means=FALSE,waicscores=FALSE,logscores=FALSE,
	dmgs=TRUE,rust=FALSE,nrust=100000,predictordata=TRUE,centering=TRUE,
	debug=FALSE,aderivs=TRUE){
#
# 1 intro
#
	if(debug)message("inside qlst_p1k3")
	stopifnot(	is.finite(x),
							!is.na(x),
							is.finite(p),
							!is.na(p),
							p>0,
							p<1,

							!is.na(kdf))
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
	ics=lst_p1k3_setics(x,t,c(0,0,0))
	opt=optim(ics,lst_p1k3_loglik,x=x,t=t,kdf=kdf,
		control=list(fnscale=-1))
	v1hat=opt$par[1]
	v2hat=opt$par[2]
	v3hat=opt$par[3]
	ml_params=c(v1hat,v2hat,v3hat)
	muhat=ml_params[1]+ml_params[2]*t
	muhat0=makemuhat0(t0,n0,t,ml_params)
	residuals=x-muhat
	if(debug)message("  ml_params=",ml_params)
#
# 4 predictordata
#
	prd=lst_p1k3_predictordata(predictordata,x,t,t0,ml_params,kdf)
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
	ml_quantiles=qlst_p1k3((1-alpha),t0,ymn=v1hat,slope=v2hat,sigma=v3hat,kdf=kdf)
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
		if(aderivs) ldd=lst_p1k3_ldda(x,t,v1hat,v2hat,v3hat,kdf=kdf)
		if(!aderivs)ldd=lst_p1k3_ldd(x,t,v1hat,d1,v2hat,d2,v3hat,fd3,kdf=kdf)
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
#
# 8 lddd
#
		if(debug)message("calculate lddd")
		if(aderivs) lddd=lst_p1k3_lddda(x,t,v1hat,v2hat,v3hat,kdf=kdf)
		if(!aderivs)lddd=lst_p1k3_lddd(x,t,v1hat,d1,v2hat,d2,v3hat,fd3,kdf=kdf)
#
# 9 mu1
#
		if(debug)message("calculate mu1")
		mu1=lst_p1k3_mu1f(alpha,t0,v1hat,d1,v2hat,d2,v3hat,fd3,kdf=kdf)
#
# 10 mu2
#
		if(debug)message("calculate mu2")
		mu2=lst_p1k3_mu2f(alpha,t0,v1hat,d1,v2hat,d2,v3hat,fd3,kdf=kdf)
#
# 11 rhp
#
		if(debug)message("  rhp")
		lambdad_rhp=c(0,0,-1/v3hat) #this is rhp
#
# 12 rhp quantiles
#
		if(debug)message("  rhp quantiles")
		fhat=dlst_p1k3(ml_quantiles,t0,ymn=v1hat,slope=v2hat,sigma=v3hat,log=FALSE,kdf=kdf)
		dq=dmgs(lddi,lddd,mu1,lambdad_rhp,mu2,dim=3)
		rh_quantiles=ml_quantiles+dq/(nx*fhat)
#
# 13 means (might as well always calculate)
#
		ml_mean=muhat0
		rh_mean=ml_mean
#
# 14 waicscores
#
		waic=lst_p1k3_waic(waicscores,x,t,v1hat,d1,v2hat,d2,v3hat,fd3,kdf=kdf,
			lddi,lddd,lambdad_rhp,aderivs)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 15 logscores
#
		logscores=lst_p1k3_logscores(logscores,x,t,d1,d2,fd3,kdf=kdf,aderivs)
		ml_oos_logscore=logscores$ml_oos_logscore
		rh_oos_logscore=logscores$rh_oos_logscore
#
# 16 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rlst_p1k3_cp(nrust,x,t=t,t0=t0,kdf=kdf,rust=TRUE,mlcp=FALSE)
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
#' @rdname lst_p1k3_cp
#' @inheritParams man
#' @export
rlst_p1k3_cp=function(n,x,t,t0=NA,n0=NA,d1=0.01,d2=0.01,fd3=0.01,
	kdf=10,rust=FALSE,mlcp=TRUE,centering=TRUE,debug=FALSE,aderivs=TRUE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),is.finite(t),!is.na(t),
#						is.finite(kdf),!is.na(kdf))
	stopifnot(is.finite(x),!is.na(x),is.finite(t),!is.na(t),
						is.finite(kdf),!is.na(kdf))

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
		q=qlst_p1k3_cp(x,t,t0=t0,n0=NA,p=runif(n),d1,d2,fd3,kdf=kdf,
			centering=centering,aderivs=aderivs)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tlst_p1k3_cp(n,x,t)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			mu=th[i,1]+t0*th[i,2]
			ru_deviates[i]=rlst(1,mu=mu,sigma=th[i,3],df=kdf)
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
#' @rdname lst_p1k3_cp
#' @inheritParams man
#' @export
dlst_p1k3_cp=function(x,t,t0=NA,n0=NA,y=x,d1=0.01,d2=0.01,fd3=0.01,
	kdf=10,rust=FALSE,nrust=1000,centering=TRUE,debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
						is.finite(t),!is.na(t),is.finite(kdf),!is.na(kdf))#

	t0=maket0(t0,n0,t)
#
# centering
#
	if(centering){
		meant=mean(t)
		t=t-meant
		t0=t0-meant
	}

	dd=dlst_p1k3sub(x=x,t=t,y=y,t0=t0,d1,d2,fd3,kdf=kdf,aderivs=aderivs)
	ru_pdf="rust not selected"
	ml_params=dd$ml_params
	if(rust){
		th=tlst_p1k3_cp(nrust,x,t,kdf)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+t0*th[ir,2]
			ru_pdf=ru_pdf+dlst(y,mu=mu,sigma=th[ir,3],df=kdf)
		}
		ru_pdf=ru_pdf/nrust
	}
#
# decentering
#
	if(centering){
		ml_params[1]=ml_params[1]-ml_params[2]*meant
	}

	op=list(	ml_params=ml_params,
				ml_pdf=dd$ml_pdf,
				cp_pdf=dd$rh_pdf,
				ru_pdf=ru_pdf,
				cp_method=rhp_dmgs_cpmethod())
	return(op)
}
#' @rdname lst_p1k3_cp
#' @inheritParams man
#' @export
plst_p1k3_cp=function(x,t,t0=NA,n0=NA,y=x,d1=0.01,d2=0.01,fd3=0.01,
	kdf=10,rust=FALSE,nrust=1000,centering=TRUE,debug=FALSE,aderivs=TRUE){

	t0=maket0(t0,n0,t)

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
						is.finite(t),!is.na(t),is.finite(kdf),!is.na(kdf))

#
# centering
#
	if(centering){
		meant=mean(t)
		t=t-meant
		t0=t0-meant
	}


	dd=dlst_p1k3sub(x=x,t=t,y=y,t0=t0,d1,d2,fd3,kdf=kdf,aderivs=aderivs)
	ru_cdf="rust not selected"
	ml_params=dd$ml_params
	if(rust){
		th=tlst_p1k3_cp(nrust,x,t,kdf)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+t0*th[ir,2]
			ru_cdf=ru_cdf+plst(y,mu=mu,sigma=th[ir,3],df=kdf)
		}
		ru_cdf=ru_cdf/nrust
	}

#
# decentering
#
	if(centering){
		ml_params[1]=ml_params[1]-ml_params[2]*meant
	}

	op=list(	ml_params=ml_params,
				ml_cdf=dd$ml_cdf,
				cp_cdf=dd$rh_cdf,
				ru_cdf=ru_cdf,
				cp_method=rhp_dmgs_cpmethod())
	return(op)
}
#' @rdname lst_p1k3_cp
#' @inheritParams man
#' @export
tlst_p1k3_cp=function(n,x,t,d1=0.01,d2=0.01,fd3=0.01,kdf=10,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),is.finite(t),!is.na(t),
#						is.finite(kdf),!is.na(kdf))
	stopifnot(is.finite(x),!is.na(x),is.finite(t),!is.na(t),
						is.finite(kdf),!is.na(kdf))

#
# centering
#
	meant=mean(t)
	t=t-meant

	th=ru(lst_p1k3_logf,x=x,t=t,kdf=kdf,n=n,d=3,init=c(0,0,1))
	theta_samples=th$sim_vals

#
# decentering
#
	theta_samples[,1]=theta_samples[,1]-theta_samples[,2]*meant

	list(theta_samples=theta_samples)

}
