#' GEV Distribution with Known Shape with a Predictor, Predictions Based on a Calibrating Prior
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
#' The GEV distribution with known shape with a predictor has distribution function
#' \deqn{F(x;a,b,\sigma)=\exp{(-t(x;\mu(a,b),\sigma))}}
#' where
#' \deqn{t(x;a,b,\sigma) =
#'    \begin{cases}
#'      {\left[1+\xi\left(\frac{x-\mu(a,b)}{\sigma}\right)\right]}^{-1/\xi} & \text{if $\xi \ne 0$}\\
#'      \exp{(-\frac{x-\mu(a,b)}{\sigma})} & \text{if $\xi=0$}
#'    \end{cases}}
#' where
#' \eqn{x} is the random variable,
#' \eqn{\mu=a+bt} is the location parameter,
#' \eqn{\sigma>0} is the shape parameter and
#' \eqn{\xi} is known (hence the \code{k3} in the name).
#'
#' The calibrating prior we use is given by
#' \deqn{\pi(\mu,\sigma) \propto \frac{1}{\sigma}}
#' as given in Jewson et al. (2025).
#'
#' @example man/examples/example_074_gev_p1k3.R
#'
#' @name gev_p1k3_cp
NULL
#' @rdname gev_p1k3_cp
#' @inheritParams man
#' @export
#'
qgev_p1k3_cp=function(x,t,t0=NA,n0=NA,p=seq(0.1,0.9,0.1),
	fdalpha=0.01,kshape=0,
	means=FALSE,waicscores=FALSE,
	pdf=FALSE,dmgs=TRUE,rust=FALSE,nrust=100000,predictordata=TRUE,
	centering=TRUE,debug=FALSE){
#
# 1 intro
#
	debug=TRUE
	debug=FALSE
	if(debug)message("inside qgev_p1k3")
	stopifnot(	is.finite(x),
							!is.na(x),
							is.finite(p),
							!is.na(p),
							p>0,
							p<1,

							!is.na(kshape))
	alpha=1-p
	nx=length(x)
	nalpha=length(alpha)
	t0=maket0(t0,n0,t)
	if(pdf){
		dalpha=pmin(fdalpha*alpha,fdalpha*(1-alpha))
		alpham=alpha-dalpha
		alphap=alpha+dalpha
	}
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
	v3start=100000000 #use a very large value to avoid missing some of the range
	opt1=optim(c(v1start,v2start,v3start),gev_p1k3_loglik,x=x,t=t,
		kshape=kshape,control=list(fnscale=-1))
	v1hat=opt1$par[1]
	v2hat=opt1$par[2]
	v3hat=opt1$par[3]
	ml_params=c(v1hat,v2hat,v3hat)
	muhat=ml_params[1]+ml_params[2]*t
	if(debug)message("  ml_params=",ml_params)
	if(kshape<=(-1)){revert2ml=TRUE}else{revert2ml=FALSE}
#
# 4 predictordata
#
	prd=gev_p1k3_predictordata(predictordata,x,t,t0,ml_params,kshape)
	predictedparameter=prd$predictedparameter
	adjustedx=prd$adjustedx
#
# 5 aic
#
	if(debug)message("calc aic")
 	ml_value=opt1$val
	maic=make_maic(ml_value,nparams=3)
#
# 6 mle quantiles
#
	if(debug)message("calc mle quantiles")
	ml_quantiles=qgev_p1k3((1-alpha),t0,ymn=v1hat,slope=v2hat,sigma=v3hat,kshape=kshape)
#
# dmgs
#
	standard_errors="dmgs not selected"
	rh_quantiles="dmgs not selected"
	ru_quantiles="dmgs not selected"
	ml_pdf="dmgs not selected"
	cp_pdf="dmgs not selected"
	waic1="dmgs not selected"
	waic2="dmgs not selected"
	ml_mean="dmgs not selected"
	rh_mean="dmgs not selected"
	cp_mean="dmgs not selected"
	cp_method="dmgs not selected"
	if((dmgs)&&(!revert2ml)){
		if(pdf){
			ml_quantilesm=qgev((1-alpham),mu=muhat,sigma=v3hat,xi=kshape)
			ml_quantilesp=qgev((1-alphap),mu=muhat,sigma=v3hat,xi=kshape)
			fhatm=dgev(ml_quantilesm,mu=muhat,sigma=v3hat,xi=kshape)
			fhatp=dgev(ml_quantilesp,mu=muhat,sigma=v3hat,xi=kshape)
		}
#
# 7 lddi
#
		if(debug)message("calc ldd")
		ldd=gev_p1k3_ldda(x,t,v1hat,v2hat,v3hat,kshape=kshape)
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
#
# 8 lddd
#
		if(debug)message("calculate lddd")
		lddd=gev_p1k3_lddda(x,t,v1hat,v2hat,v3hat,kshape=kshape)
#
# 9 mu1
#
		if(debug)message("calculate mu1")
		mu1=gev_p1k3_mu1fa(alpha,t0,v1hat,v2hat,v3hat,kshape=kshape)

		if(pdf){
			mu1m=gev_p1k3_mu1fa(alpham,t0,v1hat,v2hat,v3hat,kshape=kshape)
			mu1p=gev_p1k3_mu1fa(alphap,t0,v1hat,v2hat,v3hat,kshape=kshape)
		}
#
# 10 mu2
#
		if(debug)message("calculate mu2")
		mu2=gev_p1k3_mu2fa(alpha,t0,v1hat,v2hat,v3hat,kshape=kshape)
		if(pdf){
			mu2m=gev_p1k3_mu2fa(alpham,t0,v1hat,v2hat,v3hat,kshape=kshape)
			mu2p=gev_p1k3_mu2fa(alphap,t0,v1hat,v2hat,v3hat,kshape=kshape)
		}
#
# 11 rhp
#
		if(debug)message("  rhp")
		lambdad_cp=c(0,0,-1/v3hat) #this is rhp
#
# 12 rhp quantiles
#
		if(debug)message("  rhp quantiles")
		fhat=dgev_p1k3(ml_quantiles,t0,ymn=v1hat,slope=v2hat,sigma=v3hat,log=FALSE,kshape=kshape)
		dq=dmgs(lddi,lddd,mu1,lambdad_cp,mu2,dim=3)
		rh_quantiles=ml_quantiles+dq/(nx*fhat)
		if(pdf){
			dqm=dmgs(lddi,lddd,mu1m,lambdad_cp,mu2m,dim=3)
			dqp=dmgs(lddi,lddd,mu1p,lambdad_cp,mu2p,dim=3)
			quantilesm=ml_quantilesm+dqm/(nx*fhatm)
			quantilesp=ml_quantilesp+dqp/(nx*fhatp)
			ml_pdf=fhat
			rh_pdf=-(alphap-alpham)/(quantilesp-quantilesm)
		} else{
			ml_pdf=fhat
			rh_pdf="pdf not selected"
		}
#
# 13 means
#
		means=gev_p1k3_means(means,t0,ml_params,kshape,nx)
		ml_mean				=means$ml_mean
		rh_mean				=means$cp_mean
#
# 14 waicscores
#
		waic=gev_p1k3_waic(waicscores,x,t,v1hat,v2hat,v3hat,kshape=kshape,
			lddi,lddd,lambdad_cp)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 16 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rgev_p1k3_cp(nrust,x,t=t,t0=t0,kshape=kshape,rust=TRUE,mlcp=FALSE)
			ru_quantiles=makeq(rustsim$ru_deviates,p)
		}
	} else {
		rh_quantiles=ml_quantiles
		ru_quantiles=ml_quantiles
		rh_pdf=ml_pdf
		rh_mean=ml_mean
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
				revert2ml=revert2ml,
				ml_quantiles=ml_quantiles,
				cp_quantiles=rh_quantiles,
				ru_quantiles=ru_quantiles,
				ml_pdf=ml_pdf,
				cp_pdf=rh_pdf,
				maic=maic,
				waic1=waic1,
				waic2=waic2,
				ml_mean=ml_mean,
				cp_mean=rh_mean,
				cp_method=rhp_dmgs_cpmethod())

}
#' @rdname gev_p1k3_cp
#' @inheritParams man
#' @export
rgev_p1k3_cp=function(n,x,t,t0=NA,n0=NA,
	kshape=0,rust=FALSE,mlcp=TRUE,centering=TRUE,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),is.finite(t),!is.na(t),
#						is.finite(kshape),!is.na(kshape))
	stopifnot(is.finite(x),!is.na(x),is.finite(t),!is.na(t),
						is.finite(kshape),!is.na(kshape))

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
		q=qgev_p1k3_cp(x,t,t0=t0,n0=NA,p=runif(n),kshape=kshape,
			centering=centering)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tgev_p1k3_cp(n,x,t)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			mu=th[i,1]+t0*th[i,2]
			ru_deviates[i]=rgev(1,mu=mu,sigma=th[i,3],xi=kshape)
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
#' @rdname gev_p1k3_cp
#' @inheritParams man
#' @export
dgev_p1k3_cp=function(x,t,t0=NA,n0=NA,y=x,
	kshape=0,rust=FALSE,nrust=1000,centering=TRUE,
	debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
						is.finite(t),!is.na(t),is.finite(kshape),!is.na(kshape))

	t0=maket0(t0,n0,t)
#
# centering
#
  if(centering){
    meant=mean(t)
    t=t-meant
    t0=t0-meant
  }

	dd=dgev_p1k3sub(x=x,t=t,y=y,t0=t0,kshape=kshape)
	ru_pdf="rust not selected"
	ml_params=dd$ml_params
	if(kshape<=(-1)){revert2ml=TRUE}else{revert2ml=FALSE}

	if(rust&&(!revert2ml)){
		th=tgev_p1k3_cp(nrust,x,t,kshape)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+t0*th[ir,2]
			ru_pdf=ru_pdf+dgev(y,mu=mu,sigma=th[ir,3],xi=kshape)
		}
		ru_pdf=ru_pdf/nrust
	} else {
		ru_pdf=dd$ml_pdf
	}
#
# decentering
#
 if(centering)ml_params[1]=ml_params[1]-ml_params[2]*meant

	op=list(
				ml_params=ml_params,
				ml_pdf=dd$ml_pdf,
				revert2ml=revert2ml,
				ru_pdf=ru_pdf,
				cp_method=nopdfcdfmsg())
	return(op)
}
#' @rdname gev_p1k3_cp
#' @inheritParams man
#' @export
pgev_p1k3_cp=function(x,t,t0=NA,n0=NA,y=x,
	kshape=0,rust=FALSE,nrust=1000,centering=TRUE,debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
						is.finite(t),!is.na(t),is.finite(kshape),!is.na(kshape))

	t0=maket0(t0,n0,t)

#
# centering
#
  if(centering){
    meant=mean(t)
    t=t-meant
    t0=t0-meant
  }

	dd=dgev_p1k3sub(x=x,t=t,y=y,t0=t0,kshape=kshape)
	ru_cdf="rust not selected"
	ml_params=dd$ml_params
	if(kshape<=(-1)){revert2ml=TRUE}else{revert2ml=FALSE}

	if(rust&&(!revert2ml)){
		th=tgev_p1k3_cp(nrust,x,t,kshape)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+t0*th[ir,2]
			ru_cdf=ru_cdf+pgev(y,mu=mu,sigma=th[ir,3],xi=kshape)
		}
		ru_cdf=ru_cdf/nrust
	} else {
		ru_pdf=dd$ml_pdf
	}
#
# decentering
#
 if(centering)ml_params[1]=ml_params[1]-ml_params[2]*meant

	op=list(
				ml_params=ml_params,
				ml_cdf=dd$ml_cdf,
				revert2ml=revert2ml,
				ru_cdf=ru_cdf,
				cp_method=nopdfcdfmsg())

	return(op)
}
#' @rdname gev_p1k3_cp
#' @inheritParams man
#' @export
tgev_p1k3_cp=function(n,x,t,kshape=0,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),is.finite(t),!is.na(t),
#						is.finite(kshape),!is.na(kshape))
	stopifnot(is.finite(x),!is.na(x),is.finite(t),!is.na(t),
						is.finite(kshape),!is.na(kshape))

#
# centering
#
  meant=mean(t)
  t=t-meant

	th=ru(gev_p1k3_logf,x=x,t=t,kshape=kshape,n=n,d=3,init=c(0,0,1))
  theta_samples=th$sim_vals
#
# decentering
#
  theta_samples[,1]=theta_samples[,1]-theta_samples[,2]*meant

	list(theta_samples=theta_samples)

}
