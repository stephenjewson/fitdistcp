#' Generalized Extreme Value Distribution, Predictions Based on a Calibrating Prior
#'
#' @inherit man description author references seealso return
#' @inheritParams man
#'
#' @inheritSection man Optional Return Values
#' @inheritSection man Optional Return Values (EVT models only)
#' @inheritSection man Optional Return Values (some EVT models only)
# #' @inheritSection man Details (homogeneous models)
#' @inheritSection man Details (non-homogeneous models)
# #' @inheritSection man Details (analytic integration)
#' @inheritSection man Details (DMGS integration)
#' @inheritSection man Details (RUST)
#'
#' @section Details of the Model:
#' The GEV distribution has distribution function
#' \deqn{F(x;\mu,\sigma,\xi)=\exp{(-t(x;\mu,\sigma,\xi))}}
#' where
#' \deqn{t(x;\mu,\sigma,\xi) =
#'    \begin{cases}
#'      {\left[1+\xi\left(\frac{x-\mu}{\sigma}\right)\right]}^{-1/\xi} & \text{if $\xi \ne 0$}\\
#'      \exp{\left(-\frac{x-\mu}{\sigma}\right)} & \text{if $\xi=0$}
#'    \end{cases}}
#' where
#' \eqn{x} is the random variable and
#' \eqn{\mu,\sigma>0,\xi} are the parameters.
#'
#' The calibrating prior we use is given by
#' \deqn{\pi(\mu,\sigma,\xi) \propto \frac{1}{\sigma}}
#' as given in Jewson et al. (2025).
#'
#' The code will stop with an error if the
#' input data gives a maximum likelihood
#' value for the shape parameter that lies outside the range \code{(minxi,maxxi)},
#' since outside this range there may be numerical problems.
#' Such values seldom occur
#' in real observed data for maxima.
#'
#' @example man/examples/example_110_gev.R
#'
#' @name gev_cp
NULL
#' @rdname gev_cp
#' @inheritParams man
#' @export
#'
qgev_cp=function(x,p=seq(0.1,0.9,0.1),ics=c(0,0,0),
	fdalpha=0.01,
	minxi=-1,maxxi=1,
	means=FALSE,waicscores=FALSE,extramodels=FALSE,
	pdf=FALSE,customprior=0,
	dmgs=TRUE,rust=FALSE,nrust=100000,
	pwm=FALSE,debug=FALSE){
#
# 1 intro
#
	stopifnot(is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1,length(ics)==3,fdalpha<1)
	alpha=1-p
	nx=length(x)
	nalpha=length(alpha)
	if(pdf){
		dalpha=pmin(fdalpha*alpha,fdalpha*(1-alpha))
		alpham=alpha-dalpha
		alphap=alpha+dalpha
	}
#
# 2 ml param estimate
#
	ics=gev_setics(x,ics)
	opt1=optim(ics,gev_loglik,x=x,control=list(fnscale=-1))
	v1hat=opt1$par[1]
	v2hat=opt1$par[2]
	v3hat=opt1$par[3]
	ml_params=c(v1hat,v2hat,v3hat)
#	gev_checkmle(ml_params,minxi,maxxi)
	if(debug)message("	v1hat,v2hat,v3hat=",v1hat,v2hat,v3hat,"")
	pw_params="pwm not selected"
	if(pwm)pw_params=gev_pwm_params(x)
	if(abs(v3hat)>=1){revert2ml=TRUE}else{revert2ml=FALSE}
# 3 aic
#
	ml_value=opt1$val
	maic=make_maic(ml_value,nparams=3)
#
# 4 calc ml quantiles and densities (vectorized over alpha)
#
	ml_quantiles=qgev((1-alpha),mu=v1hat,sigma=v2hat,xi=v3hat)
	if(v3hat<0){
		ml_max=v1hat-v2hat/v3hat
	} else {
		ml_max=Inf
	}
	fhat=dgev(ml_quantiles,mu=v1hat,sigma=v2hat,xi=v3hat)
	pw_quantiles="pwm not selected"
	if(pwm)pw_quantiles=qgev((1-alpha),mu=pw_params[1],sigma=pw_params[2],xi=pw_params[3])
#
# dmgs
#
	standard_errors="dmgs not selected"
	rh_flat_quantiles="dmgs not selected"
	cp_quantiles="dmgs not selected"
	ru_quantiles="dmgs not selected"
	custom_quantiles="dmgs not selected"
	ml_pdf="dmgs not selected"
	rh_flat_pdf="dmgs not selected"
	cp_pdf="dmgs not selected"
	waic1="dmgs not selected"
	waic2="dmgs not selected"
	ml_mean="dmgs not selected"
	rh_mean="dmgs not selected"
	rh_flat_mean="dmgs not selected"
	cp_mean="dmgs not selected"
	cp_method="dmgs not selected"
	custom_mean="dmgs not selected"
#
# 5 alpha pdf stuff
#
	if((dmgs)&&(!revert2ml)){
		if(debug)message("  ml_quantiles=",ml_quantiles,"")
		if(pdf){
			ml_quantilesm=qgev((1-alpham),mu=v1hat,sigma=v2hat,xi=v3hat)
			ml_quantilesp=qgev((1-alphap),mu=v1hat,sigma=v2hat,xi=v3hat)
			fhatm=dgev(ml_quantilesm,mu=v1hat,sigma=v2hat,xi=v3hat)
			fhatp=dgev(ml_quantilesp,mu=v1hat,sigma=v2hat,xi=v3hat)
		}
#
# 7 ldd
#
		if(debug)message("  calculate ldd")
		ldd=gev_ldda(x,v1hat,v2hat,v3hat)
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
#
# 8 lddd
#
		if(debug)message("  calculate lddd")
		lddd=gev_lddda(x,v1hat,v2hat,v3hat)
#
# 9 mu1
#
		mu1=gev_mu1fa(alpha,v1hat,v2hat,v3hat)
		if(pdf){
			mu1m=gev_mu1fa(alpham,v1hat,v2hat,v3hat)
			mu1p=gev_mu1fa(alphap,v1hat,v2hat,v3hat)
		}
#
# 10 mu2
#
		mu2=gev_mu2fa(alpha,v1hat,v2hat,v3hat)
		if(pdf){
			mu2m=gev_mu2fa(alpham,v1hat,v2hat,v3hat)
			mu2p=gev_mu2fa(alphap,v1hat,v2hat,v3hat)
		}
#
# 13 model 4: rh_flat with flat prior on shape
#
		lambdad_rh_flat=c(0,-1/v2hat,0)
		dq=dmgs(lddi,lddd,mu1,lambdad_rh_flat,mu2,dim=3)
		rh_flat_quantiles=ml_quantiles+dq/(nx*fhat)
		if(pdf){
			dqm=dmgs(lddi,lddd,mu1m,lambdad_rh_flat,mu2m,dim=3)
			dqp=dmgs(lddi,lddd,mu1p,lambdad_rh_flat,mu2p,dim=3)
			quantilesm=ml_quantilesm+dqm/(nx*fhatm)
			quantilesp=ml_quantilesp+dqp/(nx*fhatp)
			ml_pdf=fhat
			rh_flat_pdf=-(alphap-alpham)/(quantilesp-quantilesm)
		} else{
			ml_pdf=fhat
			rh_flat_pdf="pdf not selected"
		}
#
# 15 model 6: custom prior on shape parameter
#
		if(extramodels|means){
			lambdad_custom=c(0,-1/v2hat,customprior)
			dq=dmgs(lddi,lddd,mu1,lambdad_custom,mu2,dim=3)
			custom_quantiles=ml_quantiles+dq/(nx*fhat)
		} else {
			custom_quantiles="extramodels not selected"
		}
#
# 16 means
#
		means=gev_means(means,ml_params,lddi,lddd,
										lambdad_rh_flat,lambdad_custom,
										nx,dim=3)
		ml_mean				=means$ml_mean
		rh_flat_mean	=means$rh_flat_mean
		custom_mean		=means$custom_mean
#
# 17 waicscores
#
		waic=gev_waic(waicscores,x,v1hat,v2hat,v3hat,lddi,lddd,
			lambdad_rh_flat)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 19 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rgev_cp(nrust,x,rust=TRUE,mlcp=FALSE)
			ru_quantiles=makeq(rustsim$ru_deviates,p)
		}
#end of if(dmgs)
	} else {
		rh_flat_quantiles=ml_quantiles
	  ru_quantiles=ml_quantiles
	  pw_quantiles=ml_quantiles
	  custom_quantiles=ml_quantiles
	  rh_flat_pdf=ml_pdf
	  rh_flat_mean=ml_mean
	  custom_mean=ml_mean
	}

	list(	ml_params=ml_params,
				pw_params=pw_params,
				ml_value=ml_value,
#				ldd=ldd,
#				lddi=lddi,
				standard_errors=standard_errors,
				ml_quantiles=ml_quantiles,
				ml_max=ml_max,
				revert2ml=revert2ml,
				cp_quantiles=rh_flat_quantiles,
				ru_quantiles=ru_quantiles,
				pw_quantiles=pw_quantiles,
				custom_quantiles=custom_quantiles,
				ml_pdf=ml_pdf,
				cp_pdf=rh_flat_pdf,
				maic=maic,
				waic1=waic1,
				waic2=waic2,
				ml_mean=ml_mean,
				cp_mean=rh_flat_mean,
				custom_mean=custom_mean,
				cp_method=crhpflat_dmgs_cpmethod())

}
#' @rdname gev_cp
#' @inheritParams man
#' @export
rgev_cp=function(n,x,ics=c(0,0,0),
	minxi=-1,maxxi=1,
	extramodels=FALSE,rust=FALSE,mlcp=TRUE,debug=FALSE){

	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),length(ics)==3)

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qgev_cp(x,runif(n),ics=ics,extramodels=extramodels)
		ml_params=q$ml_params
#		gev_checkmle(ml_params,minxi,maxxi)
		ml_deviates=q$ml_quantiles
		ru_deviates=q$ru_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tgev_cp(n,x)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			ru_deviates[i]=rgev(1,mu=th[i,1],sigma=th[i,2],xi=th[i,3])
		}
	}
	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
			 cp_method=crhpflat_dmgs_cpmethod())

	return(op)

}
#' @rdname gev_cp
#' @inheritParams man
#' @export
dgev_cp=function(x,y=x,ics=c(0,0,0),
	minxi=-1,maxxi=1,
	extramodels=FALSE,
	rust=FALSE,nrust=1000,
	boot=FALSE,nboot=1000,
	debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),length(ics)==3)

	ics=gev_setics(x,ics)
	opt1=optim(ics,gev_loglik,x=x,control=list(fnscale=-1))
	v1hat=opt1$par[1]
	v2hat=opt1$par[2]
	v3hat=opt1$par[3]
	if(v3hat<=(-1)){revert2ml=TRUE}else{revert2ml=FALSE}
	ml_params=c(v1hat,v2hat,v3hat)
#	gev_checkmle(ml_params,minxi,maxxi)
	dd=dgevsub(x=x,y=y,ics=ics,minxi,maxxi)
	ru_pdf="rust not selected"
	bs_pdf="boot not selected"

	if(rust&&(!revert2ml)){
		th=tgev_cp("rust",nrust,x)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_pdf=ru_pdf+dgev(y,mu=th[ir,1],sigma=th[ir,2],xi=th[ir,3])
		}
		ru_pdf=ru_pdf/nrust
	} else {
		ru_pdf=dd$ml_pdf
	}

	if(boot){
		th=tgev_cp("boot",nrust,x)$theta_samples
		bs_pdf=numeric(length(y))
		for (ir in 1:nrust){
			bs_pdf=bs_pdf+dgev(y,mu=th[ir,1],sigma=th[ir,2],xi=th[ir,3])
		}
		bs_pdf=bs_pdf/nboot
	}

		op=list(
					ml_params=dd$ml_params,
					ml_pdf=dd$ml_pdf,
					revert2ml=revert2ml,
					ru_pdf=ru_pdf,
					bs_pdf=bs_pdf,
					cp_method=nopdfcdfmsg())
	return(op)
}
#' @rdname gev_cp
#' @inheritParams man
#' @export
pgev_cp=function(x,y=x,ics=c(0,0,0),
	minxi=-1,maxxi=1,
	extramodels=FALSE,
	rust=FALSE,nrust=1000,
	boot=FALSE,nboot=1000,
	debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),length(ics)==3)

	ics=gev_setics(x,ics)
	opt1=optim(ics,gev_loglik,x=x,control=list(fnscale=-1))
	v1hat=opt1$par[1]
	v2hat=opt1$par[2]
	v3hat=opt1$par[3]
	if(v3hat<=(-1)){revert2ml=TRUE}else{revert2ml=FALSE}
	ml_params=c(v1hat,v2hat,v3hat)
#	gev_checkmle(ml_params,minxi,maxxi)
	dd=dgevsub(x=x,y=y,ics=ics,minxi,maxxi)
	ru_cdf="rust not selected"
	bs_cdf="rust not selected"

	if(rust&&(!revert2ml)){
		th=tgev_cp("rust",nrust,x)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_cdf=ru_cdf+pgev(y,mu=th[ir,1],sigma=th[ir,2],xi=th[ir,3])
		}
		ru_cdf=ru_cdf/nrust
	} else {
		ru_pdf=dd$ml_pdf
	}


	if(boot){
		th=tgev_cp("boot",nrust,x)$theta_samples
		bs_cdf=numeric(length(y))
		for (ir in 1:nrust){
			bs_cdf=bs_cdf+pgev(y,mu=th[ir,1],sigma=th[ir,2],xi=th[ir,3])
		}
		bs_cdf=bs_cdf/nboot
	}


	op=list(
					ml_params=dd$ml_params,
					ml_cdf=dd$ml_cdf,
					revert2ml=revert2ml,
					ru_cdf=ru_cdf,
					bs_cdf=bs_cdf,
					cp_method=nopdfcdfmsg())
	return(op)
}
#' @rdname gev_cp
#' @inheritParams man
#' @export
tgev_cp=function(method,n,x,ics=c(0,0,0),extramodels=FALSE,debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),length(ics)==3)

	if(method=="rust"){
		ics=gev_setics(x,ics)
		th=ru(gev_logf,x=x,n=n,d=3,init=ics)
	} else if (method=="boot"){
		th=bgev(x=x,n=n)
	} else{
		message("tgev method not valid so stopping.\n")
		stop()
	}

	list(theta_samples=th$sim_vals)

}
