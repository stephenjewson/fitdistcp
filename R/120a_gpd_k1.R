#' Generalized Pareto Distribution with Known Location Parameter, Predictions Based on a Calibrating Prior
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
#' @details
#' The GP distribution has exceedcance distribution function
#' \deqn{S(x;\mu,\sigma,\xi) =
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
#' @example man/examples/example_120_gpd_k1.R
#'
#' @name gpd_k1_cp
NULL
#' @rdname gpd_k1_cp
#' @inheritParams man
#' @export
#'
qgpd_k1_cp=function(x,p=seq(0.1,0.9,0.1),kloc=0,ics=c(0,0),
	fdalpha=0.01,customprior=0,
	minxi=-1,maxxi=2.0,
	means=FALSE,waicscores=FALSE,extramodels=FALSE,
	pdf=FALSE,dmgs=TRUE,rust=FALSE,nrust=100000,debug=FALSE){
#
# 1 intro
#
	stopifnot(is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1,length(ics)==2,!x<0)
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
	ics=gpd_k1_setics(x,ics)
	opt1=optim(ics,gpd_k1_loglik,x=x,kloc=kloc,control=list(fnscale=-1))
	v1hat=opt1$par[1]
	v2hat=opt1$par[2]
	ml_params=c(v1hat,v2hat)
#	gpd_k1_checkmle(ml_params,kloc,minxi,maxxi)
	if(debug)message("	v1hat,v2hat=",v1hat,v2hat)
	if((abs(v2hat)>=1)||(v2hat>30)){revert2ml=TRUE}else{revert2ml=FALSE}
#
# 3 aic
#
	ml_value=opt1$val
	maic=make_maic(ml_value,nparams=2)
#
# 4 calc ml quantiles and densities (vectorized over alpha)
#
	ml_quantiles=extraDistr::qgpd((1-alpha),mu=kloc,sigma=v1hat,xi=v2hat)
	if(v2hat<0){
		ml_max=kloc-v1hat/v2hat
	} else {
		ml_max=Inf
	}
	fhat=extraDistr::dgpd(ml_quantiles,mu=kloc,sigma=v1hat,xi=v2hat)
#
# dmgs
#
	standard_errors="dmgs not selected"
	cp_quantiles="dmgs not selected"
	rh_flat_quantiles="dmgs not selected"
	lp_quantiles="dmgs not selected"
	lp2_quantiles="dmgs not selected"
	dpi_quantiles="dmgs not selected"
	ru_quantiles="dmgs not selected"
	ml_pdf="dmgs not selected"
	cp_pdf="dmgs not selected"
	rh_flat_pdf="dmgs not selected"
	waic1="dmgs not selected"
	waic2="dmgs not selected"
	ml_mean="dmgs not selected"
	rh_flat_mean="dmgs not selected"
	cp_mean="dmgs not selected"
	cp_method="dmgs not selected"

	if((dmgs)&&(!revert2ml)){
#
# 5 alpha pdf stuff
#
		if(pdf){
			ml_quantilesm=extraDistr::qgpd((1-alpham),mu=kloc,sigma=v1hat,xi=v2hat)
			ml_quantilesp=extraDistr::qgpd((1-alphap),mu=kloc,sigma=v1hat,xi=v2hat)
			fhatm=extraDistr::dgpd(ml_quantilesm,mu=kloc,sigma=v1hat,xi=v2hat)
			fhatp=extraDistr::dgpd(ml_quantilesp,mu=kloc,sigma=v1hat,xi=v2hat)
		}
		if(debug)message("  ml_quantiles=",ml_quantiles)
#
# 7 ldd (two versions)
#
		if(debug)message("  calculate ldd")
		ldd=gpd_k1_ldda(x,v1hat,v2hat,kloc)
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)

		if(debug)message("  ldd=",ldd)
#
# 8 lddd (two versions)
#
		if(debug)message("  calculate lddd")
		lddd=gpd_k1_lddda(x,v1hat,v2hat,kloc)
#
# 9 mu1 (two versions)
#
		mu1=gpd_k1_mu1fa(alpha,v1hat,v2hat,kloc)

		if(pdf){
			mu1m=gpd_k1_mu1fa(alpham,v1hat,v2hat,kloc)
			mu1p=gpd_k1_mu1fa(alphap,v1hat,v2hat,kloc)
		}
#
# 10 mu2 (two versions)
#
		mu2=gpd_k1_mu2fa(alpha,v1hat,v2hat,kloc)

		if(pdf){
			mu2m=gpd_k1_mu2fa(alpham,v1hat,v2hat,kloc)
			mu2p=gpd_k1_mu2fa(alphap,v1hat,v2hat,kloc)
		}
#
# 13 model 4: rh_Flat with flat prior on shape (needs to use 3d version of Bayesian code)
#
		lambdad_rh_flat=c(-1/v1hat,0)
		dq=dmgs(lddi,lddd,mu1,lambdad_rh_flat,mu2,dim=2)
		rh_flat_quantiles=ml_quantiles+dq/(nx*fhat)
		if(pdf){
			dqm=dmgs(lddi,lddd,mu1m,lambdad_rh_flat,mu2m,dim=2)
			dqp=dmgs(lddi,lddd,mu1p,lambdad_rh_flat,mu2p,dim=2)
			quantilesm=ml_quantilesm+dqm/(nx*fhatm)
			quantilesp=ml_quantilesp+dqp/(nx*fhatp)
			ml_pdf=fhat
			rh_flat_pdf=-(alphap-alpham)/(quantilesp-quantilesm)
		} else{
			ml_pdf=fhat
			rh_flat_pdf="pdf not selected"
		}
#
# 15 model 6: Laplace's method
#
		if(extramodels|means){
			lambdad_lp=c(0,0)
			lddd_lp=array(0,c(2,2,2))
			dq=dmgs(lddi,lddd_lp,mu1,lambdad_lp,mu2,dim=2)
			lp_quantiles=ml_quantiles+dq/(nx*fhat)
		} else {
			lp_quantiles="extramodels not selected"
		}
#
# 16 model 7: Laplace's method, but with 1/sigma prior
#
		if(extramodels|means){
			lambdad_lp2=c(-1/v1hat,0)
			lddd_lp2=array(0,c(2,2,2))
			dq=dmgs(lddi,lddd_lp2,mu1,lambdad_lp2,mu2,dim=2)
			lp2_quantiles=ml_quantiles+dq/(nx*fhat)
		} else {
			lp2_quantiles="extramodels not selected"
		}
#
# 17 model 8: user defined xi gradient of log prior
#
		if(extramodels|means){
			lambdad_dpi=c(-1/v1hat,customprior)
			dq=dmgs(lddi,lddd,mu1,lambdad_dpi,mu2,dim=2)
			dpi_quantiles=ml_quantiles+dq/(nx*fhat)
		} else {
			dpi_quantiles="extramodels not selected"
		}
#
# 18 means
#
		means=gpd_k1_means(means,ml_params,lddi,
											lambdad_rh_flat,nx,dim=2,kloc)
		rh_flat_mean	=means$rh_flat_mean
#
# 19 waicscores
#
		waic=gpd_k1_waic(waicscores,x,v1hat,v2hat,kloc,lddi,lddd,
			lambdad_rh_flat)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 21 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rgpd_k1_cp(nrust,x,kloc,rust=TRUE,mlcp=FALSE)
			ru_quantiles=makeq(rustsim$ru_deviates,p)
		}
	} else {
		rh_flat_quantiles=ml_quantiles
		ru_quantiles=ml_quantiles
		lp_quantiles=ml_quantiles
		lp2_quantiles=ml_quantiles
		dpi_quantiles=ml_quantiles
		rh_flat_pdf=ml_pdf
		rh_flat_mean=ml_mean
	} #end of if(dmgs)

	list(	ml_params=ml_params,
				ml_value=ml_value,
#				ldd=ldd,
#				lddi=lddi,
				standard_errors=standard_errors,
				ml_quantiles=ml_quantiles,
				ml_max=ml_max,
				revert2ml=revert2ml,
				cp_quantiles=rh_flat_quantiles,
				ru_quantiles=ru_quantiles,
				lp_quantiles=lp_quantiles,
				lp2_quantiles=lp2_quantiles,
				dpi_quantiles=dpi_quantiles,
				ml_pdf=ml_pdf,
				rh_flat_pdf=rh_flat_pdf,
				maic=maic,
				waic1=waic1,
				waic2=waic2,
				ml_mean=ml_mean,
				cp_mean=rh_flat_mean,
				cp_method=crhpflat_dmgs_cpmethod())

}
#' @rdname gpd_k1_cp
#' @inheritParams man
#' @export
rgpd_k1_cp=function(n,x,kloc=0,ics=c(0,0),
	minxi=-1,maxxi=2.0,
	extramodels=FALSE,rust=FALSE,mlcp=TRUE,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),length(ics)==2,!x<0)
	stopifnot(is.finite(x),!is.na(x),length(ics)==2,!x<0)

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	ru_deviates="rust not selected"
	cp_deviates="rust not selected"

	if(mlcp){
		q=qgpd_k1_cp(x,runif(n),kloc=kloc,ics=ics,
			extramodels=extramodels)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		ru_deviates=q$ru_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tgpd_k1_cp(n,x)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			ru_deviates[i]=rgpd(1,mu=kloc,sigma=th[i,1],xi=th[i,2])
		}
	}
	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
			 cp_method=crhpflat_dmgs_cpmethod())

	return(op)

}
#' @rdname gpd_k1_cp
#' @inheritParams man
#' @export
dgpd_k1_cp=function(x,y=x,kloc=0,ics=c(0,0),
	customprior=0,
	minxi=-1,maxxi=2.0,extramodels=FALSE,
	rust=FALSE,nrust=1000,debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),length(ics)==2,!x<0,!y<0)

	ics=gpd_k1_setics(x,ics)
	opt1=optim(ics,gpd_k1_loglik,x=x,kloc=kloc,control=list(fnscale=-1))
	v1hat=opt1$par[1]
	v2hat=opt1$par[2]
	if(v2hat<=(-1)){revert2ml=TRUE}else{revert2ml=FALSE}
	ml_params=c(v1hat,v2hat)
#	gpd_k1_checkmle(ml_params,kloc,minxi,maxxi)
	dd=dgpdsub(x=x,y=y,ics=ics,kloc,customprior,
		minxi,maxxi,extramodels=extramodels)
	ru_pdf="rust not selected"

	if(rust&&(!revert2ml)){
		th=tgpd_k1_cp(nrust,x)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_pdf=ru_pdf+dgpd(y,mu=kloc,sigma=th[ir,1],xi=th[ir,2])
		}
		ru_pdf=ru_pdf/nrust
	} else {
		ru_pdf=dd$ml_pdf
	}
	op=list(
					ml_params=dd$ml_params,
					ml_pdf=dd$ml_pdf,
					revert2ml=revert2ml,
					ru_pdf=ru_pdf,
					cp_method=nopdfcdfmsg())

	return(op)
}
#' @rdname gpd_k1_cp
#' @inheritParams man
#' @export
pgpd_k1_cp=function(x,y=x,kloc=0,ics=c(0,0),customprior=0,
	minxi=-1,maxxi=2.0,extramodels=FALSE,
	rust=FALSE,nrust=1000,debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),length(ics)==2,!x<0,!y<0)

	ics=gpd_k1_setics(x,ics)
	opt1=optim(ics,gpd_k1_loglik,x=x,kloc=kloc,control=list(fnscale=-1))
	v1hat=opt1$par[1]
	v2hat=opt1$par[2]
	if(v2hat<=(-1)){revert2ml=TRUE}else{revert2ml=FALSE}
	ml_params=c(v1hat,v2hat)
#	gpd_k1_checkmle(ml_params,kloc,minxi,maxxi)
	dd=dgpdsub(x=x,y=y,ics=ics,kloc,customprior,
		minxi,maxxi,extramodels=extramodels)
	ru_cdf="rust not selected"

	if(rust&&(!revert2ml)){
		th=tgpd_k1_cp(nrust,x)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_cdf=ru_cdf+pgpd(y,mu=kloc,sigma=th[ir,1],xi=th[ir,2])
		}
		ru_cdf=ru_cdf/nrust
	} else {
		ru_pdf=dd$ml_pdf
	}

	op=list(
					ml_params=dd$ml_params,
					ml_cdf=dd$ml_cdf,
					revert2ml=revert2ml,
					ru_cdf=ru_cdf,
					cp_method=nopdfcdfmsg())
	return(op)
}
#' @rdname gpd_k1_cp
#' @inheritParams man
#' @export
tgpd_k1_cp=function(n,x,kloc=0,ics=c(0,0),extramodels=FALSE,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),length(ics)==2,!x<0)
	stopifnot(is.finite(x),!is.na(x),length(ics)==2,!x<0)

	ics=gpd_k1_setics(x,ics)
	t=ru(gpd_k1_logf,x=x,kloc=kloc,n=n,d=2,init=ics)

	list(theta_samples=t$sim_vals)

}
