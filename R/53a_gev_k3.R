#' Generalized Extreme Value Distribution with Known Shape, Predictions Based on a Calibrating Prior
#'
#' @inherit man description author references seealso
#' @inheritParams man
#'
#' @inheritSection man Default Return Values
#' @inheritSection man Optional Return Values
#' @inheritSection man Optional Return Values (EVT models only)
# #' @inheritSection man Optional Return Values (non-RHP models only)
#' @inheritSection man Details (homogeneous models)
# #' @inheritSection man Details (non-homogeneous models)
# #' @inheritSection man Details (analytic integration)
#' @inheritSection man Details (DMGS integration)
#' @inheritSection man Details (RUST)
#'
#' @section Details of the Model:
#' The GEV distribution with known shape has distribution function
#' \deqn{F(x;\mu,\sigma)=\exp{(-t(x;\mu,\sigma))}}
#' where
#' \deqn{t(x;\mu,\sigma) =
#'    \begin{cases}
#'      {\left[1+\xi\left(\frac{x-\mu}{\sigma}\right)\right]}^{-1/\xi} & \text{if $\xi \ne 0$}\\
#'      \exp{(-\frac{x-\mu}{\sigma})} & \text{if $\xi=0$}
#'    \end{cases}}
#' where
#' \eqn{x} is the random variable,
#' \eqn{\mu,\sigma>0} are the parameters
#' and \eqn{\xi} is known (hence the \code{k3} in the name).
#'
#' The calibrating prior we use is given by
#' \deqn{\pi(\mu,\sigma) \propto \frac{1}{\sigma}}
#' as given in Jewson et al. (2024).
#'
#' @example man/examples/example_53_gev_k3.R
#'
#' @name gev_k3_cp
NULL
#' @rdname gev_k3_cp
#' @inheritParams man
#' @export
#'
qgev_k3_cp=function(x,p=seq(0.1,0.9,0.1),d1=0.01,fd2=0.01,fdalpha=0.01,
	kshape=0,means=FALSE,waicscores=FALSE,pdf=FALSE,
	dmgs=TRUE,rust=FALSE,nrust=100000,
	debug=FALSE,aderivs=TRUE){
#
# 1 intro
#
#	stopifnot(is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1)
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
	if(debug)cat("2 calc ml param estimate")
	v1start=mean(x)
	v2start=sd(x)
	opt=optim(c(v1start,v2start),gev_k3_loglik,x=x,kshape=kshape,
		control=list(fnscale=-1))
	v1hat=opt$par[1]
	v2hat=opt$par[2]
	ml_params=c(v1hat,v2hat)
	if(debug)cat("  v1hat,v2hat=",v1hat,v2hat,"//")
#
# 3 aic
#
	ml_value=opt$val
	maic=make_maic(ml_value,nparams=2)
#
# 4 ml quantiles (vectorized over alpha)
#
	ml_quantiles=qgev((1-alpha),mu=v1hat,sigma=v2hat,xi=kshape)
	fhat=dgev(ml_quantiles,mu=v1hat,sigma=v2hat,xi=kshape)
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
	if(dmgs){
		if(pdf){
			ml_quantilesm=qgev((1-alpham),mu=v1hat,sigma=v2hat,xi=kshape)
			ml_quantilesp=qgev((1-alphap),mu=v1hat,sigma=v2hat,xi=kshape)
			fhatm=dgev(ml_quantilesm,mu=v1hat,sigma=v2hat,xi=kshape)
			fhatp=dgev(ml_quantilesp,mu=v1hat,sigma=v2hat,xi=kshape)
		}
#
# 5 lddi
#
		if(debug)cat("  calculate ldd,lddi\n")
		if(aderivs)	ldd=gev_k3_ldda(x,v1hat,v2hat,kshape)
		if(!aderivs)ldd=gev_k3_ldd(x,v1hat,d1,v2hat,fd2,kshape)
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
#
# 6 lddd
#
		if(debug)cat("  calculate lddd\n")
		if(aderivs)	lddd=gev_k3_lddda(x,v1hat,v2hat,kshape)
		if(!aderivs)lddd=gev_k3_lddd(x,v1hat,d1,v2hat,fd2,kshape)
#
# 7 mu1
#
		if(debug)cat("  calculate mu1\n")
		if(aderivs) mu1=gev_k3_mu1fa(alpha,v1hat,v2hat,kshape)
		if(!aderivs)mu1=gev_k3_mu1f(alpha,v1hat,d1,v2hat,fd2,kshape)

		if(pdf){
			if(aderivs){
				mu1m=gev_k3_mu1fa(alpham,v1hat,v2hat,kshape)
				mu1p=gev_k3_mu1fa(alphap,v1hat,v2hat,kshape)
			} else {
				mu1m=gev_k3_mu1f(alpham,v1hat,d1,v2hat,fd2,kshape)
				mu1p=gev_k3_mu1f(alphap,v1hat,d1,v2hat,fd2,kshape)
			}
		}
#
# 8 mu2
#
		if(debug)cat("  calculate mu2\n")
		if(aderivs) mu2=gev_k3_mu2fa(alpha,v1hat,v2hat,kshape)
		if(!aderivs)mu2=gev_k3_mu2f(alpha,v1hat,d1,v2hat,fd2,kshape)

		if(pdf){
			if(aderivs){
				mu2m=gev_k3_mu2fa(alpham,v1hat,v2hat,kshape)
				mu2p=gev_k3_mu2fa(alphap,v1hat,v2hat,kshape)
			} else {
				mu2m=gev_k3_mu2f(alpham,v1hat,d1,v2hat,fd2,kshape)
				mu2p=gev_k3_mu2f(alphap,v1hat,d1,v2hat,fd2,kshape)
			}
		}
#
# 9 rhp
#
		lambdad_rhp=c(0,-1/v2hat)
#
# 10 fhat, dq and quantiles
#
		if(debug)cat("  fhat, dq and quantiles\n")
		dq=dmgs(lddi,lddd,mu1,lambdad_rhp,mu2,dim=2)
		rh_quantiles=ml_quantiles+dq/(nx*fhat)
		if(pdf){
			dqm=dmgs(lddi,lddd,mu1m,lambdad_rhp,mu2m,dim=2)
			dqp=dmgs(lddi,lddd,mu1p,lambdad_rhp,mu2p,dim=2)
			quantilesm=ml_quantilesm+dqm/(nx*fhatm)
			quantilesp=ml_quantilesp+dqp/(nx*fhatp)
			ml_pdf=fhat
			rh_pdf=-(alphap-alpham)/(quantilesp-quantilesm)
		} else{
			ml_pdf=fhat
			rh_pdf="pdf not selected"
		}
#
# 11 means
#
		means=gev_k3_means(means,ml_params,lddi,lddd,lambdad_rhp,nx,dim=2,kshape=kshape)
		ml_mean=means$ml_mean
		rh_mean=means$rh_mean
#
# 12 waicscores
#
		waic=gev_k3_waic(waicscores,x,v1hat,d1,v2hat,fd2,kshape,lddi,lddd,
			lambdad_rhp,aderivs)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 14 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rgev_k3_cp(n=nrust,x=x,kshape=kshape,rust=TRUE,mlcp=FALSE)
			ru_quantiles=makeq(rustsim$ru_deviates,p)
		}
	} #end of if(dmgs)

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
				ml_pdf=ml_pdf,
				cp_pdf=rh_pdf,
				maic=maic,
				waic1=waic1,
				waic2=waic2,
				ml_mean=ml_mean,
				cp_mean=rh_mean,
				cp_method=rhp_dmgs_cpmethod())

}
#' @rdname gev_k3_cp
#' @inheritParams man
#' @export
rgev_k3_cp=function(n,x,d1=0.01,fd2=0.01,kshape=0,rust=FALSE,mlcp=TRUE,
	debug=FALSE,aderivs=TRUE){

# this next line was creating the crazy error on install and I don't know
#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x))
#	stopifnot(is.finite(x),!is.na(x))

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qgev_k3_cp(x,runif(n),d1=d1,fd2=fd2,kshape=kshape,aderivs=aderivs)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tgev_k3_cp(n,x)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			ru_deviates[i]=rgev(1,mu=th[i,1],sigma=th[i,2],xi=kshape)
		}
	}

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
				cp_method=rhp_dmgs_cpmethod())

	return(op)

}
#' @rdname gev_k3_cp
#' @inheritParams man
#' @export
dgev_k3_cp=function(x,y=x,d1=0.01,fd2=0.01,kshape=0,rust=FALSE,nrust=1000,
	debug=FALSE,aderivs=TRUE){

#	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y))

	dd=dgev_k3sub(x=x,y=y,d1,fd2,kshape=kshape,aderivs=aderivs)
	ru_pdf="rust not selected"

	if(rust){
		th=tgev_k3_cp(nrust,x,kshape)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_pdf=ru_pdf+dgev(y,mu=th[ir,1],sigma=th[ir,2],xi=kshape)
		}
		ru_pdf=ru_pdf/nrust
	}
		op=list(
				ml_params=dd$ml_params,
				ml_pdf=dd$ml_pdf,
				ru_pdf=ru_pdf,
				cp_method=nopdfcdfmsg())
	return(op)

}
#' @rdname gev_k3_cp
#' @inheritParams man
#' @export
pgev_k3_cp=function(x,y=x,d1=0.01,fd2=0.01,kshape=0,rust=FALSE,nrust=1000,
	debug=FALSE,aderivs=TRUE){

#	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y))

	dd=dgev_k3sub(x=x,y=y,d1,fd2,kshape=kshape,aderivs=aderivs)
	ru_cdf="rust not selected"

	if(rust){
		th=tgev_k3_cp(nrust,x,kshape)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			ru_cdf=ru_cdf+pgev(y,mu=th[ir,1],sigma=th[ir,2],xi=kshape)
		}
		ru_cdf=ru_cdf/nrust
	}
		op=list(
				ml_params=dd$ml_params,
				ml_cdf=dd$ml_cdf,
				ru_cdf=ru_cdf,
				cp_method=nopdfcdfmsg())
	return(op)
}
#' @rdname gev_k3_cp
#' @inheritParams man
#' @export
tgev_k3_cp=function(n,x,d1=0.01,fd2=0.01,kshape=0,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x))
#	stopifnot(is.finite(x),!is.na(x))

	t=ru(gev_k3_logf,x=x,kshape=kshape,n=n,d=2,init=c(mean(x),sd(x)))

	list(theta_samples=t$sim_vals)

}
