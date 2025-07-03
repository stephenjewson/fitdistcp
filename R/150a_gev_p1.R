#' Generalized Extreme Value Distribution with a Predictor, Predictions Based on a Calibrating Prior
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
#' The GEV distribution with a predictor has distribution function
#' \deqn{F(x;a,b,\sigma,\xi)=\exp{(-t(x;\mu(a,b),\sigma,\xi))}}
#' where
#' \deqn{t(x;\mu(a,b),\sigma,\xi) =
#'    \begin{cases}
#'      {\left[1+\xi\left(\frac{x-\mu(a,b)}{\sigma}\right)\right]}^{-1/\xi} & \text{if $\xi \ne 0$}\\
#'      \exp{\left(-\frac{x-\mu(a,b)}{\sigma}\right)} & \text{if $\xi=0$}
#'    \end{cases}}
#' where
#' \eqn{x} is the random variable,
#' \eqn{\mu=a+bt} is the location parameter,
#' modelled as a function of parameters \eqn{a,b} and predictor \eqn{t},
#' and \eqn{\sigma>0,\xi} are the scale and shape parameters.
#'
#' The calibrating prior we use is given by
#' \deqn{\pi(a,b,\sigma,\xi) \propto \frac{1}{\sigma}}
#' as given in Jewson et al. (2025).
#'
#' The code will stop with an error if the
#' input data gives a maximum likelihood
#' value for the shape parameter that lies outside the range \code{(minxi,maxxi)},
#' since outside this range there may be numerical problems.
#' Such values seldom occur
#' in real observed data for maxima.
#'
#' @example man/examples/example_150_gev_p1.R
#'
#' @name gev_p1_cp
NULL
#' @rdname gev_p1_cp
#' @inheritParams man
#' @export
#'

qgev_p1_cp=function(x,t,t0=NA,n0=NA,p=seq(0.1,0.9,0.1),ics=c(0,0,0,0),
	d1=0.01,d2=0.01,fd3=0.01,d4=0.01,fdalpha=0.01,
	minxi=-0.45,maxxi=0.45,
	means=FALSE,waicscores=FALSE,extramodels=FALSE,
	pdf=FALSE,dmgs=TRUE,rust=FALSE,nrust=100000,predictordata=TRUE,
	centering=TRUE,debug=FALSE,aderivs=TRUE){

	stopifnot(	is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1,
							length(ics)==4)
#
# 1 intro
#
	alpha=1-p
	nx=length(x)
	nalpha=length(alpha)
	t0=maket0(t0,n0,t)
	if(debug)message(" t0=",t0)
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
	if(debug)message(" ml param estimate")
	ics=gev_p1_setics(x,t,ics)
	opt1=optim(ics,gev_p1_loglik,x=x,t=t,control=list(fnscale=-1)) #this one uses the evd routine for dgev
	v1hat=opt1$par[1]
	v2hat=opt1$par[2]
	v3hat=opt1$par[3]
	v4hat=opt1$par[4]
	ml_params=c(v1hat,v2hat,v3hat,v4hat)
#	gev_p1_checkmle(ml_params,minxi,maxxi)
	if(debug)message(" ml_params=",ml_params)
	if(abs(v4hat)>=1){revert2ml=TRUE}else{revert2ml=FALSE}
#
# 4 predictordata
#
	prd=gev_p1_predictordata(predictordata,x,t,t0,ml_params)
	predictedparameter=prd$predictedparameter
	adjustedx=prd$adjustedx
#
# 5 aic
#
	if(debug)message(" aic")
	ml_value=opt1$val
	maic=make_maic(ml_value,nparams=4)
#
# 6 calc ml quantiles and density
#
	if(debug)message(" ml_quantiles")
	ml_quantiles=qgev_p1((1-alpha),t0,ymn=v1hat,slope=v2hat,sigma=v3hat,xi=v4hat)
	if(v4hat<0){
		ml_max=(v1hat+v2hat*t0)-v3hat/v4hat
	} else {
		ml_max=Inf
	}
	fhat=dgev_p1(ml_quantiles,t0,ymn=v1hat,slope=v2hat,sigma=v3hat,xi=v4hat,log=FALSE)
	if(debug)message(" ml_quantiles=",ml_quantiles)
	if(debug)message(" fhat=",fhat)
#
# dmgs
#
	standard_errors="dmgs not selected"
	rh_flat_quantiles="dmgs not selected"
	flat_quantiles="dmgs not selected"
	rh_ml_quantiles="dmgs not selected"
	cp_quantiles="dmgs not selected"
	jp_quantiles="dmgs not selected"
	ru_quantiles="dmgs not selected"
	rh_flat_pdf="dmgs not selected"
	ml_pdf="dmgs not selected"
	cp_pdf="dmgs not selected"
	waic1="dmgs not selected"
	waic2="dmgs not selected"
	ml_mean="dmgs not selected"
	cp_mean="dmgs not selected"
	rh_flat_mean="dmgs not selected"
	flat_mean="dmgs not selected"
	rh_ml_mean="dmgs not selected"
	jp_mean="dmgs not selected"
	cp_method="dmgs not selected"
	if((dmgs)&&(!revert2ml)){
#
# 7 alpha pdf stuff
# -for now, I only make the pdf. I could add cdf too I suppose. Somehow.
		if(pdf){
			ml_quantilesm=qgev_p1((1-alpham),t0,ymn=v1hat,slope=v2hat,sigma=v3hat,xi=v4hat)
			ml_quantilesp=qgev_p1((1-alphap),t0,ymn=v1hat,slope=v2hat,sigma=v3hat,xi=v4hat)
			fhatm=dgev_p1(ml_quantilesm,t0,ymn=v1hat,slope=v2hat,sigma=v3hat,xi=v4hat,log=FALSE)
			fhatp=dgev_p1(ml_quantilesp,t0,ymn=v1hat,slope=v2hat,sigma=v3hat,xi=v4hat,log=FALSE)
		}
#
# 8 ldd
#
		if(aderivs) ldd=gev_p1_ldda(x,t,v1hat,v2hat,v3hat,v4hat)
		if(!aderivs)ldd=gev_p1_ldd(x,t,v1hat,d1,v2hat,d2,v3hat,fd3,v4hat,d4)
		if(debug)message(" ldd=",ldd)
		if(debug)message(" det(ldd)=",det(ldd))
		lddi=solve(ldd)

		standard_errors=make_se(nx,lddi)
		if(extramodels|means){
			if(aderivs) ldd_k4=gev_p1k3_ldda(x,t,v1hat,v2hat,v3hat,v4hat)
			if(!aderivs)ldd_k4=gev_p1k3_ldd(x,t,v1hat,d1,v2hat,d2,v3hat,fd3,v4hat)
			lddi_k4=solve(ldd_k4)
		}
#
# 9 information matrix and related (for Jeffreys prior)
# -because of difficulty of calculating expected information, I just use observed information
		if(debug)message(" call gev.infomat")
		if(extramodels|means){
			gg=-ldd
			ggi=solve(gg)
			detg=det(gg)
			ggd=gev_p1_ggd(x,t,v1hat,d1,v2hat,d2,v3hat,fd3,v4hat,d4)	#uses ldd! so not expected at all, just observed
		}
#
# 10 calculate lddd (two versions)
#
		if(debug)message(" lddd")
		if(aderivs) lddd=gev_p1_lddda(x,t,v1hat,v2hat,v3hat,v4hat)
		if(!aderivs)lddd=gev_p1_lddd(x,t,v1hat,d1,v2hat,d2,v3hat,fd3,v4hat,d4)
		if(extramodels|means){
			if(aderivs) lddd_k4=gev_p1k3_lddda(x,t,v1hat,v2hat,v3hat,v4hat)
			if(!aderivs)lddd_k4=gev_p1k3_lddd(x,t,v1hat,d1,v2hat,d2,v3hat,fd3,v4hat)
		}
#
# 11 mu1 (two versions)
#
		if(debug)message(" calculate mu1")
		if(aderivs) mu1=gev_p1_mu1fa(alpha,t0,v1hat,v2hat,v3hat,v4hat)
		if(!aderivs)mu1=gev_p1_mu1f(alpha,t0,v1hat,d1,v2hat,d2,v3hat,fd3,v4hat,d4)

		if(extramodels|means)mu1_k4=gev_p1k3_mu1f(alpha,t0,v1hat,d1,v2hat,d2,v3hat,fd3,v4hat)
		if(pdf){
			if(aderivs){
				mu1m=gev_p1_mu1fa(alpham,t0,v1hat,v2hat,v3hat,v4hat)
				mu1p=gev_p1_mu1fa(alphap,t0,v1hat,v2hat,v3hat,v4hat)
			} else{
				mu1m=gev_p1_mu1f(alpham,t0,v1hat,d1,v2hat,d2,v3hat,fd3,v4hat,d4)
				mu1p=gev_p1_mu1f(alphap,t0,v1hat,d1,v2hat,d2,v3hat,fd3,v4hat,d4)
			}
		}
#
# 12 mu2 (two versions)
#
		if(debug)message(" calculate mu2")
		if(aderivs) mu2=gev_p1_mu2fa(alpha,t0,v1hat,v2hat,v3hat,v4hat)
		if(!aderivs)mu2=gev_p1_mu2f(alpha,t0,v1hat,d1,v2hat,d2,v3hat,fd3,v4hat,d4)

		if(extramodels|means)mu2_k4=gev_p1k3_mu2f(alpha,t0,v1hat,d1,v2hat,d2,v3hat,fd3,v4hat)
		if(pdf){
			if(aderivs){
				mu2m=gev_p1_mu2fa(alpham,t0,v1hat,v2hat,v3hat,v4hat)
				mu2p=gev_p1_mu2fa(alphap,t0,v1hat,v2hat,v3hat,v4hat)
			} else {
				mu2m=gev_p1_mu2f(alpham,t0,v1hat,d1,v2hat,d2,v3hat,fd3,v4hat,d4)
				mu2p=gev_p1_mu2f(alphap,t0,v1hat,d1,v2hat,d2,v3hat,fd3,v4hat,d4)
			}
		}
#
# 13 model 2: flat prior
#
		if(extramodels|means){
			lambdad_flat=matrix(0,4)
			dq=dmgs(lddi,lddd,mu1,lambdad_flat,mu2,dim=4)
			flat_quantiles=ml_quantiles+dq/(nx*fhat)
		} else {
			flat_quantiles="extramodels not selected"
		}
#
# 14 model 3: rh_ML (needs to use 3d version of Bayesian code, and fixedshape ldd,lddd,mu1,mu2)
#
		if(extramodels|means){
			lambdad_rh_mle=c(0,0,-1/v3hat)
			dq=dmgs(lddi_k4,lddd_k4,mu1_k4,lambdad_rh_mle,mu2_k4,dim=3)
			rh_ml_quantiles=ml_quantiles+dq/(nx*fhat)
		} else {
			rh_ml_quantiles="extramodels not selected"
		}
#
# 15 model 4: rh_Flat with flat prior on shape (needs to use 4d version of Bayesian code)
#
		lambdad_rh_flat=c(0,0,-1/v3hat,0)
		dq=dmgs(lddi,lddd,mu1,lambdad_rh_flat,mu2,dim=4)
		rh_flat_quantiles=ml_quantiles+dq/(nx*fhat)
		if(pdf){
			dqm=dmgs(lddi,lddd,mu1m,lambdad_rh_flat,mu2m,dim=4)
			dqp=dmgs(lddi,lddd,mu1p,lambdad_rh_flat,mu2p,dim=4)
			quantilesm=ml_quantilesm+dqm/(nx*fhatm)
			quantilesp=ml_quantilesp+dqp/(nx*fhatp)
			ml_pdf=fhat
			rh_flat_pdf=-(alphap-alpham)/(quantilesp-quantilesm)
		} else{
			ml_pdf=fhat
			rh_flat_pdf="pdf not selected"
		}
#
# 16 model 5:
#
		if(extramodels|means){
			lambdad_jp=jpf4p(ggd,detg,ggi) #this is jp
			dq=dmgs(lddi,lddd,mu1,lambdad_jp,mu2,dim=4)
			jp_quantiles=ml_quantiles+dq/(nx*fhat)
		} else {
			jp_quantiles="extramodels not selected"
		}
#
# 17 means
#
		means=gev_p1_means(means,t0,ml_params,lddi,lddi_k4,lddd,lddd_k4,
											lambdad_flat,lambdad_rh_mle,lambdad_rh_flat,lambdad_jp,nx,dim=4)
		ml_mean				=means$ml_mean
		flat_mean			=means$flat_mean
		rh_ml_mean		=means$rh_ml_mean
		rh_flat_mean	=means$rh_flat_mean
		jp_mean				=means$jp_mean
#
# 18 waicscores
#
		waic=gev_p1_waic(waicscores,x,t,v1hat,d1,v2hat,d2,v3hat,fd3,v4hat,d4,
			lddi,lddd,lambdad_rh_flat,aderivs)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 20 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rgev_p1_cp(nrust,x,t=t,t0=t0,rust=TRUE,mlcp=FALSE)
			ru_quantiles=makeq(rustsim$ru_deviates,p)
		}
	} else {
		flat_quantiles=ml_quantiles
		rh_ml_quantiles=ml_quantiles
		rh_flat_quantiles=ml_quantiles
		ru_quantiles=ml_quantiles
		jp_quantiles=ml_quantiles
		rh_flat_pdf=ml_pdf
		flat_mean=ml_mean
		rh_ml_mean=ml_mean
		rh_flat_mean=ml_mean
		jp_mean=ml_mean
	} #end of if(dmgs)
#
# 21 decentering
#
  if(centering){
    ml_params[1]=ml_params[1]-ml_params[2]*meant
    if(predictordata)predictedparameter=predictedparameter-ml_params[2]*meant
  }

	list(	ml_params=ml_params,
				ml_value=ml_value,
				predictedparameter=predictedparameter,
				adjustedx=adjustedx,
#				expinfmat=expinfmat,
#				expinfmati=expinfmati,
				standard_errors=standard_errors,
				revert2ml=revert2ml,
				ml_quantiles=ml_quantiles,
				ml_max=ml_max,
				flat_quantiles=flat_quantiles,
				rh_ml_quantiles=rh_ml_quantiles,
				cp_quantiles=rh_flat_quantiles,
				ru_quantiles=ru_quantiles,
				jp_quantiles=jp_quantiles,
				ml_pdf=ml_pdf,
				cp_pdf=rh_flat_pdf,
				maic=maic,
				waic1=waic1,
				waic2=waic2,
				ml_mean=ml_mean,
				flat_mean=flat_mean,
				rh_ml_mean=rh_ml_mean,
				cp_mean=rh_flat_mean,
				jp_mean=jp_mean,
				cp_method=crhpflat_dmgs_cpmethod())

}
#' @rdname gev_p1_cp
#' @inheritParams man
#' @export
rgev_p1_cp=function(n,x,t,t0=NA,n0=NA,ics=c(0,0,0,0),
	d1=0.01,d2=0.01,fd3=0.01,d4=0.01,
		minxi=-0.45,maxxi=0.45,
		extramodels=FALSE,rust=FALSE,mlcp=TRUE,debug=FALSE,aderivs=TRUE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),is.finite(t),!is.na(t),
#						length(ics)==4)
	stopifnot(is.finite(x),!is.na(x),is.finite(t),!is.na(t),
						length(ics)==4)

	t0=maket0(t0=t0,n0=n0,t=t)

#
# centering
#
  meant=mean(t)
  t=t-meant
  t0=t0-meant

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	flat_deviates="mlcp not selected"
	rh_ml_deviates="mlcp not selected"
	jp_deviates="mlcp not selected"
	ru_deviates="rust not selected"
	cp_deviates="rust not selected"

	if(mlcp){
		q=qgev_p1_cp(x,t=t,t0=t0,n0=NA,p=runif(n),ics=ics,d1=d1,d2=d2,fd3=fd3,d4=d4,
			extramodels=extramodels,aderivs=aderivs)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		flat_deviates=q$flat_quantiles
		rh_ml_deviates=q$rh_ml_quantiles
		ru_deviates=q$ru_quantiles
		cp_deviates=q$cp_quantiles
		jp_deviates=q$jp_quantiles
	}

	if(rust){
		th=tgev_p1_cp(n,x,t)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			mu=th[i,1]+t0*th[i,2]
			ru_deviates[i]=rgev(1,mu=mu,sigma=th[i,3],xi=th[i,4])
		}
	}

#
# decentering
#
  if(mlcp)ml_params[1]=ml_params[1]-ml_params[2]*meant

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 flat_deviates=flat_deviates,
			 rh_ml_deviates=rh_ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
			 jp_deviates=jp_deviates,
				cp_method=crhpflat_dmgs_cpmethod())

	return(op)

}
#' @rdname gev_p1_cp
#' @inheritParams man
#' @export
dgev_p1_cp=function(x,t,t0=NA,n0=NA,y=x,ics=c(0,0,0,0),
	d1=0.01,d2=0.01,fd3=0.01,d4=0.01,
	minxi=-0.45,maxxi=0.45,extramodels=FALSE,
	rust=FALSE,nrust=1000,centering=TRUE,debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
						is.finite(t),!is.na(t),length(ics)==4)

	t0=maket0(t0,n0,t)

#
# centering
#
  if(centering){
    meant=mean(t)
    t=t-meant
    t0=t0-meant
  }

	ics=gev_p1_setics(x,t,ics)
	opt1=optim(ics,gev_p1_loglik,x=x,t=t,control=list(fnscale=-1)) #this one uses the evd routine for dgev
	v1hat=opt1$par[1]
	v2hat=opt1$par[2]
	v3hat=opt1$par[3]
	v4hat=opt1$par[4]
	if(v4hat<=(-1)){revert2ml=TRUE}else{revert2ml=FALSE}
	ml_params=c(v1hat,v2hat,v3hat,v4hat)
#	gev_p1_checkmle(ml_params,minxi,maxxi)
	dd=dgev_p1sub(x=x,t=t,y=y,t0=t0,ics=ics,d1,d2,fd3,d4,
		minxi,maxxi,extramodels=extramodels,aderivs=aderivs)
	ru_pdf="rust not selected"
		ml_params=dd$ml_params

	if(rust&&(!revert2ml)){
		th=tgev_p1_cp(nrust,x,t)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+t0*th[ir,2]
			ru_pdf=ru_pdf+dgev(y,mu=mu,sigma=th[ir,3],xi=th[ir,4])
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
#' @rdname gev_p1_cp
#' @inheritParams man
#' @export
pgev_p1_cp=function(x,t,t0=NA,n0=NA,y=x,ics=c(0,0,0,0),
	d1=0.01,d2=0.01,fd3=0.01,d4=0.01,
	minxi=-0.45,maxxi=0.45,extramodels=FALSE,
	rust=FALSE,nrust=1000,centering=TRUE,debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
						is.finite(t),!is.na(t),length(ics)==4)

	t0=maket0(t0,n0,t)
#
# centering
#
  if(centering){
    meant=mean(t)
    t=t-meant
    t0=t0-meant
  }

	ics=gev_p1_setics(x,t,ics)
	opt1=optim(ics,gev_p1_loglik,x=x,t=t,control=list(fnscale=-1)) #this one uses the evd routine for dgev
	v1hat=opt1$par[1]
	v2hat=opt1$par[2]
	v3hat=opt1$par[3]
	v4hat=opt1$par[4]
	if(v4hat<=(-1)){revert2ml=TRUE}else{revert2ml=FALSE}
	ml_params=c(v1hat,v2hat,v3hat,v4hat)
#	gev_p1_checkmle(ml_params,minxi,maxxi)
	dd=dgev_p1sub(x=x,t=t,y=y,t0=t0,ics=ics,d1,d2,fd3,d4,
		minxi,maxxi,extramodels=extramodels,aderivs=aderivs)
	ru_cdf="rust not selected"
		ml_params=dd$ml_params

	if(rust&&(!revert2ml)){
		th=tgev_p1_cp(nrust,x,t)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+t0*th[ir,2]
			ru_cdf=ru_cdf+pgev(y,mu=mu,sigma=th[ir,3],xi=th[ir,4])
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
#' @rdname gev_p1_cp
#' @inheritParams man
#' @export
tgev_p1_cp=function(n,x,t,ics=c(0,0,0,0),
	d1=0.01,d2=0.01,fd3=0.01,d4=0.01,
	extramodels=FALSE,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),is.finite(t),!is.na(t),
#						length(ics)==4)
	stopifnot(is.finite(x),!is.na(x),is.finite(t),!is.na(t),
						length(ics)==4)

#
# centering
#
  meant=mean(t)
  t=t-meant

	ics=gev_p1_setics(x,t,ics)
#	th=ru(gev_p1_logf,x=x,t=t,n=n,d=4,init=c(0,0,1,0))
	th=ru(gev_p1_logf,x=x,t=t,n=n,d=4,init=ics)
  theta_samples=th$sim_vals
#
# decentering
#
  theta_samples[,1]=theta_samples[,1]-theta_samples[,2]*meant

	list(theta_samples=theta_samples)

}
