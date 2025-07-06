#' Generalized Extreme Value Distribution with Three Predictors, Predictions based on a Calibrating Prior
#'
#' @inherit man description author references seealso return
#' @inheritParams man
#'
#' @inheritSection man Optional Return Values
#' @inheritSection man Optional Return Values (EVT models only)
# #' @inheritSection man Details (homogeneous models)
#' @inheritSection man Details (non-homogeneous models)
# #' @inheritSection man Details (analytic integration)
#' @inheritSection man Details (DMGS integration)
#' @inheritSection man Details (RUST)
#'
#' @section Details of the Model:
#' The GEV distribution with three predictors has distribution function
#' \deqn{F(x;a_1,b_1,a_2,b_2,a_3,b_3)=\exp{(-t(x;\mu(a_1,b_1),\sigma(a_2,b_2),\xi(a_3,b_3))}}
#' where
#' \deqn{t(x;\mu(a_1,b_1),\sigma(a_2,b_2),\xi(a_3,b_3)) =
#'    \begin{cases}
#'      {\left[1+\xi(a_3,b_3)\left(\frac{x-\mu(a_1,b_1)}{\sigma(a_2,b_2)}\right)\right]}
#'      ^{-1/\xi(a_3,b_3)} & \text{if $\xi(a_3,b_3) \ne 0$}\\
#'      \exp{\left(-\frac{x-\mu(a_1,b_1)}{\sigma(a_2,b_2)}\right)} & \text{if $\xi(a_3,b_3)=0$}
#'    \end{cases}}
#' where
#' \eqn{x} is the random variable,
#' \eqn{\mu=a_1+b_1t_1} is the location parameter,
#' modelled as a function of parameters \eqn{a_1,b_1} and predictor \eqn{t_1},
#' \eqn{\sigma=e^{a_2+b_2t_2}} is the scale parameter,
#' modelled as a function of parameters \eqn{a_2,b_2} and predictor \eqn{t_2},
#' and
#' \eqn{\xi=a_3+b_3t_3} is the shape parameter,
#' modelled as a function of parameters \eqn{a_3,b_3} and predictor \eqn{t_3}.
#'
#' The calibrating prior we use is given by
#' \deqn{\pi(a_1,b_1,a_2,b_2,a_3,b_3) \propto 1}
#' as given in Jewson et al. (2025).
#'
#' The code will stop with an error if the
#' input data gives a maximum likelihood
#' value for the shape parameter that lies outside the range \code{(minxi,maxxi)},
#' since outside this range there may be numerical problems.
#' Such values seldom occur
#' in real observed data for maxima.
#'
#' @example man/examples/example_152_gev_p123.R
#'
#' @name gev_p123_cp
NULL
#' @rdname gev_p123_cp
#' @inheritParams man
#' @export
#'

qgev_p123_cp=function(x,t1,t2,t3,t01=NA,t02=NA,t03=NA,n01=NA,n02=NA,n03=NA,p=seq(0.1,0.9,0.1),ics=c(0,0,0,0,0,0),
	fdalpha=0.01,
	minxi=-1,maxxi=1,
	means=FALSE,waicscores=FALSE,extramodels=FALSE,
	pdf=FALSE,dmgs=TRUE,rust=FALSE,nrust=100000,
	centering=TRUE,debug=FALSE){

	stopifnot(	is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1,
							length(t1)==length(x),length(t2)==length(x),length(t3)==length(x),
							is.finite(t1),is.finite(t2),is.finite(t3),
							!is.na(t1),!is.na(t2),!is.na(t3),
							length(ics)==6)
#
# 1 intro
#
	alpha=1-p
	nx=length(x)
	nalpha=length(alpha)
	t01=maket0(t01,n01,t1)
	t02=maket0(t02,n02,t2)
	t03=maket0(t03,n03,t3)
	if(pdf){
		dalpha=pmin(fdalpha*alpha,fdalpha*(1-alpha))
		alpham=alpha-dalpha
		alphap=alpha+dalpha
	}
#
# 2 centering
#
  if(centering){
	  meant1=mean(t1)
	  meant2=mean(t2)
	  meant3=mean(t3)
  	t1=t1-meant1
  	t2=t2-meant2
  	t3=t3-meant3
  	t01=t01-meant1
  	t02=t02-meant2
  	t03=t03-meant3
		if(debug)message(" meant1=",meant1)
   }

#
# 3 ml param estimate
#
	if(debug)message(" ml param estimate")
	ics=gev_p123_setics(x,t1,t2,t3,ics)
	opt1=optim(ics,gev_p123_loglik,x=x,t1=t1,t2=t2,t3=t3,control=list(fnscale=-1))
	v1h=opt1$par[1]
	v2h=opt1$par[2]
	v3h=opt1$par[3]
	v4h=opt1$par[4]
	v5h=opt1$par[5]
	v6h=opt1$par[6]
	ml_params=c(v1h,v2h,v3h,v4h,v5h,v6h)
#	gev_p123_checkmle(ml_params,minxi,maxxi,t1,t2,t3)
# if at any point xi goes below -1, then revert for the whole thing
	revert2ml=calc_revert2ml(v5h,v6h,t3)
	if(debug)message(" ml_params=",ml_params)
#
# 4 predictordata
#
	prd=gev_p123_predictordata(x,t1,t2,t3,t01,t02,t03,ml_params)
	predictedparameter=prd$predictedparameter
	adjustedx=prd$adjustedx
#
# 5 aic
#
	if(debug)message(" aic")
	ml_value=opt1$val
	maic=make_maic(ml_value,nparams=6)
#
# 6 calc ml quantiles and density
#
	if(debug)message(" ml_quantiles")
	ml_quantiles=qgev_p123((1-alpha),t01,t02,t03,ymn=v1h,slope=v2h,sigma1=v3h,sigma2=v4h,xi1=v5h,xi2=v6h)
	fhat=dgev_p123(ml_quantiles,t01,t02,t03,ymn=v1h,slope=v2h,sigma1=v3h,sigma2=v4h,xi1=v5h,xi2=v6h,log=FALSE)
	if(debug)message(" ml_quantiles=",ml_quantiles)
	if(debug)message(" fhat=",fhat)
#
# dmgs
#
	standard_errors="dmgs not selected"
	cp_quantiles="dmgs not selected"
	ru_quantiles="dmgs not selected"
	ml_pdf="dmgs not selected"
	cp_pdf="dmgs not selected"
	rh_flat_pdf="dmgs not selected"
	waic1="dmgs not selected"
	waic2="dmgs not selected"
	ml_mean="dmgs not selected"
	cp_mean="dmgs not selected"
	rh_flat_mean="dmgs not selected"
	if((dmgs)&&(!revert2ml)){
#
# 7 alpha pdf stuff
#
		if(pdf){
			ml_quantilesm=qgev_p123((1-alpham),t01,t02,t03,ymn=v1h,slope=v2h,sigma1=v3h,sigma2=v4h,xi1=v5h,xi2=v6h)
			ml_quantilesp=qgev_p123((1-alphap),t01,t02,t03,ymn=v1h,slope=v2h,sigma1=v3h,sigma2=v4h,xi1=v5h,xi2=v6h)
			fhatm=dgev_p123(ml_quantilesm,t01,t02,t03,ymn=v1h,slope=v2h,sigma1=v3h,sigma2=v4h,xi1=v5h,xi2=v6h,log=FALSE)
			fhatp=dgev_p123(ml_quantilesp,t01,t02,t03,ymn=v1h,slope=v2h,sigma1=v3h,sigma2=v4h,xi1=v5h,xi2=v6h,log=FALSE)
		}
#
# 8 ldd (two versions)
#
		ldd=gev_p123_ldda(x,t1,t2,t3,v1h,v2h,v3h,v4h,v5h,v6h)

		if(debug)message(" ldd=",ldd)
		if(debug)message(" det(ldd)=",det(ldd))
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
#
# 9 (jp taken away for this one)
#

#
# 10 calculate lddd (two versions)
#
		if(debug)message(" lddd")
		lddd=gev_p123_lddda(x,t1,t2,t3,v1h,v2h,v3h,v4h,v5h,v6h)
#
# 11 mu1 (two versions)
#
		if(debug)message(" calculate mu1")

		mu1=gev_p123_mu1fa(alpha,t01,t02,t03,v1h,v2h,v3h,v4h,v5h,v6h)

		if(pdf){
			mu1m=gev_p123_mu1fa(alpham,t01,t02,t03,v1h,v2h,v3h,v4h,v5h,v6h)
			mu1p=gev_p123_mu1fa(alphap,t01,t02,t03,v1h,v2h,v3h,v4h,v5h,v6h)
		}
#
# 12 mu2 (two versions)
#
		if(debug)message(" calculate mu2")

		mu2=gev_p123_mu2fa(alpha,t01,t02,t03,v1h,v2h,v3h,v4h,v5h,v6h)

		if(extramodels|means){
		}
		if(pdf){
			if(debug)message(" alpha pdf option")
			mu2m=gev_p123_mu2fa(alpham,t01,t02,t03,v1h,v2h,v3h,v4h,v5h,v6h)
			mu2p=gev_p123_mu2fa(alphap,t01,t02,t03,v1h,v2h,v3h,v4h,v5h,v6h)
		}
#
# 13 model 1: cp=flat prior
#
		if(debug)message(" lddi=",lddi)
		if(debug)message(" lddd=",lddd)
		if(debug)message(" mu1=",mu1)
		if(debug)message(" mu2=",mu2)
		lambdad_cp=matrix(0,6)
		dq=dmgs(lddi,lddd,mu1,lambdad_cp,mu2,dim=6)
		if(debug)message(" dq=",dq)
		rh_flat_quantiles=ml_quantiles+dq/(nx*fhat)
		if(debug)message(" rh_flat_quantiles=",rh_flat_quantiles)
#
# 14 model 4: rh_Flat with flat prior on shape (needs to use 4d version of Bayesian code)
#
		if(extramodels|means){
			lambdad_crhp_mle=matrix(0,5)
			dq=dmgs(lddi,lddd,mu1,lambdad_crhp_mle,mu2,dim=5)
			crhp_mle_quantiles=ml_quantiles+dq/(nx*fhat)
		} else {
			crhp_mle_quantiles="extramodels not selected"
		}
#
# 15 alpha pdf
#
		if(pdf){
			lambdad_rh_flat=matrix(0,6)
			dqm=dmgs(lddi,lddd,mu1m,lambdad_rh_flat,mu2m,dim=6)
			dqp=dmgs(lddi,lddd,mu1p,lambdad_rh_flat,mu2p,dim=6)
			quantilesm=ml_quantilesm+dqm/(nx*fhatm)
			quantilesp=ml_quantilesp+dqp/(nx*fhatp)
			ml_pdf=fhat
			rh_flat_pdf=-(alphap-alpham)/(quantilesp-quantilesm)
		} else{
			ml_pdf=fhat
			rh_flat_pdf="pdf not selected"
		}
#
# 16 means
#
		means=gev_p123_means(means,t01,t02,t03,ml_params,nx)
		ml_mean				=means$ml_mean
		rh_flat_mean	=means$cp_mean
		crhp_mle_mean	=means$flat_mean
#
# 17 waicscores
#
		waic=gev_p123_waic(waicscores,x,t1,t2,t3,v1h,v2h,v3h,v4h,v5h,v6h,
			lddi,lddd,lambdad_cp)
		waic1=waic$waic1
		waic2=waic$waic2

		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rgev_p123_cp(nrust,x,t1=t1,t2=t2,t3=t3,t01=t01,t02=t02,t03=t03,rust=TRUE,mlcp=FALSE,debug=debug)
			ru_quantiles=makeq(rustsim$ru_deviates,p)
		}
	} else {
		rh_flat_quantiles=ml_quantiles
		ru_quantiles=ml_quantiles
		rh_flat_pdf=ml_pdf
		rh_flat_mean=ml_mean

	} #end of if(dmgs)
#
# 19 decentering
#
  if(centering){
  	if(debug)message(" qgev:ml_params,meant1=",ml_params,meant1)
    ml_params[1]=ml_params[1]-ml_params[2]*meant1
    ml_params[3]=ml_params[3]-ml_params[4]*meant2
   	if(debug)message("predictedparameter=",predictedparameter)
  }

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
				cp_quantiles=rh_flat_quantiles,
				ru_quantiles=ru_quantiles,
				ml_pdf=ml_pdf,
				cp_pdf=rh_flat_pdf,
				maic=maic,
				waic1=waic1,
				waic2=waic2,
				ml_mean=ml_mean,
				cp_mean=rh_flat_mean,
				cp_method=crhpflat_dmgs_cpmethod())

}
#' @rdname gev_p123_cp
#' @inheritParams man
#' @export
rgev_p123_cp=function(n,x,t1,t2,t3,t01=NA,t02=NA,t03=NA,n01=NA,n02=NA,n03=NA,ics=c(0,0,0,0,0,0),
		minxi=-1,maxxi=1,
		extramodels=FALSE,rust=FALSE,mlcp=TRUE,centering=TRUE,
	debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),
							length(t1)==length(x),length(t2)==length(x),length(t3)==length(x),
							is.finite(t1),is.finite(t2),is.finite(t3),
							!is.na(t1),!is.na(t2),!is.na(t3),
						length(ics)==6)

	t01=maket0(t01,n01,t1)
	t02=maket0(t02,n02,t2)
	t03=maket0(t03,n03,t3)

#
# centering
#
	if(centering){
	  meant1=mean(t1)
	  meant2=mean(t2)
	  meant3=mean(t3)
  	t1=t1-meant1
  	t2=t2-meant2
  	t3=t3-meant3
  	t01=t01-meant1
  	t02=t02-meant2
  	t03=t03-meant3
	}

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qgev_p123_cp(x,t1=t1,t2=t2,t3=t3,t01=t01,t02=t02,t03=t03,
			n01=NA,n02=NA,n03=NA,p=runif(n),ics=ics,
			extramodels=extramodels)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		ru_deviates=q$ru_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tgev_p123_cp(n,x,t1,t2,t3)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			mu=th[i,1]+t01*th[i,2]
			sigma=exp(th[i,3]+t02*th[i,4])
			xi=th[i,5]+t03*th[i,6]
			ru_deviates[i]=rgev(1,mu=mu,sigma=sigma,xi=xi)
		}
	}

#
# decentering
#
 	if(debug)message(" rgev:ml_params,meant1=",ml_params,meant1)
	if(mlcp&centering){
	  ml_params[1]=ml_params[1]-ml_params[2]*meant1
		ml_params[3]=ml_params[3]-ml_params[4]*meant2
		ml_params[5]=ml_params[5]-ml_params[6]*meant3
	}

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
				cp_method=crhpflat_dmgs_cpmethod())

	return(op)

}
#' @rdname gev_p123_cp
#' @inheritParams man
#' @export
dgev_p123_cp=function(x,t1,t2,t3,t01=NA,t02=NA,t03=NA,n01=NA,n02=NA,n03=NA,y=x,ics=c(0,0,0,0,0,0),
	minxi=-1,maxxi=1,extramodels=FALSE,
	rust=FALSE,nrust=10,centering=TRUE,debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
							length(t1)==length(x),length(t2)==length(x),length(t3)==length(x),
							is.finite(t1),is.finite(t2),is.finite(t3),
							!is.na(t1),!is.na(t2),!is.na(t3),
		length(ics)==6)

	if(debug)message(" maket0")
	t01=maket0(t01,n01,t1)
	t02=maket0(t02,n02,t2)
	t03=maket0(t03,n03,t3)

#
# centering
#
  if(centering){
		if(debug)message(" centering")
	  meant1=mean(t1)
	  meant2=mean(t2)
	  meant3=mean(t3)
  	t1=t1-meant1
  	t2=t2-meant2
  	t3=t3-meant3
  	t01=t01-meant1
  	t02=t02-meant2
  	t03=t03-meant3
  }

	if(debug)message(" ics and optim")
	ics=gev_p123_setics(x,t1,t2,t3,ics)
	opt1=optim(ics,gev_p123_loglik,x=x,t1=t1,t2=t2,t3=t3,control=list(fnscale=-1)) #this one uses the evd routine for dgev
	v1h=opt1$par[1]
	v2h=opt1$par[2]
	v3h=opt1$par[3]
	v4h=opt1$par[4]
	v5h=opt1$par[5]
	v6h=opt1$par[6]
	revert2ml=calc_revert2ml(v5h,v6h,t3)
	ml_params=c(v1h,v2h,v3h,v4h,v5h,v6h)
#	gev_p123_checkmle(ml_params,minxi,maxxi,t1,t2,t3)
	if(debug)message(" call sub")
	dd=dgev_p123sub(x=x,t1=t1,t2=t2,t3=t3,y=y,t01=t01,t02=t02,t03=t03,ics=ics,
		extramodels,debug=debug)
	ru_pdf="rust not selected"
		ml_params=dd$ml_params

	if(rust&&(!revert2ml)){
		if(debug)message(" rust")
		th=tgev_p123_cp(nrust,x=x,t1=t1,t2=t2,t3=t3,debug=debug)$theta_samples
		if(debug)message(" tgev call done")
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+t01*th[ir,2]
			sigma=exp(th[ir,3]+t02*th[ir,4])
			xi=th[ir,5]+t03*th[ir,6]
			ru_pdf=ru_pdf+dgev(y,mu=mu,sigma=sigma,xi=xi)
		}
		ru_pdf=ru_pdf/nrust
	} else {
		ru_pdf=dd$ml_pdf
	}
#
# decentering
#
  if(centering){
		if(debug)message(" decentering")
	 	if(debug)message(" dgev:ml_params,meant1=",ml_params,meant1)
    ml_params[1]=ml_params[1]-ml_params[2]*meant1
    ml_params[3]=ml_params[3]-ml_params[4]*meant2
    ml_params[5]=ml_params[5]-ml_params[6]*meant3
  }

		op=list(
					ml_params=ml_params,
					ml_pdf=dd$ml_pdf,
					revert2ml=revert2ml,
#					cp_pdf=dd$cp_pdf,
#					ru_pdf=ru_pdf,
#					cp_method=crhpflat_dmgs_cpmethod())
					ru_pdf=ru_pdf,
					cp_method=nopdfcdfmsg())

	return(op)

}
#' @rdname gev_p123_cp
#' @inheritParams man
#' @export
pgev_p123_cp=function(x,t1,t2,t3,t01=NA,t02=NA,t03=NA,n01=NA,n02=NA,n03=NA,y=x,ics=c(0,0,0,0,0,0),
	minxi=-1,maxxi=1,extramodels=FALSE,
	rust=FALSE,nrust=1000,centering=TRUE,debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
							length(t1)==length(x),length(t2)==length(x),length(t3)==length(x),
							is.finite(t1),is.finite(t2),is.finite(t3),
							!is.na(t1),!is.na(t2),!is.na(t3),
		length(ics)==6)

	t01=maket0(t01,n01,t1)
	t02=maket0(t02,n02,t2)
	t03=maket0(t03,n03,t3)
#
# centering
#
  if(centering){
	  meant1=mean(t1)
	  meant2=mean(t2)
	  meant3=mean(t3)
  	t1=t1-meant1
  	t2=t2-meant2
  	t3=t3-meant3
  	t01=t01-meant1
  	t02=t02-meant2
  	t03=t03-meant3
  }

	ics=gev_p123_setics(x,t1,t2,t3,ics)
	opt1=optim(ics,gev_p123_loglik,x=x,t1=t1,t2=t2,t3=t3,control=list(fnscale=-1)) #this one uses the evd routine for dgev
	v1h=opt1$par[1]
	v2h=opt1$par[2]
	v3h=opt1$par[3]
	v4h=opt1$par[4]
	v5h=opt1$par[5]
	v6h=opt1$par[6]
	revert2ml=calc_revert2ml(v5h,v6h,t3)
	ml_params=c(v1h,v2h,v3h,v4h,v5h,v6h)
#	gev_p123_checkmle(ml_params,minxi,maxxi,t1,t2,t3)
	dd=dgev_p123sub(x=x,t1=t1,t2=t2,t3=t3,y=y,t01=t01,t02=t02,t03=t03,ics=ics,
		extramodels,debug=debug)
	ru_cdf="rust not selected"
		ml_params=dd$ml_params

	if(rust&&(!revert2ml)){
		th=tgev_p123_cp(nrust,x,t1,t2,t3)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+t01*th[ir,2]
			sigma=exp(th[ir,3]+t02*th[ir,4])
			xi=th[ir,5]+t03*th[ir,6]
			ru_cdf=ru_cdf+pgev(y,mu=mu,sigma=exp(sigma),xi=xi)
		}
		ru_cdf=ru_cdf/nrust
	} else {
		ru_pdf=dd$ml_pdf
	}
#
# decentering
#
  if(centering){
	 	if(debug)message(" pgev:ml_params,meant1=",ml_params,meant1)
    ml_params[1]=ml_params[1]-ml_params[2]*meant1
    ml_params[3]=ml_params[3]-ml_params[4]*meant2
    ml_params[5]=ml_params[5]-ml_params[6]*meant3
  }

	op=list(
					ml_params=ml_params,
					ml_cdf=dd$ml_cdf,
					revert2ml=revert2ml,
#					cp_cdf=dd$cp_cdf,
#					ru_cdf=ru_cdf,
#					cp_method=crhpflat_dmgs_cpmethod())
					ru_cdf=ru_cdf,
					cp_method=nopdfcdfmsg())
	return(op)
}
#' @rdname gev_p123_cp
#' @inheritParams man
#' @export
tgev_p123_cp=function(n,x,t1,t2,t3,ics=c(0,0,0,0,0,0),
	extramodels=FALSE,debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),
							length(t1)==length(x),length(t2)==length(x),length(t3)==length(x),
							is.finite(t1),is.finite(t2),is.finite(t3),
							!is.na(t1),!is.na(t2),!is.na(t3),
							length(ics)==6)

	warning("Sorry...rust doesn't work for gev_p123. I'm working on it.")
	stop()

#
# centering
#
	  meant1=mean(t1)
	  meant2=mean(t2)
	  meant3=mean(t3)
  	t1=t1-meant1
  	t2=t2-meant2
  	t3=t3-meant3

	if(debug)message("sums=",sum(x),sum(t1),sum(t2),n)
	th=ru(gev_p123_logf,x=x,t1=t1,t2=t2,t3=t3,n=n,d=6,init=c(0,0,0,0,0,0))
	if(debug)message(" back from rust")
  theta_samples=th$sim_vals
#
# decentering
#
  theta_samples[,1]=theta_samples[,1]-theta_samples[,2]*meant1
  theta_samples[,3]=theta_samples[,3]-theta_samples[,4]*meant2
  theta_samples[,5]=theta_samples[,5]-theta_samples[,6]*meant3

	list(theta_samples=theta_samples)

}
