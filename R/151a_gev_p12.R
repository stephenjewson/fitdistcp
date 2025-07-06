#' Generalized Extreme Value Distribution with Two Predictors, Predictions based on a Calibrating Prior
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
#' The GEV distribution with two predictors has distribution function
#' \deqn{F(x;a_1,b_1,a_2,b_2,\xi)=\exp{(-t(x;\mu(a_1,b_1),\sigma(a_2,b_2),\xi))}}
#' where
#' \deqn{t(x;\mu(a_1,b_1),\sigma(a_2,b_2),\xi) =
#'    \begin{cases}
#'      {\left[1+\xi\left(\frac{x-\mu(a_1,b_1)}{\sigma(a_2,b_2)}\right)\right]}^{-1/\xi} & \text{if $\xi \ne 0$}\\
#'      \exp{\left(-\frac{x-\mu(a_1,b_1)}{\sigma(a_2,b_2)}\right)} & \text{if $\xi=0$}
#'    \end{cases}}
#' where
#' \eqn{x} is the random variable,
#' \eqn{\mu=a_1+b_1t_1} is the location parameter,
#' modelled as a function of parameters \eqn{a_1,b_1} and predictor \eqn{t_1},
#' \eqn{\sigma=e^{a_2+b_2t_2}} is the scale parameter,
#' modelled as a function of parameters \eqn{a_2,b_2} and predictor \eqn{t_2},
#' and \eqn{\xi} is the shape parameter.
#'
#' The calibrating prior we use is given by
#' \deqn{\pi(a_1,b_1,a_2,b_2,\xi) \propto 1}
#' as given in Jewson et al. (2025).
#'
#' The code will stop with an error if the
#' input data gives a maximum likelihood
#' value for the shape parameter that lies outside the range \code{(minxi,maxxi)},
#' since outside this range there may be numerical problems.
#' Such values seldom occur
#' in real observed data for maxima.
#'
#' @example man/examples/example_151_gev_p12.R
#'
#' @name gev_p12_cp
NULL
#' @rdname gev_p12_cp
#' @inheritParams man
#' @export
#'

qgev_p12_cp=function(x,t1,t2,t01=NA,t02=NA,n01=NA,n02=NA,p=seq(0.1,0.9,0.1),ics=c(0,0,0,0,0),
	fdalpha=0.01,
	minxi=-1,maxxi=1,
	means=FALSE,waicscores=FALSE,extramodels=FALSE,
	pdf=FALSE,dmgs=TRUE,rust=FALSE,nrust=100000,predictordata=TRUE,
	centering=TRUE,debug=FALSE){

	stopifnot(	is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1,
							length(t1)==length(x),length(t2)==length(x),
							length(ics)==5)
#
# 1 intro
#
	alpha=1-p
	nx=length(x)
	nalpha=length(alpha)
	t01=maket0(t01,n01,t1)
	t02=maket0(t02,n02,t2)
	if(debug)message(" t01=",t01)
	if(debug)message(" t02=",t02)
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
  	t1=t1-meant1
  	t2=t2-meant2
  	t01=t01-meant1
  	t02=t02-meant2
   }
#
# 3 ml param estimate
#
	if(debug)message(" ml param estimate")
#
# in a small number of cases, direct maxlik fails. I don't really know why, but it gives nonsense
# so I'm going to try using gev_p1 maxlik first, to determine initial conditions.
#
	ics=gev_p1_setics(x,t1,ics)
	opt1=optim(ics,gev_p1_loglik,x=x,t=t1,control=list(fnscale=-1))
	ics[1]=opt1$par[1]
	ics[2]=opt1$par[2]
	ics[3]=log(opt1$par[3])
	ics[4]=0
#	ics[5]=opt1$par[4]
	ics[5]=0 #to avoid an initial error occurring sometimes
#		gev_p12_setics(x,t1,t2,ics)
	opt1=optim(ics,gev_p12_loglik,x=x,t1=t1,t2=t2,control=list(fnscale=-1)) #this one uses the evd routine for dgev
	v1hat=opt1$par[1]
	v2hat=opt1$par[2]
	v3hat=opt1$par[3]
	v4hat=opt1$par[4]
	v5hat=opt1$par[5]
	ml_params=c(v1hat,v2hat,v3hat,v4hat,v5hat)
#	gev_p12_checkmle(ml_params,minxi,maxxi)
	if(debug)message(" ml_params=",ml_params)
	if((abs(v5hat)>=1)){revert2ml=TRUE}else{revert2ml=FALSE}
# I'm having some numerical problems with ldd in reliability testing...only in gev_p12...for nx=10 and xi=0.4
# maybe limiting v5hat to +1 in this way will help
# for samples of nx=10, and xi=0.4 this will be triggered very frequently, but so be it
#
# 4 predictordata
#
	prd=gev_p12_predictordata(predictordata,x,t1,t2,t01,t02,ml_params)
	predictedparameter=prd$predictedparameter
	adjustedx=prd$adjustedx
#
# 5 aic
#
	if(debug)message(" aic")
	ml_value=opt1$val
	maic=make_maic(ml_value,nparams=5)
#
# 6 calc ml quantiles and density
#
	if(debug)message(" ml_quantiles")
	ml_quantiles=qgev_p12((1-alpha),t01,t02,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat,xi=v5hat)
	if(v5hat<0){
		ml_max=(v1hat+v2hat*t01)-exp((v3hat+v4hat*t02))/v5hat
	} else {
		ml_max=Inf
	}
	fhat=dgev_p12(ml_quantiles,t01,t02,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat,xi=v5hat,log=FALSE)
	if(debug)message(" 1: ml_quantiles=",ml_quantiles)
  if(debug)message(" 1: fhat=",fhat)
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
			ml_quantilesm=qgev_p12((1-alpham),t01,t02,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat,xi=v5hat)
			ml_quantilesp=qgev_p12((1-alphap),t01,t02,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat,xi=v5hat)
			fhatm=dgev_p12(ml_quantilesm,t01,t02,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat,xi=v5hat,log=FALSE)
			fhatp=dgev_p12(ml_quantilesp,t01,t02,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat,xi=v5hat,log=FALSE)
		}
#
# 8 ldd
#
		if(debug)message("calc ldd")
		ldd=gev_p12_ldda(x,t1,t2,v1hat,v2hat,v3hat,v4hat,v5hat)
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
#
# 10 calculate lddd
#
		if(debug)message(" calc lddd")
		lddd=gev_p12_lddda(x,t1,t2,v1hat,v2hat,v3hat,v4hat,v5hat)
#
# 11 mu1
#
		if(debug)message(" calculate mu1")

		mu1=gev_p12_mu1fa(alpha,t01,t02,v1hat,v2hat,v3hat,v4hat,v5hat)
		if(pdf){
			mu1m=gev_p12_mu1fa(alpham,t01,t02,v1hat,v2hat,v3hat,v4hat,v5hat)
			mu1p=gev_p12_mu1fa(alphap,t01,t02,v1hat,v2hat,v3hat,v4hat,v5hat)
		}
#
# 12 mu2
#
		if(debug)message(" calculate mu2")
		mu2=gev_p12_mu2fa(alpha,t01,t02,v1hat,v2hat,v3hat,v4hat,v5hat)
		if(pdf){
			if(debug)message(" alpha pdf option")
			mu2m=gev_p12_mu2fa(alpham,t01,t02,v1hat,v2hat,v3hat,v4hat,v5hat)
			mu2p=gev_p12_mu2fa(alphap,t01,t02,v1hat,v2hat,v3hat,v4hat,v5hat)
		}
#
# 13 rh_flat model
#
		if(debug)message("call dmgs")
		lambdad_cp=matrix(0,5)
		dq=dmgs(lddi,lddd,mu1,lambdad_cp,mu2,dim=5)
		if(debug)message("make cp quantiles")
		rh_flat_quantiles=ml_quantiles+dq/(nx*fhat)
#
# 15 alpha pdf
#
		if(debug)message("step 15")
		if(pdf){
			lambdad_crhp_mle=matrix(0,4)
			dqm=dmgs(lddi,lddd,mu1m,lambdad_cp,mu2m,dim=5)
			dqp=dmgs(lddi,lddd,mu1p,lambdad_cp,mu2p,dim=5)
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
		if(debug)message("step 16")
		means=gev_p12_means(means,t01,t02,ml_params,nx)
		ml_mean				=means$ml_mean
		rh_flat_mean	=means$cp_mean
#
# 17 waicscores
#
		if(debug)message("step 17")
		waic=gev_p12_waic(waicscores,x,t1,t2,v1hat,v2hat,v3hat,v4hat,v5hat,
			lddi,lddd,lambdad_cp)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 19 rust
#
		if(debug)message("step 19")
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rgev_p12_cp(nrust,x,t1=t1,t2=t2,t01=t01,t02=t02,rust=TRUE,mlcp=FALSE,debug=debug)
			ru_quantiles=makeq(rustsim$ru_deviates,p)
		}
	} else {
		rh_flat_quantiles=ml_quantiles
		ru_quantiles=ml_quantiles
		rh_flat_pdf=ml_pdf
		rh_flat_mean=ml_mean
	} #end of if(dmgs)
#
# 20 decentering
#
	if(debug)message("step 20")
  if(centering){
  	if(debug)message(" qgev:ml_params,meant1=",ml_params,meant1)
    ml_params[1]=ml_params[1]-ml_params[2]*meant1
    ml_params[3]=ml_params[3]-ml_params[4]*meant2
    if(predictordata){
#    	if(debug)message("predictedparameter=",predictedparameter)
    	predictedparameter=predictedparameter-ml_params[2]*meant1
    }
  }

	if(debug)message("step 21")
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
				ml_max=ml_max,
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
#' @rdname gev_p12_cp
#' @inheritParams man
#' @export
rgev_p12_cp=function(n,x,t1,t2,t01=NA,t02=NA,n01=NA,n02=NA,ics=c(0,0,0,0,0),
		minxi=-1,maxxi=1,
		extramodels=FALSE,rust=FALSE,mlcp=TRUE,centering=TRUE,
	debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),is.finite(t1),is.finite(t2),!is.na(t1),!is.na(t2),
#						length(ics)==4)
	stopifnot(is.finite(x),!is.na(x),
							length(t1)==length(x),length(t2)==length(x),
		is.finite(t1),is.finite(t2),!is.na(t1),!is.na(t2),
						length(ics)==5)

	t01=maket0(t01,n01,t1)
	t02=maket0(t02,n02,t2)

#
# centering
#
	if(centering){
	  meant1=mean(t1)
	  meant2=mean(t2)
  	t1=t1-meant1
  	t2=t2-meant2
  	t01=t01-meant1
  	t02=t02-meant2
	}

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qgev_p12_cp(x,t1=t1,t2=t2,t01=t01,t02=t02,n01=NA,n02=NA,p=runif(n),ics=ics,
			extramodels=extramodels,centering=centering)
		ml_params=q$ml_params
		if(debug)message(" inside rgev_p12_cp: ml_params=",ml_params)
		ml_deviates=q$ml_quantiles
		ru_deviates=q$ru_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tgev_p12_cp(n,x,t1,t2)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			mu=th[i,1]+t01*th[i,2]
			sigma=exp(th[i,3]+t02*th[i,4])
			ru_deviates[i]=rgev(1,mu=mu,sigma=sigma,xi=th[i,5])
		}
	}

#
# decentering
#
 	if(debug)message(" rgev:ml_params,meant1=",ml_params,meant1)
	if(mlcp&centering){
	  ml_params[1]=ml_params[1]-ml_params[2]*meant1
		ml_params[3]=ml_params[3]-ml_params[4]*meant2
	}

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
			 cp_method=crhpflat_dmgs_cpmethod())

	return(op)

}
#' @rdname gev_p12_cp
#' @inheritParams man
#' @export
dgev_p12_cp=function(x,t1,t2,t01=NA,t02=NA,n01=NA,n02=NA,y=x,ics=c(0,0,0,0,0),
	minxi=-1,maxxi=1,extramodels=FALSE,
	rust=FALSE,nrust=10,centering=TRUE,debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
							length(t1)==length(x),length(t2)==length(x),
						is.finite(t1),is.finite(t2),!is.na(t1),!is.na(t2),
		length(ics)==5)

	if(debug)message(" maket0")
	t01=maket0(t01,n01,t1)
	t02=maket0(t02,n02,t2)

#
# centering
#
  if(centering){
		if(debug)message(" centering")
	  meant1=mean(t1)
	  meant2=mean(t2)
  	t1=t1-meant1
  	t2=t2-meant2
  	t01=t01-meant1
  	t02=t02-meant2
  }

	if(debug)message(" ics and optim")
	ics=gev_p12_setics(x,t1,t2,ics)
	opt1=optim(ics,gev_p12_loglik,x=x,t1=t1,t2=t2,control=list(fnscale=-1))
	v1hat=opt1$par[1]
	v2hat=opt1$par[2]
	v3hat=opt1$par[3]
	v4hat=opt1$par[4]
	v5hat=opt1$par[5]
	if(v5hat<=(-1)){revert2ml=TRUE}else{revert2ml=FALSE}
	ml_params=c(v1hat,v2hat,v3hat,v4hat,v5hat)
#	gev_p12_checkmle(ml_params,minxi,maxxi)
	if(debug)message(" call sub")
	dd=dgev_p12sub(x=x,t1=t1,t2=t2,y=y,t01=t01,t02=t02,ics=ics,
		minxi,maxxi,extramodels=extramodels,debug=debug)
	ru_pdf="rust not selected"
		ml_params=dd$ml_params

	if(rust&&(!revert2ml)){
		if(debug)message(" rust")
		th=tgev_p12_cp(nrust,x=x,t1=t1,t2=t2,debug=debug)$theta_samples
		if(debug)message(" tgev call done")
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+t01*th[ir,2]
			sigma=exp(th[ir,3]+t02*th[ir,4])
			ru_pdf=ru_pdf+dgev(y,mu=mu,sigma=sigma,xi=th[ir,5])
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
  }

		op=list(
					ml_params=ml_params,
					ml_pdf=dd$ml_pdf,
					revert2ml=revert2ml,
					ru_pdf=ru_pdf,
					cp_method=nopdfcdfmsg())

	return(op)

}
#' @rdname gev_p12_cp
#' @inheritParams man
#' @export
pgev_p12_cp=function(x,t1,t2,t01=NA,t02=NA,n01=NA,n02=NA,y=x,ics=c(0,0,0,0,0),
	minxi=-1,maxxi=1,extramodels=FALSE,
	rust=FALSE,nrust=1000,centering=TRUE,debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
							length(t1)==length(x),length(t2)==length(x),
						is.finite(t1),is.finite(t2),!is.na(t1),!is.na(t2),
		length(ics)==5)

	t01=maket0(t01,n01,t1)
	t02=maket0(t02,n02,t2)
#
# centering
#
  if(centering){
	  meant1=mean(t1)
	  meant2=mean(t2)
  	t1=t1-meant1
  	t2=t2-meant2
  	t01=t01-meant1
  	t02=t02-meant2
  }

	ics=gev_p12_setics(x,t1,t2,ics)
	opt1=optim(ics,gev_p12_loglik,x=x,t1=t1,t2=t2,control=list(fnscale=-1)) #this one uses the evd routine for dgev
	v1hat=opt1$par[1]
	v2hat=opt1$par[2]
	v3hat=opt1$par[3]
	v4hat=opt1$par[4]
	v5hat=opt1$par[5]
	if(v5hat<=(-1)){revert2ml=TRUE}else{revert2ml=FALSE}
	ml_params=c(v1hat,v2hat,v3hat,v4hat,v5hat)
#	gev_p12_checkmle(ml_params,minxi,maxxi)
	dd=dgev_p12sub(x=x,t1=t1,t2=t2,y=y,t01=t01,t02=t02,ics=ics,
		minxi,maxxi,extramodels=extramodels,debug=debug)
	ru_cdf="rust not selected"
		ml_params=dd$ml_params

	if(rust&&(!revert2ml)){
		th=tgev_p12_cp(nrust,x,t1,t2)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+t01*th[ir,2]
			sigma=exp(th[ir,3]+t02*th[ir,4])
			ru_cdf=ru_cdf+pgev(y,mu=mu,sigma=exp(sigma),xi=th[ir,5])
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
  }

	op=list(
					ml_params=ml_params,
					ml_cdf=dd$ml_cdf,
					revert2ml=revert2ml,
					ru_cdf=ru_cdf,
					cp_method=nopdfcdfmsg())
	return(op)
}
#' @rdname gev_p12_cp
#' @inheritParams man
#' @export
tgev_p12_cp=function(n,x,t1,t2,ics=c(0,0,0,0,0),
	extramodels=FALSE,debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),
							length(t1)==length(x),length(t2)==length(x),
							is.finite(t1),is.finite(t2),!is.na(t1),!is.na(t2),
							length(ics)==5)

#
# centering
#
	  meant1=mean(t1)
	  meant2=mean(t2)
  	t1=t1-meant1
  	t2=t2-meant2

	if(debug)message("sums=",sum(x),sum(t1),sum(t2),n)
	th=ru(gev_p12_logf,x=x,t1=t1,t2=t2,n=n,d=5,init=c(0,0,0,0,0))
	if(debug)message(" back from rust")
  theta_samples=th$sim_vals
#
# decentering
#
  theta_samples[,1]=theta_samples[,1]-theta_samples[,2]*meant1
  theta_samples[,3]=theta_samples[,3]-theta_samples[,4]*meant2

	list(theta_samples=theta_samples)

}
