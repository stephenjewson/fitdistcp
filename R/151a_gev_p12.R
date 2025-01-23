#' Generalized Extreme Value Distribution with Two Predictors, Predictions based on a Calibrating Prior
#'
#' @inherit man description author references seealso
#' @inheritParams man
#'
#' @inheritSection man Default Return Values
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
#' as given in Jewson et al. (2024).
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
	d1=0.01,d2=0.01,d3=0.01,d4=0.01,d5=0.01,fdalpha=0.01,
	minxi=-0.45,maxxi=0.45,
	means=FALSE,waicscores=FALSE,extramodels=FALSE,
	pdf=FALSE,dmgs=TRUE,rust=FALSE,nrust=100000,predictordata=TRUE,
	centering=TRUE,debug=FALSE,aderivs=TRUE){

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
	if(debug)cat(" t01=",t01,"\n")
	if(debug)cat(" t02=",t02,"\n")
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
	if(debug)cat(" ml param estimate\n")
	ics=gev_p12_setics(x,t1,t2,ics)
	opt1=optim(ics,gev_p12_loglik,x=x,t1=t1,t2=t2,control=list(fnscale=-1)) #this one uses the evd routine for dgev
	v1hat=opt1$par[1]
	v2hat=opt1$par[2]
	v3hat=opt1$par[3]
	v4hat=opt1$par[4]
	v5hat=opt1$par[5]
	ml_params=c(v1hat,v2hat,v3hat,v4hat,v5hat)
	gev_p12_checkmle(ml_params,minxi,maxxi)
#	cat(" ml_params=",ml_params,"\n")
	if(debug)cat(" ml_params=",ml_params,"\n")
#
# 4 predictordata
#
	prd=gev_p12_predictordata(predictordata,x,t1,t2,t01,t02,ml_params)
	predictedparameter=prd$predictedparameter
	adjustedx=prd$adjustedx
#
# 5 aic
#
	if(debug)cat(" aic\n")
	ml_value=opt1$val
	maic=make_maic(ml_value,nparams=5)
#
# 6 calc ml quantiles and density
#
	if(debug)cat(" ml_quantiles\n")
	ml_quantiles=qgev_p12((1-alpha),t01,t02,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat,xi=v5hat)
	if(v5hat<0){
		ml_max=(v1hat+v2hat*t01)-exp((v3hat+v4hat*t02))/v5hat
	} else {
		ml_max=Inf
	}
	fhat=dgev_p12(ml_quantiles,t01,t02,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat,xi=v5hat,log=FALSE)
	if(debug)cat(" 1: ml_quantiles=",ml_quantiles,"\n")
  if(debug)cat(" 1: fhat=",fhat,"\n")
#	cat("****************************************************\n")
#	cat(" 1: ml_quantiles=",ml_quantiles,"\n")
# cat(" 1: fhat=",fhat,"\n")
#
# dmgs
#
	standard_errors="dmgs not selected"
	cp_quantiles="dmgs not selected"
	ru_quantiles="dmgs not selected"
	ml_pdf="dmgs not selected"
	cp_pdf="dmgs not selected"
	waic1="dmgs not selected"
	waic2="dmgs not selected"
	ml_mean="dmgs not selected"
	cp_mean="dmgs not selected"
	if(dmgs){
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
# 8 ldd (two versions)
#
		if(debug)cat("calc ldd\n")
		if(aderivs) ldd=gev_p12_ldda(x,t1,t2,v1hat,v2hat,v3hat,v4hat,v5hat)
		if(!aderivs)ldd=gev_p12_ldd(x,t1,t2,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4,v5hat,d5)
#		if(debug)cat(" ldd=",ldd,"\n")
#		if(debug)cat(" det(ldd)=",det(ldd),"\n")
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
		if(extramodels|means)ldd_k5=gev_p12k3_ldd(x,t1,t2,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4,v5hat)
		if(extramodels|means)lddi_k5=solve(ldd_k5)
#
# 9 information matrix and related (for Jeffreys prior)
# -because of difficulty of calculating expected information, I just use observed information
		if(debug)cat(" call gev.infomat\n")
		if(extramodels|means){
			gg=-ldd
			ggi=solve(gg)
			detg=det(gg)
			ggd=gev_p12_ggd(x,t1,t2,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4,v5hat,d5)	#uses ldd! so not expected at all, just observed
		}
#
# 10 calculate lddd (two versions)
#
		if(debug)cat(" calc lddd\n")
		if(aderivs) lddd=gev_p12_lddda(x,t1,t2,v1hat,v2hat,v3hat,v4hat,v5hat)
		if(!aderivs)lddd=gev_p12_lddd(x,t1,t2,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4,v5hat,d5)
		if(extramodels|means)lddd_k5=gev_p12k3_lddd(x,t1,t2,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d5,v5hat)
#
# 11 mu1 (two versions)
#
		if(debug)cat(" calculate mu1\n")

		if(aderivs) mu1=gev_p12_mu1fa(alpha,t01,t02,v1hat,v2hat,v3hat,v4hat,v5hat)
		if(!aderivs)mu1=gev_p12_mu1f(alpha,t01,t02,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4,v5hat,d5)

		if(extramodels|means)mu1_k5=gev_p12k3_mu1f(alpha,t01,t02,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d5,v5hat)
		if(pdf){
			if(aderivs){
				mu1m=gev_p12_mu1fa(alpham,t01,t02,v1hat,v2hat,v3hat,v4hat,v5hat)
				mu1p=gev_p12_mu1fa(alphap,t01,t02,v1hat,v2hat,v3hat,v4hat,v5hat)
			} else {
				mu1m=gev_p12_mu1f(alpham,t01,t02,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4,v5hat,d5)
				mu1p=gev_p12_mu1f(alphap,t01,t02,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4,v5hat,d5)
			}
		}
#
# 12 mu2 (two versions)
#
		if(debug)cat(" calculate mu2\n")
		mu2=gev_p12_mu2f(alpha,t01,t02,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4,v5hat,d5)
		if(extramodels|means){
			if(debug)cat(" calculate mu2_k5\n")
			mu2_k5=gev_p12k3_mu2f(alpha,t01,t02,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4,v5hat)
		}
		if(pdf){
			if(debug)cat(" alpha pdf option\n")
			mu2m=gev_p12_mu2f(alpham,t01,t02,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4,v5hat,d5)
			mu2p=gev_p12_mu2f(alphap,t01,t02,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4,v5hat,d5)
		}
#
# compare ldd methods
#
#		ldda=gev_p12_ldda(x,t1,t2,v1hat,v2hat,v3hat,v4hat,v5hat)
#		cat("ldd aderivs=",sum(ldd),"\n")
#		lddn=gev_p12_ldd(x,t1,t2,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4,v5hat,d5)
#		cat("ldd nderivs=",sum(ldd),"\n")
#		diff=lddn-ldda
#		cat("diff=",sum(diff),"\n")
#
# compare lddd methods
#
#		lddda=gev_p12_lddda(x,t1,t2,v1hat,v2hat,v3hat,v4hat,v5hat)
#		cat("lddd aderivs=",sum(lddda),"\n")
#		ldddn=gev_p12_lddd(x,t1,t2,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4,v5hat,d5)
#		cat("lddd nderivs=",sum(lddd),"\n")
#		diff=ldddn-lddda
#		cat("diff=",sum(diff),"\n")
#		stop()
##
# 13 model 1: cp=flat prior
#
#		if(debug)cat(" lddi=",lddi,"\n")
#		if(debug)cat(" lddd=",lddd,"\n")
#		if(debug)cat(" mu1=",mu1,"\n")
#		if(debug)cat(" mu2=",mu2,"\n")
		if(debug)cat("call dmgs\n")
		lambdad_cp=matrix(0,5)
		dq=dmgs(lddi,lddd,mu1,lambdad_cp,mu2,dim=5)
#		if(debug)cat(" dq=",dq,"\n")
		if(debug)cat("make cp quantiles\n")
		cp_quantiles=ml_quantiles+dq/(nx*fhat)
#		cat("lddi=",lddi,"\n")
#		cat("lddd=",lddd,"\n")
#		cat("mu1=",mu1,"\n")
#		cat("mu2=",mu2,"\n")
#		cat("dq=",dq,"\n")
#		cat("fhat=",fhat,"\n")
#		cat("cp_quantiles=",cp_quantiles,"\n")
#		if(sum(is.na(cp_quantiles))>0)stop()

#		if(debug)cat(" cp_quantiles=",cp_quantiles,"\n")
#		cat("ml_quantiles=",ml_quantiles,"\n")
#		cat("flat_quantiles=",flat_quantiles,"\n")
#
# 14 model 4: rh_Flat with flat prior on shape (needs to use 4d version of Bayesian code)
#
		if(debug)cat("step 14\n")
		if(extramodels|means){
			lambdad_crhp_mle=matrix(0,4)
			dq=dmgs(lddi,lddd,mu1,lambdad_crhp_mle,mu2,dim=4)
			crhp_mle_quantiles=ml_quantiles+dq/(nx*fhat)
		} else {
			crhp_mle_quantiles="extramodels not selected"
		}
#
# 15 alpha pdf
#
		if(debug)cat("step 15\n")
		if(pdf){
			lambdad_crhp_mle=matrix(0,4)
			dqm=dmgs(lddi,lddd,mu1m,lambdad_crhp_mle,mu2m,dim=4)
			dqp=dmgs(lddi,lddd,mu1p,lambdad_crhp_mle,mu2p,dim=4)
			quantilesm=ml_quantilesm+dqm/(nx*fhatm)
			quantilesp=ml_quantilesp+dqp/(nx*fhatp)
			ml_pdf=fhat
			cp_pdf=-(alphap-alpham)/(quantilesp-quantilesm)
		} else{
			ml_pdf=fhat
			cp_pdf="pdf not selected"
		}
#
# 16 means
#
		if(debug)cat("step 16\n")
		means=gev_p12_means(means,t01,t02,ml_params,nx)
		ml_mean				=means$ml_mean
		cp_mean				=means$cp_mean
		crhp_mle_mean	=means$flat_mean
#
# 17 waicscores
#
		if(debug)cat("step 17\n")
		waic=gev_p12_waic(waicscores,x,t1,t2,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4,v5hat,d5,
			lddi,lddd,lambdad_cp,aderivs)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 19 rust
#
		if(debug)cat("step 19\n")
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rgev_p12_cp(nrust,x,t1=t1,t2=t2,t01=t01,t02=t02,rust=TRUE,mlcp=FALSE,debug=debug)
			ru_quantiles=makeq(rustsim$ru_deviates,p)
		}
	} #end of if(dmgs)
#
# 20 decentering
#
	if(debug)cat("step 20\n")
  if(centering){
  	if(debug)cat(" qgev:ml_params,meant1=",ml_params,meant1,"\n")
    ml_params[1]=ml_params[1]-ml_params[2]*meant1
    ml_params[3]=ml_params[3]-ml_params[4]*meant2
    if(predictordata){
#    	if(debug)cat("predictedparameter=",predictedparameter,"\n")
    	predictedparameter=predictedparameter-ml_params[2]*meant1
    }
  }

	if(debug)cat("step 21\n")
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
				ml_max=ml_max,
				cp_quantiles=cp_quantiles,
				ru_quantiles=ru_quantiles,
				ml_pdf=ml_pdf,
				cp_pdf=cp_pdf,
				maic=maic,
				waic1=waic1,
				waic2=waic2,
				ml_mean=ml_mean,
				cp_mean=cp_mean,
				cp_method=crhpflat_dmgs_cpmethod())

}
#' @rdname gev_p12_cp
#' @inheritParams man
#' @export
rgev_p12_cp=function(n,x,t1,t2,t01=NA,t02=NA,n01=NA,n02=NA,ics=c(0,0,0,0,0),
	d1=0.01,d2=0.01,d3=0.01,d4=0.01,d5=0.01,
		minxi=-0.45,maxxi=0.45,
		extramodels=FALSE,rust=FALSE,mlcp=TRUE,centering=TRUE,
	debug=FALSE,aderivs=TRUE){

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
		q=qgev_p12_cp(x,t1=t1,t2=t2,t01=t01,t02=t02,n01=NA,n02=NA,p=runif(n),ics=ics,d1=d1,d2=d2,d3=d3,d4=d4,d5=d5,
			extramodels=extramodels,centering=centering,aderivs=aderivs)
		ml_params=q$ml_params
		if(debug)cat(" inside rgev_p12_cp: ml_params=",ml_params,"\n")
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
 	if(debug)cat(" rgev:ml_params,meant1=",ml_params,meant1,"\n")
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
	d1=0.01,d2=0.01,d3=0.01,d4=0.01,d5=0.01,
	minxi=-0.45,maxxi=0.45,extramodels=FALSE,
	rust=FALSE,nrust=10,centering=TRUE,debug=FALSE,aderivs=TRUE){

#	cat("inside dgev_p12_cp\n")
#	cat("debug=",debug,"\n")

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
							length(t1)==length(x),length(t2)==length(x),
						is.finite(t1),is.finite(t2),!is.na(t1),!is.na(t2),
		length(ics)==5)

	if(debug)cat(" maket0\n")
	t01=maket0(t01,n01,t1)
	t02=maket0(t02,n02,t2)

#
# centering
#
  if(centering){
		if(debug)cat(" centering\n")
	  meant1=mean(t1)
	  meant2=mean(t2)
  	t1=t1-meant1
  	t2=t2-meant2
  	t01=t01-meant1
  	t02=t02-meant2
  }

	if(debug)cat(" ics and optim\n")
	ics=gev_p12_setics(x,t1,t2,ics)
	opt1=optim(ics,gev_p12_loglik,x=x,t1=t1,t2=t2,control=list(fnscale=-1))
	v1hat=opt1$par[1]
	v2hat=opt1$par[2]
	v3hat=opt1$par[3]
	v4hat=opt1$par[4]
	v5hat=opt1$par[5]
	ml_params=c(v1hat,v2hat,v3hat,v4hat,v5hat)
	gev_p12_checkmle(ml_params,minxi,maxxi)
	if(debug)cat(" call sub\n")
	dd=dgev_p12sub(x=x,t1=t1,t2=t2,y=y,t01=t01,t02=t02,ics=ics,d1,d2,d3,d4,d5,
		minxi,maxxi,extramodels=extramodels,debug=debug,aderivs=aderivs)
	ru_pdf="rust not selected"
		ml_params=dd$ml_params

	if(rust){
		if(debug)cat(" rust\n")
		th=tgev_p12_cp(nrust,x=x,t1=t1,t2=t2,debug=debug)$theta_samples
		if(debug)cat(" tgev call done\n")
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+t01*th[ir,2]
			sigma=exp(th[ir,3]+t02*th[ir,4])
			ru_pdf=ru_pdf+dgev(y,mu=mu,sigma=sigma,xi=th[ir,5])
		}
		ru_pdf=ru_pdf/nrust
	}
#
# decentering
#
  if(centering){
		if(debug)cat(" decentering\n")
	 	if(debug)cat(" dgev:ml_params,meant1=",ml_params,meant1,"\n")
    ml_params[1]=ml_params[1]-ml_params[2]*meant1
    ml_params[3]=ml_params[3]-ml_params[4]*meant2
  }

		op=list(
					ml_params=ml_params,
					ml_pdf=dd$ml_pdf,
					ru_pdf=ru_pdf,
					cp_method=nopdfcdfmsg())

	return(op)

}
#' @rdname gev_p12_cp
#' @inheritParams man
#' @export
pgev_p12_cp=function(x,t1,t2,t01=NA,t02=NA,n01=NA,n02=NA,y=x,ics=c(0,0,0,0,0),
	d1=0.01,d2=0.01,d3=0.01,d4=0.01,d5=0.01,
	minxi=-0.45,maxxi=0.45,extramodels=FALSE,
	rust=FALSE,nrust=1000,centering=TRUE,debug=FALSE,aderivs=TRUE){

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
	ml_params=c(v1hat,v2hat,v3hat,v4hat,v5hat)
	gev_p12_checkmle(ml_params,minxi,maxxi)
	dd=dgev_p12sub(x=x,t1=t1,t2=t2,y=y,t01=t01,t02=t02,ics=ics,d1,d2,d3,d4,d5,
		minxi,maxxi,extramodels=extramodels,debug=debug,aderivs=aderivs)
	ru_cdf="rust not selected"
		ml_params=dd$ml_params

	if(rust){
		th=tgev_p12_cp(nrust,x,t1,t2)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+t01*th[ir,2]
			sigma=exp(th[ir,3]+t02*th[ir,4])
			ru_cdf=ru_cdf+pgev(y,mu=mu,sigma=exp(sigma),xi=th[ir,5])
		}
		ru_cdf=ru_cdf/nrust
	}
#
# decentering
#
  if(centering){
	 	if(debug)cat(" pgev:ml_params,meant1=",ml_params,meant1,"\n")
    ml_params[1]=ml_params[1]-ml_params[2]*meant1
    ml_params[3]=ml_params[3]-ml_params[4]*meant2
  }

	op=list(
					ml_params=ml_params,
					ml_cdf=dd$ml_cdf,
					ru_cdf=ru_cdf,
					cp_method=nopdfcdfmsg())
	return(op)
}
#' @rdname gev_p12_cp
#' @inheritParams man
#' @export
tgev_p12_cp=function(n,x,t1,t2,ics=c(0,0,0,0,0),
	d1=0.01,d2=0.01,d3=0.01,d4=0.01,d5=0.01,
	minxi=-0.45,maxxi=0.45,
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

	if(debug)cat("sums=",sum(x),sum(t1),sum(t2),n,"\n")
	th=ru(gev_p12_logf,x=x,t1=t1,t2=t2,n=n,d=5,init=c(0,0,0,0,0))
	if(debug)cat(" back from rust\n")
  theta_samples=th$sim_vals
#
# decentering
#
  theta_samples[,1]=theta_samples[,1]-theta_samples[,2]*meant1
  theta_samples[,3]=theta_samples[,3]-theta_samples[,4]*meant2

	list(theta_samples=theta_samples)

}
