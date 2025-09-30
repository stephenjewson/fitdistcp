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
#' The GEV distribution with a 2d predictor has distribution function
#' \deqn{F(x;a,b_1,b_2,\sigma,\xi)=\exp{(-t(x;\mu(a,b_1,b_2),\sigma,\xi))}}
#' where
#' \deqn{t(x;\mu(a,b_1,b_2),\sigma,\xi) =
#'    \begin{cases}
#'      {\left[1+\xi\left(\frac{x-\mu(a,b_1,b_2)}{\sigma}\right)\right]}^{-1/\xi} & \text{if $\xi \ne 0$}\\
#'      \exp{\left(-\frac{x-\mu(a,b_1,b_2)}{\sigma}\right)} & \text{if $\xi=0$}
#'    \end{cases}}
#' where
#' \eqn{x} is the random variable,
#' \eqn{\mu=a+b_1t_1+b_2t_2} is the location parameter,
#' modelled as a function of parameters \eqn{a,b_1,b_2} and predictor \eqn{t_1,t_2},
#' and \eqn{\sigma>0,\xi} are the scale and shape parameters.
#'
#' The calibrating prior we use is given by
#' \deqn{\pi(a,b_1,b_2,\sigma,\xi) \propto \frac{1}{\sigma}}
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
#' @name gev_p1n_cp
NULL
#' @rdname gev_p1n_cp
#' @inheritParams man
#' @export
#'
qgev_p1n_cp=function(x,t,t0=NA,n0=NA,p=seq(0.1,0.9,0.1),
	fdalpha=0.01,
	minxi=-1,maxxi=1,
	means=FALSE,waicscores=FALSE,extramodels=FALSE,
	pdf=FALSE,dmgs=TRUE,rust=FALSE,nrust=100000,predictordata=TRUE,
	centering=TRUE,debug=FALSE){

# notes
# -I've removed the option to pass in ics, because I couldn't make that depend on nt
# -routines made using auto deriv need individual parameters
# -and need to be called differently depending on nt
# -all other routines use the vhat vector, and are called the same

	stopifnot(	is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1)
#
# 1 intro
#
	t=ifvectorthenmatrix(t)
	nt=findnt(t)
	alpha=1-p
	nx=length(x)
	nalpha=length(alpha)
	for (i in 1:nt){
		t0[i]=maket0(t0[i],n0[i],t[,i])
	}
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
		meant=matrix(0,nt)
   	for (i in 1:nt){
		 meant[i]=mean(t[,i])
	   t[,i]=t[,i]-meant[i]
	   t0[i]=t0[i]-meant[i]
  	}
  }
#
# 3 ml param estimate
#
	if(debug)message(" ml param estimate")
	ics=gev_p1n_setics(x,t)
	opt1=optim(ics,gev_p1n_loglik,x=x,t=t,control=list(fnscale=-1)) #this one uses the evd routine for dgev
	vhat=opt1$par
	ml_params=vhat
#	gev_p1n_checkmle(ml_params,minxi,maxxi)
	if(debug)message(" ml_params=",ml_params)
	if(abs(vhat[3+nt])>=1){revert2ml=TRUE}else{revert2ml=FALSE}
#
# 4 predictordata
#
	prd=gev_p1n_predictordata(predictordata,x,t,t0,ml_params)
	predictedparameter=prd$predictedparameter
	adjustedx=prd$adjustedx
#
# 5 aic
#
	if(debug)message(" aic")
	ml_value=opt1$val
	maic=make_maic(ml_value,nparams=(3+nt))
#
# 6 calc ml quantiles and density
#
	if(debug)message(" ml_quantiles")
	ml_quantiles=qgev_p1n((1-alpha),t0,params=vhat)
	if(vhat[3+nt]<0){
		mu=vhat[1]+makebetat0(nt,vhat,t0)
		ml_max=mu-vhat[nt+2]/vhat[nt+3]
	} else {
		ml_max=Inf
	}
	fhat=dgev_p1n(ml_quantiles,t0,params=vhat,log=FALSE)
	if(debug)message(" ml_quantiles=",ml_quantiles)
	if(debug)message(" fhat=",fhat)
#
# dmgs
#
	standard_errors="dmgs not selected"
	rh_flat_quantiles="dmgs not selected"
	cp_quantiles="dmgs not selected"
	ru_quantiles="dmgs not selected"
	rh_flat_pdf="dmgs not selected"
	ml_pdf="dmgs not selected"
	cp_pdf="dmgs not selected"
	waic1="dmgs not selected"
	waic2="dmgs not selected"
	ml_mean="dmgs not selected"
	cp_mean="dmgs not selected"
	rh_flat_mean="dmgs not selected"
	cp_method="dmgs not selected"
	if((dmgs)&&(!revert2ml)){
#
# 7 alpha pdf stuff
# -for now, I only make the pdf. I could add cdf too I suppose. Somehow.
		if(pdf){

			ml_quantilesm=qgev_p1n((1-alpham),t0,vhat)
			ml_quantilesp=qgev_p1n((1-alphap),t0,vhat)
			fhatm=dgev_p1n(ml_quantilesm,t0,params=vhat,log=FALSE)
			fhatp=dgev_p1n(ml_quantilesp,t0,params=vhat,log=FALSE)
		}
#
# 8 ldd
#
		if(nt==1)ldd=gev_p1a_ldda(x,t[,1],						vhat[1],vhat[2],vhat[3],vhat[4])
		if(nt==2)ldd=gev_p1b_ldda(x,t[,1],t[,2],			vhat[1],vhat[2],vhat[3],vhat[4],vhat[5])
		if(nt==3)ldd=gev_p1c_ldda(x,t[,1],t[,2],t[,3],vhat[1],vhat[2],vhat[3],vhat[4],vhat[5],vhat[6])
		if(debug)message(" ldd=",ldd)
		if(debug)message(" det(ldd)=",det(ldd))
		lddi=solve(ldd)

		standard_errors=make_se(nx,lddi)
#
# 10 calculate lddd
#
		if(debug)message(" lddd")
		if(nt==1)lddd=gev_p1a_lddda(x,t[,1],						vhat[1],vhat[2],vhat[3],vhat[4])
		if(nt==2)lddd=gev_p1b_lddda(x,t[,1],t[,2],			vhat[1],vhat[2],vhat[3],vhat[4],vhat[5])
		if(nt==3)lddd=gev_p1c_lddda(x,t[,1],t[,2],t[,3],vhat[1],vhat[2],vhat[3],vhat[4],vhat[5],vhat[6])
#
# 11 mu1
#
		if(debug)message(" calculate mu1")
		if(nt==1)mu1=gev_p1a_mu1fa(alpha,t0[1],							vhat[1],vhat[2],vhat[3],vhat[4])
		if(nt==2)mu1=gev_p1b_mu1fa(alpha,t0[1],t0[2],				vhat[1],vhat[2],vhat[3],vhat[4],vhat[5])
		if(nt==3)mu1=gev_p1c_mu1fa(alpha,t0[1],t0[2],t0[3],	vhat[1],vhat[2],vhat[3],vhat[4],vhat[5],vhat[6])

		if(pdf){

			if(nt==1)mu1m=gev_p1a_mu1fa(alpham,t0[1],							vhat[1],vhat[2],vhat[3],vhat[4])
			if(nt==2)mu1m=gev_p1b_mu1fa(alpham,t0[1],t0[2],				vhat[1],vhat[2],vhat[3],vhat[4],vhat[5])
			if(nt==3)mu1m=gev_p1c_mu1fa(alpham,t0[1],t0[2],t0[3],	vhat[1],vhat[2],vhat[3],vhat[4],vhat[5],vhat[6])

			if(nt==1)mu1p=gev_p1a_mu1fa(alphap,t0[1],							vhat[1],vhat[2],vhat[3],vhat[4])
			if(nt==2)mu1p=gev_p1b_mu1fa(alphap,t0[1],t0[2],				vhat[1],vhat[2],vhat[3],vhat[4],vhat[5])
			if(nt==3)mu1p=gev_p1c_mu1fa(alphap,t0[1],t0[2],t0[3],	vhat[1],vhat[2],vhat[3],vhat[4],vhat[5],vhat[6])
		}
#
# 12 mu2
#
		if(debug)message(" calculate mu2")
		if(nt==1)mu2=gev_p1a_mu2fa(alpha,t0[1],							vhat[1],vhat[2],vhat[3],vhat[4])
		if(nt==2)mu2=gev_p1b_mu2fa(alpha,t0[1],t0[2],				vhat[1],vhat[2],vhat[3],vhat[4],vhat[5])
		if(nt==3)mu2=gev_p1c_mu2fa(alpha,t0[1],t0[2],t0[3],	vhat[1],vhat[2],vhat[3],vhat[4],vhat[5],vhat[6])

		if(pdf){

			if(nt==1)mu2m=gev_p1a_mu2fa(alpham,t0[1],							vhat[1],vhat[2],vhat[3],vhat[4])
			if(nt==2)mu2m=gev_p1b_mu2fa(alpham,t0[1],t0[2],				vhat[1],vhat[2],vhat[3],vhat[4],vhat[5])
			if(nt==3)mu2m=gev_p1c_mu2fa(alpham,t0[1],t0[2],t0[3],	vhat[1],vhat[2],vhat[3],vhat[4],vhat[5],vhat[6])

			if(nt==1)mu2p=gev_p1a_mu2fa(alphap,t0[1],							vhat[1],vhat[2],vhat[3],vhat[4])
			if(nt==2)mu2p=gev_p1b_mu2fa(alphap,t0[1],t0[2],				vhat[1],vhat[2],vhat[3],vhat[4],vhat[5])
			if(nt==3)mu2p=gev_p1c_mu2fa(alphap,t0[1],t0[2],t0[3],	vhat[1],vhat[2],vhat[3],vhat[4],vhat[5],vhat[6])
		}
#
# 15 model 4: rh_Flat with flat prior on shape (needs to use 4d version of Bayesian code)
#
		if(nt==1)lambdad_rh_flat=c(0,0,			-1/vhat[3],0)
		if(nt==2)lambdad_rh_flat=c(0,0,0,		-1/vhat[4],0)
		if(nt==3)lambdad_rh_flat=c(0,0,0,0,	-1/vhat[5],0)
		dq=dmgs(lddi,lddd,mu1,lambdad_rh_flat,mu2,dim=(nt+3))
		rh_flat_quantiles=ml_quantiles+dq/(nx*fhat)
		if(pdf){
			dqm=dmgs(lddi,lddd,mu1m,lambdad_rh_flat,mu2m,dim=(nt+3))
			dqp=dmgs(lddi,lddd,mu1p,lambdad_rh_flat,mu2p,dim=(nt+3))
			quantilesm=ml_quantilesm+dqm/(nx*fhatm)
			quantilesp=ml_quantilesp+dqp/(nx*fhatp)
			ml_pdf=fhat
			rh_flat_pdf=-(alphap-alpham)/(quantilesp-quantilesm)
		} else{
			ml_pdf=fhat
			rh_flat_pdf="pdf not selected"
		}
#
# 17 means
#
		means=gev_p1n_means(means,t0,ml_params,lddi,lddd,
											lambdad_rh_flat,nx,dim=(nt+3))
		ml_mean				=means$ml_mean
		rh_flat_mean	=means$rh_flat_mean
#
# 18 waicscores
#
		waic=gev_p1n_waic(waicscores,x,t,vhat,lddi,lddd,lambdad_rh_flat)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 20 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rgev_p1n_cp(nrust,x,t=t,t0=t0,rust=TRUE,mlcp=FALSE)
			ru_quantiles=makeq(rustsim$ru_deviates,p)
		}
	} else {
		rh_flat_quantiles=ml_quantiles
		ru_quantiles=ml_quantiles
		rh_flat_pdf=ml_pdf
		rh_flat_mean=ml_mean
	} #end of if(dmgs)
#
# 21 decentering
#
  if(centering){
	    ml_params[1]=ml_params[1]-makebetat0(nt,ml_params,meant)
	   if(predictordata)predictedparameter=predictedparameter-makebetat0(nt,ml_params,meant)
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
#' @rdname gev_p1n_cp
#' @inheritParams man
#' @export
rgev_p1n_cp=function(n,x,t,t0=NA,n0=NA,
		minxi=-1,maxxi=1,
		extramodels=FALSE,rust=FALSE,mlcp=TRUE,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),is.finite(t),!is.na(t),
#						length(ics)==5)
	stopifnot(is.finite(x),!is.na(x),is.finite(t),!is.na(t))

	t=ifvectorthenmatrix(t)
	nt=findnt(t)

	for (i in 1:nt){
		t0[i]=maket0(t0=t0[i],n0=n0[i],t=t[,i])
	}

#
# centering
#
	meant=matrix(0,nt)
	for (i in 1:nt){
	  meant[i]=mean(t[,i])
		t[,i]=t[,i]-meant[i]
		t0[i]=t0[i]-meant[i]
	}

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	ru_deviates="rust not selected"
	cp_deviates="rust not selected"

	if(mlcp){
		q=qgev_p1n_cp(x,t=t,t0=t0,n0=NA,p=runif(n),extramodels=extramodels)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		ru_deviates=q$ru_quantiles
		cp_deviates=q$cp_quantiles
	}

	if(rust){
		th=tgev_p1n_cp(n,x,t)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			mu=th[i,1]+makebetat0(nt,th[i,],t0)
			ru_deviates[i]=rgev(1,mu=mu,sigma=th[i,(nt+2)],xi=th[i,(nt+3)])
		}
	}

#
# decentering
#
  if(mlcp)ml_params[1]=ml_params[1]-makebetat0(nt,ml_params,meant)

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
			 cp_method=crhpflat_dmgs_cpmethod())

	return(op)

}
#' @rdname gev_p1n_cp
#' @inheritParams man
#' @export
dgev_p1n_cp=function(x,t,t0=NA,n0=NA,y=x,
	minxi=-1,maxxi=1,extramodels=FALSE,
	rust=FALSE,nrust=1000,centering=TRUE,debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
						is.finite(t),!is.na(t))

	t=ifvectorthenmatrix(t)
	nt=findnt(t)

	for (i in 1:nt){
		t0[i]=maket0(t0[i],n0[i],t[,i])
	}

#
# centering
#
  if(centering){
		meant=matrix(0,nt)
  	for (i in 1:nt){
     meant[i]=mean(t[,i])
     t[,i]=t[,i]-meant[i]
     t0[i]=t0[i]-meant[i]
  	}
  }

	ics=gev_p1n_setics(x,t)
	opt1=optim(ics,gev_p1n_loglik,x=x,t=t,control=list(fnscale=-1)) #this one uses the evd routine for dgev
	vhat=opt1$par
	if(vhat[nt+3]<=(-1)){revert2ml=TRUE}else{revert2ml=FALSE}
	ml_params=vhat
#	gev_p1n_checkmle(ml_params,minxi,maxxi)
	dd=dgev_p1nsub(x=x,t=t,y=y,t0=t0,ics=ics,minxi,maxxi,extramodels=extramodels)
	ru_pdf="rust not selected"
		ml_params=dd$ml_params

	if(rust&&(!revert2ml)){
		th=tgev_p1n_cp(nrust,x,t)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+makebetat0(nt,th[ir,],t0)
			sigma=th[ir,(nt+2)]
			xi=th[ir,(nt+3)]
			dpdf=extraDistr::dgev(y,mu=mu,sigma=sigma,xi=sigma)
			ru_pdf=ru_pdf+dpdf
#			if(is.na(mean(ru_pdf)))stop()
		}
		ru_pdf=ru_pdf/nrust
	} else {
		ru_pdf=dd$ml_pdf
	}
#
# decentering
#
	 if(centering)ml_params[1]=ml_params[1]-makebetat0(nt,ml_params,meant)

		op=list(
					ml_params=ml_params,
					ml_pdf=dd$ml_pdf,
					revert2ml=revert2ml,
					ru_pdf=ru_pdf,
					cp_method=nopdfcdfmsg())

	return(op)

}
#' @rdname gev_p1n_cp
#' @inheritParams man
#' @export
pgev_p1n_cp=function(x,t,t0=NA,n0=NA,y=x,
	minxi=-1,maxxi=1,extramodels=FALSE,
	rust=FALSE,nrust=1000,centering=TRUE,debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
						is.finite(t),!is.na(t))

	t=ifvectorthenmatrix(t)
	nt=findnt(t)

	for (i in 1:nt){
		t0[i]=maket0(t0[i],n0[i],t[,i])
	}
#
# centering
#
  if(centering){
		meant=matrix(0,nt)
  	for (i in 1:nt){
	    meant[i]=mean(t[,i])
	    t[,i]=t[,i]-meant[i]
	  	t0[i]=t0[i]-meant[i]
  	}
  }

	ics=gev_p1n_setics(x,t)
	opt1=optim(ics,gev_p1n_loglik,x=x,t=t,control=list(fnscale=-1)) #this one uses the evd routine for dgev
	vhat=opt1$par
	if(vhat[nt+3]<=(-1)){revert2ml=TRUE}else{revert2ml=FALSE}
	ml_params=vhat
#	gev_p1n_checkmle(ml_params,minxi,maxxi)
	dd=dgev_p1nsub(x=x,t=t,y=y,t0=t0,ics=ics,minxi,maxxi,extramodels=extramodels)
	ru_cdf="rust not selected"
		ml_params=dd$ml_params

	if(rust&&(!revert2ml)){
		th=tgev_p1n_cp(nrust,x,t)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+makebetat0(nt,th[ir,],t0)
			sigma=th[ir,(nt+2)]
			ru_cdf=ru_cdf+pgev(y,mu=mu,sigma=sigma,xi=th[ir,(nt+3)])
		}
		ru_cdf=ru_cdf/nrust
	} else {
		ru_pdf=dd$ml_pdf
	}
#
# decentering
#
	 if(centering)ml_params[1]=ml_params[1]-makebetat0(nt,ml_params,meant)

		op=list(
					ml_params=ml_params,
					ml_cdf=dd$ml_cdf,
					revert2ml=revert2ml,
					ru_cdf=ru_cdf,
					cp_method=nopdfcdfmsg())
	return(op)
}
#' @rdname gev_p1n_cp
#' @inheritParams man
#' @export
tgev_p1n_cp=function(n,x,t,
	extramodels=FALSE,debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),is.finite(t),!is.na(t),
#						length(ics)==5)
	stopifnot(is.finite(x),!is.na(x),is.finite(t),!is.na(t))

	t=ifvectorthenmatrix(t)
	nt=findnt(t)
#
# centering
#
	meant=matrix(0,nt)
	for (i in 1:nt){
	 meant[i]=mean(t[,i])
	 t[,i]=t[,i]-meant[i]
	}

	ics=gev_p1n_setics(x,t)
#	ics=c(ics) #this converts the matrix to a vector (but now it's a vector anyway)
	th=ru(gev_p1n_logf,x=x,t=t,n=n,d=(nt+3),init=ics)
  theta_samples=th$sim_vals
#
# decentering
#
  for (i in 1:n){
	  theta_samples[i,1]=theta_samples[i,1]-makebetat0(nt,theta_samples[i,],meant)
	}

	list(theta_samples=theta_samples)

}
