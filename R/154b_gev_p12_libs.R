#
# NOTE THAT PARAMETER 4 HAS A DEFAULT VALUE OF ZERO AND RANGE FROM -INF TO INF
# BECAUSE ITS INSIDE AN EXPONENTIAL
#
#' rgev for gev_p12 but with maxlik xi within bounds
#' @return Vector
#' @inheritParams manf
rgev_p12_minmax=function(nx,mu=0,sigma=1,xi=0,t1,t2,minxi=-0.45,maxxi=0.45,centering=TRUE){
	xihat=-9999
  if(centering){
  	t1=t1-mean(t1)
  	t2=t2-mean(t2)
  }
	while((xihat<minxi)||(xihat>maxxi)){ #0.46 also works...0.47 doesn't
		xx=extraDistr::rgev(nx,mu=mu,sigma=sigma,xi=xi)
		ics=gev_p12_setics(xx,t1,t2,c(0,0,0,0,0))
		opt1=optim(ics,gev_p12_loglik,x=xx,t1=t1,t2=t2,control=list(fnscale=-1))
		xihat=opt1$par[5]
	}
	return(xx)
}
#' Waic
#' @inherit manwaic return
#' @inheritParams manf
gev_p12_waic=function(waicscores,x,t1,t2,v1hat,v2hat,v3hat,v4hat,v5hat,
	lddi,lddd,lambdad){
		if(waicscores){
			f1f=gev_p12_f1fw(x,t1,t2,v1hat,v2hat,v3hat,v4hat,v5hat)
			f2f=gev_p12_f2fw(x,t1,t2,v1hat,v2hat,v3hat,v4hat,v5hat)
			fhatx=dgev_p12(x,t1,t2,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat,xi=v5hat,log=FALSE)
			waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad,f2f,dim=5)
			waic1=waic$waic1
			waic2=waic$waic2
		}else{
			waic1="waicscores not selected"
			waic2="waicscores not selected"
		}
		list(waic1=waic1,waic2=waic2)
}
#' Predicted Parameter and Generalized Residuals
#' @inherit manpredictor return
#' @inheritParams manf
gev_p12_predictordata=function(predictordata,x,t1,t2,t01,t02,params){
	if(predictordata){
#
# calculate the probabilities of the data using the fited model
#
# t is always a matrix
# t0 is always a vector
		a=params[1]
		b=params[2]
		sc1=params[3]
		sc2=params[4]
		sh=params[5]
		mu=a+b*t1
		sigma=exp(sc1+sc2*t2)
		px=extraDistr::pgev(x,mu=mu,sigma=sigma,xi=sh)
#
# calculate the quantiles for those probabilities at t01,t02
#
		mu0=a+b*t01
		sigma0=exp(sc1+sc2*t02)
		qx=extraDistr::qgev(px,mu=mu0,sigma=sigma0,xi=sh)
	} else{
		predictedparameter="predictordata not selected"
		adjustedx="predictordata not selected"
	}


	list(predictedparameter=mu,adjustedx=qx)
}#' Logf for RUST
#' @inherit manlogf return
#' @inheritParams manf
gev_p12_logf=function(params,x,t1,t2){
	a=params[1]
	b=params[2]
	sc1=params[3]
	sc2=params[4]
	sh=params[5]
#	if(is.vector(t)){
#		mu=a+b*t1
#		sigma=exp(sc1+sc2*t2)
#	} else {
#		mu=a+b*t1
#		sigma=exp(sc1+sc2*t2)
#	}
	mu=a+b*t1
	sigma=exp(sc1+sc2*t2)
	logf=sum(extraDistr::dgev(x,mu=mu,sigma=sigma,xi=sh,log=TRUE))
	return(logf)
}#' Set initial conditions
#' @return Vector
#' @inheritParams manf
gev_p12_setics=function(x,t1,t2,ics){
# t is always a matrix
	nx=length(x)
	if((ics[1]==0)&&(ics[2]==0)&&(ics[3]==0)&&(ics[4]==0)&&(ics[5]==0)){
		lm=lm(x~t1)
		ics[1]=lm$coefficients[1]
		ics[2]=lm$coefficients[2]
		xhat=ics[1]+ics[2]*t1
		ics[3]=0
		ics[4]=0 #zero because it's inside an exponential
		ics[5]=0
	}
	return(ics)
}
#'  observed log-likelihood function
#' @inherit manloglik return
#' @inheritParams manf
gev_p12_loglik=function(vv,x,t1,t2){
# t is always a matrix
	n=length(x)
	mu=vv[1]+vv[2]*t1 #so mean is a vector, just like x
	sigma=exp(vv[3]+vv[4]*t2)
	loglik=sum(extraDistr::dgev(x,mu=mu,sigma=sigma,xi=vv[5],log=TRUE))
	return(loglik)
}
#' Check MLE
#' @inherit mancheckmle return
#' @inheritParams manf
gev_p12_checkmle=function(ml_params,minxi=-1,maxxi=1){
# currently not used, because instead I revert2ml
	v1hat=ml_params[1]
	v2hat=ml_params[2]
	v3hat=ml_params[3]
	v4hat=ml_params[4]
	v5hat=ml_params[5]
	if(is.na(v1hat))stop()
	if(is.na(v2hat))stop()
	if(is.na(v3hat))stop()
	if(is.na(v4hat))stop()
	if(is.na(v5hat))stop()
	if(v5hat<minxi){warning("\n***v5hat=",v5hat,"=> execution halted because maxlik shape parameter <",minxi,"***");stop()}
	if(v5hat>maxxi){warning("\n***v5hat=",v5hat,"=> execution halted because maxlik shape parameter >",maxxi,"***");stop()}
}
#' GEVD-with-p1: Quantile function
#' @inherit manvector return
#' @inheritParams manf
qgev_p12=function(p,t1,t2,ymn,slope,sigma1,sigma2,xi){
# t is sometimes a vector, sometimes a matrix

#	if(is.vector(t)){
#		mu=ymn+slope*t1
#		sigma=exp(sigma1+sigma2*t2)
#	} else {
#		mu=ymn+slope*t1
#		sigma=exp(sigma1+sigma2*t2)
#	}
	mu=ymn+slope*t1
	sigma=exp(sigma1+sigma2*t2)

	return(extraDistr::qgev(p,mu=mu,sigma=sigma,xi=xi))

}
#' GEVD-with-p1: Density function
#' @inherit manvector return
#' @inheritParams manf
dgev_p12=function(x,t1,t2,ymn,slope,sigma1,sigma2,xi,log=FALSE){
# t is sometimes a vector, sometimes a matrix

#	if(is.vector(t)){
#		mu=ymn+slope*t1
#		sigma=exp(sigma1+sigma2*t2)
#	} else {
#		mu=ymn+slope*t1
#		sigma=exp(sigma1+sigma2*t2)
#	}
		mu=ymn+slope*t1
		sigma=exp(sigma1+sigma2*t2)

	return(extraDistr::dgev(x,mu=mu,sigma=sigma,xi=xi,log=log))

}
#' GEVD-with-p1: Distribution function
#' @inherit manvector return
#' @inheritParams manf
pgev_p12=function(y,t1,t2,ymn,slope,sigma1,sigma2,xi){
# t is sometimes a vector, sometimes a matrix

#	if(is.vector(t)){
#		mu=ymn+slope*t1
#		sigma=exp(sigma1+sigma2*t2)
#	} else {
#		mu=ymn+slope*t1
#		sigma=exp(sigma1+sigma2*t2)
#	}
		mu=ymn+slope*t1
		sigma=exp(sigma1+sigma2*t2)

	return(extraDistr::pgev(y,mu=mu,sigma=sigma,xi=xi))

}
#' Analytical expressions for Predictive Means
#' RHP mean based on the expectation of DMGS equation 2.1
#' @inherit manmeans return
#' @inheritParams manf
gev_p12_means=function(means,t01,t02,ml_params,nx){

	if(means){
# intro
		eulerconstant=0.57721566490153286060651209008240243104215933593992
		ymn=ml_params[1]
		slope=ml_params[2]
		sigma1=ml_params[3]
		sigma2=ml_params[4]
		sigma=exp(sigma1+sigma2*t02)
		xi=ml_params[5]

		if(ml_params[5]==0){
# xi=0 case
			ml_mean=ymn+slope*t01+sigma*eulerconstant
		} else{
# xi!=0 case
			g0=gamma(1-xi)
			g1=g0*digamma(1-xi)
			g2=(trigamma(1-xi)*g0*g0+g1*g1)/g0
			ml_mean=ymn+slope*t01+sigma*(g0-1)/xi
		}
# return
		crhp_mle_mean="haven't worked it out yet"
		pu_mean=Inf
	}else{
		pu_mean="means not selected"
		ml_mean="means not selected"
		crhp_mle_mean="means not selected"
	}
	list(ml_mean=ml_mean,crhp_mle_mean=crhp_mle_mean,pu_mean=pu_mean)
}
#' Densities for 5 predictions
#' @inherit mandsub return
#' @inheritParams manf
dgev_p12sub=function(x,t1,t2,y,t01,t02,ics,
	minxi,maxxi,debug,extramodels=FALSE){

		nx=length(x)

		ics=gev_p12_setics(x,t1,t2,ics)
		opt=optim(ics,gev_p12_loglik,x=x,t1=t1,t2=t2,control=list(fnscale=-1))
		v1hat=opt$par[1]
		v2hat=opt$par[2]
		v3hat=opt$par[3]
		v4hat=opt$par[4]
		v5hat=min(maxxi,max(minxi,opt$par[5]))
		ml_params=c(v1hat,v2hat,v3hat,v4hat,v5hat)

# now that I've dropped dmgs d and p, I don't think I need this any more
#		muhat0=v1hat+v2hat*t01
#		sghat0=exp(v3hat+v4hat*t02)
#		y=fixgevrange(y,muhat0,sghat0,v5hat)

# mle
		ml_pdf=dgev_p12(y,t01,t02,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat,xi=v5hat)
		ml_cdf=pgev_p12(y,t01,t02,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat,xi=v5hat)

# return
		list(ml_params=ml_params,
					ml_pdf=ml_pdf,
					ml_cdf=ml_cdf)
}

