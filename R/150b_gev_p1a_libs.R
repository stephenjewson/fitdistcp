#' rgev for gev_p1 but with maxlik xi within bounds
#' @return Vector
#' @inheritParams manf
rgev_p1_minmax=function(nx,mu=0,sigma=1,xi=0,tt,minxi=-0.45,maxxi=0.45,centering=TRUE){
	xihat=-9999
  if(centering)tt=tt-mean(tt)
	while((xihat<minxi)||(xihat>maxxi)){ #0.46 also works...0.47 doesn't
		xx=extraDistr::rgev(nx,mu=mu,sigma=sigma,xi=xi)
		ics=gev_p1_setics(xx,tt,c(0,0,0,0))
		opt1=optim(ics,gev_p1_loglik,x=xx,t=tt,control=list(fnscale=-1))
		xihat=opt1$par[4]
	}
	return(xx)
}
#' Waic
#' @inherit manwaic return
#' @inheritParams manf
gev_p1_waic=function(waicscores,x,t0,v1hat,v2hat,v3hat,v4hat,
	lddi,lddd,lambdad){
		if(waicscores){
			f1f=gev_p1_f1fw(x,t0,v1hat,v2hat,v3hat,v4hat)
			f2f=gev_p1_f2fw(x,t0,v1hat,v2hat,v3hat,v4hat)
			fhatx=dgev_p1(x,t0,ymn=v1hat,slope=v2hat,sigma=v3hat,xi=v4hat,log=FALSE)
			waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad,f2f,dim=4)
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
gev_p1_predictordata=function(predictordata,x,t,t0,params){
	if(predictordata){
#
# calculate the probabilities of the data using the fited model
#
		a=params[1]
		b=params[2]
		sc=params[3]
		sh=params[4]
		mu=a+b*t
		px=extraDistr::pgev(x,mu=mu,sigma=sc,xi=sh)
#
# calculate the quantiles for those probabilities at t0
#
		mu0=a+b*t0
		qx=extraDistr::qgev(px,mu=mu0,sigma=sc,xi=sh)
	} else{
		predictedparameter="predictordata not selected"
		adjustedx="predictordata not selected"
	}

	list(predictedparameter=mu,adjustedx=qx)
}#' Logf for RUST
#' @inherit manlogf return
#' @inheritParams manf
gev_p1_logf=function(params,x,t){
#	a=params[1]
#	b=params[2]
#	sc=params[3]
#	sh=params[4]
#	mu=a+b*t
#	if(sc>0){
#		logf=sum(extraDistr::dgev(x,mu=mu,sigma=sc,xi=sh,log=TRUE))-log(sc)
#	}else{
#		logf=-Inf
#	}
	a=params[1]
	b=params[2]
	sc=pmax(params[3],sqrt(.Machine$double.eps))
	sh=params[4]
	mu=a+b*t
	logf=sum(extraDistr::dgev(x,mu=mu,sigma=sc,xi=sh,log=TRUE))-log(sc)
	return(logf)
}
#' Set initial conditions
#' @return Vector
#' @inheritParams manf
gev_p1_setics=function(x,t,ics){
	nx=length(x)
	if((ics[1]==0)&&(ics[2]==0)&&(ics[3]==0)&&(ics[4]==0)){
		lm=lm(x~t)
		ics[1]=lm$coefficients[1]
		ics[2]=lm$coefficients[2]
		xhat=ics[1]+ics[2]*t
		ics[3]=sqrt((sum((x-xhat)^2))/nx)
		ics[4]=0
	}
	return(ics)
}
#'  observed log-likelihood function
#' @inherit manloglik return
#' @inheritParams manf
gev_p1_loglik=function(vv,x,t){
	n=length(x)
	mu=vv[1]+vv[2]*t #so mean is a vector, just like x
	loglik=sum(extraDistr::dgev(x,mu=mu,sigma=max(vv[3],.Machine$double.eps),xi=vv[4],log=TRUE))
	return(loglik)
}
#' Check MLE
#' @return No return value (just a message to the screen).
#' @inheritParams manf
gev_p1_checkmle=function(ml_params,minxi=-1,maxxi=1){
# currently not used, because instead I revert2ml
	v1hat=ml_params[1]
	v2hat=ml_params[2]
	v3hat=ml_params[3]
	v4hat=ml_params[4]
	if(is.na(v1hat))stop()
	if(is.na(v2hat))stop()
	if(is.na(v3hat))stop()
	if(is.na(v4hat))stop()
	if(v4hat<minxi){warning("\n***v4hat=",v4hat,"=> execution halted because maxlik shape parameter <",minxi,"***\n");stop()}
	if(v4hat>maxxi){warning("\n***v4hat=",v4hat,"=> execution halted because maxlik shape parameter >",maxxi,"***\n");stop()}
}
#' GEVD-with-p1: Quantile function
#' @inherit manvector return
#' @inheritParams manf
qgev_p1=function(p,t0,ymn,slope,sigma,xi){

	return(extraDistr::qgev(p,mu=(ymn+slope*t0),sigma=sigma,xi=xi))

}
#' GEVD-with-p1: Density function
#' @inherit manvector return
#' @inheritParams manf
dgev_p1=function(x,t0,ymn,slope,sigma,xi,log=FALSE){

	mu=(ymn+slope*t0)
	return(extraDistr::dgev(x,mu=mu,sigma=sigma,xi=xi,log=log))

}
#' GEVD-with-p1: Distribution function
#' @inherit manvector return
#' @inheritParams manf
pgev_p1=function(y,t0,ymn,slope,sigma,xi){

	return(extraDistr::pgev(y,mu=(ymn+slope*t0),sigma=sigma,xi=xi))

}
#' Analytical expressions for Predictive Means
#' RHP mean based on the expectation of DMGS equation 2.1
#' @inherit manmeans return
#' @inheritParams manf
gev_p1_means=function(means,t0,ml_params,lddi,lddd,
									lambdad_rh_flat,nx,dim=4){

	if(means){
# intro
		eulerconstant=0.57721566490153286060651209008240243104215933593992
		ymn=ml_params[1]
		slope=ml_params[2]
		sigma=ml_params[3]
		xi=ml_params[4]

# set up derivative arrays
		meand1=array(0,c(4,1))
		meand2=array(0,c(4,4,1)) #but all zero for gumbel

		if(ml_params[4]==0){
# xi=0 case

# mle
			ml_mean=ymn+slope*t0+sigma*eulerconstant
# calculate first derivative array for bayesian xi=0 cases
			meand1[1,1]=1
			meand1[2,1]=t0
			meand1[3,1]=eulerconstant
			meand1[4,1]=0
# meand2 is all zero as initialized

		} else{
# non-gumbel case

			g0=gamma(1-xi)
			g1=g0*digamma(1-xi)
			g2=(trigamma(1-xi)*g0*g0+g1*g1)/g0
# mle
			ml_mean=ymn+slope*t0+sigma*(g0-1)/xi
# calculate first derivative array for bayesian xi!=0 cases
			meand1[1,1]=1
			meand1[2,1]=t0
			meand1[3,1]=(g0-1)/xi
			meand1[4,1]=(1-g0-xi*g1)*sigma/(xi*xi)
# calculate second derivative array (only has 1 non-zero term!)
			meand2[3,4,1]=(1-g0-xi*g1)/(xi*xi)
			meand2[4,3,1]=meand2[3,4,1]
			meand2[4,4,1]=sigma*(-2+2*g0+2*xi*g1+xi*xi*g2)/(xi*xi*xi)
		}
# I've realized now that when I integrate over xi, the mean in Inf

		rh_flat_mean	=Inf


	}else{
		ml_mean="means not selected"
		rh_flat_mean="means not selected"
	}

# return
	list(ml_mean=ml_mean,rh_flat_mean=rh_flat_mean)
}
#' Densities for 5 predictions
#' @inherit mandsub return
#' @inheritParams manf
dgev_p1sub=function(x,t,y,t0,ics,minxi,maxxi,extramodels=FALSE){

		nx=length(x)

		ics=gev_p1_setics(x,t,ics)
		opt=optim(ics,gev_p1_loglik,x=x,t=t,control=list(fnscale=-1))
		v1hat=opt$par[1]
		v2hat=opt$par[2]
		v3hat=opt$par[3]
		v4hat=min(maxxi,max(minxi,opt$par[4]))
		ml_params=c(v1hat,v2hat,v3hat,v4hat)

# now that I've dropped dmgs, not sure I need this anymore
#		muhat0=v1hat+v2hat*t0
#		y=fixgevrange(y,muhat0,v3hat,v4hat)

# mle
		ml_pdf=dgev_p1(y,t0,ymn=v1hat,slope=v2hat,sigma=v3hat,xi=v4hat)
		ml_cdf=pgev_p1(y,t0,ymn=v1hat,slope=v2hat,sigma=v3hat,xi=v4hat)

# return
		list(
					ml_params=ml_params,
					ml_pdf=ml_pdf,
					ml_cdf=ml_cdf)
}

