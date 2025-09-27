#' rgev for gev_p1b but with maxlik xi within bounds
#' @return Vector
#' @inheritParams manf
rgev_p1b_minmax=function(nx,mu=0,sigma=1,xi=0,tt,minxi=-0.45,maxxi=0.45,centering=TRUE){
	xihat=-9999
  if(centering){
  	tt[,1]=tt[,1]-mean(tt[,1])
  	tt[,2]=tt[,2]-mean(tt[,2])
  }
	while((xihat<minxi)||(xihat>maxxi)){ #0.46 also works...0.47 doesn't
		xx=extraDistr::rgev(nx,mu=mu,sigma=sigma,xi=xi)
		ics=gev_p1b_setics(xx,tt,c(0,0,0,0,0))
		opt1=optim(ics,gev_p1b_loglik,x=xx,t=tt,control=list(fnscale=-1))
		xihat=opt1$par[5]
	}
	return(xx)
}
#' Waic
#' @inherit manwaic return
#' @inheritParams manf
gev_p1b_waic=function(waicscores,x,t0,v1hat,v2hat,v3hat,v4hat,v5hat,
	lddi,lddd,lambdad){
		if(waicscores){
			f1f=gev_p1b_f1fw(x,t0[1],t0[2],v1hat,v2hat,v3hat,v4hat,v5hat)
			f2f=gev_p1b_f2fw(x,t0[1],t0[2],v1hat,v2hat,v3hat,v4hat,v5hat)
			fhatx=dgev_p1b(x,t0,ymn=v1hat,slope1=v2hat,slope2=v3hat,sigma=v4hat,xi=v5hat,log=FALSE)
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
gev_p1b_predictordata=function(predictordata,x,t,t0,params){
	if(predictordata){
#
# calculate the probabilities of the data using the fitted model
#
		a=params[1]
		b1=params[2]
		b2=params[3]
		sc=params[4]
		sh=params[5]
		mu=a+b1*t[,1]+b2*t[,2]
		px=extraDistr::pgev(x,mu=mu,sigma=sc,xi=sh)
#
# calculate the quantiles for those probabilities at t0
#
		mu0=a+b1*t0[1]+b2*t0[2]
		qx=extraDistr::qgev(px,mu=mu0,sigma=sc,xi=sh)
	} else{
		predictedparameter="predictordata not selected"
		adjustedx="predictordata not selected"
	}

	list(predictedparameter=mu,adjustedx=qx)
}#' Logf for RUST
#' @inherit manlogf return
#' @inheritParams manf
gev_p1b_logf=function(params,x,t){
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
	b1=params[2]
	b2=params[3]
	sc=pmax(params[4],sqrt(.Machine$double.eps))
	sh=params[5]
	mu=a+b1*t[1,]+b2*t[,2]
	logf=sum(extraDistr::dgev(x,mu=mu,sigma=sc,xi=sh,log=TRUE))-log(sc)
	return(logf)
}
#' Set initial conditions
#' @return Vector
#' @inheritParams manf
gev_p1b_setics=function(x,t,ics){
	nx=length(x)
	if((ics[1]==0)&&(ics[2]==0)&&(ics[3]==0)&&(ics[4]==0)&&(ics[4]==0)){
		lm=lm(x~t)
		ics[1]=lm$coefficients[1]
		ics[2]=lm$coefficients[2]
		ics[3]=lm$coefficients[3]
		xhat=ics[1]+ics[2]*t[,1]+ics[3]*t[,2]
		ics[4]=sqrt((sum((x-xhat)^2))/nx)
		ics[5]=0
	}
	return(ics)
}
#'  observed log-likelihood function
#' @inherit manloglik return
#' @inheritParams manf
gev_p1b_loglik=function(vv,x,t){
	n=length(x)
	mu=vv[1]+vv[2]*t[,1]+vv[3]*t[,2] #so mean is a vector, just like x
	loglik=sum(extraDistr::dgev(x,mu=mu,sigma=max(vv[4],.Machine$double.eps),xi=vv[5],log=TRUE))
	return(loglik)
}
#' Check MLE
#' @return No return value (just a message to the screen).
#' @inheritParams manf
gev_p1b_checkmle=function(ml_params,minxi=-1,maxxi=1){
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
	if(v5hat<minxi){warning("\n***v5hat=",v5hat,"=> execution halted because maxlik shape parameter <",minxi,"***\n");stop()}
	if(v5hat>maxxi){warning("\n***v5hat=",v5hat,"=> execution halted because maxlik shape parameter >",maxxi,"***\n");stop()}
}
#' GEVD-with-p1: Quantile function
#' @inherit manvector return
#' @inheritParams manf
qgev_p1b=function(p,t0,ymn,slope1,slope2,sigma,xi){

	return(extraDistr::qgev(p,mu=(ymn+slope1*t0[1]+slope2*t0[2]),sigma=sigma,xi=xi))

}
#' GEVD-with-p1: Density function
#' @inherit manvector return
#' @inheritParams manf
dgev_p1b=function(x,t0,ymn,slope1,slope2,sigma,xi,log=FALSE){

	mu=(ymn+slope1*t0[1]+slope2*t0[2])
	return(extraDistr::dgev(x,mu=mu,sigma=sigma,xi=xi,log=log))

}
#' GEVD-with-p1: Distribution function
#' @inherit manvector return
#' @inheritParams manf
pgev_p1b=function(y,t0,ymn,slope1,slope2,sigma,xi){

	return(extraDistr::pgev(y,mu=(ymn+slope1*t0[1]+slope2*t0[2]),sigma=sigma,xi=xi))

}
#' Analytical expressions for Predictive Means
#' RHP mean based on the expectation of DMGS equation 2.1
#' @inherit manmeans return
#' @inheritParams manf
gev_p1b_means=function(means,t0,ml_params,lddi,lddd,
									lambdad_rh_flat,nx,dim=5){

	if(means){
# intro
		eulerconstant=0.57721566490153286060651209008240243104215933593992
		ymn=ml_params[1]
		slope1=ml_params[2]
		slope2=ml_params[3]
		sigma=ml_params[4]
		xi=ml_params[5]

# set up derivative arrays
		meand1=array(0,c(5,1))
		meand2=array(0,c(5,5,1)) #but all zero for gumbel

		if(ml_params[5]==0){
# xi=0 case

# mle
			ml_mean=ymn+slope1*t0[1]+slope2*t0[2]+sigma*eulerconstant
# calculate first derivative array for bayesian xi=0 cases
			meand1[1,1]=1
			meand1[2,1]=t0[1]
			meand1[3,1]=t0[2]
			meand1[4,1]=eulerconstant
			meand1[5,1]=0
# meand2 is all zero as initialized

		} else{
# non-gumbel case

			g0=gamma(1-xi)
			g1=g0*digamma(1-xi)
			g2=(trigamma(1-xi)*g0*g0+g1*g1)/g0
# mle
			ml_mean=ymn+slope1*t0[1]+slope2*t0[2]+sigma*(g0-1)/xi
# calculate first derivative array for bayesian xi!=0 cases
			meand1[1,1]=1
			meand1[2,1]=t0[1]
			meand1[3,1]=t0[2]
			meand1[4,1]=(g0-1)/xi
			meand1[5,1]=(1-g0-xi*g1)*sigma/(xi*xi)
# calculate second derivative array (only has 1 non-zero term!)
			meand2[4,5,1]=(1-g0-xi*g1)/(xi*xi)
			meand2[5,4,1]=meand2[4,5,1]
			meand2[5,5,1]=sigma*(-2+2*g0+2*xi*g1+xi*xi*g2)/(xi*xi*xi)
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
dgev_p1bsub=function(x,t,y,t0,ics,minxi,maxxi,extramodels=FALSE){

		nx=length(x)

		ics=gev_p1b_setics(x,t,ics)
		opt=optim(ics,gev_p1b_loglik,x=x,t=t,control=list(fnscale=-1))
		v1hat=opt$par[1]
		v2hat=opt$par[2]
		v3hat=opt$par[3]
		v4hat=opt$par[4]
		v5hat=min(maxxi,max(minxi,opt$par[5]))
		ml_params=c(v1hat,v2hat,v3hat,v4hat,v5hat)

# now that I've dropped dmgs, not sure I need this anymore
#		muhat0=v1hat+v2hat*t0
#		y=fixgevrange(y,muhat0,v3hat,v4hat,v5hat)

# mle
		ml_pdf=dgev_p1b(y,t0,ymn=v1hat,slope1=v2hat,slope2=v3hat,sigma=v4hat,xi=v5hat)
		ml_cdf=pgev_p1b(y,t0,ymn=v1hat,slope1=v2hat,slope2=v3hat,sigma=v4hat,xi=v5hat)

# return
		list(
					ml_params=ml_params,
					ml_pdf=ml_pdf,
					ml_cdf=ml_cdf)
}

