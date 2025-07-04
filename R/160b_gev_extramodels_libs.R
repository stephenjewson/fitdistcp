#' rgev but with maxlik xi guaranteed within bounds
#' @return Vector
#' @inheritParams manf
rgev_minmax=function(nx,mu,sigma,xi,minxi=-0.45,maxxi=0.45){
	xihat=-9999
	while((xihat<minxi)||(xihat>maxxi)){
		xx=extraDistr::rgev(nx,mu=mu,sigma=sigma,xi=xi)
		ics=gev_setics(xx,c(0,0,0))
		opt1=optim(ics,gev_loglik,x=xx,control=list(fnscale=-1))
		xihat=opt1$par[3]
	}
	return(xx)
}
#' Waic
#' @inherit manwaic return
#' @inheritParams manf
gev_waic=function(waicscores,x,v1hat,d1,v2hat,fd2,v3hat,d3,lddi,lddd,
	lambdad){
		if(waicscores){

			f1f=gev_f1fa(x,v1hat,v2hat,v3hat)
			f2f=gev_f2fa(x,v1hat,v2hat,v3hat)

			fhatx=extraDistr::dgev(x,mu=v1hat,sigma=v2hat,xi=v3hat)
			waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad,f2f,dim=3)
			waic1=waic$waic1
			waic2=waic$waic2
		}else{
			waic1="waicscores not selected"
			waic2="waicscores not selected"
		}
		list(waic1=waic1,waic2=waic2)
}
#' Logf for RUST
#' @inherit manlogf return
#' @inheritParams manf
gev_logf=function(params,x){
	mu=params[1]
	sc=pmax(params[2],.Machine$double.eps)
	sh=params[3]
	logf=sum(extraDistr::dgev(x,mu=mu,sigma=sc,xi=sh,log=TRUE))-log(sc)
	return(logf)
}
#' Set initial conditions
#' @return Vector
#' @inheritParams manf
gev_setics=function(x,ics){
	if((ics[1]==0)&&(ics[2]==0)&&(ics[3]==0)){
		ics[1]=mean(x)
		ics[2]=sd(x)
		ics[3]=0
	}
	return(ics)
}
#' PWM parameter estimation
#' @return Vector
#' @inheritParams manf
gev_pwm_params=function(x){
		pw_params=matrix(0,3)
		g=fExtremes::gevFit(x,type="pwm")
		temp=g@fit$par.ests
# reorder to match the mu, sigma, xi ordering that I use
		pw_params[1]=temp[2]
		pw_params[2]=temp[3]
		pw_params[3]=temp[1]
return(pw_params)
}
#'  log-likelihood function
#' @inherit manloglik return
#' @inheritParams manf
gev_loglik=function(vv,x){
	n=length(x)
	loglik=sum(extraDistr::dgev(x,mu=vv[1],sigma=max(vv[2],.Machine$double.eps),xi=vv[3],log=TRUE))
	return(loglik)
}
#' Check MLE
#' @return No return value (just a message to the screen).
#' @inheritParams manf
gev_checkmle=function(ml_params,minxi,maxxi){
	v1hat=ml_params[1]
	v2hat=ml_params[2]
	v3hat=ml_params[3]
	if(is.na(v1hat))stop()
	if(is.na(v2hat))stop()
	if(is.na(v3hat))stop()
	if(v3hat<minxi){warning("\n***v3hat=",v3hat,"=> execution halted because maxlik shape parameter <",minxi,"***\n");stop()}
	if(v3hat>maxxi){warning("\n***v3hat=",v3hat,"=> execution halted because maxlik shape parameter >",maxxi,"***\n");stop()}
}
#' Derivative of expected information matrix, based on MEV routine gev.infomat
#' @inherit manldd return
#' @inheritParams manf
gev_ggd_mev=function(v1,d1,v2,fd2,v3,d3){
	ggd=array(0,c(3,3,3))
  d2=fd2*v2

	ggm1=gev.infomat(c(v1-d1,v2,v3),dat=c(1),method=c("exp"))
	ggp1=gev.infomat(c(v1+d1,v2,v3),dat=c(1),method=c("exp"))
	ggd[1,,]=(ggp1-ggm1)/(2*d1)

	ggm2=gev.infomat(c(v1,v2-d2,v3),dat=c(1),method=c("exp"))
	ggp2=gev.infomat(c(v1,v2+d2,v3),dat=c(1),method=c("exp"))
	ggd[2,,]=(ggp2-ggm2)/(2*d2)

	ggm3=gev.infomat(c(v1,v2,v3-d3),dat=c(1),method=c("exp"))
	ggp3=gev.infomat(c(v1,v2,v3+d3),dat=c(1),method=c("exp"))
	ggd[3,,]=(ggp3-ggm3)/(2*d3)

  return(ggd)
}
#' Derivative of inverse expected information matrix, based on MEV routine gev.infomat
#' @inherit manlddd return
#' @inheritParams manf
gev_ggid_mev=function(v1,d1,v2,fd2,v3,d3){
	ggid=array(0,c(3,3,3))
  d2=fd2*v2

	ggm1=gev.infomat(c(v1-d1,v2,v3),dat=c(1),method=c("exp"))
	ggp1=gev.infomat(c(v1+d1,v2,v3),dat=c(1),method=c("exp"))
	ggim1=solve(ggm1)
	ggip1=solve(ggp1)
	ggid[1,,]=(ggip1-ggim1)/(2*d1)

	ggm2=gev.infomat(c(v1,v2-d2,v3),dat=c(1),method=c("exp"))
	ggp2=gev.infomat(c(v1,v2+d2,v3),dat=c(1),method=c("exp"))
	ggim2=solve(ggm2)
	ggip2=solve(ggp2)
	ggid[2,,]=(ggip2-ggim2)/(2*d2)

	ggm3=gev.infomat(c(v1,v2,v3-d3),dat=c(1),method=c("exp"))
	ggp3=gev.infomat(c(v1,v2,v3+d3),dat=c(1),method=c("exp"))
	ggim3=solve(ggm3)
	ggip3=solve(ggp3)
	ggid[3,,]=(ggip3-ggim3)/(2*d3)

  return(ggid)
}
#' Analytical Expressions for Predictive Means
#' RHP mean based on the expectation of DMGS equation 2.1
#' @inherit manmeans return
#' @inheritParams manf
gev_means=function(means,ml_params,lddi,lddi_k3,lddd,lddd_k3,
									lambdad_flat,lambdad_rh_mle,
									lambdad_rh_flat,lambdad_jp,lambdad_custom,
									nx,dim=3){

	if(means){
# intro
		eulerconstant=0.57721566490153286060651209008240243104215933593992
		mu=ml_params[1]
		sigma=ml_params[2]
		xi=ml_params[3]

# set up derivative arrays
		meand1=array(0,c(3,1))
		meand2=array(0,c(3,3,1)) #but all zero for gumbel
		meand1_k3=array(0,c(2,1))
		meand2_k3=array(0,c(2,2,1)) #but all zero for gumbel

		if(ml_params[3]==0){
# xi=0 case

# mle
			ml_mean=mu+sigma*eulerconstant
# calculate first derivative array for bayesian xi=0 cases
			meand1[1,1]=1
			meand1[2,1]=eulerconstant
# and copy for the 2d case
			meand1_k3[1,1]=meand1[1,1]
			meand1_k3[2,1]=meand1[2,1]
# meand2 is all zero as initialized

		} else{
# non-gumbel case

			g0=gamma(1-xi)
			g1=g0*digamma(1-xi)
			g2=(trigamma(1-xi)*g0*g0+g1*g1)/g0
# mle
			ml_mean=mu+sigma*(g0-1)/xi
# calculate first derivative array for bayesian xi!=0 cases
			meand1[1,1]=1
			meand1[2,1]=(g0-1)/xi
			meand1[3,1]=(1-g0-xi*g1)*sigma/(xi*xi)
# and copy for the 2d case
			meand1_k3[1,1]=meand1[1,1]
			meand1_k3[2,1]=meand1[2,1]
# calculate second derivative array (only has 1 non-zero term!)
			meand2[2,3,1]=(1-g0-xi*g1)/(xi*xi)
			meand2[3,2,1]=meand2[2,3,1]
			meand2[3,3,1]=sigma*(-2+2*g0+2*xi*g1+xi*xi*g2)/(xi*xi*xi)
		}
# I've realized now that when we integrate over xi, the mean is Inf

	flat_mean			=Inf
	rh_ml_mean		=ml_mean+dmgs(lddi_k3,lddd_k3,meand1_k3,lambdad_rh_mle,meand2_k3,dim=2)/nx
	rh_flat_mean	=Inf
	jp_mean				=Inf
	custom_mean		=Inf

	}else{
		flat_mean="means not selected"
		ml_mean="means not selected"
		rh_ml_mean="means not selected"
		rh_flat_mean="means not selected"
		jp_mean="means not selected"
		custom_mean="means not selected"
	}

# return
	list(ml_mean=ml_mean,flat_mean=flat_mean,rh_ml_mean=rh_ml_mean,
				rh_flat_mean=rh_flat_mean,jp_mean=jp_mean,custom_mean=custom_mean)
}
#' Densities for 5 predictions
#' @inherit mandsub return
#' @inheritParams manf
dgevsub=function(x,y,ics,d1=0.01,fd2=0.01,d3=0.01,customprior,
	minxi,maxxi,extramodels=FALSE){

		nx=length(x)

		ics=gev_setics(x,ics)
		opt=optim(ics,gev_loglik,x=x,control=list(fnscale=-1))
		v1hat=opt$par[1]
		v2hat=opt$par[2]
		v3hat=opt$par[3]
		ml_params=c(v1hat,v2hat,v3hat)

# mle
		ml_pdf=extraDistr::dgev(y,mu=v1hat,sigma=v2hat,xi=v3hat)
		ml_cdf=extraDistr::pgev(y,mu=v1hat,sigma=v2hat,xi=v3hat)

# return
		list(
					ml_params=ml_params,
					ml_pdf=ml_pdf,
					ml_cdf=ml_cdf)
}

