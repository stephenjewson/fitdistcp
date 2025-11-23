#' GEV_p1n n=1 example data
#' @returns A list containing data to run an example
#' @param  iseed The random seed
gev_p1n_n1_exampledata=function(iseed){
	set.seed(iseed)
	nx		=100
	t			=c(1:nx)/nx
#	t=t-mean(t)
	n0		=nx
	t0		=t[nx]
	params=c(0,1,1,0.1)
	x			=rgev_p1n_minmax(nx,mu=params[1]+params[2]*t,sigma=params[3],xi=params[4],t,minxi=-0.45,maxxi=0.45)
	y			=params[1]+params[2]*t #shouldn't be sorted because the x are ordered
#	x=x-mean(x)
	return(list(nx=nx,x=x,t=t,n0=n0,t0=t0,params=params,y=y))
}
#' GEV_p1n n=2 example data
#' @returns A list containing data to run an example
#' @param  iseed The random seed
gev_p1n_n2_exampledata=function(iseed){
	set.seed(iseed)
	nx=100
	t=matrix(0,nx,2)
	t[,1]=c(1:nx)/nx
	t[,2]=sample(t[,1])
#	t[,1]=t[,1]-mean(t[,1])
#	t[,2]=t[,2]-mean(t[,2])
	n0=c(nx,nx)
	t0=t[nx,]
#	params=c(0,1,1,0.01,0.1)
	params=c(0,1,0,1,0.1)
	x=rgev_p1n_minmax(nx,
		mu=params[1]+params[2]*t[,1]+params[3]*t[,2],
		sigma=params[4],xi=params[5],t,minxi=-0.45,maxxi=0.45)
#	x=x-mean(x)
	y=params[1]+params[2]*t[,1]+params[3]*t[,2] #shouldn't be sorted because the x are ordered
	return(list(nx=nx,x=x,t=t,n0=n0,t0=t0,params=params,y=y))
}
#' rgev for gev_p1n but with maxlik xi within bounds
#' @return Vector
#' @inheritParams manf
rgev_p1n_minmax=function(nx,mu=0,sigma=1,xi=0,tt,minxi=-0.45,maxxi=0.45,centering=TRUE){
	xihat=-9999
	tt=ifvectorthenmatrix(tt)
	nt=findnt(tt)
  if(centering){
		for (i in 1:nt){
	  	tt[,i]=tt[,i]-mean(tt[,i])
		}
  }
	while((xihat<minxi)||(xihat>maxxi)){ #0.46 also works...0.47 doesn't
		xx=extraDistr::rgev(nx,mu=mu,sigma=sigma,xi=xi)
		ics=gev_p1n_setics(xx,tt)
#	stop()
		opt1=optim(ics,gev_p1n_loglik,x=xx,t=tt,control=list(fnscale=-1))
		xihat=opt1$par[nt+3]
	}
	return(xx)
}
#' Waic
#' @inherit manwaic return
#' @inheritParams manf
gev_p1n_waic=function(waicscores,x,t0,vhat,
	lddi,lddd,lambdad){

		nt=findnt(t0)
		if(waicscores){

			if(nt==1)f1f=gev_p1a_f1fw(x,t0[1],							vhat[1],vhat[2],vhat[3],vhat[4])
			if(nt==2)f1f=gev_p1b_f1fw(x,t0[1],t0[2],				vhat[1],vhat[2],vhat[3],vhat[4],vhat[5])
			if(nt==3)f1f=gev_p1c_f1fw(x,t0[1],t0[2],t0[3],	vhat[1],vhat[2],vhat[3],vhat[4],vhat[5],vhat[6])

			if(nt==1)f2f=gev_p1a_f2fw(x,t0[1],							vhat[1],vhat[2],vhat[3],vhat[4])
			if(nt==2)f2f=gev_p1b_f2fw(x,t0[1],t0[2],				vhat[1],vhat[2],vhat[3],vhat[4],vhat[5])
			if(nt==3)f2f=gev_p1c_f2fw(x,t0[1],t0[2],t0[3],	vhat[1],vhat[2],vhat[3],vhat[4],vhat[5],vhat[6])

			fhatx=dgev_p1n(x,t0,params=vhat,log=FALSE)
			waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad,f2f,dim=(nt+3))
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
gev_p1n_predictordata=function(predictordata,x,t,t0,params){

	if(predictordata){
#
# calculate the probabilities of the data using the fitted model
#
		nt=findnt(t)
		a=params[1]
		sc=params[nt+2]
		sh=params[nt+3]
		mu=a+makebetatm(nt,params,t)
		px=extraDistr::pgev(x,mu=mu,sigma=sc,xi=sh)
#
# calculate the quantiles for those probabilities at t0
#
		mu0=a+makebetat0(nt,params,t0)
		qx=extraDistr::qgev(px,mu=mu0,sigma=sc,xi=sh)
	} else{
		predictedparameter="predictordata not selected"
		adjustedx="predictordata not selected"
	}

	list(predictedparameter=mu,adjustedx=qx)
}#' Logf for RUST
#' @inherit manlogf return
#' @inheritParams manf
gev_p1n_logf=function(params,x,t){
#	t=matrix(t)
	nt=findnt(t)
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
#	sc=pmax(params[nt+2],sqrt(.Machine$double.eps))
	mu=params[1]+makebetatm(nt,params,t)
	sigma=max(params[nt+2],sqrt(.Machine$double.eps))
	xi=params[nt+3]
	logf=sum(extraDistr::dgev(x,mu=mu,sigma=sigma,xi=xi,log=TRUE))-log(sigma)
	return(logf)
}
#' Set initial conditions
#' @return Vector
#' @inheritParams manf
gev_p1n_setics=function(x,t){
	nx=length(x)
	nt=findnt(t)
#	if((ics[1]==0)&&(ics[2]==0)&&(ics[3]==0)&&(ics[4]==0)&&(ics[4]==0)){
		lm=lm(x~t)
		ics=numeric(nt+3)
		ics[1:(nt+1)]=lm$coefficients[1:(nt+1)]
		xhat=ics[1]+makebetatm(nt,ics,t)
		ics[2+nt]=sqrt((sum((x-xhat)^2))/nx)
		ics[3+nt]=0 #set xi to 0
#	}
	return(ics)
}
#'  observed log-likelihood function
#' @inherit manloglik return
#' @inheritParams manf
gev_p1n_loglik=function(vv,x,t){
	n=length(x)
	t=ifvectorthenmatrix(t)
	nt=findnt(t)
	mu=vv[1]+makebetatm(nt,vv,t)
	loglik=sum(extraDistr::dgev(x,mu=mu,sigma=max(vv[nt+2],.Machine$double.eps),xi=vv[nt+3],log=TRUE))
	return(loglik)
}
#' Check MLE
#' @inherit mancheckmle return
#' @inheritParams manf
gev_p1n_checkmle=function(ml_params,minxi=-1,maxxi=1){
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
qgev_p1n=function(p,t0,params){

	nt=findnt(t0)
	mu=params[1]+makebetat0(nt,params,t0)
	sigma=params[nt+2]
	xi=params[nt+3]
	return(extraDistr::qgev(p,mu=mu,sigma=sigma,xi=xi))

}
#' GEVD-with-p1: Density function
#' @inherit manvector return
#' @inheritParams manf
dgev_p1n=function(x,t0,params,log=FALSE){

	nt=findnt(t0)
	mu=params[1]+makebetat0(nt,params,t0)
	return(extraDistr::dgev(x,mu=mu,sigma=params[nt+2],xi=params[nt+3],log=log))

}
#' GEVD-with-p1: Distribution function
#' @inherit manvector return
#' @inheritParams manf
pgev_p1n=function(y,t0,params){

	nt=findnt(t0)
	mu=params[1]+makebetat0(nt,params,t0)
	return(extraDistr::pgev(y,mu=mu,sigma=params[nt+2],xi=params[nt+3]))

}
#' Analytical expressions for Predictive Means
#' RHP mean based on the expectation of DMGS equation 2.1
#' @inherit manmeans return
#' @inheritParams manf
gev_p1n_means=function(means,t0,ml_params,lddi,lddd,
									lambdad_rh_flat,nx,dim=(nt+3)){

	nt=findnt(t0)
	if(means){
# intro
		eulerconstant=0.57721566490153286060651209008240243104215933593992
#		ymn=ml_params[1]
#		slope1=ml_params[2]
#		slope2=ml_params[3]
		sigma=ml_params[(nt+2)]
		xi=ml_params[(nt+3)]

# set up derivative arrays
		meand1=array(0,c((nt+3),1))
		meand2=array(0,c((nt+3),(nt+3),1)) #but all zero for gumbel

		if(ml_params[nt+3]==0){
# xi=0 case

# mle
			ml_mean=ml_params[1]+makebetat0(nt,ml_params,t0)+sigma*eulerconstant
# calculate first derivative array for bayesian xi=0 cases
			meand1[1,1]=1
			for (i in 1:nt){
				meand1[(i+1),1]=t0[i]
			}
			meand1[(nt+2),1]=eulerconstant
			meand1[(nt+3),1]=0
# meand2 is all zero as initialized

		} else{
# non-gumbel case

			g0=gamma(1-xi)
			g1=g0*digamma(1-xi)
			g2=(trigamma(1-xi)*g0*g0+g1*g1)/g0
# mle
			ml_mean=ml_params[1]+makebetat0(nt,ml_params,t0)+sigma*(g0-1)/xi
# calculate first derivative array for bayesian xi!=0 cases
			meand1[1,1]=1
			for (i in 1:nt){
				meand1[(i+1),1]=t0[i]
			}
			meand1[(nt+2),1]=(g0-1)/xi
			meand1[(nt+3),1]=(1-g0-xi*g1)*sigma/(xi*xi)
# calculate second derivative array (only has 1 non-zero term!)
			meand2[(nt+2),(nt+3),1]=(1-g0-xi*g1)/(xi*xi)
			meand2[(nt+3),(nt+2),1]=meand2[(nt+2),(nt+3),1]
			meand2[(nt+3),(nt+3),1]=sigma*(-2+2*g0+2*xi*g1+xi*xi*g2)/(xi*xi*xi)
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
dgev_p1nsub=function(x,t,y,t0,ics,minxi,maxxi,extramodels=FALSE){

		nt=findnt(t)
		nx=length(x)

		ics=gev_p1n_setics(x,t)
		opt=optim(ics,gev_p1n_loglik,x=x,t=t,control=list(fnscale=-1))
		vhat=opt$par
		vhat[(nt+3)]=min(maxxi,max(minxi,opt$par[(nt+3)]))
		ml_params=vhat

# now that I've dropped dmgs, not sure I need this anymore
#		muhat0=v1hat+v2hat*t0
#		y=fixgevrange(y,muhat0,v3hat,v4hat,v5hat)

# mle
		ml_pdf=dgev_p1n(y,t0,params=vhat)
		ml_cdf=pgev_p1n(y,t0,params=vhat)

# return
		list(
					ml_params=ml_params,
					ml_pdf=ml_pdf,
					ml_cdf=ml_cdf)
}

