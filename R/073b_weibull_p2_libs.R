#' Waic
#' @inherit manwaic return
#' @inheritParams manf
weibull_p2_waic=function(waicscores,x,t,v1hat,v2hat,v3hat,
	lddi,lddd,lambdad){
		if(waicscores){
			f1f=weibull_p2_f1fw(x,t,v1hat,v2hat,v3hat)
			f2f=weibull_p2_f2fw(x,t,v1hat,v2hat,v3hat)
			fhatx=dweibull_p2(x,t,shape=v1hat,ymn=v2hat,slope=v3hat,log=FALSE)
			waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad,f2f,dim=3)
			waic1=waic$waic1
			waic2=waic$waic2
		}else{
			waic1="extras not selected"
			waic2="extras not selected"
		}
	list(waic1=waic1,waic2=waic2)
}
#' Predicted Parameter and Generalized Residuals
#' @inherit manpredictor return
#' @inheritParams manf
weibull_p2_predictordata=function(predictordata,x,t,t0,params){
	if(predictordata){
#
# calculate the probabilities of the data using the fited model
#
		sh=min(20,params[1])
		a=params[2]
		b=params[3]
		sc=exp(a+b*t)
		px=pweibull(x,shape=sh,scale=sc)
#
# calculate the quantiles for those probabilities at t0
#
		sc0=exp(a+b*t0)
		qx=qweibull(px,shape=sh,scale=sc0)
	} else{
		predictedparameter="predictordata not selected"
		adjustedx="predictordata not selected"
	}

	list(predictedparameter=sc,adjustedx=qx)
}
#' Logf for RUST
#' @inherit manlogf return
#' @inheritParams manf
weibull_p2_logf=function(params,x,t){
#	sh=min(50,params[1]) #high values of shape give NaNs in dweibull
#	a=params[2]
#	b=params[3]
#	sc=exp(a+b*t)
#	if((min(sh)>sqrt(.Machine$double.eps))&(min(sc)>sqrt(.Machine$double.eps))){
## if sc or s get too small, dweibull crashes, it seems
#	logf=sum(dweibull(x,shape=sh,scale=sc,log=TRUE))-log(sh)
#	}else{
#		logf=-Inf
#	}
	sh=pmax(min(20,params[1]),sqrt(.Machine$double.eps)) #high values of shape give NaNs in dweibull
	a=params[2]
	b=params[3]
	sc=pmax(exp(a+b*t),sqrt(.Machine$double.eps))
	logf=sum(dweibull(x,shape=sh,scale=sc,log=TRUE))-log(sh)
	if(is.na(logf)){
		message("dweibull is giving NaNs again...let's have a look why.")
		message("logf=",logf)
		message("a=",a)
		message("b=",b)
		message("sh=",sh)
		message("sc=",sc)
		for(i in 1:length(x)){
			message("x,d=",x[i],dweibull(x[i],shape=sh,scale=sc,log=TRUE))
		}
	}
	return(logf)
}
#'  observed log-likelihood function
#' @inherit manloglik return
#' @inheritParams	manf
weibull_p2_loglik=function(vv,x,t){
	sh=pmax(min(20,vv[1]),.Machine$double.eps)
	sc=pmax(exp(vv[2]+vv[3]*t),.Machine$double.eps)
	loglik=sum(dweibull(x,shape=sh,scale=sc,log=TRUE))
	if(is.na(loglik))message("\n B:",loglik,sh,sc)
	return(loglik)
}
#' Weibull-with-p1 quantile function
#' @inherit manvector return
#' @inheritParams	manf
qweibull_p2=function(p,t0,shape,ymn,slope){

	sc=exp(ymn+slope*t0)
	return(qweibull(p,shape=shape,scale=sc))

}
#' Weibull-with-p1 density function
#' @inherit manvector return
#' @inheritParams	manf
dweibull_p2=function(x,t0,shape,ymn,slope,log=FALSE){

	shape=pmax(min(20,shape),.Machine$double.eps)
	sc=pmax(exp(ymn+slope*t0),.Machine$double.eps)
	return(dweibull(x,shape=shape,scale=sc,log=log))

}
#' Weibull-with-p1 distribution function
#' @inherit manvector return
#' @inheritParams	manf
pweibull_p2=function(x,t0,shape,ymn,slope){

	sc=exp(ymn+slope*t0)
	return(pweibull(x,shape=shape,scale=sc))

}
#' weibull distribution: RHP mean
#' @inherit manmeans return
#' @inheritParams	manf
weibull_p2_means=function(means,t0,ml_params,lddi,lddd,lambdad_rhp,nx,dim){

	if(means){
# intro
		v1=ml_params[1]
		v2=ml_params[2]
		v3=ml_params[3]

# ml mean
		mu_hat=v2+v3*t0
		arg=1+(1/v1)
		ml_mean=exp(mu_hat)*gamma(arg)

# rhp mean
# starts with derivatives of the log of the mean, to make the algebra easier
#
		lmeand2=array(0,c(3,1))
		lmeand2[1,1]=1
		lmeand2[2,1]=t0
		lmeand2[3,1]=(1/(v1*v1))*digamma(arg)
		meand2=array(0,c(3,1))
		meand2[1,1]=ml_mean*lmeand2[1,1]
		meand2[2,1]=ml_mean*lmeand2[2,1]
		meand2[3,1]=ml_mean*lmeand2[3,1]
#
		lmeand3=array(0,c(3,3,1))
		lmeand3[1,1,1]=0 #alpha-alpha
		lmeand3[1,2,1]=0 #alpha-beta
		lmeand3[1,3,1]=0 #alpha-k
		lmeand3[2,1,1]=0 #beta-alpha
		lmeand3[2,2,1]=0 #beta-beta
		lmeand3[2,2,1]=0 #beta-k
		lmeand3[3,1,1]=0 #alpha-k
		lmeand3[3,2,1]=0 #beta-k
		lmeand3[3,3,1]=2*digamma(arg)/(v1*v1*v1)+trigamma(arg)/(v1*v1*v1*v1)
#
		meand3=array(0,c(3,3,1))
		meand3[1,1,1]=lmeand2[1,1]/ml_mean								#alpha-alpha
		meand3[1,2,1]=lmeand2[1,1]*lmeand2[2,1]/ml_mean	#alpha-beta
		meand3[1,3,1]=lmeand2[1,1]*lmeand2[3,1]/ml_mean	#alpha-k
		meand3[2,1,1]=meand3[1,2,1] 										#beta-alpha
		meand3[2,2,1]=lmeand2[2,1]*lmeand2[2,1]/ml_mean	#beta-beta
		meand3[2,3,1]=lmeand2[2,1]*lmeand2[3,1]/ml_mean	#beta-k
		meand3[3,1,1]=meand3[1,3,1] 										#alpha-k
		meand3[3,2,1]=ml_mean*lmeand3[3,3,1]+lmeand2[3,1]*lmeand2[3,1]/ml_mean #beta-k
#
		dmean=dmgs(lddi,lddd,meand2,lambdad_rhp,meand3,dim=3)
		rh_mean=ml_mean+dmean/nx
	}else{
		ml_mean="means not selected"
		rh_mean="means not selected"
	}


	list(ml_mean=ml_mean,rh_mean=rh_mean)

}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inherit manlogscores return
#' @inheritParams	manf
weibull_p2_logscores=function(logscores,x,t){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){

			x1=x[-i]
			t1=t[-i]

			dd=dweibull_p2sub(x1,t1,x[i],t[i])

			ml_pdf=dd$ml_pdf
			ml_oos_logscore=ml_oos_logscore+log(ml_pdf)

			rh_pdf=dd$rh_pdf
			rh_oos_logscore=rh_oos_logscore+log(rh_pdf)
		}
	}else{
		ml_oos_logscore="extras not selected"
		rh_oos_logscore="extras not selected"
	}
	list(ml_oos_logscore=ml_oos_logscore,rh_oos_logscore=rh_oos_logscore)
}
#' Densities from MLE and RHP
#' @inherit mandsub return
#' @inheritParams	manf
dweibull_p2sub=function(x,t,y,t0){

		nx=length(x)

		opt1=optim(c(1,0,0),weibull_p2_loglik,x=x,t=t,
			control=list(fnscale=-1))
		v1hat=opt1$par[1]
		v2hat=opt1$par[2]
		v3hat=opt1$par[3]
		ml_params=c(v1hat,v2hat,v3hat)

# ml
		sc=exp(v2hat+v3hat*t0)
		ml_pdf=dweibull(y,shape=v1hat,scale=sc)
		ml_cdf=pweibull(y,shape=v1hat,scale=sc)

# rhp
		ldd=weibull_p2_ldda(x,t,v1hat,v2hat,v3hat)
		lddi=solve(ldd)
		lddd=weibull_p2_lddda(x,t,v1hat,v2hat,v3hat)

		f1=weibull_p2_f1fa(y,t0,v1hat,v2hat,v3hat)
		f2=weibull_p2_f2fa(y,t0,v1hat,v2hat,v3hat)

		p1=weibull_p2_p1fa(y,t0,v1hat,v2hat,v3hat)
		p2=weibull_p2_p2fa(y,t0,v1hat,v2hat,v3hat)

		lambdad_rhp=c(-1/v1hat,0,0)
		df=dmgs(lddi,lddd,f1,lambdad_rhp,f2,dim=3)
		dp=dmgs(lddi,lddd,p1,lambdad_rhp,p2,dim=3)
		rh_pdf=pmax(ml_pdf+df/nx,0)
		rh_cdf=pmin(pmax(ml_cdf+dp/nx,0),1)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					rh_pdf=rh_pdf,
					ml_cdf=ml_cdf,
					rh_cdf=rh_cdf)
}
