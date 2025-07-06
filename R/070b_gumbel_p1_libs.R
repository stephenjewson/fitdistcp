#' Waic
#' @inherit manwaic return
#' @inheritParams manf
gumbel_p1_waic=function(waicscores,x,t,v1hat,v2hat,v3hat,
	lddi,lddd,lambdad){
		if(waicscores){
			f1f=gumbel_p1_f1fw(x,t,v1hat,v2hat,v3hat)
			f2f=gumbel_p1_f2fw(x,t,v1hat,v2hat,v3hat)
			fhatx=dgumbel_p1(x,t,ymn=v1hat,slope=v2hat,sigma=v3hat,log=FALSE)
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
gumbel_p1_predictordata=function(predictordata,x,t,t0,params){
	if(predictordata){
#
# calculate the probabilities of the data using the fited model
#
		a=params[1]
		b=params[2]
		s=params[3]
		mu=a+b*t
		px=pgumbel(x,mu=mu,sigma=s)
#
# calculate the quantiles for those probabilities at t0
#
		mu0=a+b*t0
		qx=qgumbel(px,mu=mu0,sigma=s)
	} else{
		predictedparameter="predictordata not selected"
		adjustedx="predictordata not selected"
	}

	list(predictedparameter=mu,adjustedx=qx)
}
#' Logf for RUST
#' @inherit manlogf return
#' @inheritParams manf
gumbel_p1_logf=function(params,x,t){
#	a=params[1]
#	b=params[2]
#	s=params[3]
#	mu=a+b*t
#	if(s>0){
#		logf=sum(dgumbel(x,mu=mu,sigma=s,log=TRUE))-log(s)
#	}else{
#		logf=-Inf
#	}
	a=params[1]
	b=params[2]
	s=pmax(params[3],.Machine$double.eps)
	mu=a+b*t
	logf=sum(dgumbel(x,mu=mu,sigma=s,log=TRUE))-log(s)
	return(logf)
}
#'  observed log-likelihood function
#' @inherit manloglik return
#' @inheritParams	manf
gumbel_p1_loglik=function(vv,x,t){
	mu=vv[1]+vv[2]*t
	loglik=sum(dgumbel(x,mu=mu,sigma=max(vv[3],.Machine$double.eps),log=TRUE))
	return(loglik)
}
#' Gumbel-with-p1 quantile function
#' @inherit manvector return
#' @inheritParams	manf
qgumbel_p1=function(p,t0,ymn,slope,sigma){

	return(qgumbel(p,mu=(ymn+slope*t0),sigma=sigma))

}
#' Gumbel-with-p1 density function
#' @inherit manvector return
#' @inheritParams	manf
dgumbel_p1=function(x,t0,ymn,slope,sigma,log=FALSE){

	return(dgumbel(x,mu=(ymn+slope*t0),sigma=sigma,log=log))

}
#' Gumbel-with-p1 distribution function
#' @inherit manvector return
#' @inheritParams	manf
pgumbel_p1=function(x,t0,ymn,slope,sigma){

	return(pgumbel(x,mu=(ymn+slope*t0),sigma=sigma))

}
#' Gumbel distribution: RHP mean
#' @inherit manmeans return
#' @inheritParams	manf
gumbel_p1_means=function(means,t0,ml_params,lddi,lddd,lambdad_rhp,nx,dim){

	if(means){
# intro
		eulerconstant=0.57721566490153286060651209008240243104215933593992
		v1=ml_params[1]
		v2=ml_params[2]
		v3=ml_params[3]

# ml mean
		eulerconstant=0.57721566490153286060651209008240243104215933593992
		mu_hat=v1+v2*t0
		ml_mean=mu_hat+v3*eulerconstant

# rhp mean
		meand1=array(0,c(3,1))
		meand1[1,1]=1
		meand1[2,1]=t0
		meand1[3,1]=eulerconstant
		meand2=array(0,c(3,3,1)) #but all zero for gumbel
		dmean=dmgs(lddi,lddd,meand1,lambdad_rhp,meand2,dim=3)
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
gumbel_p1_logscores=function(logscores,x,t){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){

			x1=x[-i]
			t1=t[-i]

			dd=dgumbel_p1sub(x1,t1,x[i],t[i])

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
dgumbel_p1sub=function(x,t,y,t0){

		nx=length(x)

		lm=lm(x~t)
		v1start=lm$coefficients[1]
		v2start=lm$coefficients[2]
		xhat=v1start+v2start*t
		v3start=sqrt((sum((x-xhat)^2))/(nx-1))
		opt1=optim(c(v1start,v2start,v3start),gumbel_p1_loglik,x=x,t=t,
			control=list(fnscale=-1))
		v1hat=opt1$par[1]
		v2hat=opt1$par[2]
		v3hat=opt1$par[3]
		ml_params=c(v1hat,v2hat,v3hat)

# ml
		muhat=v1hat+v2hat*t0
		ml_pdf=dgumbel(y,mu=muhat,sigma=v3hat)
		ml_cdf=pgumbel(y,mu=muhat,sigma=v3hat)

# rhp
		ldd=gumbel_p1_ldda(x,t,v1hat,v2hat,v3hat)
		lddi=solve(ldd)
		lddd=gumbel_p1_lddda(x,t,v1hat,v2hat,v3hat)

		f1=gumbel_p1_f1fa(y,t0,v1hat,v2hat,v3hat)
		f2=gumbel_p1_f2fa(y,t0,v1hat,v2hat,v3hat)

		p1=gumbel_p1_p1fa(y,t0,v1hat,v2hat,v3hat)
		p2=gumbel_p1_p2fa(y,t0,v1hat,v2hat,v3hat)

		lambdad_rhp=c(0,0,-1/v3hat)
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
