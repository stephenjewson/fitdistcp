#' Waic
#' @inherit manwaic return
#' @inheritParams manf
frechet_p2k1_waic=function(waicscores,x,t,v1hat,v2hat,v3hat,kloc,
	lddi,lddd,lambdad){
		if(waicscores){
			f1f=frechet_p2k1_f1fw(x,t,v1hat,v2hat,v3hat,kloc=kloc)
			f2f=frechet_p2k1_f2fw(x,t,v1hat,v2hat,v3hat,kloc=kloc)
			fhatx=dfrechet_p2k1(x,t,ymn=v1hat,slope=v2hat,lambda=v3hat,log=FALSE,kloc=kloc)
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
frechet_p2k1_predictordata=function(predictordata,x,t,t0,params,kloc){
	if(predictordata){
#
# calculate the probabilities of the data using the fited model
#
		a=params[1]
		b=params[2]
		s=params[3]
		sg=exp(a+b*t)
		px=pfrechet(x,mu=kloc,sigma=sg,lambda=s)
#
# calculate the quantiles for those probabilities at t0
#
		sg0=exp(a+b*t0)
		qx=qfrechet(px,mu=kloc,sigma=sg0,lambda=s)
	} else{
		predictedparameter="predictordata not selected"
		adjustedx="predictordata not selected"
	}

	list(predictedparameter=sg,adjustedx=qx)
}
#' Logf for RUST
#' @inherit manlogf return
#' @inheritParams manf
frechet_p2k1_logf=function(params,x,t,kloc){
#	a=params[1]
#	b=params[2]
#	s=params[3]
#	sg=exp(a+b*t)
#	if(s>0){
#		logf=sum(dfrechet(x,mu=kloc,sigma=sg,lambda=s,log=TRUE))-log(s)
#	}else{
#		logf=-Inf
#	}
	a=params[1]
	b=params[2]
	s=pmax(params[3],.Machine$double.eps)
	sg=exp(a+b*t)
	logf=sum(dfrechet(x,mu=kloc,sigma=sg,lambda=s,log=TRUE))-log(s)
	return(logf)
}
#'  observed log-likelihood function
#' @inherit manloglik return
#' @inheritParams	manf
frechet_p2k1_loglik=function(vv,x,t,kloc){
	mu=vv[1]+vv[2]*t
	loglik=sum(dfrechet(x,mu=kloc,sigma=exp(mu),lambda=max(vv[3],.Machine$double.eps),log=TRUE))
	return(loglik)
}
#' Frechet_k1-with-p2 quantile function
#' @inherit manvector return
#' @inheritParams	manf
qfrechet_p2k1=function(p,t0,ymn,slope,lambda,kloc){

	mu=(ymn+slope*t0)
	return(qfrechet(p,mu=kloc,sigma=exp(mu),lambda=lambda))

}
#' Frechet_k1-with-p2 density function
#' @inherit manvector return
#' @inheritParams	manf
dfrechet_p2k1=function(x,t0,ymn,slope,lambda,log=FALSE,kloc){

	mu=(ymn+slope*t0)
	return(dfrechet(x,mu=kloc,sigma=exp(mu),lambda=lambda,log=log))

}
#' Frechet_k1-with-p2 distribution function
#' @inherit manvector return
#' @inheritParams	manf
pfrechet_p2k1=function(x,t0,ymn,slope,lambda,kloc){

	mu=(ymn+slope*t0)
	return(pfrechet(x,mu=kloc,sigma=exp(mu),lambda=lambda))

}
#' frechet_k1 distribution: RHP mean
#' @inherit manmeans return
#' @inheritParams	manf
frechet_p2k1_means=function(means,t0,ml_params,lddi,lddd,lambdad_rhp,nx,dim,kloc){

	if(means){
# intro
		v1=ml_params[1]
		v2=ml_params[2]
		v3=ml_params[3]

# ml mean
		mu_hat=v1+v2*t0
		if(mu_hat>1){
			ml_mean=mu_hat*gamma(1-(1/v3))
		} else {
			ml_mean=Inf
		}

# rhp mean
		rh_mean=Inf

	}else{
		ml_mean="means not selected"
		rh_mean="means not selected"
	}

	list(ml_mean=ml_mean,rh_mean=rh_mean)

}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inherit manlogscores return
#' @inheritParams	manf
frechet_p2k1_logscores=function(logscores,x,t,kloc){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){

			x1=x[-i]
			t1=t[-i]

			dd=dfrechet_p2k1sub(x1,t1,x[i],t[i],kloc)

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
dfrechet_p2k1sub=function(x,t,y,t0,kloc){

		nx=length(x)

		lm=lm(x~t)
		v1start=lm$coefficients[1]
		v2start=lm$coefficients[2]
		xhat=v1start+v2start*t
		v3start=sqrt((sum((x-xhat)^2))/(nx-1))
		v1start=kloc
		v2start=0
		v3start=2
		opt1=optim(c(v1start,v2start,v3start),frechet_p2k1_loglik,x=x,t=t,
			kloc=kloc,control=list(fnscale=-1))
		v1hat=opt1$par[1]
		v2hat=opt1$par[2]
		v3hat=opt1$par[3]
		ml_params=c(v1hat,v2hat,v3hat)

# ml
		muhat=v1hat+v2hat*t0
		ml_pdf=dfrechet(y,mu=kloc,sigma=exp(muhat),lambda=v3hat)
		ml_cdf=pfrechet(y,mu=kloc,sigma=exp(muhat),lambda=v3hat)

# rhp
		ldd=frechet_p2k1_ldda(x,t,v1hat,v2hat,v3hat,kloc=kloc)
		lddi=solve(ldd)
		lddd=frechet_p2k1_lddda(x,t,v1hat,v2hat,v3hat,kloc=kloc)

		f1=frechet_p2k1_f1fa(y,t0,v1hat,v2hat,v3hat,kloc=kloc)
		f2=frechet_p2k1_f2fa(y,t0,v1hat,v2hat,v3hat,kloc=kloc)

		p1=frechet_p2k1_p1fa(y,t0,v1hat,v2hat,v3hat,kloc=kloc)
		p2=frechet_p2k1_p2fa(y,t0,v1hat,v2hat,v3hat,kloc=kloc)

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
