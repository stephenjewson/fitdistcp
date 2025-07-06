#' Waic
#' @inherit manwaic return
#' @inheritParams manf
exp_p1_waic=function(waicscores,x,t,v1hat,v2hat,lddi,lddd,lambdad){
		if(waicscores){

			f1f=exp_p1_f1fw(x,t,v1hat,v2hat)
			f2f=exp_p1_f2fw(x,t,v1hat,v2hat)

			fhatx=dexp_p1(x,t,ymn=v1hat,slope=v2hat,log=FALSE)
			waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad,f2f,dim=2)
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
exp_p1_predictordata=function(predictordata,x,t,t0,params){
	if(predictordata){
#
# calculate the probabilities of the data using the fited model
#
		a=params[1]
		b=params[2]
		rr=1/exp(a+b*t)
		px=pexp(x,rate=rr)
#
# calculate the quantiles for those probabilities at t0
#
		rr0=1/exp(a+b*t0)
		qx=qexp(px,rate=rr0)
	} else{
		predictedparameter="predictordata not selected"
		adjustedx="predictordata not selected"
	}

	list(predictedparameter=rr,adjustedx=qx)
}
#' Logf for RUST
#' @inherit manlogf return
#' @inheritParams manf
exp_p1_logf=function(params,x,t){
#	a=params[1]
#	b=params[2]
#	r=1/exp(a+b*t)
#	logf=sum(dexp(x,rate=r,log=TRUE))
	a=params[1]
	b=params[2]
	r=1/exp(a+b*t) #rate can be zero, that's ok...it's 1/sigma
	logf=sum(dexp(x,rate=r,log=TRUE))
	return(logf)
}
#'  observed log-likelihood function
#' @inherit manloglik return
#' @inheritParams	manf
exp_p1_loglik=function(vv,x,t){
	mu=vv[1]+vv[2]*t
# the 1/ parametrisation is beter because there is an expression for the mean
# -really? but the only difference is  signs
# in the ML model
# so this is the version in which sigma is loglinear, not rate
# for standardisation, beter to have sigma as loglinear, not rate. That's the best argument.
	loglik=sum(dexp(x,rate=1/exp(mu),log=TRUE))
#	loglik=-sum(dexp(x,exp(mu),log=TRUE))
	return(loglik)
}
#' -with-p1 quantile function
#' @inherit manvector return
#' @inheritParams	manf
qexp_p1=function(p,t0,ymn,slope){

	mu=(ymn+slope*t0)
	return(qexp(p,rate=(1/exp(mu))))
#	return(qexp(p,rate=exp(mu)))

}
#' Exponential-with-p1 density function
#' @inherit manvector return
#' @inheritParams	manf
dexp_p1=function(x,t0,ymn,slope,log=FALSE){

	mu=(ymn+slope*t0)
	return(dexp(x,rate=(1/exp(mu)),log=log))
#	return(dexp(x,rate=exp(mu),log=log))

}
#' Exponential-with-p1 distribution function
#' @inherit manvector return
#' @inheritParams	manf
pexp_p1=function(x,t0,ymn,slope){

	mu=(ymn+slope*t0)
	return(pexp(x,rate=(1/exp(mu))))
#	return(pexp(x,rate=exp(mu)))

}
#' exp distribution: RHP means
#' @inherit manmeans return
#' @inheritParams	manf
exp_p1_means=function(means,t0,ml_params,lddi,lddd,lambdad_rhp,nx,dim=2){

	if(means){
# intro
		v1=ml_params[1]
		v2=ml_params[2]

# ml mean
		mu_hat=v1+v2*t0
		ml_mean=exp(mu_hat)

# rhp mean
		meand1=array(0,c(2,1))
		meand1[1,1]=ml_mean
		meand1[2,1]=t0*ml_mean
		meand2=array(0,c(2,2,1))
		meand2[1,1,1]=ml_mean
		meand2[1,2,1]=t0*ml_mean
		meand2[2,1,1]=t0*ml_mean
		meand2[2,2,1]=t0*t0*ml_mean
		dmean=dmgs(lddi,lddd,meand1,lambdad_rhp,meand2,dim=2)
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
exp_p1_logscores=function(logscores,x,t){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){

			x1=x[-i]
			t1=t[-i]

			dd=dexp_p1sub(x1,t1,x[i],t[i])

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
dexp_p1sub=function(x,t,y,t0){

		nx=length(x)

		lm=lm(x~t)
		v1start=lm$coefficients[1]
		v2start=lm$coefficients[2]
		v1start=0
		v2start=0
#		xhat=v1start+v2start*t
		opt1=optim(c(v1start,v2start),exp_p1_loglik,x=x,t=t,control=list(fnscale=-1))
		v1hat=opt1$par[1]
		v2hat=opt1$par[2]
		ml_params=c(v1hat,v2hat)

# ml
		muhat=v1hat+v2hat*t0
		ml_pdf=dexp(y,rate=1/exp(muhat))
		ml_cdf=pexp(y,rate=1/exp(muhat))
#		ml_pdf=dexp(y,rate=exp(muhat))
#		ml_cdf=pexp(y,rate=exp(muhat))

# rhp
		ldd=exp_p1_ldda(x,t,v1hat,v2hat)
		lddi=solve(ldd)
		lddd=exp_p1_lddda(x,t,v1hat,v2hat)

		f1=exp_p1_f1fa(y,t0,v1hat,v2hat)
		f2=exp_p1_f2fa(y,t0,v1hat,v2hat)

		p1=exp_p1_p1fa(y,t0,v1hat,v2hat)
		p2=exp_p1_p2fa(y,t0,v1hat,v2hat)

		lambdad_rhp=c(0,0)
		df=dmgs(lddi,lddd,f1,lambdad_rhp,f2,dim=2)
		dp=dmgs(lddi,lddd,p1,lambdad_rhp,p2,dim=2)
		rh_pdf=pmax(ml_pdf+df/nx,0)
		rh_cdf=pmin(pmax(ml_cdf+dp/nx,0),1)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					rh_pdf=rh_pdf,
					ml_cdf=ml_cdf,
					rh_cdf=rh_cdf)
}
