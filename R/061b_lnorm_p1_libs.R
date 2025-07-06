#' Waic
#' @inherit manwaic return
#' @inheritParams manf
lnorm_p1_waic=function(waicscores,x,t,v1hat,v2hat,v3hat){
# norm and lnorm waic code works differently from other cases
# normally lddi and lddd have already been calculated, and I just pass them in
# but in this case they haven't been calculated
# so I need to calculate them here
# for which I need t, in addition to t0
# confused about whether I should pass in t,t0 or ta, ta0 though
	if(waicscores){
		f1f=lnorm_p1_f1fw(x,t,v1hat,v2hat,v3hat)
		f2f=lnorm_p1_f2fw(x,t,v1hat,v2hat,v3hat)

		ldd=lnorm_p1_ldda(x,t,v1hat,v2hat,v3hat)
		lddi=solve(ldd)
		lddd=lnorm_p1_lddda(x,t,v1hat,v2hat,v3hat)

		fhatx=dlnorm_p1(x,t,ymn=v1hat,slope=v2hat,sigma=v3hat,log=FALSE)
		lambdad_rhp=c(0,0,-1/v3hat)
		waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad_rhp,f2f,dim=3)
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
lnorm_p1_predictordata=function(x,t,t0,params){
#
# calculate the probabilities of the data using the fited model
#
		a=params[1]
		b=params[2]
		s=params[3]
		mu=a+b*t
		px=plnorm(x,meanlog=mu,sdlog=s)
#
# calculate the quantiles for those probabilities at t0
#
		mu0=a+b*t0
		qx=qlnorm(px,meanlog=mu0,sdlog=s)

	list(predictedparameter=mu,adjustedx=qx)
}
#' Logf for RUST
#' @inherit manlogf return
#' @inheritParams manf
lnorm_p1_logf=function(params,x,t){
#	a=params[1]
#	b=params[2]
#	s=params[3]
#	mu=a+b*t
#	if(s>0){
#		logf=sum(dlnorm(x,meanlog=mu,sdlog=s,log=TRUE))-log(s)
#	}else{
#		logf=-Inf
#	}
	a=params[1]
	b=params[2]
	s=pmax(params[3],.Machine$double.eps)
	mu=a+b*t
	logf=sum(dlnorm(x,meanlog=mu,sdlog=s,log=TRUE))-log(s)
	return(logf)
}
#' Log-normal-with-p1  observed log-likelihood function
#' @inherit manloglik return
#' @inheritParams manf
lnorm_p1_loglik=function(vv,x,t){
	mu=vv[1]+vv[2]*t #so mean is a vector, just like x
	loglik=sum(dlnorm(x,meanlog=mu,sdlog=max(vv[3],0),log=TRUE))
	return(loglik)
}
#' Normal-with-p1 quantile function
#' @inherit manvector return
#' @inheritParams manf
qlnorm_p1=function(p,t0,ymn,slope,sigma){

	return(qlnorm(p,meanlog=(ymn+slope*t0),sdlog=sigma))

}
#' Normal-with-p1 density function
#' @inherit manvector return
#' @inheritParams manf
dlnorm_p1=function(x,t0,ymn,slope,sigma,log=FALSE){

	return(dlnorm(x,meanlog=(ymn+slope*t0),sdlog=sigma,log=log))

}
#' Normal-with-p1 distribution function
#' @inherit manvector return
#' @inheritParams manf
plnorm_p1=function(x,t0,ymn,slope,sigma){

	return(plnorm(x,meanlog=(ymn+slope*t0),sdlog=sigma))

}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inherit manlogscores return
#' @inheritParams manf
lnorm_p1_logscores=function(logscores,x,t){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){

			x1=x[-i]
			t1=t[-i]

			dd=dlnorm_p1sub(x1,t1,x[i],t[i])

			ml_params=dd$ml_params

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
#' @inheritParams manf
dlnorm_p1sub=function(x,t,y,t0,debug=FALSE){
# y are the future values of x
# z is log x

		nx=length(x)

		meant=mean(t)
		ta=t-mean(t)
		ta0=t0-meant

		z=log(x)
		logy=log(y)

		ml_params=norm_p1_mlparams(z,t) #this one really is norm not lnorm
		if(debug)message(" inside dlnorm_p1, ml_params=",ml_params)
		v1hat=ml_params[1]
		v2hat=ml_params[2]
		v3hat=ml_params[3]
		muhat0=v1hat+v2hat*t0
# ml
		ml_pdf=dlnorm(y,meanlog=muhat0,sdlog=v3hat)
		ml_cdf=plnorm(y,meanlog=muhat0,sdlog=v3hat)

# rhp
		rh_pdf=dnorm_p1_formula(logy,ta,ta0,nx,muhat0,v3hat)/y #using the formula from norm_p1
		rh_pdf=pmax(rh_pdf,0)

		rh_cdf=pnorm_p1_formula(logy,ta,ta0,nx,muhat0,v3hat) #using the formula from norm_p1
		rh_cdf=pmin(pmax(rh_cdf,0),1)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					rh_pdf=rh_pdf,
					ml_cdf=ml_cdf,
					rh_cdf=rh_cdf)
}

