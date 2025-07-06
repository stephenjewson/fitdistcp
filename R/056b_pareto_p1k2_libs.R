#' Waic
#' @inherit manwaic return
#' @inheritParams manf
pareto_p1k2_waic=function(waicscores,x,t,v1hat,v2hat,kscale,lddi,lddd,
	lambdad){
		if(waicscores){

			f1f=pareto_p1k2_f1fw(x,t,v1hat,v2hat,kscale)
			f2f=pareto_p1k2_f2fw(x,t,v1hat,v2hat,kscale)

			fhatx=dpareto_p1k2(x,t,ymn=v1hat,slope=v2hat,kscale=kscale,log=FALSE)
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
pareto_p1k2_predictordata=function(predictordata,x,t,t0,params,kscale){
	if(predictordata){
#
# calculate the probabilities of the data using the fited model
#
		a=params[1]
		b=params[2]
		sh=1/exp(a+b*t)
		px=ppareto(x,a=sh,b=kscale)
#
# calculate the quantiles for those probabilities at t0
#
		sh0=1/exp(a+b*t0)
		qx=qpareto(px,a=sh0,b=kscale)
	} else{
		predictedparameter="predictordata not selected"
		adjustedx="predictordata not selected"
	}

	list(predictedparameter=sh,adjustedx=qx)
}
#' Logf for RUST
#' @inherit manlogf return
#' @inheritParams manf
pareto_p1k2_logf=function(params,x,t,kscale){
	a=params[1]
	b=params[2]
	sh=pmax(1/exp(a+b*t),.Machine$double.eps)
	logf=sum(dpareto(x,a=sh,b=kscale,log=TRUE))
	return(logf)
}
#'  observed log-likelihood function
#' @inherit manloglik return
#' @inheritParams	manf
pareto_p1k2_loglik=function(vv,x,t,kscale){
	mu=vv[1]+vv[2]*t
	loglik=sum(dpareto(x,a=1/exp(mu),b=kscale,log=TRUE))
	return(loglik)
}
#' pareto_k1-with-p2 quantile function
#' @inherit manvector return
#' @inheritParams	manf
qpareto_p1k2=function(p,t0,ymn,slope,kscale){

	mu=(ymn+slope*t0)
	return(qpareto(p,a=1/exp(mu),b=kscale))

}
#' pareto_k1-with-p2 density function
#' @inherit manvector return
#' @inheritParams	manf
dpareto_p1k2=function(x,t0,ymn,slope,kscale,log=FALSE){

	mu=(ymn+slope*t0)
	return(dpareto(x,a=1/exp(mu),b=kscale,log=log))

}
#' pareto_k1-with-p2 distribution function
#' @inherit manvector return
#' @inheritParams	manf
ppareto_p1k2=function(x,t0,ymn,slope,kscale){

	mu=(ymn+slope*t0)
	return(ppareto(x,a=1/exp(mu),b=kscale))

}
#' pareto_k1 distribution: RHP mean
#' @inherit manmeans return
#' @inheritParams	manf
pareto_p1k2_means=function(means,t0,ml_params,lddi,lddd,lambdad_rhp,nx,dim=2,kscale){

	if(means){
# intro
		v1=ml_params[1]
		v2=ml_params[2]

# ml mean
		mu_hat=v1+v2*t0
		if(mu_hat>1){
			ml_mean=mu_hat*kscale/(mu_hat-1)
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
pareto_p1k2_logscores=function(logscores,x,t,kscale,debug){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){

			x1=x[-i]
			t1=t[-i]

			dd=dpareto_p1k2sub(x1,t1,x[i],t[i],kscale=kscale,debug=debug)

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
dpareto_p1k2sub=function(x,t,y,t0,kscale,debug=FALSE){


		if(debug)message("inside pareto_p1k2sub")
		nx=length(x)

		lm=lm(x~t)
		v1start=lm$coefficients[1]
		v2start=lm$coefficients[2]
		v1start=0
		v2start=1
		xhat=v1start+v2start*t
		opt1=optim(c(v1start,v2start),pareto_p1k2_loglik,x=x,t=t,
			kscale=kscale,control=list(fnscale=-1))
		v1hat=opt1$par[1]
		v2hat=opt1$par[2]
		ml_params=c(v1hat,v2hat)

# ml
		muhat=v1hat+v2hat*t0
		ml_pdf=dpareto(y,a=1/exp(muhat),b=kscale)
		ml_cdf=ppareto(y,a=1/exp(muhat),b=kscale)

# rhp
		if(debug)message("calc ldd")
		ldd=pareto_p1k2_ldda(x,t,v1hat,v2hat,kscale)
		lddi=solve(ldd)

		if(debug)message("calc lddd")
		lddd=pareto_p1k2_lddda(x,t,v1hat,v2hat,kscale)

		if(debug)message("calc f1f")
		f1=pareto_p1k2_f1fa(y,t0,v1hat,v2hat,kscale)
		if(debug)message("calc f2f")
		f2=pareto_p1k2_f2fa(y,t0,v1hat,v2hat,kscale)

		if(debug)message("calc p1f")
		p1=pareto_p1k2_p1fa(y,t0,v1hat,v2hat,kscale)
		if(debug)message("calc p2f")
		p2=pareto_p1k2_p2fa(y,t0,v1hat,v2hat,kscale)

		if(debug)message("call dmgs")
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
