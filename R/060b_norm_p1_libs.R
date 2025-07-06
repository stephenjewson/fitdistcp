#' Linear regression formula, densities
#' @inherit manvector return
#' @inheritParams manf
pnorm_p1_formula=function(y,ta,ta0,nx,muhat0,v3hat){
	top=ta0*ta0
	sx=sum(ta*ta)

# convert maxlik sg to sum^2/(n-2), which is used in the standard formula
# noting that maxlik is equivalent to using 1/(n-1)
	sg1=v3hat*sqrt((nx-1)/(nx-2))
	sg2=sg1*sqrt((1+1/nx+top/sx))

	yd=(y-muhat0)/sg2
	rh_cdf=(pt(yd,df=nx-2))
	return(rh_cdf)
}
#' Linear regression formula, densities
#' @inherit manvector return
#' @inheritParams manf
dnorm_p1_formula=function(y,ta,ta0,nx,muhat0,v3hat){
	top=ta0*ta0
	sx=sum(ta*ta)

# convert maxlik sg to sum^2/(n-2), which is used in the standard formula
	sg1=v3hat*sqrt((nx-1)/(nx-2))
	sg2=sg1*sqrt((1+1/nx+top/sx))

	yd=(y-muhat0)/sg2
	rh_pdf=(dt(yd,df=nx-2))/sg2
	return(rh_pdf)
}
#' Linear regression formula, quantiles
#' @inherit manvector return
#' @inheritParams manf
qnorm_p1_formula=function(alpha,ta,ta0,nx,muhat0,v3hat){
	top=ta0*ta0
	sx=sum(ta*ta)

# convert maxlik sg to sum^2/(n-2), which is used in the standard formula
	sg1=v3hat*sqrt((nx-1)/(nx-2))
	sg2=sg1*sqrt((1+1/nx+top/sx))

	temp=qt((1-alpha),df=nx-2)
	rh_quantiles=muhat0+temp*sg2
	return(rh_quantiles)
}
#' Waic
#' @inherit manwaic return
#' @inheritParams manf
norm_p1_waic=function(waicscores,x,t,v1hat,v2hat,v3hat){
	if(waicscores){
		f1f=norm_p1_f1fw(x,t,v1hat,v2hat,v3hat)
		f2f=norm_p1_f2fw(x,t,v1hat,v2hat,v3hat)

		ldd=norm_p1_ldda(x,t,v1hat,v2hat,v3hat)
		lddi=solve(ldd)
		lddd=norm_p1_lddda(x,t,v1hat,v2hat,v3hat)

		fhatx=dnorm_p1(x,t,ymn=v1hat,slope=v2hat,sigma=v3hat,log=FALSE)
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
norm_p1_predictordata=function(x,t,t0,params){
#
# calculate the probabilities of the data using the fited model
#
# note that t may be centred
# -but the params take that into account
# -so the mu is not affected by centering
		a=params[1]
		b=params[2]
		s=params[3]
		mu=a+b*t
		px=pnorm(x,mean=mu,sd=s)
#
# calculate the quantiles for those probabilities at t0
#
		mu0=a+b*t0
		qx=qnorm(px,mean=mu0,sd=s)

	list(predictedparameter=mu,adjustedx=qx)
}
#' Logf for RUST
#' @inherit manlogf return
#' @inheritParams manf
norm_p1_logf=function(params,x,t){
	a=params[1]
	b=params[2]
	s=pmax(params[3],.Machine$double.eps)
	mu=a+b*t
	logf=sum(dnorm(x,mean=mu,sd=s,log=TRUE))
	return(logf)
}
#' Normal-with-p1  observed log-likelihood function
#' @inherit manloglik return
#' @inheritParams manf
norm_p1_loglik=function(vv,x,t){
	n=length(x)
	mean=vv[1]+vv[2]*t #so mean is a vector, just like x
	loglik=sum(dnorm(x,mean=mean,sd=max(vv[3],0),log=TRUE))
	return(loglik)
}
#' Maximum likelihood estimator
#' @inherit manvector return
#' @inheritParams manf
norm_p1_mlparams=function(x,t){
	mlparams=matrix(0,3)
	reg=lm(x~t)
	mlparams[1]=reg$coefficients[1]
	mlparams[2]=reg$coefficients[2]
	mlparams[3]=sd(reg$residuals)
	return(mlparams)
}
#' Normal-with-p1 quantile function
#' @inherit manvector return
#' @inheritParams manf
qnorm_p1=function(p,t0,ymn,slope,sigma){

	return(qnorm(p,mean=(ymn+slope*t0),sd=sigma))

}
#' Normal-with-p1 density function
#' @inherit manvector return
#' @inheritParams manf
dnorm_p1=function(x,t0,ymn,slope,sigma,log=FALSE){

	return(dnorm(x,mean=(ymn+slope*t0),sd=sigma,log=log))

}
#' Normal-with-p1 distribution function
#' @inherit manvector return
#' @inheritParams manf
pnorm_p1=function(x,t0,ymn,slope,sigma){

	return(pnorm(x,mean=(ymn+slope*t0),sd=sigma))

}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inherit manlogscores return
#' @inheritParams manf
norm_p1_logscores=function(logscores,x,t){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){
			x1=x[-i]
			t1=t[-i]

			dd=dnorm_p1sub(x1,t1,x[i],t[i])
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
dnorm_p1sub=function(x,t,y,t0){

		nx=length(x)

		meant=mean(t)
		ta=t-meant
		ta0=t0-meant
# we have to centre this here again because this is used in a cross-validation loop
# note that t0 should never be centred

		ml_params=norm_p1_mlparams(x,t)
		v1hat=ml_params[1]
		v2hat=ml_params[2]
		v3hat=ml_params[3]
		muhat0=v1hat+v2hat*t0
# ml
		ml_pdf=dnorm(y,mean=muhat0,sd=v3hat)
		ml_cdf=pnorm(y,mean=muhat0,sd=v3hat)
# rhp

		rh_pdf=dnorm_p1_formula(y,ta,ta0,nx,muhat0,v3hat)
		rh_pdf=pmax(rh_pdf,0)

		rh_cdf=pnorm_p1_formula(y,ta,ta0,nx,muhat0,v3hat)
		rh_cdf=pmin(pmax(rh_cdf,0),1)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					rh_pdf=rh_pdf,
					ml_cdf=ml_cdf,
					rh_cdf=rh_cdf)
}

