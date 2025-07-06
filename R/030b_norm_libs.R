#' Waic
#' @inherit manwaic return
#' @inheritParams manf
norm_waic=function(waicscores,x,v1hat,v2hat){
	if(waicscores){

		f1f=norm_f1fa(x,v1hat,v2hat)
		f2f=norm_f2fa(x,v1hat,v2hat)
		ldd=norm_ldda(x,v1hat,v2hat)
		lddi=solve(ldd)
		lddd=norm_lddda(x,v1hat,v2hat)
		fhatx=dnorm(x,mean=v1hat,sd=v2hat)
		lambdad_rhp=c(0,-1/v2hat)
		waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad_rhp,f2f,dim=2)
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
norm_logf=function(params,x){
	m=params[1]
	s=pmax(params[2],.Machine$double.eps)
	logf=sum(dnorm(x,mean=m,sd=s,log=TRUE))-log(s)
	return(logf)
}
#' Maximum likelihood estimator
#' @inherit manloglik return
#' @inheritParams manf
norm_ml_params=function(x){
	mlparams=matrix(0,2)
	nx=length(x)
	mlparams[1]=mean(x)
	mlparams[2]=sqrt((nx-1)/nx)*sd(x)
	return(mlparams)
}
#' Method of moments estimator
#' @return Vector
#' @inheritParams manf
norm_unbiasedv_params=function(x){
	params=matrix(0,2)
	nx=length(x)
	params[1]=mean(x)
	params[2]=sd(x)
	return(params)
}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inherit manlogscores return
#' @inheritParams manf
norm_logscores=function(logscores,x){

	if(logscores){

		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){
			x1=x[-i]

			dd=dnormsub(x1,x[i])
			ml_params=dd$ml_params

			ml_pdf=dd$ml_pdf
			ml_oos_logscore=ml_oos_logscore+log(ml_pdf)

			rh_pdf=dd$rh_pdf
			rh_oos_logscore=rh_oos_logscore+log(rh_pdf)
		}
	}else{
		ml_oos_logscore="logscores not selected"
		rh_oos_logscore="logscores not selected"
	}
	list(ml_oos_logscore=ml_oos_logscore,rh_oos_logscore=rh_oos_logscore)
}
#' Densities from MLE and RHP
#' @inherit mandsub return
#' @inheritParams manf
dnormsub=function(x,y){

		nx=length(x)

# ml
		ml_params=norm_ml_params(x)
		ml_pdf=dnorm(y,mean=ml_params[1],sd=ml_params[2])
		ml_cdf=pnorm(y,mean=ml_params[1],sd=ml_params[2])

# rhp pdf
# -first, convert sigma from maxlik to unbiased
		sgu=ml_params[2]*sqrt(nx/(nx-1))
# then, convert sigma to predictive sigma
		sg1=sgu*sqrt((nx+1)/nx)

		yd=(y-ml_params[1])/sg1
		rh_pdf=(dt(yd,df=nx-1))/sg1
		rh_pdf=pmax(rh_pdf,0)

# rhp cdf
		rh_cdf=pt(yd,df=nx-1)
		rh_cdf=pmin(pmax(rh_cdf,0),1)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					rh_pdf=rh_pdf,
					ml_cdf=ml_cdf,
					rh_cdf=rh_cdf)
}
