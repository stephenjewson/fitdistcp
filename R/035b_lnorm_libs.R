#' Waic for RUST
#' @inherit manwaic return
#' @inheritParams manf
lnorm_waic=function(waicscores,x,v1hat,v2hat){
	if(waicscores){

		f1f=lnorm_f1fa(x,v1hat,v2hat)
		f2f=lnorm_f2fa(x,v1hat,v2hat)

		ldd=lnorm_ldda(x,v1hat,v2hat)
		lddi=solve(ldd)
		lddd=lnorm_lddda(x,v1hat,v2hat)
		fhatx=dlnorm(x,meanlog=v1hat,sdlog=v2hat)
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
lnorm_logf=function(params,x){
	m=params[1]
	s=pmax(params[2],.Machine$double.eps)
	logf=sum(dlnorm(x,meanlog=m,sdlog=s,log=TRUE))-log(s)
	return(logf)
}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inherit manlogscores return
#' @inheritParams manf
lnorm_logscores=function(logscores,x){

	if(logscores){
		y=log(x)
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){

			x1=x[-i]

			dd=dlnormsub(x1,x[i])

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
dlnormsub=function(x,y){

# great potential for confusion here
# x and y are both lognormal
# logx and logy are both normal

		nx=length(x)
		logx=log(x)
		logy=log(y)
# ml
		ml_params=norm_ml_params(logx) #this really should be norm not lnorm
		ml_pdf=dlnorm(y,meanlog=ml_params[1],sdlog=ml_params[2])
		ml_cdf=plnorm(y,meanlog=ml_params[1],sdlog=ml_params[2])

# rhp pdf
# -first, convert sigma from maxlik to unbiased
		sgu=ml_params[2]*sqrt(nx/(nx-1))
# then, convert sigma to predictive sigma
		sg1=sgu*sqrt((nx+1)/nx)

		logyd=(logy-ml_params[1])/sg1
		rh_pdf=(dt(logyd,df=nx-1))/(y*sg1) #note the extra exp term here, to convert to loglnormal predictive density
		rh_pdf=pmax(rh_pdf,0)

# rhp cdf
		rh_cdf=pt(logyd,df=nx-1)
		rh_cdf=pmin(pmax(rh_cdf,0),1)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					rh_pdf=rh_pdf,
					ml_cdf=ml_cdf,
					rh_cdf=rh_cdf)
}
