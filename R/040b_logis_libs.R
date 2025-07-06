#' Waic
#' @inherit manwaic return
#' @inheritParams manf
logis_waic=function(waicscores,x,v1hat,v2hat,lddi,lddd,lambdad){
		if(waicscores){

			f1f=logis_f1fa(x,v1hat,v2hat)
			f2f=logis_f2fa(x,v1hat,v2hat)

			fhatx=dlogis(x,location=v1hat,scale=v2hat)
			waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad,f2f,dim=2)
			waic1=waic$waic1
			waic2=waic$waic2
		}else{
			waic1="extras not selected"
			waic2="extras not selected"
		}
		list(waic1=waic1,waic2=waic2)
}
#' Logf for RUST
#' @inherit manlogf return
#' @inheritParams manf
logis_logf=function(params,x){
#	l=params[1]
#	s=params[2]
#	if(s>0){
#		logf=sum(dlogis(x,location=l,scale=s,log=TRUE))-log(s)
#	}else{
#		logf=-Inf
#	}
	l=params[1]
	s=pmax(params[2],.Machine$double.eps)
	logf=sum(dlogis(x,location=l,scale=s,log=TRUE))-log(s)
	return(logf)
}
#'  log-likelihood function
#' @inherit manloglik return
#' @inheritParams manf
logis_loglik=function(vv,x){
	n=length(x)
	loglik=sum(dlogis(x,location=vv[1],scale=max(vv[2],.Machine$double.eps),log=TRUE))
	return(loglik)
}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inherit manlogscores return
#' @inheritParams manf
logis_logscores=function(logscores,x){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){
			x1=x[-i]

			dd=dlogis2sub(x1,x[i])

			ml_params=dd$ml_params

			ml_pdf=dd$ml_pdf
			ml_oos_logscore=ml_oos_logscore+log(ml_pdf)

			rh_pdf=dd$rh_pdf
			rh_oos_logscore=rh_oos_logscore+log(max(rh_pdf,.Machine$double.eps))
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
dlogis2sub=function(x,y){

		nx=length(x)

		v1start=mean(x)
		v2start=sd(x)
		opt=optim(c(v1start,v2start),logis_loglik,x=x,control=list(fnscale=-1))
		v1hat=opt$par[1]
		v2hat=opt$par[2]
		ml_params=c(v1hat,v2hat)

# mle
		ml_pdf=dlogis(y,location=v1hat,scale=v2hat)
		ml_cdf=plogis(y,location=v1hat,scale=v2hat)

# rhp
		ldd=logis_ldda(x,v1hat,v2hat)
		lddi=solve(ldd)
		lddd=logis_lddda(x,v1hat,v2hat)

		f1=logis_f1fa(y,v1hat,v2hat)

		f2=logis_f2fa(y,v1hat,v2hat)

		p1=logis_p1fa(y,v1hat,v2hat)

		p2=logis_p2fa(y,v1hat,v2hat)

		lambdad_rhp=c(0,-1/v2hat)
		df1=dmgs(lddi,lddd,f1,lambdad_rhp,f2,dim=2)
		dp1=dmgs(lddi,lddd,p1,lambdad_rhp,p2,dim=2)
		rh_pdf=pmax(ml_pdf+df1/nx,0)
		rh_cdf=pmin(pmax(ml_cdf+dp1/nx,0),1)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					rh_pdf=rh_pdf,
					ml_cdf=ml_cdf,
					rh_cdf=rh_cdf)
}
