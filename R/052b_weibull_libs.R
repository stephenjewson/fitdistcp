#' Waic for RUST
#' @inherit manwaic return
#' @inheritParams manf
weibull_waic=function(waicscores,x,v1hat,v2hat,lddi,lddd,lambdad){
		if(waicscores){

			f1f=weibull_f1fa(x,v1hat,v2hat)

			f2f=weibull_f2fa(x,v1hat,v2hat)

			fhatx=dweibull(x,shape=v1hat,scale=v2hat)
			waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad,f2f,dim=2)
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
weibull_logf=function(params,x){
	sh=pmax(min(20,params[1]),sqrt(.Machine$double.eps))
	sc=pmax(params[2],sqrt(.Machine$double.eps))
	logf=sum(dweibull(x,shape=sh,scale=sc,log=TRUE))-log(sh)-log(sc)
	return(logf)
}
#'  log-likelihood function
#' @inherit manloglik return
#' @inheritParams manf
weibull_loglik=function(vv,x){
	n=length(x)
	loglik=sum(dweibull(x,shape=max(min(20,vv[1]),sqrt(.Machine$double.eps)),scale=max(vv[2],sqrt(.Machine$double.eps)),log=TRUE))
	return(loglik)
}
#' MLE and RHP predictive means
#' @inherit manmeans return
#' @inheritParams manf
weibull_means=function(means,ml_params,lddi,lddd,lambdad_rhp,nx,dim=2){
# v1 is shape
# v2 is scale
	if(means){
		v1=ml_params[1]
		v2=ml_params[2]
		f1=1+(1/v1)
		iv12=1/(v1*v1)
		iv13=1/(v1*v1*v1)

# ml mean
		ml_mean=v2*gamma(f1)

# rhp mean

		meand1=array(0,c(2,1))
		meand1[1,1]=-v2*gamma(f1)*digamma(f1)*iv12
		meand1[2,1]=gamma(f1)

		meand2=array(0,c(2,2,1))
		meand2[1,1,1]=v2*(gamma(f1)*trigamma(f1)*iv12+2*gamma(f1)*digamma(f1)*iv13)
		meand2[1,2,1]=-gamma(f1)*digamma(f1)*iv12
		meand2[2,1,1]=meand2[1,2,1]
		meand2[2,2,1]=0

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
#' @inheritParams manf
weibull_logscores=function(logscores,x){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){
			x1=x[-i]

			dd=dweibullsub(x1,x[i])

			ml_params=dd$ml_params

			ml_pdf=dd$ml_pdf
			ml_oos_logscore=ml_oos_logscore+log(ml_pdf)

			rh_pdf=dd$rh_pdf
			rh_oos_logscore=rh_oos_logscore+log(max(rh_pdf,.Machine$double.eps))
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
dweibullsub=function(x,y){

		nx=length(x)

		opt=optim(c(1,1),weibull_loglik,x=x,control=list(fnscale=-1))
		v1hat=opt$par[1]
		v2hat=opt$par[2]
		ml_params=c(v1hat,v2hat)

# mle
		ml_pdf=dweibull(y,shape=v1hat,scale=v2hat)
		ml_cdf=pweibull(y,shape=v1hat,scale=v2hat)

# rhp
		ldd=weibull_ldda(x,v1hat,v2hat)
		lddi=solve(ldd)
		lddd=weibull_lddda(x,v1hat,v2hat)

		f1=weibull_f1fa(y,v1hat,v2hat)
		f2=weibull_f2fa(y,v1hat,v2hat)

		p1=weibull_p1fa(y,v1hat,v2hat)
		p2=weibull_p2fa(y,v1hat,v2hat)

		lambdad_rhp=c(-1/v1hat,-1/v2hat)
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
