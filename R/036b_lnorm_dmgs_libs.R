#'  Waic
#' @inherit manwaic return
#' @inheritParams manf
lnorm_dmgs_waic=function(waicscores,x,v1hat,v2hat,lddi,lddd,lambdad){
		if(waicscores){

			f1f=lnorm_f1fa(x,v1hat,v2hat)
			f2f=lnorm_f2fa(x,v1hat,v2hat)

			fhatx=dlnorm(x,meanlog=v1hat,sdlog=v2hat)
			waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad,f2f,dim=2)
			waic1=waic$waic1
			waic2=waic$waic2
		}else{
			waic1="waicscores not selected"
			waic2="waicscores not selected"
		}
		list(waic1=waic1,waic2=waic2)
}
#'  log-likelihood function
#' @inherit manlogf return
#' @inheritParams manf
lnorm_dmgs_loglik=function(vv,x){
	loglik=sum(dlnorm(x,meanlog=vv[1],sdlog=max(vv[2],.Machine$double.eps),log=TRUE))
	return(loglik)
}
#' MLE and RHP predictive means
#' @inherit manmeans return
#' @inheritParams manf
lnorm_dmgs_means=function(means,ml_params,lddi,lddd,lambdad_rhp,nx,dim=2){

	if(means){
# intro
	  v1=ml_params[1]
	 v2=ml_params[2]

# ml mean
		ml_mean=exp(v1+0.5*v2*v2)

# rhp mean
		meand1=array(0,c(2,1))
		meand1[1,1]=ml_mean
		meand1[2,1]=v2*ml_mean
		meand2=array(0,c(2,2,1)) #but all zero for lnorm_dmgs
		meand2[1,1,1]=ml_mean
		meand2[1,2,1]=v2*ml_mean
		meand2[2,1,1]=meand2[1,2,1]
		meand2[2,2,1]=v2*v2*ml_mean
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
lnorm_dmgs_logscores=function(logscores,x){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){

			x1=x[-i]

			dd=dlnorm_dmgssub(x1,x[i])

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
dlnorm_dmgssub=function(x,y){

		nx=length(x)

		v1start=mean(x)
		v2start=sd(x)
		opt=optim(c(v1start,v2start),lnorm_dmgs_loglik,x=x,control=list(fnscale=-1))
		v1hat=opt$par[1]
		v2hat=opt$par[2]
		ml_params=c(v1hat,v2hat)

# mle
		ml_pdf=dlnorm(y,meanlog=v1hat,sdlog=v2hat)
		ml_cdf=plnorm(y,meanlog=v1hat,sdlog=v2hat)

# rhp
		ldd=lnorm_ldda(x,v1hat,v2hat)
		lddi=solve(ldd)
		lddd=lnorm_lddda(x,v1hat,v2hat)

		f1=lnorm_f1fa(y,v1hat,v2hat)
		f2=lnorm_f2fa(y,v1hat,v2hat)

		p1=lnorm_p1fa(y,v1hat,v2hat)
		p2=lnorm_p2fa(y,v1hat,v2hat)

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

