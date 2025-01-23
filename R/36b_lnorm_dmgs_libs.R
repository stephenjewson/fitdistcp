#'  Waic
#' @inheritParams manf
lnorm_dmgs_waic=function(waicscores,x,v1hat,d1,v2hat,fd2,lddi,lddd,lambdad,aderivs){
		if(waicscores){

			if(aderivs) f1f=lnorm_f1fa(x,v1hat,v2hat)
			if(!aderivs)f1f=lnorm_f1f(x,v1hat,d1,v2hat,fd2)

			if(aderivs) f2f=lnorm_f2fa(x,v1hat,v2hat)
			if(!aderivs)f2f=lnorm_f2f(x,v1hat,d1,v2hat,fd2)

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
#' @inheritParams manf
lnorm_dmgs_loglik=function(vv,x){
	loglik=sum(dlnorm(x,meanlog=vv[1],sdlog=max(vv[2],.Machine$double.eps),log=TRUE))
	return(loglik)
}
#' One component of the second derivative of the expected log-likelihood
#' @inheritParams manf
lnorm_dmgs_gg11=function(alpha,v1,d1,v2,fd2){
  x=qlnorm((1-alpha),meanlog=v1,sdlog=v2)
  d2=v2*fd2
	l0=dlnorm(x,meanlog=v1,   sdlog=v2,log=TRUE)
	lm=dlnorm(x,meanlog=v1-d1,sdlog=v2,log=TRUE)
	lp=dlnorm(x,meanlog=v1+d1,sdlog=v2,log=TRUE)
	d2ld2=(lp-2*l0+lm)/(d1*d1)
	integrand=d2ld2
	return(integrand)
}
#' One component of the second derivative of the expected log-likelihood
#' @inheritParams manf
lnorm_dmgs_gg12=function(alpha,v1,d1,v2,fd2){
  x=qlnorm((1-alpha),meanlog=v1,sdlog=v2)
  d2=v2*fd2
	lmm=dlnorm(x,meanlog=v1-d1,sdlog=v2-d2,log=TRUE)
	lpm=dlnorm(x,meanlog=v1+d1,sdlog=v2-d2,log=TRUE)
	lmp=dlnorm(x,meanlog=v1-d1,sdlog=v2+d2,log=TRUE)
	lpp=dlnorm(x,meanlog=v1+d1,sdlog=v2+d2,log=TRUE)
	d2ld12=(lpp-lmp-lpm+lmm)/(4*d1*d2)
	integrand=d2ld12
	return(integrand)
}
#' One component of the second derivative of the expected log-likelihood
#' @inheritParams manf
lnorm_dmgs_gg22=function(alpha,v1,d1,v2,fd2){
  x=qlnorm((1-alpha),meanlog=v1,sdlog=v2)
  d2=v2*fd2
	l0=dlnorm(x,meanlog=v1,sdlog=v2,log=TRUE)
	lm=dlnorm(x,meanlog=v1,sdlog=v2-d2,log=TRUE)
	lp=dlnorm(x,meanlog=v1,sdlog=v2+d2,log=TRUE)
	d2ld2=(lp-2*l0+lm)/(d2*d2)
	integrand=d2ld2
	return(integrand)
}
#' DMGS equation 3.3, p1 term
#' @inheritParams manf
lnorm_dmgs_p1f=function(y,v1,d1,v2,fd2){
	d2=fd2*v2
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v2 stuff
	v2m1=v2-1*d2
	v200=v2+0*d2
	v2p1=v2+1*d2
# v1 derivatives
	F1m1=plnorm(y,meanlog=v1m1,sdlog=v200)
	F1p1=plnorm(y,meanlog=v1p1,sdlog=v200)
# v2 derivatives
	F2m1=plnorm(y,meanlog=v100,sdlog=v2m1)
	F2p1=plnorm(y,meanlog=v100,sdlog=v2p1)
	p1=matrix(0,2,length(y))
	p1[1,]=(F1p1-F1m1)/(2*d1)
	p1[2,]=(F2p1-F2m1)/(2*d2)
	return(p1)
}
#' DMGS equation 3.3, mu1 term
#' @inheritParams manf
lnorm_dmgs_mu1f=function(alpha,v1,d1,v2,fd2){
	q00=qlnorm((1-alpha),meanlog=v1,sdlog=v2)
	d2=fd2*v2
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v2 stuff
	v2m1=v2-1*d2
	v200=v2+0*d2
	v2p1=v2+1*d2
# v1 derivatives
	F1m1=plnorm(q00,meanlog=v1m1,sdlog=v200)
	F1p1=plnorm(q00,meanlog=v1p1,sdlog=v200)
# v2 derivatives
	F2m1=plnorm(q00,meanlog=v100,sdlog=v2m1)
	F2p1=plnorm(q00,meanlog=v100,sdlog=v2p1)
	mu1=matrix(0,2,length(alpha))
	mu1[1,]=-(F1p1-F1m1)/(2*d1)
	mu1[2,]=-(F2p1-F2m1)/(2*d2)
	return(mu1)
}
#' DMGS equation 3.3, p2 term
#' @inheritParams manf
lnorm_dmgs_p2f=function(y,v1,d1,v2,fd2){
	d2=fd2*v2
# v1 stuff
	v1m2=v1-2*d1
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
	v1p2=v1+2*d1
# v2 stuff
	v2m2=v2-2*d2
	v2m1=v2-1*d2
	v200=v2+0*d2
	v2p1=v2+1*d2
	v2p2=v2+2*d2
	F1m2=plnorm(y,meanlog=v1m2,sdlog=v200)
	F1m1=plnorm(y,meanlog=v1m1,sdlog=v200)
	F100=plnorm(y,meanlog=v100,sdlog=v200)
	F1p1=plnorm(y,meanlog=v1p1,sdlog=v200)
	F1p2=plnorm(y,meanlog=v1p2,sdlog=v200)
# v2 derivative
	F2m2=plnorm(y,meanlog=v100,sdlog=v2m2)
	F2m1=plnorm(y,meanlog=v100,sdlog=v2m1)
	F200=plnorm(y,meanlog=v100,sdlog=v200)
	F2p1=plnorm(y,meanlog=v100,sdlog=v2p1)
	F2p2=plnorm(y,meanlog=v100,sdlog=v2p2)
# cross derivative
	Fcm1m1=plnorm(y,meanlog=v1m1,sdlog=v2m1)
	Fcm1p1=plnorm(y,meanlog=v1m1,sdlog=v2p1)
	Fcp1m1=plnorm(y,meanlog=v1p1,sdlog=v2m1)
	Fcp1p1=plnorm(y,meanlog=v1p1,sdlog=v2p1)
	p2=array(0,c(2,2,length(y)))
	p2[1,1,]=(F1p1-2*F100+F1m1)/(d1*d1)
	p2[2,2,]=(F2p1-2*F200+F2m1)/(d2*d2)
	p2[1,2,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2)
	# copy
	p2[2,1,]=p2[1,2,]
	return(p2)
}
#' DMGS equation 3.3, mu2 term
#' @inheritParams manf
lnorm_dmgs_mu2f=function(alpha,v1,d1,v2,fd2){
	q00=qlnorm((1-alpha),meanlog=v1,sdlog=v2)
	d2=fd2*v2
# v1 stuff
	v1m2=v1-2*d1
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
	v1p2=v1+2*d1
# v2 stuff
	v2m2=v2-2*d2
	v2m1=v2-1*d2
	v200=v2+0*d2
	v2p1=v2+1*d2
	v2p2=v2+2*d2
	mu2=array(0,c(2,2,length(alpha)))
	F1m2=plnorm(q00,meanlog=v1m2,sdlog=v200)
	F1m1=plnorm(q00,meanlog=v1m1,sdlog=v200)
	F100=plnorm(q00,meanlog=v100,sdlog=v200)
	F1p1=plnorm(q00,meanlog=v1p1,sdlog=v200)
	F1p2=plnorm(q00,meanlog=v1p2,sdlog=v200)
# v2 derivative
	F2m2=plnorm(q00,meanlog=v100,sdlog=v2m2)
	F2m1=plnorm(q00,meanlog=v100,sdlog=v2m1)
	F200=plnorm(q00,meanlog=v100,sdlog=v200)
	F2p1=plnorm(q00,meanlog=v100,sdlog=v2p1)
	F2p2=plnorm(q00,meanlog=v100,sdlog=v2p2)
# cross derivative
	Fcm1m1=plnorm(q00,meanlog=v1m1,sdlog=v2m1)
	Fcm1p1=plnorm(q00,meanlog=v1m1,sdlog=v2p1)
	Fcp1m1=plnorm(q00,meanlog=v1p1,sdlog=v2m1)
	Fcp1p1=plnorm(q00,meanlog=v1p1,sdlog=v2p1)
	mu2[1,1,]=-(F1p1-2*F100+F1m1)/(d1*d1)
	mu2[2,2,]=-(F2p1-2*F200+F2m1)/(d2*d2)
	mu2[1,2,]=-(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2)
	# copy
	mu2[2,1,]=mu2[1,2,]
	return(mu2)
}
#' MLE and RHP predictive means
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
#' @inheritParams manf
lnorm_dmgs_logscores=function(logscores,x,d1=0.01,fd2=0.01){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){

			x1=x[-i]

			dd=dlnorm_dmgssub(x1,x[i],d1,fd2)

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
#' @inheritParams manf
dlnorm_dmgssub=function(x,y,d1=0.01,fd2=0.01,aderivs=TRUE){

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
		if(aderivs)	ldd=lnorm_ldda(x,v1hat,v2hat)
		if(!aderivs)ldd=lnorm_ldd(x,v1hat,d1,v2hat,fd2)
		lddi=solve(ldd)
		if(aderivs)	lddd=lnorm_lddda(x,v1hat,v2hat)
		if(!aderivs)lddd=lnorm_lddd(x,v1hat,d1,v2hat,fd2)

		if(aderivs) f1=lnorm_f1fa(y,v1hat,v2hat)
		if(!aderivs)f1=lnorm_f1f(y,v1hat,d1,v2hat,fd2)

		if(aderivs) f2=lnorm_f2fa(y,v1hat,v2hat)
		if(!aderivs)f2=lnorm_f2f(y,v1hat,d1,v2hat,fd2)

		if(aderivs) p1=lnorm_p1fa(y,v1hat,v2hat)
		if(!aderivs)p1=lnorm_dmgs_p1f(y,v1hat,d1,v2hat,fd2)

		if(aderivs) p2=lnorm_p2fa(y,v1hat,v2hat)
		if(!aderivs)p2=lnorm_dmgs_p2f(y,v1hat,d1,v2hat,fd2)

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

