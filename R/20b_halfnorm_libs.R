#' Waic
#' @inheritParams manf
halfnorm_waic=function(waicscores,x,v1hat,fd1,lddi,lddd,lambdad,aderivs){
		if(waicscores){

			if(aderivs) f1f=halfnorm_f1fa(x,v1hat)
			if(!aderivs)f1f=halfnorm_f1f(x,v1hat,fd1)

			if(aderivs) f2f=halfnorm_f2fa(x,v1hat)
			if(!aderivs)f2f=halfnorm_f2f(x,v1hat,fd1)

			fhatx=dhalfnorm(x,theta=v1hat)
			waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad,f2f,dim=1)
			waic1=waic$waic1
			waic2=waic$waic2
		}else{
			waic1="waicscores not selected"
			waic2="waicscores not selected"
		}
	list(waic1=waic1,waic2=waic2)
}
#' Logf for RUST
#' @inheritParams manf
halfnorm_logf=function(params,x){
	sc=pmax(params[1],.Machine$double.eps)
	logf=sum(dhalfnorm(x,theta=sc,log=TRUE))-log(sc)
	return(logf)
}
#' Log-likelihood function
#' @inheritParams manf
halfnorm_loglik=function(vv,x){
	loglik=sum(dhalfnorm(x,theta=max(vv[1],.Machine$double.eps),log=TRUE))
	return(loglik)
}
#' The second derivative of the normalized log-likelihood
#' @inheritParams manf
halfnorm_ldd=function(x,v1,fd1){
	nx=length(x)
	d1=fd1*v1
	v1m1=v1-1*d1
	v100=v1
	v1p1=v1+1*d1
	lm1=sum(log(dhalfnorm(x,theta=v1m1)))/nx
	l00=sum(log(dhalfnorm(x,theta=v100)))/nx
	lp1=sum(log(dhalfnorm(x,theta=v1p1)))/nx
	ldd=matrix(0,1,1)
	ldd[1,1]=(lp1-2*l00+lm1)/(d1*d1)
	return(ldd)
}
#' Second derivative of the expected log-likelihood
#' @inheritParams manf
halfnorm_gg11=function(alpha,v1,fd1){
  x=qhalfnorm((1-alpha),theta=v1)
	d1=fd1*v1
	lm=dhalfnorm(x,theta=v1-d1,log=TRUE)
	l0=dhalfnorm(x,theta=v1,log=TRUE)
	lp=dhalfnorm(x,theta=v1+d1,log=TRUE)
	gg=(lp-2*l0+lm)/(d1*d1)
	return(gg)
}
#' Expected information matrix
#' @inheritParams manf
halfnorm_gg=function(v1,fd1){
	expinfmat=matrix(0,1,1)
	expinfmat[1,1]=-quad(halfnorm_gg11,xa=0,xb=1,v1=v1,fd1=fd1)
 return(expinfmat)
}
#' Third derivative of the normalized log-likelihood
#' @inheritParams manf
halfnorm_l111=function(x,v1,fd1){
	nx=length(x)
	d1=fd1*v1
	v1m2=v1-2*d1
	v1m1=v1-1*d1
	v1p1=v1+1*d1
	v1p2=v1+2*d1
	lm2=sum(dhalfnorm(x,theta=v1m2,log=TRUE))/nx
	lm1=sum(dhalfnorm(x,theta=v1m1,log=TRUE))/nx
	lp1=sum(dhalfnorm(x,theta=v1p1,log=TRUE))/nx
	lp2=sum(dhalfnorm(x,theta=v1p2,log=TRUE))/nx
	dld111=(lp2-2*lp1+2*lm1-lm2)/(2*d1*d1*d1)
	return(dld111)
}
#' Third derivative tensor of the log-likelihood
#' @inheritParams manf
halfnorm_lddd=function(x,v1,fd1){
	nx=length(x)
	lddd=array(0,c(1,1,1))
	lddd[1,1,1]=halfnorm_l111(x,v1,fd1)
	return(lddd)
}
#' DMGS equation 2.1, f1 term
#' @inheritParams manf
halfnorm_f1f=function(y,v1,fd1){
	d1=fd1*v1
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v1 derivatives
	F1m1=dhalfnorm(y,theta=v1m1)
	F1p1=dhalfnorm(y,theta=v1p1)
	f1=matrix(0,1,length(y))
	f1[1,]=(F1p1-F1m1)/(2*d1)
	return(f1)
}
#' DMGS equation 2.1, p1 term
#' @inheritParams manf
halfnorm_p1f=function(y,v1,fd1){
	d1=fd1*v1
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v1 derivatives
	F1m1=phalfnorm(y,theta=v1m1)
	F1p1=phalfnorm(y,theta=v1p1)
	p1=matrix(0,1,length(y))
	p1[1,]=(F1p1-F1m1)/(2*d1)
	return(p1)
}
#' DMGS equation 3.3, mu1 term
#' @inheritParams manf
halfnorm_mu1f=function(alpha,v1,fd1){
	q00=qhalfnorm((1-alpha),theta=v1)
	d1=fd1*v1
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v1 derivatives
	F1m1=phalfnorm(q00,theta=v1m1)
	F1p1=phalfnorm(q00,theta=v1p1)
	mu1=matrix(0,1,length(alpha))
	mu1[1,]=-(F1p1-F1m1)/(2*d1)
	return(mu1)
}
#' DMGS equation 2.1, f2 term
#' @inheritParams manf
halfnorm_f2f=function(y,v1,fd1){
	d1=fd1*v1
# v1 stuff
	v1m2=v1-2*d1
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
	v1p2=v1+2*d1
	F1m2=dhalfnorm(y,theta=v1m2)
	F1m1=dhalfnorm(y,theta=v1m1)
	F100=dhalfnorm(y,theta=v100)
	F1p1=dhalfnorm(y,theta=v1p1)
	F1p2=dhalfnorm(y,theta=v1p2)
	f2=array(0,c(1,1,length(y)))
	f2[1,1,]=(F1p1-2*F100+F1m1)/(d1*d1)
	return(f2)
}
#' DMGS equation 2.1, p2 term
#' @inheritParams manf
halfnorm_p2f=function(y,v1,fd1){
	d1=fd1*v1
# v1 stuff
	v1m2=v1-2*d1
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
	v1p2=v1+2*d1
	F1m2=phalfnorm(y,theta=v1m2)
	F1m1=phalfnorm(y,theta=v1m1)
	F100=phalfnorm(y,theta=v100)
	F1p1=phalfnorm(y,theta=v1p1)
	F1p2=phalfnorm(y,theta=v1p2)
	p2=array(0,c(1,1,length(y)))
	p2[1,1,]=(F1p1-2*F100+F1m1)/(d1*d1)
	return(p2)
}
#' DMGS equation 3.3, mu2 term
#' @inheritParams manf
halfnorm_mu2f=function(alpha,v1,fd1){
	q00=qhalfnorm((1-alpha),theta=v1)
	d1=fd1*v1
# v1 stuff
	v1m2=v1-2*d1
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
	v1p2=v1+2*d1
	mu2=array(0,c(1,1,length(alpha)))
	F1m2=phalfnorm(q00,theta=v1m2)
	F1m1=phalfnorm(q00,theta=v1m1)
	F100=phalfnorm(q00,theta=v100)
	F1p1=phalfnorm(q00,theta=v1p1)
	F1p2=phalfnorm(q00,theta=v1p2)
	mu2[1,1,]=-(F1p1-2*F100+F1m1)/(d1*d1)
	return(mu2)
}
#' MLE and RHP predictive means
#' RHP mean based on the expectation of DMGS equation 2.1
#' @inheritParams manf
halfnorm_means=function(means,ml_params,lddi,lddd,lambdad_rhp,nx,dim=1){

	if(means){
# ml mean (remembering that the parameter is theta=1/sigma)
		cc=sqrt(2)/sqrt(3.1415926535)
#		ml_mean=cc/ml_params[1]
		theta=ml_params[1]
		theta2=theta*theta
		theta3=theta2*theta
		ml_mean=1/theta

# rhp mean
		meand1=array(-1/theta2,c(1,1))
		meand2=array(1/theta3,c(1,1,1))
		dmean=dmgs(lddi,lddd,meand1,lambdad_rhp,meand2,dim=1)
		rh_mean=ml_mean+dmean/nx

	} else{
		ml_mean="means not selected"
		rh_mean="means not selected"
	}

	list(ml_mean=ml_mean,rh_mean=rh_mean)

}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inheritParams manf
halfnorm_logscores=function(logscores,x,fd1=0.01,aderivs=TRUE){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){
			x1=x[-i]

			dd=dhalfnormsub(x1,x[i],fd1,aderivs)

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
dhalfnormsub=function(x,y,fd1=0.01,aderivs=TRUE){

		nx=length(x)

		v1start=sqrt((sum(x*x))/nx)
		opt=optim(c(v1start),halfnorm_loglik,x=x,method="Brent",
			lower=.Machine$double.eps,upper=999999999999,control=list(fnscale=-1))
		v1hat=opt$par[1]
		ml_params=c(v1hat)
# ml
		ml_pdf=dhalfnorm(y,theta=v1hat)
		ml_cdf=phalfnorm(y,theta=v1hat)

# rhp
		if(aderivs)	ldd=halfnorm_ldda(x,v1hat)
		if(!aderivs)ldd=halfnorm_ldd(x,v1hat,fd1)
		lddi=1/ldd
		if(aderivs)lddd=halfnorm_lddda(x,v1hat)
		if(!aderivs)lddd=halfnorm_lddd(x,v1hat,fd1)

		if(aderivs) f1=halfnorm_f1fa(y,v1hat)
		if(!aderivs)f1=halfnorm_f1f(y,v1hat,fd1)

		if(aderivs) f2=halfnorm_f2fa(y,v1hat)
		if(!aderivs)f2=halfnorm_f2f(y,v1hat,fd1)

		p1=halfnorm_p1f(y,v1hat,fd1)
		p2=halfnorm_p2f(y,v1hat,fd1)
		lambdad_rhp=matrix(-1/v1hat)
		df1=dmgs(lddi,lddd,f1,lambdad_rhp,f2,dim=1)
		dp1=dmgs(lddi,lddd,p1,lambdad_rhp,p2,dim=1)
		rh_pdf=pmax(ml_pdf+df1/nx,0)
		rh_cdf=pmin(pmax(ml_cdf+dp1/nx,0),1)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					rh_pdf=rh_pdf,
					ml_cdf=ml_cdf,
					rh_cdf=rh_cdf)
}

