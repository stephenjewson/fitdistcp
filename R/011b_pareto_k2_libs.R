# extraDistr is very confusing:
# a = shape parameter (they call scale), that I'm varying here
# b = scale parameter (they call location), that I'm calling kscale
#
#' Waic
#' @inherit manwaic return
#' @inheritParams manf
pareto_k2_waic=function(waicscores,x,v1hat,fd1,kscale,aderivs){
	if(waicscores){
		if(aderivs) f1f=pareto_k2_f1fa(x,v1hat,kscale=kscale)
		if(!aderivs)f1f=pareto_k2_f1f(x,v1hat,fd1,kscale=kscale)

		if(aderivs) f2f=pareto_k2_f2fa(x,v1hat,kscale=kscale)
		if(!aderivs)f2f=pareto_k2_f2f(x,v1hat,fd1,kscale=kscale)
		if(aderivs)	ldd=pareto_k2_ldda(x,v1hat,kscale=kscale)
		if(!aderivs)ldd=pareto_k2_ldd(x,v1hat,fd1,kscale=kscale)
		lddi=solve(ldd)
		if(aderivs)	lddd=pareto_k2_lddda(x,v1hat,kscale=kscale)
		if(!aderivs)lddd=pareto_k2_lddd(x,v1hat,fd1,kscale=kscale)
		fhatx=extraDistr::dpareto(x,a=v1hat,b=kscale)
		lambdad_rhp=-1/v1hat
		waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad_rhp,f2f,dim=1)
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
pareto_k2_logf=function(params,x,kscale){
	sh=pmax(params[1],.Machine$double.eps)
	logf=sum(extraDistr::dpareto(x,a=sh,b=kscale,log=TRUE))-log(sh)
	return(logf)
}
#' Maximum likelihood estimator
#' @inherit manloglik return
#' @inheritParams manf
pareto_k2_ml_params=function(x,kscale){
	nx=length(x)
#	lnx=log(x)
	lnx=log(x/kscale)
	slnx=sum(lnx)
	mlparams=nx/slnx
	return(mlparams)
}
#' The second derivative of the normalized log-likelihood
#' @inherit manldd return
#' @inheritParams manf
pareto_k2_ldd=function(x,v1,fd1,kscale){
	nx=length(x)
	d1=fd1*v1
	v1m1=v1-1*d1
	v100=v1
	v1p1=v1+1*d1
	lm1=sum(log(extraDistr::dpareto(x,b=kscale,a=v1m1)))/nx
	l00=sum(log(extraDistr::dpareto(x,b=kscale,a=v100)))/nx
	lp1=sum(log(extraDistr::dpareto(x,b=kscale,a=v1p1)))/nx
	ldd=matrix(0,1,1)
	ldd[1,1]=(lp1-2*l00+lm1)/(d1*d1)
	return(ldd)
}
#' Third derivative of the normalized log-likelihood
#' @inherit manlnnn return
#' @inheritParams manf
pareto_k2_l111=function(x,v1,fd1,kscale){
	nx=length(x)
	d1=fd1*v1
	v1m2=v1-2*d1
	v1m1=v1-1*d1
	v1p1=v1+1*d1
	v1p2=v1+2*d1
	lm2=sum(extraDistr::dpareto(x,b=kscale,a=v1m2,log=TRUE))/nx
	lm1=sum(extraDistr::dpareto(x,b=kscale,a=v1m1,log=TRUE))/nx
	lp1=sum(extraDistr::dpareto(x,b=kscale,a=v1p1,log=TRUE))/nx
	lp2=sum(extraDistr::dpareto(x,b=kscale,a=v1p2,log=TRUE))/nx
	dld111=(lp2-2*lp1+2*lm1-lm2)/(2*d1*d1*d1)
	return(dld111)
}
#' Third derivative tensor of the log-likelihood
#' @inherit manlddd return
#' @inheritParams manf
pareto_k2_lddd=function(x,v1,fd1,kscale){
	nx=length(x)
	lddd=array(0,c(1,1,1))
	lddd[1,1,1]=pareto_k2_l111(x,v1,fd1,kscale)
	return(lddd)
}
#' DMGS equation 2.1, f1 term
#' @inherit man1f return
#' @inheritParams manf
pareto_k2_f1f=function(y,v1,fd1,kscale){
	d1=fd1*v1
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v1 derivatives
	F1m1=extraDistr::dpareto(y,b=kscale,a=v1m1)
	F1p1=extraDistr::dpareto(y,b=kscale,a=v1p1)
	f1=matrix(0,1,length(y))
	f1[1,]=(F1p1-F1m1)/(2*d1)
	return(f1)
}
#' DMGS equation 2.1, f2 term
#' @inherit man2f return
#' @inheritParams manf
pareto_k2_f2f=function(y,v1,fd1,kscale){
	d1=fd1*v1
# v1 stuff
	v1m2=v1-2*d1
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
	v1p2=v1+2*d1
	F1m2=extraDistr::dpareto(y,b=kscale,a=v1m2)
	F1m1=extraDistr::dpareto(y,b=kscale,a=v1m1)
	F100=extraDistr::dpareto(y,b=kscale,a=v100)
	F1p1=extraDistr::dpareto(y,b=kscale,a=v1p1)
	F1p2=extraDistr::dpareto(y,b=kscale,a=v1p2)
	f2=array(0,c(1,1,length(y)))
	f2[1,1,]=(F1p1-2*F100+F1m1)/(d1*d1)
	return(f2)
}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inherit manlogscores return
#' @inheritParams manf
pareto_k2_logscores=function(logscores,x,kscale){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){
			x1=x[-i]
			dd=dpareto_k2_sub(x1,x[i],kscale)
			ml_params1=dd$ml_params
# ml
			ml_pdf=dd$ml_pdf
			ml_oos_logscore=ml_oos_logscore+log(ml_pdf)
# rhp
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
dpareto_k2_sub=function(x,y,kscale,aderivs=TRUE){

		nx=length(x)

# ml
		ml_params=pareto_k2_ml_params(x,kscale)
		ml_pdf=extraDistr::dpareto(y,a=ml_params,b=kscale)
		ml_cdf=extraDistr::ppareto(y,a=ml_params,b=kscale)

# rhp
		slnx=sum(log(x/kscale))
		bigx=slnx
		top=nx*(bigx**nx)
		bigz=log(y/kscale)+slnx
		bot=y*((bigz**(nx+1)))
		rh_pdf=top/bot
		rh_pdf=pmax(rh_pdf,0)

		rh_cdf=1-(bigx**nx)/(bigz**nx)
		rh_cdf=pmin(pmax(rh_cdf,0),1)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					rh_pdf=rh_pdf,
					ml_cdf=ml_cdf,
					rh_cdf=rh_cdf)
}

