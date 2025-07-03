#' Waicscores
#' @inherit manwaic return
#' @inheritParams manf
exp_waic=function(waicscores,x,v1hat,fd1,aderivs){
	if(waicscores){
		if(aderivs) f1f=exp_f1fa(x,v1hat)
		if(!aderivs)f1f=exp_f1f(x,v1hat,fd1)

		if(aderivs) f2f=exp_f2fa(x,v1hat)
		if(!aderivs)f2f=exp_f2f(x,v1hat,fd1)

		if(aderivs)	ldd=exp_ldda(x,v1hat)
		if(!aderivs)ldd=exp_ldd(x,v1hat,fd1)
		lddi=solve(ldd)
		if(aderivs)	lddd=exp_lddda(x,v1hat)
		if(!aderivs)lddd=exp_lddd(x,v1hat,fd1)
		fhatx=dexp(x,rate=v1hat)
		lambdad_rhp=-1/v1hat
		waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad_rhp,f2f,dim=1)
		waic1=waic$waic1
		waic2=waic$waic2
	}else{
		waic1="waicscores not selected"
		waic2="waicscores not selected"
	}
	return(list(waic1=waic1,waic2=waic2))
}
#' Logf for RUST
#' @inherit manlogf return
#' @inheritParams manf
exp_logf=function(params,x){
	l=pmax(params[1],.Machine$double.eps)
	logf=sum(dexp(x,rate=l,log=TRUE))-log(l)
	return(logf)
}
#' The second derivative of the normalized log-likelihood
#' @inherit manldd return
#' @inheritParams manf
exp_ldd=function(x,v1,fd1){
	nx=length(x)
	d1=fd1*v1
	v1m1=v1-1*d1
	v100=v1
	v1p1=v1+1*d1
	lm1=sum(log(dexp(x,rate=v1m1)))/nx
	l00=sum(log(dexp(x,rate=v100)))/nx
	lp1=sum(log(dexp(x,rate=v1p1)))/nx
	ldd=matrix(0,1,1)
	ldd[1,1]=(lp1-2*l00+lm1)/(d1*d1)
	return(ldd)
}
#' Third derivative of the normalized log-likelihood
#' @inherit manlnnn return
#' @inheritParams manf
exp_l111=function(x,v1,fd1){
	nx=length(x)
	d1=fd1*v1
	v1m2=v1-2*d1
	v1m1=v1-1*d1
	v1p1=v1+1*d1
	v1p2=v1+2*d1
	lm2=sum(dexp(x,rate=v1m2,log=TRUE))/nx
	lm1=sum(dexp(x,rate=v1m1,log=TRUE))/nx
	lp1=sum(dexp(x,rate=v1p1,log=TRUE))/nx
	lp2=sum(dexp(x,rate=v1p2,log=TRUE))/nx
	dld111=(lp2-2*lp1+2*lm1-lm2)/(2*d1*d1*d1)
	return(dld111)
}
#' Third derivative tensor of the log-likelihood
#' @inherit manlddd return
#' @inheritParams manf
exp_lddd=function(x,v1,fd1){
	nx=length(x)
	lddd=array(0,c(1,1,1))
	lddd[1,1,1]=exp_l111(x,v1,fd1)
	return(lddd)
}
#' DMGS equation 2.1, f1 term
#' @inherit man1f return
#' @inheritParams manf
exp_f1f=function(y,v1,fd1){
	d1=fd1*v1
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v1 derivatives
	F1m1=dexp(y,rate=v1m1)
	F1p1=dexp(y,rate=v1p1)
	f1=matrix(0,1,length(y))
	f1[1,]=(F1p1-F1m1)/(2*d1)
	return(f1)
}
#' DMGS equation 2.1, f2 term
#' @inherit man2f return
#' @inheritParams manf
exp_f2f=function(y,v1,fd1){
	d1=fd1*v1
# v1 stuff
	v1m2=v1-2*d1
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
	v1p2=v1+2*d1
	F1m2=dexp(y,rate=v1m2)
	F1m1=dexp(y,rate=v1m1)
	F100=dexp(y,rate=v100)
	F1p1=dexp(y,rate=v1p1)
	F1p2=dexp(y,rate=v1p2)
	f2=array(0,c(1,1,length(y)))
	f2[1,1,]=(F1p1-2*F100+F1m1)/(d1*d1)
	return(f2)
}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inherit manlogscores return
#' @inheritParams manf
exp_logscores=function(logscores,x){

	if(logscores){
# I could put the logs inside dexpsub, but I'd have toa actually calculate the log for the rhp case

		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){
			x1=x[-i]
			dd=dexpsub(x1,x[i])
			ml_params1=(nx-1)/sum(x1)
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
dexpsub=function(x,y,aderivs=TRUE){

		nx=length(x)

# ml
		ml_params=nx/sum(x)
		ml_pdf=dexp(y,rate=ml_params)
		ml_cdf=pexp(y,rate=ml_params)

# rhp pdf
		sx=sum(x)
		top=sx**nx
		bot1=(sx+y)**(nx+1)
		rh_pdf=nx*top/bot1
		rh_pdf=pmax(rh_pdf,0)

# rhp cdf
		bot2=(sx+y)**nx
		rh_cdf=1-top/bot2
		rh_cdf=pmin(pmax(rh_cdf,0),1)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					rh_pdf=rh_pdf,
					ml_cdf=ml_cdf,
					rh_cdf=rh_cdf)
}
