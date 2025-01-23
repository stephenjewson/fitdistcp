#' Waic
#' @inheritParams manf
exp_p1_waic=function(waicscores,x,t,v1hat,d1,v2hat,d2,lddi,lddd,lambdad,aderivs){
		if(waicscores){

			if(aderivs) f1f=exp_p1_f1fa(x,t,v1hat,v2hat)
			if(!aderivs)f1f=exp_p1_f1f(x,t,v1hat,d1,v2hat,d2)

			if(aderivs) f2f=exp_p1_f2fa(x,t,v1hat,v2hat)
			if(!aderivs)f2f=exp_p1_f2f(x,t,v1hat,d1,v2hat,d2)

			fhatx=dexp_p1(x,t,ymn=v1hat,slope=v2hat,log=FALSE)
			waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad,f2f,dim=2)
			waic1=waic$waic1
			waic2=waic$waic2
		}else{
			waic1="extras not selected"
			waic2="extras not selected"
		}
		list(waic1=waic1,waic2=waic2)
}
#' Predicted Parameter and Generalized Residuals
#' @inheritParams manf
exp_p1_predictordata=function(predictordata,x,t,t0,params){
	if(predictordata){
#
# calculate the probabilities of the data using the fited model
#
		a=params[1]
		b=params[2]
		rr=1/exp(a+b*t)
		px=pexp(x,rate=rr)
#
# calculate the quantiles for those probabilities at t0
#
		rr0=1/exp(a+b*t0)
		qx=qexp(px,rate=rr0)
	} else{
		predictedparameter="predictordata not selected"
		adjustedx="predictordata not selected"
	}

	list(predictedparameter=rr,adjustedx=qx)
}
#' Logf for RUST
#' @inheritParams manf
exp_p1_logf=function(params,x,t){
#	a=params[1]
#	b=params[2]
#	r=1/exp(a+b*t)
#	logf=sum(dexp(x,rate=r,log=TRUE))
	a=params[1]
	b=params[2]
	r=1/exp(a+b*t) #rate can be zero, that's ok...it's 1/sigma
	logf=sum(dexp(x,rate=r,log=TRUE))
	return(logf)
}
#'  observed log-likelihood function
#' @inheritParams	manf
exp_p1_loglik=function(vv,x,t){
	mu=vv[1]+vv[2]*t
# the 1/ parametrisation is beter because there is an expression for the mean
# -really? but the only difference is  signs
# in the ML model
# so this is the version in which sigma is loglinear, not rate
# for standardisation, beter to have sigma as loglinear, not rate. That's the best argument.
	loglik=sum(dexp(x,rate=1/exp(mu),log=TRUE))
#	loglik=-sum(dexp(x,exp(mu),log=TRUE))
	return(loglik)
}
#' -with-p1 quantile function
#' @inheritParams	manf
qexp_p1=function(p,t0,ymn,slope){

	mu=(ymn+slope*t0)
	return(qexp(p,rate=(1/exp(mu))))
#	return(qexp(p,rate=exp(mu)))

}
#' Exponential-with-p1 density function
#' @inheritParams	manf
dexp_p1=function(x,t0,ymn,slope,log=FALSE){

	mu=(ymn+slope*t0)
	return(dexp(x,rate=(1/exp(mu)),log=log))
#	return(dexp(x,rate=exp(mu),log=log))

}
#' Exponential-with-p1 distribution function
#' @inheritParams	manf
pexp_p1=function(x,t0,ymn,slope){

	mu=(ymn+slope*t0)
	return(pexp(x,rate=(1/exp(mu))))
#	return(pexp(x,rate=exp(mu)))

}
#' One component of the second derivative of the normalized log-likelihood
#' @inheritParams manf
exp_p1_lmn=function(x,t,v1,d1,v2,d2,mm,nn){
	net3=matrix(0,3,2)
	net4=matrix(0,4,2)
	lmn=matrix(0,4)
	dd=c(d1,d2)
	vv=c(v1,v2)
	vvd=matrix(0,2)
	nx=length(x)
# different
	if(mm!=nn){
		net4[,mm]=c(-1,-1,1,1)
		net4[,nn]=c(-1,1,-1,1)
		for (i in 1:4){
			for (j in 1:2){
				vvd[j]=vv[j]+net4[i,j]*dd[j]
			}
			lmn[i]=sum(dexp_p1(x,t,ymn=vvd[1],slope=vvd[2],log=TRUE))/nx
		}
		dld=(lmn[1]-lmn[2]-lmn[3]+lmn[4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:3){
			for (j in 1:2){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[i]=sum(dexp_p1(x,t,ymn=vvd[1],slope=vvd[2],log=TRUE))/nx
		}
		dld=(lmn[1]-2*lmn[2]+lmn[3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' Second derivative matrix of the normalized log-likelihood
#' @inheritParams	manf
exp_p1_ldd=function(x,t,v1,d1,v2,d2){
	ldd=matrix(0,2,2)
	for (i in 1:2){
		for (j in i:2){
			ldd[i,j]=exp_p1_lmn(x,t,v1,d1,v2,d2,i,j)
		}
	}
	for (i in 2:1){
		for (j in 1:(i-1)){
			ldd[i,j]=ldd[j,i]
		}
	}
	return(ldd)
}
#' One component of the second derivative of the normalized log-likelihood
#' @inheritParams manf
exp_p1_lmnp=function(x,t,v1,d1,v2,d2,mm,nn,rr){
	net4=matrix(0,4,2)
	net6=matrix(0,6,2)
	net8=matrix(0,8,2)
	lmn=matrix(0,8)
	dd=c(d1,d2)
	vv=c(v1,v2)
	vvd=matrix(0,2)
	nx=length(x)
# all diff
	if ((mm!=nn)&(nn!=rr)&(rr!=mm)){
		net8[,mm]=c(-1,1,-1,1,-1,1,-1,1)
		net8[,nn]=c(-1,-1,1,1,-1,-1,1,1)
		net8[,rr]=c(-1,-1,-1,-1,1,1,1,1)
		for (i in 1:8){
			for (j in 1:3){
				vvd[j]=vv[j]+net8[i,j]*dd[j]
			}
			lmn[i]=sum(dexp_p1(x,t,ymn=vvd[1],slope=vvd[2],log=TRUE))/nx
		}
		dld1=(lmn[2]-lmn[1])/(2*dd[mm])
		dld2=(lmn[4]-lmn[3])/(2*dd[mm])
		dld21=(dld2-dld1)/(2*dd[nn])
		dld3=(lmn[6]-lmn[5])/(2*dd[mm])
		dld4=(lmn[8]-lmn[7])/(2*dd[mm])
		dld43=(dld4-dld3)/(2*dd[nn])
		dld=(dld43-dld21)/(2*dd[rr])
# all 3 the same
	} else if ((mm==nn)&(nn==rr)){
		net4[,mm]=c(-2,-1,1,2)
		for (i in 1:4){
			for (j in 1:2){
				vvd[j]=vv[j]+net4[i,j]*dd[j]
			}
			lmn[i]=sum(dexp_p1(x,t,ymn=vvd[1],slope=vvd[2],log=TRUE))/nx
		}
		dld=(-lmn[1]+2*lmn[2]-2*lmn[3]+lmn[4])/(2*dd[mm]*dd[mm]*dd[mm])
	} else {
# 2 the same
# mm is the repeated one, nn is the other one
		if(mm==nn){m2=mm;n2=rr}
		if(mm==rr){m2=mm;n2=nn}
		if(nn==rr){m2=nn;n2=mm}
		net6[,m2]=c(-1,0,1,-1,0,1)
		net6[,n2]=c(-1,-1,-1,1,1,1)
		for (i in 1:6){
			for (j in 1:2){
				vvd[j]=vv[j]+net6[i,j]*dd[j]
			}
			lmn[i]=sum(dexp_p1(x,t,ymn=vvd[1],slope=vvd[2],log=TRUE))/nx
		}
		dld1=(lmn[3]-2*lmn[2]+lmn[1])/(dd[m2]*dd[m2])
		dld2=(lmn[6]-2*lmn[5]+lmn[4])/(dd[m2]*dd[m2])
		dld=(dld2-dld1)/(2*dd[n2])
	}
	return(dld)
}
#' Third derivative tensor of the normalized log-likelihood
#' @inheritParams	manf
exp_p1_lddd=function(x,t,v1,d1,v2,d2){
# calculate the unique values
	lddd=array(0,c(2,2,2))
	for (i in 1:2){
		for (j in i:2){
			for (k in j:2){
				lddd[i,j,k]=exp_p1_lmnp(x,t,v1,d1,v2,d2,i,j,k)
			}
		}
	}
# steves dumb algorithm for filling in the non-unique values
	for (i in 1:2){
		for (j in 1:2){
			for (k in 1:2){
				a=c(i,j,k)
				b=sort(a)
				lddd[a[1],a[2],a[3]]=lddd[b[1],b[2],b[3]]
			}
		}
	}
	return(lddd)
}
#' DMGS equation 2.1, f1 term
#' @inheritParams	manf
exp_p1_f1f=function(y,t0,v1,d1,v2,d2){
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v2 stuff
	v2m1=v2-1*d2
	v200=v2+0*d2
	v2p1=v2+1*d2
# v1 derivatives
	F1m1=dexp_p1(y,t0,ymn=v1m1,slope=v200)
	F1p1=dexp_p1(y,t0,ymn=v1p1,slope=v200)
# v2 derivatives
	F2m1=dexp_p1(y,t0,ymn=v100,slope=v2m1)
	F2p1=dexp_p1(y,t0,ymn=v100,slope=v2p1)
	f1=matrix(0,2,length(y))
	f1[1,]=(F1p1-F1m1)/(2*d1)
	f1[2,]=(F2p1-F2m1)/(2*d2)
	return(f1)
}
#' DMGS equation 2.1, p1 term
#' @inheritParams	manf
exp_p1_p1f=function(y,t0,v1,d1,v2,d2){
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v2 stuff
	v2m1=v2-1*d2
	v200=v2+0*d2
	v2p1=v2+1*d2
# v1 derivatives
	F1m1=pexp_p1(y,t0,ymn=v1m1,slope=v200)
	F1p1=pexp_p1(y,t0,ymn=v1p1,slope=v200)
# v2 derivatives
	F2m1=pexp_p1(y,t0,ymn=v100,slope=v2m1)
	F2p1=pexp_p1(y,t0,ymn=v100,slope=v2p1)
	p1=matrix(0,2,length(y))
	p1[1,]=(F1p1-F1m1)/(2*d1)
	p1[2,]=(F2p1-F2m1)/(2*d2)
	return(p1)
}
#' DMGS equation 3.3, mu1 term
#' @inheritParams	manf
exp_p1_mu1f=function(alpha,t0,v1,d1,v2,d2){
	q00=qexp_p1((1-alpha),t0,ymn=v1,slope=v2)
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v2 stuff
	v2m1=v2-1*d2
	v200=v2+0*d2
	v2p1=v2+1*d2
# v1 derivatives
	F1m1=pexp_p1(q00,t0,ymn=v1m1,slope=v200)
	F1p1=pexp_p1(q00,t0,ymn=v1p1,slope=v200)
# v2 derivatives
	F2m1=pexp_p1(q00,t0,ymn=v100,slope=v2m1)
	F2p1=pexp_p1(q00,t0,ymn=v100,slope=v2p1)
	mu1=matrix(0,2,length(alpha))
	mu1[1,]=-(F1p1-F1m1)/(2*d1)
	mu1[2,]=-(F2p1-F2m1)/(2*d2)
	return(mu1)
}
#' DMGS equation 2.1, f2 term
#' @inheritParams	manf
exp_p1_f2f=function(y,t0,v1,d1,v2,d2){
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
	f2=array(0,c(2,2,length(y)))
# v1
	F1m2=dexp_p1(y,t0,ymn=v1m2,slope=v2)
	F1m1=dexp_p1(y,t0,ymn=v1m1,slope=v2)
	F100=dexp_p1(y,t0,ymn=v100,slope=v2)
	F1p1=dexp_p1(y,t0,ymn=v1p1,slope=v2)
	F1p2=dexp_p1(y,t0,ymn=v1p2,slope=v2)
	f2[1,1,]=(F1p1-2*F100+F1m1)/(d1*d1)
# v2 derivative
	F2m2=dexp_p1(y,t0,ymn=v1,slope=v2m2)
	F2m1=dexp_p1(y,t0,ymn=v1,slope=v2m1)
	F200=dexp_p1(y,t0,ymn=v1,slope=v200)
	F2p1=dexp_p1(y,t0,ymn=v1,slope=v2p1)
	F2p2=dexp_p1(y,t0,ymn=v1,slope=v2p2)
	f2[2,2,]=(F2p1-2*F200+F2m1)/(d2*d2)
# cross derivative12
	Fcm1m1=dexp_p1(y,t0,ymn=v1m1,slope=v2m1)
	Fcm1p1=dexp_p1(y,t0,ymn=v1m1,slope=v2p1)
	Fcp1m1=dexp_p1(y,t0,ymn=v1p1,slope=v2m1)
	Fcp1p1=dexp_p1(y,t0,ymn=v1p1,slope=v2p1)
	f2[1,2,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2)
	f2[2,1,]=f2[1,2,]
	return(f2)
}
#' DMGS equation 2.1, p2 term
#' @inheritParams	manf
exp_p1_p2f=function(y,t0,v1,d1,v2,d2){
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
	p2=array(0,c(2,2,length(y)))
# v1
	F1m2=pexp_p1(y,t0,ymn=v1m2,slope=v2)
	F1m1=pexp_p1(y,t0,ymn=v1m1,slope=v2)
	F100=pexp_p1(y,t0,ymn=v100,slope=v2)
	F1p1=pexp_p1(y,t0,ymn=v1p1,slope=v2)
	F1p2=pexp_p1(y,t0,ymn=v1p2,slope=v2)
	p2[1,1,]=(F1p1-2*F100+F1m1)/(d1*d1)
# v2 derivative
	F2m2=pexp_p1(y,t0,ymn=v1,slope=v2m2)
	F2m1=pexp_p1(y,t0,ymn=v1,slope=v2m1)
	F200=pexp_p1(y,t0,ymn=v1,slope=v200)
	F2p1=pexp_p1(y,t0,ymn=v1,slope=v2p1)
	F2p2=pexp_p1(y,t0,ymn=v1,slope=v2p2)
	p2[2,2,]=(F2p1-2*F200+F2m1)/(d2*d2)
# cross derivative12
	Fcm1m1=pexp_p1(y,t0,ymn=v1m1,slope=v2m1)
	Fcm1p1=pexp_p1(y,t0,ymn=v1m1,slope=v2p1)
	Fcp1m1=pexp_p1(y,t0,ymn=v1p1,slope=v2m1)
	Fcp1p1=pexp_p1(y,t0,ymn=v1p1,slope=v2p1)
	p2[1,2,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2)
	p2[2,1,]=p2[1,2,]
	return(p2)
}
#' DMGS equation 3.3, mu2 term
#' @inheritParams	manf
exp_p1_mu2f=function(alpha,t0,v1,d1,v2,d2){
	q00=qexp_p1((1-alpha),t0,ymn=v1,slope=v2)
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
# v1
	F1m2=pexp_p1(q00,t0,ymn=v1m2,slope=v2)
	F1m1=pexp_p1(q00,t0,ymn=v1m1,slope=v2)
	F100=pexp_p1(q00,t0,ymn=v100,slope=v2)
	F1p1=pexp_p1(q00,t0,ymn=v1p1,slope=v2)
	F1p2=pexp_p1(q00,t0,ymn=v1p2,slope=v2)
	mu2[1,1,]=-(F1p1-2*F100+F1m1)/(d1*d1)
# v2 derivative
	F2m2=pexp_p1(q00,t0,ymn=v1,slope=v2m2)
	F2m1=pexp_p1(q00,t0,ymn=v1,slope=v2m1)
	F200=pexp_p1(q00,t0,ymn=v1,slope=v200)
	F2p1=pexp_p1(q00,t0,ymn=v1,slope=v2p1)
	F2p2=pexp_p1(q00,t0,ymn=v1,slope=v2p2)
	mu2[2,2,]=-(F2p1-2*F200+F2m1)/(d2*d2)
# cross derivative12
	Fcm1m1=pexp_p1(q00,t0,ymn=v1m1,slope=v2m1)
	Fcm1p1=pexp_p1(q00,t0,ymn=v1m1,slope=v2p1)
	Fcp1m1=pexp_p1(q00,t0,ymn=v1p1,slope=v2m1)
	Fcp1p1=pexp_p1(q00,t0,ymn=v1p1,slope=v2p1)
	mu2[1,2,]=-(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2)
	mu2[2,1,]=mu2[1,2,]
	return(mu2)
}
#' exp distribution: RHP mean
#' @inheritParams	manf
exp_p1_means=function(means,t0,ml_params,lddi,lddd,lambdad_rhp,nx,dim=2){

	if(means){
# intro
		v1=ml_params[1]
		v2=ml_params[2]

# ml mean
		mu_hat=v1+v2*t0
		ml_mean=exp(mu_hat)

# rhp mean
		meand1=array(0,c(2,1))
		meand1[1,1]=ml_mean
		meand1[2,1]=t0*ml_mean
		meand2=array(0,c(2,2,1))
		meand2[1,1,1]=ml_mean
		meand2[1,2,1]=t0*ml_mean
		meand2[2,1,1]=t0*ml_mean
		meand2[2,2,1]=t0*t0*ml_mean
		dmean=dmgs(lddi,lddd,meand1,lambdad_rhp,meand2,dim=2)
		rh_mean=ml_mean+dmean/nx
	}else{
		ml_mean="means not selected"
		rh_mean="means not selected"
	}

	list(ml_mean=ml_mean,rh_mean=rh_mean)

}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inheritParams	manf
exp_p1_logscores=function(logscores,x,t,d1,d2,aderivs){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){

			x1=x[-i]
			t1=t[-i]

			dd=dexp_p1sub(x1,t1,x[i],t[i],d1,d2,aderivs)

			ml_pdf=dd$ml_pdf
			ml_oos_logscore=ml_oos_logscore+log(ml_pdf)

			rh_pdf=dd$rh_pdf
			rh_oos_logscore=rh_oos_logscore+log(rh_pdf)
#			cat("i,ml_pdf,rh_pdf=",i,ml_pdf,rh_pdf,"\n")
		}
	}else{
		ml_oos_logscore="extras not selected"
		rh_oos_logscore="extras not selected"
	}
	list(ml_oos_logscore=ml_oos_logscore,rh_oos_logscore=rh_oos_logscore)
}
#' Densities from MLE and RHP
#' @inheritParams	manf
dexp_p1sub=function(x,t,y,t0,d1,d2,aderivs=TRUE){

		nx=length(x)

		lm=lm(x~t)
		v1start=lm$coefficients[1]
		v2start=lm$coefficients[2]
		v1start=0
		v2start=0
#		xhat=v1start+v2start*t
		opt1=optim(c(v1start,v2start),exp_p1_loglik,x=x,t=t,control=list(fnscale=-1))
		v1hat=opt1$par[1]
		v2hat=opt1$par[2]
		ml_params=c(v1hat,v2hat)

# ml
		muhat=v1hat+v2hat*t0
		ml_pdf=dexp(y,rate=1/exp(muhat))
		ml_cdf=pexp(y,rate=1/exp(muhat))
#		ml_pdf=dexp(y,rate=exp(muhat))
#		ml_cdf=pexp(y,rate=exp(muhat))

# rhp
		if(aderivs)	ldd=exp_p1_ldda(x,t,v1hat,v2hat)
		if(!aderivs)ldd=exp_p1_ldd(x,t,v1hat,d1,v2hat,d2)
		lddi=solve(ldd)

		if(aderivs) lddd=exp_p1_lddda(x,t,v1hat,v2hat)
		if(!aderivs)lddd=exp_p1_lddd(x,t,v1hat,d1,v2hat,d2)

		if(aderivs) f1=exp_p1_f1fa(y,t0,v1hat,v2hat)
		if(!aderivs)f1=exp_p1_f1f(y,t0,v1hat,d1,v2hat,d2)

		if(aderivs) f2=exp_p1_f2fa(y,t0,v1hat,v2hat)
		if(!aderivs)f2=exp_p1_f2f(y,t0,v1hat,d1,v2hat,d2)

		if(aderivs) p1=exp_p1_p1fa(y,t0,v1hat,v2hat)
		if(!aderivs)p1=exp_p1_p1f(y,t0,v1hat,d1,v2hat,d2)

		if(aderivs) p2=exp_p1_p2fa(y,t0,v1hat,v2hat)
		if(!aderivs)p2=exp_p1_p2f(y,t0,v1hat,d1,v2hat,d2)

		lambdad_rhp=c(0,0)
		df=dmgs(lddi,lddd,f1,lambdad_rhp,f2,dim=2)
		dp=dmgs(lddi,lddd,p1,lambdad_rhp,p2,dim=2)
		rh_pdf=pmax(ml_pdf+df/nx,0)
		rh_cdf=pmin(pmax(ml_cdf+dp/nx,0),1)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					rh_pdf=rh_pdf,
					ml_cdf=ml_cdf,
					rh_cdf=rh_cdf)
}
