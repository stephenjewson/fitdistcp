#'  Waic
#' @inheritParams manf
norm_dmgs_waic=function(waicscores,x,v1hat,d1,v2hat,fd2,lddi,lddd,lambdad,aderivs){
		if(waicscores){

			if(aderivs) f1f=norm_f1fa(x,v1hat,v2hat)
			if(!aderivs)f1f=norm_f1f(x,v1hat,d1,v2hat,fd2)

			if(aderivs) f2f=norm_f2fa(x,v1hat,v2hat)
			if(!aderivs)f2f=norm_f2f(x,v1hat,d1,v2hat,fd2)

			fhatx=dnorm(x,mean=v1hat,sd=v2hat)
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
norm_dmgs_loglik=function(vv,x){
	loglik=sum(dnorm(x,mean=vv[1],sd=max(vv[2],.Machine$double.eps),log=TRUE))
	return(loglik)
}
#' One component of the second derivative of the normalized log-likelihood
#' @inheritParams manf
norm_lmn=function(x,v1,d1,v2,fd2,mm,nn){
	d2=fd2*v2
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
			lmn[i]=sum(dnorm(x,mean=vvd[1],sd=vvd[2],log=TRUE))/nx
		}
		dld=(lmn[1]-lmn[2]-lmn[3]+lmn[4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:3){
			for (j in 1:2){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[i]=sum(dnorm(x,mean=vvd[1],sd=vvd[2],log=TRUE))/nx
		}
		dld=(lmn[1]-2*lmn[2]+lmn[3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' One component of the second derivative of the normalized log-likelihood
#' @inheritParams manf
norm_lmnp=function(x,v1,d1,v2,fd2,mm,nn,rr){
	d2=fd2*v2
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
			lmn[i]=sum(dnorm(x,mean=vvd[1],sd=vvd[2],log=TRUE))/nx
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
			lmn[i]=sum(dnorm(x,mean=vvd[1],sd=vvd[2],log=TRUE))/nx
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
			lmn[i]=sum(dnorm(x,mean=vvd[1],sd=vvd[2],log=TRUE))/nx
		}
		dld1=(lmn[3]-2*lmn[2]+lmn[1])/(dd[m2]*dd[m2])
		dld2=(lmn[6]-2*lmn[5]+lmn[4])/(dd[m2]*dd[m2])
		dld=(dld2-dld1)/(2*dd[n2])
	}
	return(dld)
}
#' DMGS equation 3.3, p1 term
#' @inheritParams manf
norm_dmgs_p1f=function(y,v1,d1,v2,fd2){
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
	F1m1=pnorm(y,mean=v1m1,sd=v200)
	F1p1=pnorm(y,mean=v1p1,sd=v200)
# v2 derivatives
	F2m1=pnorm(y,mean=v100,sd=v2m1)
	F2p1=pnorm(y,mean=v100,sd=v2p1)
	p1=matrix(0,2,length(y))
	p1[1,]=(F1p1-F1m1)/(2*d1)
	p1[2,]=(F2p1-F2m1)/(2*d2)
	return(p1)
}
#' DMGS equation 3.3, mu1 term
#' @inheritParams manf
norm_dmgs_mu1f=function(alpha,v1,d1,v2,fd2){
	q00=qnorm((1-alpha),mean=v1,sd=v2)
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
	F1m1=pnorm(q00,mean=v1m1,sd=v200)
	F1p1=pnorm(q00,mean=v1p1,sd=v200)
# v2 derivatives
	F2m1=pnorm(q00,mean=v100,sd=v2m1)
	F2p1=pnorm(q00,mean=v100,sd=v2p1)
	mu1=matrix(0,2,length(alpha))
	mu1[1,]=-(F1p1-F1m1)/(2*d1)
	mu1[2,]=-(F2p1-F2m1)/(2*d2)
	return(mu1)
}
#' DMGS equation 3.3, p2 term
#' @inheritParams manf
norm_dmgs_p2f=function(y,v1,d1,v2,fd2){
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
	F1m2=pnorm(y,mean=v1m2,sd=v200)
	F1m1=pnorm(y,mean=v1m1,sd=v200)
	F100=pnorm(y,mean=v100,sd=v200)
	F1p1=pnorm(y,mean=v1p1,sd=v200)
	F1p2=pnorm(y,mean=v1p2,sd=v200)
# v2 derivative
	F2m2=pnorm(y,mean=v100,sd=v2m2)
	F2m1=pnorm(y,mean=v100,sd=v2m1)
	F200=pnorm(y,mean=v100,sd=v200)
	F2p1=pnorm(y,mean=v100,sd=v2p1)
	F2p2=pnorm(y,mean=v100,sd=v2p2)
# cross derivative
	Fcm1m1=pnorm(y,mean=v1m1,sd=v2m1)
	Fcm1p1=pnorm(y,mean=v1m1,sd=v2p1)
	Fcp1m1=pnorm(y,mean=v1p1,sd=v2m1)
	Fcp1p1=pnorm(y,mean=v1p1,sd=v2p1)
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
norm_dmgs_mu2f=function(alpha,v1,d1,v2,fd2){
	q00=qnorm((1-alpha),mean=v1,sd=v2)
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
	F1m2=pnorm(q00,mean=v1m2,sd=v200)
	F1m1=pnorm(q00,mean=v1m1,sd=v200)
	F100=pnorm(q00,mean=v100,sd=v200)
	F1p1=pnorm(q00,mean=v1p1,sd=v200)
	F1p2=pnorm(q00,mean=v1p2,sd=v200)
# v2 derivative
	F2m2=pnorm(q00,mean=v100,sd=v2m2)
	F2m1=pnorm(q00,mean=v100,sd=v2m1)
	F200=pnorm(q00,mean=v100,sd=v200)
	F2p1=pnorm(q00,mean=v100,sd=v2p1)
	F2p2=pnorm(q00,mean=v100,sd=v2p2)
# cross derivative
	Fcm1m1=pnorm(q00,mean=v1m1,sd=v2m1)
	Fcm1p1=pnorm(q00,mean=v1m1,sd=v2p1)
	Fcp1m1=pnorm(q00,mean=v1p1,sd=v2m1)
	Fcp1p1=pnorm(q00,mean=v1p1,sd=v2p1)
	mu2[1,1,]=-(F1p1-2*F100+F1m1)/(d1*d1)
	mu2[2,2,]=-(F2p1-2*F200+F2m1)/(d2*d2)
	mu2[1,2,]=-(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2)
	# copy
	mu2[2,1,]=mu2[1,2,]
	return(mu2)
}
#' MLE and RHP predictive means
#' @inheritParams manf
norm_dmgs_means=function(means,ml_params,lddi,lddd,lambdad_rhp,nx,dim=2){

	if(means){
# intro

# ml mean
		ml_mean=ml_params[1]

# rhp mean
		meand1=array(0,c(2,1))
		meand1[1,1]=1
		meand1[2,1]=0
		meand2=array(0,c(2,2,1)) #but all zero for norm_dmgs
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
norm_dmgs_logscores=function(logscores,x,d1=0.01,fd2=0.01){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){

			x1=x[-i]

			dd=dnorm_dmgssub(x1,x[i],d1,fd2)

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
dnorm_dmgssub=function(x,y,d1=0.01,fd2=0.01,aderivs=TRUE){

		nx=length(x)

		v1start=mean(x)
		v2start=sd(x)
		opt=optim(c(v1start,v2start),norm_dmgs_loglik,x=x,control=list(fnscale=-1))
		v1hat=opt$par[1]
		v2hat=opt$par[2]
		ml_params=c(v1hat,v2hat)

# mle
		ml_pdf=dnorm(y,mean=v1hat,sd=v2hat)
		ml_cdf=pnorm(y,mean=v1hat,sd=v2hat)

# rhp
		if(aderivs)	ldd=norm_ldda(x,v1hat,v2hat)
		if(!aderivs)ldd=norm_ldd(x,v1hat,d1,v2hat,fd2)
		lddi=solve(ldd)
		if(aderivs)	lddd=norm_lddda(x,v1hat,v2hat)
		if(!aderivs)lddd=norm_lddd(x,v1hat,d1,v2hat,fd2)

		if(aderivs) f1=norm_f1fa(y,v1hat,v2hat)
		if(!aderivs)f1=norm_f1f(y,v1hat,d1,v2hat,fd2)

		if(aderivs) f2=norm_f2fa(y,v1hat,v2hat)
		if(!aderivs)f2=norm_f2f(y,v1hat,d1,v2hat,fd2)

		if(aderivs) p1=norm_p1fa(y,v1hat,v2hat)
		if(!aderivs)p1=norm_dmgs_p1f(y,v1hat,d1,v2hat,fd2)

		if(aderivs) p2=norm_p2fa(y,v1hat,v2hat)
		if(!aderivs)p2=norm_dmgs_p2f(y,v1hat,d1,v2hat,fd2)

		lambdad_rhp=c(0,-1/v2hat)
		df1=dmgs(lddi,lddd,f1,lambdad_rhp,f2,dim=2)
		dp1=dmgs(lddi,lddd,p1,lambdad_rhp,p2,dim=2)
		rh_pdf=pmax(ml_pdf+df1/nx,0)
		rh_cdf=pmin(pmax(ml_cdf+dp1/nx,0),1)

# slipping the rhp: mu1 term
		lambdad_rhp=c(0,-1/v2hat)
		df1=dmgs(lddi,lddd,f1,lambdad_rhp,f2*0,dim=2)
		dp1=dmgs(lddi,lddd,p1,lambdad_rhp,p2*0,dim=2)
		rh_pdf1=df1/nx
		rh_cdf1=dp1/nx

# slipping the rhp: mu1 term
		lambdad_rhp=c(0,-1/v2hat)
		df1=dmgs(lddi,lddd,f1*0,lambdad_rhp,f2,dim=2)
		dp1=dmgs(lddi,lddd,p1*0,lambdad_rhp,p2,dim=2)
		rh_pdf2=df1/nx
		rh_cdf2=dp1/nx

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					rh_pdf=rh_pdf,
					ml_cdf=ml_cdf,
					rh_cdf=rh_cdf,
					rh_pdf1=rh_pdf1,
					rh_cdf1=rh_cdf1,
					rh_pdf2=rh_pdf2,
					rh_cdf2=rh_cdf2)

}

