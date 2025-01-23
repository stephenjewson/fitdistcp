#' Waic
#' @inheritParams manf
pareto_p1k2_waic=function(waicscores,x,t,v1hat,d1,v2hat,d2,kscale,lddi,lddd,lambdad){
		if(waicscores){
			f1f=pareto_p1k2_f1f(x,t,v1hat,d1,v2hat,d2,kscale)
			f2f=pareto_p1k2_f2f(x,t,v1hat,d1,v2hat,d2,kscale)
			fhatx=dpareto_p1k2(x,t,ymn=v1hat,slope=v2hat,kscale=kscale,log=FALSE)
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
pareto_p1k2_predictordata=function(predictordata,x,t,t0,params,kscale){
	if(predictordata){
#
# calculate the probabilities of the data using the fited model
#
		a=params[1]
		b=params[2]
		sh=1/exp(a+b*t)
		px=ppareto(x,a=sh,b=kscale)
#
# calculate the quantiles for those probabilities at t0
#
		sh0=1/exp(a+b*t0)
		qx=qpareto(px,a=sh0,b=kscale)
	} else{
		predictedparameter="predictordata not selected"
		adjustedx="predictordata not selected"
	}

	list(predictedparameter=sh,adjustedx=qx)
}
#' Logf for RUST
#' @inheritParams manf
pareto_p1k2_logf=function(params,x,t,kscale){
#	a=params[1]
#	b=params[2]
#	sh=1/exp(a+b*t)
#	logf=sum(dpareto(x,a=sh,b=kscale,log=TRUE))
	a=params[1]
	b=params[2]
	sh=pmax(1/exp(a+b*t),.Machine$double.eps)
	logf=sum(dpareto(x,a=sh,b=kscale,log=TRUE))
	return(logf)
}
#'  observed log-likelihood function
#' @inheritParams	manf
pareto_p1k2_loglik=function(vv,x,t,kscale){
	mu=vv[1]+vv[2]*t
	loglik=sum(dpareto(x,a=1/exp(mu),b=kscale,log=TRUE))
	return(loglik)
}
#' pareto_k1-with-p2 quantile function
#' @inheritParams	manf
qpareto_p1k2=function(p,t0,ymn,slope,kscale){

	mu=(ymn+slope*t0)
	return(qpareto(p,a=1/exp(mu),b=kscale))

}
#' pareto_k1-with-p2 density function
#' @inheritParams	manf
dpareto_p1k2=function(x,t0,ymn,slope,kscale,log=FALSE){

	mu=(ymn+slope*t0)
	return(dpareto(x,a=1/exp(mu),b=kscale,log=log))

}
#' pareto_k1-with-p2 distribution function
#' @inheritParams	manf
ppareto_p1k2=function(x,t0,ymn,slope,kscale){

	mu=(ymn+slope*t0)
	return(ppareto(x,a=1/exp(mu),b=kscale))

}
#' One component of the second derivative of the normalized log-likelihood
#' @inheritParams manf
pareto_p1k2_lmn=function(x,t,v1,d1,v2,d2,kscale,mm,nn){
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
			lmn[i]=sum(dpareto_p1k2(x,t,ymn=vvd[1],slope=vvd[2],kscale=kscale,log=TRUE))/nx
		}
		dld=(lmn[1]-lmn[2]-lmn[3]+lmn[4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:3){
			for (j in 1:2){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[i]=sum(dpareto_p1k2(x,t,ymn=vvd[1],slope=vvd[2],kscale=kscale,log=TRUE))/nx
		}
		dld=(lmn[1]-2*lmn[2]+lmn[3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' Second derivative matrix of the normalized log-likelihood
#' @inheritParams	manf
pareto_p1k2_ldd=function(x,t,v1,d1,v2,d2,kscale){
	nx=length(x)
	ldd=matrix(0,2,2)
	for (i in 1:2){
		for (j in i:2){
			ldd[i,j]=pareto_p1k2_lmn(x,t,v1,d1,v2,d2,kscale=kscale,i,j)
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
pareto_p1k2_lmnp=function(x,t,v1,d1,v2,d2,kscale,mm,nn,rr){
	d2=d2*v2
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
			lmn[i]=sum(dpareto_p1k2(x,t,ymn=vvd[1],slope=vvd[2],kscale=kscale,log=TRUE))/nx
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
			lmn[i]=sum(dpareto_p1k2(x,t,ymn=vvd[1],slope=vvd[2],kscale=kscale,log=TRUE))/nx
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
			lmn[i]=sum(dpareto_p1k2(x,t,ymn=vvd[1],slope=vvd[2],kscale=kscale,log=TRUE))/nx
		}
		dld1=(lmn[3]-2*lmn[2]+lmn[1])/(dd[m2]*dd[m2])
		dld2=(lmn[6]-2*lmn[5]+lmn[4])/(dd[m2]*dd[m2])
		dld=(dld2-dld1)/(2*dd[n2])
	}
	return(dld)
}
#' Third derivative tensor of the normalized log-likelihood
#' @inheritParams	manf
pareto_p1k2_lddd=function(x,t,v1,d1,v2,d2,kscale){
# calculate the unique values
	lddd=array(0,c(2,2,2))
	for (i in 1:2){
		for (j in i:2){
			for (k in j:2){
				lddd[i,j,k]=pareto_p1k2_lmnp(x,t,v1,d1,v2,d2,kscale=kscale,i,j,k)
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
pareto_p1k2_f1f=function(y,t0,v1,d1,v2,d2,kscale){
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v2 stuff
	v2m1=v2-1*d2
	v200=v2+0*d2
	v2p1=v2+1*d2
# v1 derivatives
	F1m1=dpareto_p1k2(y,t0,ymn=v1m1,slope=v200,kscale=kscale)
	F1p1=dpareto_p1k2(y,t0,ymn=v1p1,slope=v200,kscale=kscale)
# v2 derivatives
	F2m1=dpareto_p1k2(y,t0,ymn=v100,slope=v2m1,kscale=kscale)
	F2p1=dpareto_p1k2(y,t0,ymn=v100,slope=v2p1,kscale=kscale)
	f1=matrix(0,2,length(y))
	f1[1,]=(F1p1-F1m1)/(2*d1)
	f1[2,]=(F2p1-F2m1)/(2*d2)
	return(f1)
}
#' DMGS equation 2.1, p1 term
#' @inheritParams	manf
pareto_p1k2_p1f=function(y,t0,v1,d1,v2,d2,kscale){
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v2 stuff
	v2m1=v2-1*d2
	v200=v2+0*d2
	v2p1=v2+1*d2
# v1 derivatives
	F1m1=ppareto_p1k2(y,t0,ymn=v1m1,slope=v200,kscale=kscale)
	F1p1=ppareto_p1k2(y,t0,ymn=v1p1,slope=v200,kscale=kscale)
# v2 derivatives
	F2m1=ppareto_p1k2(y,t0,ymn=v100,slope=v2m1,kscale=kscale)
	F2p1=ppareto_p1k2(y,t0,ymn=v100,slope=v2p1,kscale=kscale)
	p1=matrix(0,2,length(y))
	p1[1,]=(F1p1-F1m1)/(2*d1)
	p1[2,]=(F2p1-F2m1)/(2*d2)
	return(p1)
}
#' DMGS equation 3.3, mu1 term
#' @inheritParams	manf
pareto_p1k2_mu1f=function(alpha,t0,v1,d1,v2,d2,kscale){
	q00=qpareto_p1k2((1-alpha),t0,ymn=v1,slope=v2,kscale=kscale)
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v2 stuff
	v2m1=v2-1*d2
	v200=v2+0*d2
	v2p1=v2+1*d2
# v1 derivatives
	F1m1=ppareto_p1k2(q00,t0,ymn=v1m1,slope=v200,kscale=kscale)
	F1p1=ppareto_p1k2(q00,t0,ymn=v1p1,slope=v200,kscale=kscale)
# v2 derivatives
	F2m1=ppareto_p1k2(q00,t0,ymn=v100,slope=v2m1,kscale=kscale)
	F2p1=ppareto_p1k2(q00,t0,ymn=v100,slope=v2p1,kscale=kscale)
	mu1=matrix(0,2,length(alpha))
	mu1[1,]=-(F1p1-F1m1)/(2*d1)
	mu1[2,]=-(F2p1-F2m1)/(2*d2)
	return(mu1)
}
#' DMGS equation 2.1, f2 term
#' @inheritParams	manf
pareto_p1k2_f2f=function(y,t0,v1,d1,v2,d2,kscale){
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
	F1m2=dpareto_p1k2(y,t0,ymn=v1m2,slope=v2,kscale=kscale)
	F1m1=dpareto_p1k2(y,t0,ymn=v1m1,slope=v2,kscale=kscale)
	F100=dpareto_p1k2(y,t0,ymn=v100,slope=v2,kscale=kscale)
	F1p1=dpareto_p1k2(y,t0,ymn=v1p1,slope=v2,kscale=kscale)
	F1p2=dpareto_p1k2(y,t0,ymn=v1p2,slope=v2,kscale=kscale)
	f2[1,1,]=(F1p1-2*F100+F1m1)/(d1*d1)
# v2 derivative
	F2m2=dpareto_p1k2(y,t0,ymn=v1,slope=v2m2,kscale=kscale)
	F2m1=dpareto_p1k2(y,t0,ymn=v1,slope=v2m1,kscale=kscale)
	F200=dpareto_p1k2(y,t0,ymn=v1,slope=v200,kscale=kscale)
	F2p1=dpareto_p1k2(y,t0,ymn=v1,slope=v2p1,kscale=kscale)
	F2p2=dpareto_p1k2(y,t0,ymn=v1,slope=v2p2,kscale=kscale)
	f2[2,2,]=(F2p1-2*F200+F2m1)/(d2*d2)
# cross derivative12
	Fcm1m1=dpareto_p1k2(y,t0,ymn=v1m1,slope=v2m1,kscale=kscale)
	Fcm1p1=dpareto_p1k2(y,t0,ymn=v1m1,slope=v2p1,kscale=kscale)
	Fcp1m1=dpareto_p1k2(y,t0,ymn=v1p1,slope=v2m1,kscale=kscale)
	Fcp1p1=dpareto_p1k2(y,t0,ymn=v1p1,slope=v2p1,kscale=kscale)
	f2[1,2,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2)
	f2[2,1,]=f2[1,2,]
	return(f2)
}
#' DMGS equation 2.1, p2 term
#' @inheritParams	manf
pareto_p1k2_p2f=function(y,t0,v1,d1,v2,d2,kscale){
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
	F1m2=ppareto_p1k2(y,t0,ymn=v1m2,slope=v2,kscale=kscale)
	F1m1=ppareto_p1k2(y,t0,ymn=v1m1,slope=v2,kscale=kscale)
	F100=ppareto_p1k2(y,t0,ymn=v100,slope=v2,kscale=kscale)
	F1p1=ppareto_p1k2(y,t0,ymn=v1p1,slope=v2,kscale=kscale)
	F1p2=ppareto_p1k2(y,t0,ymn=v1p2,slope=v2,kscale=kscale)
	p2[1,1,]=(F1p1-2*F100+F1m1)/(d1*d1)
# v2 derivative
	F2m2=ppareto_p1k2(y,t0,ymn=v1,slope=v2m2,kscale=kscale)
	F2m1=ppareto_p1k2(y,t0,ymn=v1,slope=v2m1,kscale=kscale)
	F200=ppareto_p1k2(y,t0,ymn=v1,slope=v200,kscale=kscale)
	F2p1=ppareto_p1k2(y,t0,ymn=v1,slope=v2p1,kscale=kscale)
	F2p2=ppareto_p1k2(y,t0,ymn=v1,slope=v2p2,kscale=kscale)
	p2[2,2,]=(F2p1-2*F200+F2m1)/(d2*d2)
# cross derivative12
	Fcm1m1=ppareto_p1k2(y,t0,ymn=v1m1,slope=v2m1,kscale=kscale)
	Fcm1p1=ppareto_p1k2(y,t0,ymn=v1m1,slope=v2p1,kscale=kscale)
	Fcp1m1=ppareto_p1k2(y,t0,ymn=v1p1,slope=v2m1,kscale=kscale)
	Fcp1p1=ppareto_p1k2(y,t0,ymn=v1p1,slope=v2p1,kscale=kscale)
	p2[1,2,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2)
	p2[2,1,]=p2[1,2,]
	return(p2)
}
#' DMGS equation 3.3, mu2 term
#' @inheritParams	manf
pareto_p1k2_mu2f=function(alpha,t0,v1,d1,v2,d2,kscale){
	q00=qpareto_p1k2((1-alpha),t0,ymn=v1,slope=v2,kscale=kscale)
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
	F1m2=ppareto_p1k2(q00,t0,ymn=v1m2,slope=v2,kscale=kscale)
	F1m1=ppareto_p1k2(q00,t0,ymn=v1m1,slope=v2,kscale=kscale)
	F100=ppareto_p1k2(q00,t0,ymn=v100,slope=v2,kscale=kscale)
	F1p1=ppareto_p1k2(q00,t0,ymn=v1p1,slope=v2,kscale=kscale)
	F1p2=ppareto_p1k2(q00,t0,ymn=v1p2,slope=v2,kscale=kscale)
	mu2[1,1,]=-(F1p1-2*F100+F1m1)/(d1*d1)
# v2 derivative
	F2m2=ppareto_p1k2(q00,t0,ymn=v1,slope=v2m2,kscale=kscale)
	F2m1=ppareto_p1k2(q00,t0,ymn=v1,slope=v2m1,kscale=kscale)
	F200=ppareto_p1k2(q00,t0,ymn=v1,slope=v200,kscale=kscale)
	F2p1=ppareto_p1k2(q00,t0,ymn=v1,slope=v2p1,kscale=kscale)
	F2p2=ppareto_p1k2(q00,t0,ymn=v1,slope=v2p2,kscale=kscale)
	mu2[2,2,]=-(F2p1-2*F200+F2m1)/(d2*d2)
# cross derivative12
	Fcm1m1=ppareto_p1k2(q00,t0,ymn=v1m1,slope=v2m1,kscale=kscale)
	Fcm1p1=ppareto_p1k2(q00,t0,ymn=v1m1,slope=v2p1,kscale=kscale)
	Fcp1m1=ppareto_p1k2(q00,t0,ymn=v1p1,slope=v2m1,kscale=kscale)
	Fcp1p1=ppareto_p1k2(q00,t0,ymn=v1p1,slope=v2p1,kscale=kscale)
	mu2[1,2,]=-(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2)
	mu2[2,1,]=mu2[1,2,]
	return(mu2)
}
#' pareto_k1 distribution: RHP mean
#' @inheritParams	manf
pareto_p1k2_means=function(means,t0,ml_params,lddi,lddd,lambdad_rhp,nx,dim=2,kscale){

	if(means){
# intro
		v1=ml_params[1]
		v2=ml_params[2]

# ml mean
		mu_hat=v1+v2*t0
		if(mu_hat>1){
			ml_mean=mu_hat*kscale/(mu_hat-1)
		} else {
			ml_mean=Inf
		}

# rhp mean
		rh_mean=Inf
	}else{
		ml_mean="means not selected"
		rh_mean="means not selected"
	}

	list(ml_mean=ml_mean,rh_mean=rh_mean)

}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inheritParams	manf
pareto_p1k2_logscores=function(logscores,x,t,d1,d2,kscale,aderivs,debug){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){

			x1=x[-i]
			t1=t[-i]

			dd=dpareto_p1k2sub(x1,t1,x[i],t[i],d1,d2,kscale=kscale,aderivs,debug=debug)

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
dpareto_p1k2sub=function(x,t,y,t0,d1,d2,kscale,aderivs=TRUE,debug=FALSE){


		if(debug)cat("inside pareto_p1k2sub\n")
		nx=length(x)

		lm=lm(x~t)
		v1start=lm$coefficients[1]
		v2start=lm$coefficients[2]
		v1start=0
		v2start=1
		xhat=v1start+v2start*t
		opt1=optim(c(v1start,v2start),pareto_p1k2_loglik,x=x,t=t,
			kscale=kscale,control=list(fnscale=-1))
		v1hat=opt1$par[1]
		v2hat=opt1$par[2]
		ml_params=c(v1hat,v2hat)

# ml
		muhat=v1hat+v2hat*t0
		ml_pdf=dpareto(y,a=1/exp(muhat),b=kscale)
		ml_cdf=ppareto(y,a=1/exp(muhat),b=kscale)

# rhp
		if(debug)cat("calc ldd\n")
		if(aderivs) ldd=pareto_p1k2_ldda(x,t,v1hat,v2hat,kscale)
		if(!aderivs)ldd=pareto_p1k2_ldd(x,t,v1hat,d1,v2hat,d2,kscale)
		lddi=solve(ldd)

		if(debug)cat("calc lddd\n")
		if(aderivs) lddd=pareto_p1k2_lddda(x,t,v1hat,v2hat,kscale)
		if(!aderivs)lddd=pareto_p1k2_lddd(x,t,v1hat,d1,v2hat,d2,kscale)

		if(debug)cat("calc f1f\n")
		if(aderivs) f1=pareto_p1k2_f1fa(y,t0,v1hat,v2hat,kscale)
		if(!aderivs)f1=pareto_p1k2_f1f(y,t0,v1hat,d1,v2hat,d2,kscale)

		if(debug)cat("calc f2f\n")
		if(aderivs) f2=pareto_p1k2_f2fa(y,t0,v1hat,v2hat,kscale)
		if(!aderivs)f2=pareto_p1k2_f2f(y,t0,v1hat,d1,v2hat,d2,kscale)

		if(debug)cat("calc p1f\n")
		if(aderivs) p1=pareto_p1k2_p1fa(y,t0,v1hat,v2hat,kscale)
		if(!aderivs)p1=pareto_p1k2_p1f(y,t0,v1hat,d1,v2hat,d2,kscale)

		if(debug)cat("calc p2f\n")
		if(aderivs) p2=pareto_p1k2_p2fa(y,t0,v1hat,v2hat,kscale)
		if(!aderivs)p2=pareto_p1k2_p2f(y,t0,v1hat,d1,v2hat,d2,kscale)

		if(debug)cat("call dmgs\n")
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
