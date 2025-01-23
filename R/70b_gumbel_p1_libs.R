#' Waic
#' @inheritParams manf
gumbel_p1_waic=function(waicscores,x,t,v1hat,d1,v2hat,d2,v3hat,fd3,
	lddi,lddd,lambdad,aderivs){
		if(waicscores){
			if(aderivs){
				f1f=gumbel_p1_f1fa(x,t,v1hat,v2hat,v3hat)
				f2f=gumbel_p1_f2fa(x,t,v1hat,v2hat,v3hat)
			} else {
				f1f=gumbel_p1_f1f(x,t,v1hat,d1,v2hat,d2,v3hat,fd3)
				f2f=gumbel_p1_f2f(x,t,v1hat,d1,v2hat,d2,v3hat,fd3)
			}
			fhatx=dgumbel_p1(x,t,ymn=v1hat,slope=v2hat,sigma=v3hat,log=FALSE)
			waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad,f2f,dim=3)
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
gumbel_p1_predictordata=function(predictordata,x,t,t0,params){
	if(predictordata){
#
# calculate the probabilities of the data using the fited model
#
		a=params[1]
		b=params[2]
		s=params[3]
		mu=a+b*t
		px=pgumbel(x,mu=mu,sigma=s)
#
# calculate the quantiles for those probabilities at t0
#
		mu0=a+b*t0
		qx=qgumbel(px,mu=mu0,sigma=s)
	} else{
		predictedparameter="predictordata not selected"
		adjustedx="predictordata not selected"
	}

	list(predictedparameter=mu,adjustedx=qx)
}
#' Logf for RUST
#' @inheritParams manf
gumbel_p1_logf=function(params,x,t){
#	a=params[1]
#	b=params[2]
#	s=params[3]
#	mu=a+b*t
#	if(s>0){
#		logf=sum(dgumbel(x,mu=mu,sigma=s,log=TRUE))-log(s)
#	}else{
#		logf=-Inf
#	}
	a=params[1]
	b=params[2]
	s=pmax(params[3],.Machine$double.eps)
	mu=a+b*t
	logf=sum(dgumbel(x,mu=mu,sigma=s,log=TRUE))-log(s)
	return(logf)
}
#'  observed log-likelihood function
#' @inheritParams	manf
gumbel_p1_loglik=function(vv,x,t){
	mu=vv[1]+vv[2]*t
	loglik=sum(dgumbel(x,mu=mu,sigma=max(vv[3],.Machine$double.eps),log=TRUE))
	return(loglik)
}
#' Gumbel-with-p1 quantile function
#' @inheritParams	manf
qgumbel_p1=function(p,t0,ymn,slope,sigma){

	return(qgumbel(p,mu=(ymn+slope*t0),sigma=sigma))

}
#' Gumbel-with-p1 density function
#' @inheritParams	manf
dgumbel_p1=function(x,t0,ymn,slope,sigma,log=FALSE){

	return(dgumbel(x,mu=(ymn+slope*t0),sigma=sigma,log=log))

}
#' Gumbel-with-p1 distribution function
#' @inheritParams	manf
pgumbel_p1=function(x,t0,ymn,slope,sigma){

	return(pgumbel(x,mu=(ymn+slope*t0),sigma=sigma))

}
#' One component of the second derivative of the normalized log-likelihood
#' @inheritParams manf
gumbel_p1_lmn=function(x,t,v1,d1,v2,d2,v3,fd3,mm,nn){
	d3=fd3*v3
	net3=matrix(0,3,3)
	net4=matrix(0,4,3)
	lmn=matrix(0,4)
	dd=c(d1,d2,d3)
	vv=c(v1,v2,v3)
	vvd=matrix(0,3)
	nx=length(x)
# different
	if(mm!=nn){
		net4[,mm]=c(-1,-1,1,1)
		net4[,nn]=c(-1,1,-1,1)
		for (i in 1:4){
			for (j in 1:3){
				vvd[j]=vv[j]+net4[i,j]*dd[j]
			}
			lmn[i]=sum(dgumbel_p1(x,t,ymn=vvd[1],slope=vvd[2],sigma=vvd[3],log=TRUE))/nx
		}
		dld=(lmn[1]-lmn[2]-lmn[3]+lmn[4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:3){
			for (j in 1:3){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[i]=sum(dgumbel_p1(x,t,ymn=vvd[1],slope=vvd[2],sigma=vvd[3],log=TRUE))/nx
		}
		dld=(lmn[1]-2*lmn[2]+lmn[3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' Second derivative matrix of the normalized log-likelihood
#' @inheritParams	manf
gumbel_p1_ldd=function(x,t,v1,d1,v2,d2,v3,fd3){
	ldd=matrix(0,3,3)
	for (i in 1:3){
		for (j in i:3){
			ldd[i,j]=gumbel_p1_lmn(x,t,v1,d1,v2,d2,v3,fd3,i,j)
		}
	}
	for (i in 3:2){
		for (j in 1:(i-1)){
			ldd[i,j]=ldd[j,i]
		}
	}
	return(ldd)
}
#' One component of the second derivative of the normalized log-likelihood
#' @inheritParams manf
gumbel_p1_lmnp=function(x,t,v1,d1,v2,d2,v3,fd3,mm,nn,rr){
	d3=fd3*v3
	net4=matrix(0,4,3)
	net6=matrix(0,6,3)
	net8=matrix(0,8,3)
	lmn=matrix(0,8)
	dd=c(d1,d2,d3)
	vv=c(v1,v2,v3)
	vvd=matrix(0,3)
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
			lmn[i]=sum(dgumbel_p1(x,t,ymn=vvd[1],slope=vvd[2],sigma=vvd[3],log=TRUE))/nx
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
			for (j in 1:3){
				vvd[j]=vv[j]+net4[i,j]*dd[j]
			}
			lmn[i]=sum(dgumbel_p1(x,t,ymn=vvd[1],slope=vvd[2],sigma=vvd[3],log=TRUE))/nx
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
			for (j in 1:3){
				vvd[j]=vv[j]+net6[i,j]*dd[j]
			}
			lmn[i]=sum(dgumbel_p1(x,t,ymn=vvd[1],slope=vvd[2],sigma=vvd[3],log=TRUE))/nx
		}
		dld1=(lmn[3]-2*lmn[2]+lmn[1])/(dd[m2]*dd[m2])
		dld2=(lmn[6]-2*lmn[5]+lmn[4])/(dd[m2]*dd[m2])
		dld=(dld2-dld1)/(2*dd[n2])
	}
	return(dld)
}#' Third derivative tensor of the normalized log-likelihood
#' @inheritParams	manf
gumbel_p1_lddd=function(x,t,v1,d1,v2,d2,v3,fd3){
	lddd=array(0,c(3,3,3))
	for (i in 1:3){
		for (j in i:3){
			for (k in j:3){
				lddd[i,j,k]=gumbel_p1_lmnp(x,t,v1,d1,v2,d2,v3,fd3,i,j,k)
			}
		}
	}
# steves dumb algorithm for filling in the non-unique values
	for (i in 1:3){
		for (j in 1:3){
			for (k in 1:3){
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
gumbel_p1_f1f=function(y,t0,v1,d1,v2,d2,v3,fd3){
	d3=fd3*v3
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v2 stuff
	v2m1=v2-1*d2
	v200=v2+0*d2
	v2p1=v2+1*d2
# v3 stuff
	v3m1=v3-1*d3
	v300=v3+0*d3
	v3p1=v3+1*d3
# v1 derivatives
	F1m1=dgumbel_p1(y,t0,ymn=v1m1,slope=v200,sigma=v3)
	F1p1=dgumbel_p1(y,t0,ymn=v1p1,slope=v200,sigma=v3)
# v2 derivatives
	F2m1=dgumbel_p1(y,t0,ymn=v100,slope=v2m1,sigma=v3)
	F2p1=dgumbel_p1(y,t0,ymn=v100,slope=v2p1,sigma=v3)
# v3 derivatives
	F3m1=dgumbel_p1(y,t0,ymn=v100,slope=v200,sigma=v3m1)
	F3p1=dgumbel_p1(y,t0,ymn=v100,slope=v200,sigma=v3p1)
	f1=matrix(0,3,length(y))
	f1[1,]=(F1p1-F1m1)/(2*d1)
	f1[2,]=(F2p1-F2m1)/(2*d2)
	f1[3,]=(F3p1-F3m1)/(2*d3)
	return(f1)
}
#' DMGS equation 2.1, p1 term
#' @inheritParams	manf
gumbel_p1_p1f=function(y,t0,v1,d1,v2,d2,v3,fd3){
	d3=fd3*v3
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v2 stuff
	v2m1=v2-1*d2
	v200=v2+0*d2
	v2p1=v2+1*d2
# v3 stuff
	v3m1=v3-1*d3
	v300=v3+0*d3
	v3p1=v3+1*d3
# v1 derivatives
	F1m1=pgumbel_p1(y,t0,ymn=v1m1,slope=v200,sigma=v3)
	F1p1=pgumbel_p1(y,t0,ymn=v1p1,slope=v200,sigma=v3)
# v2 derivatives
	F2m1=pgumbel_p1(y,t0,ymn=v100,slope=v2m1,sigma=v3)
	F2p1=pgumbel_p1(y,t0,ymn=v100,slope=v2p1,sigma=v3)
# v3 derivatives
	F3m1=pgumbel_p1(y,t0,ymn=v100,slope=v200,sigma=v3m1)
	F3p1=pgumbel_p1(y,t0,ymn=v100,slope=v200,sigma=v3p1)
	p1=matrix(0,3,length(y))
	p1[1,]=(F1p1-F1m1)/(2*d1)
	p1[2,]=(F2p1-F2m1)/(2*d2)
	p1[3,]=(F3p1-F3m1)/(2*d3)
	return(p1)
}
#' DMGS equation 3.3, mu1 term
#' @inheritParams	manf
gumbel_p1_mu1f=function(alpha,t0,v1,d1,v2,d2,v3,fd3){
	q00=qgumbel_p1((1-alpha),t0,ymn=v1,slope=v2,sigma=v3)
	d3=fd3*v3
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v2 stuff
	v2m1=v2-1*d2
	v200=v2+0*d2
	v2p1=v2+1*d2
# v3 stuff
	v3m1=v3-1*d3
	v300=v3+0*d3
	v3p1=v3+1*d3
# v1 derivatives
	F1m1=pgumbel_p1(q00,t0,ymn=v1m1,slope=v200,sigma=v3)
	F1p1=pgumbel_p1(q00,t0,ymn=v1p1,slope=v200,sigma=v3)
# v2 derivatives
	F2m1=pgumbel_p1(q00,t0,ymn=v100,slope=v2m1,sigma=v3)
	F2p1=pgumbel_p1(q00,t0,ymn=v100,slope=v2p1,sigma=v3)
# v3 derivatives
	F3m1=pgumbel_p1(q00,t0,ymn=v100,slope=v200,sigma=v3m1)
	F3p1=pgumbel_p1(q00,t0,ymn=v100,slope=v200,sigma=v3p1)
	mu1=matrix(0,3,length(alpha))
	mu1[1,]=-(F1p1-F1m1)/(2*d1)
	mu1[2,]=-(F2p1-F2m1)/(2*d2)
	mu1[3,]=-(F3p1-F3m1)/(2*d3)
	return(mu1)
}
#' DMGS equation 2.1, f2 term
#' @inheritParams	manf
gumbel_p1_f2f=function(y,t0,v1,d1,v2,d2,v3,fd3){
	d3=fd3*v3
# new method
	dd=c(d1,d2,d3)
	vv=c(v1,v2,v3)
	f2=array(0,c(3,3,length(y)))
	for (i in 1:3){
		for (j in 1:3){
			if(i==j){
				vvm=vv
				vv0=vv
				vvp=vv
				vvm[i]=vv[i]-dd[i]
				vvp[i]=vv[i]+dd[i]
				Fm1=dgumbel_p1(y,t0,ymn=vvm[1],slope=vvm[2],sigma=vvm[3])
				F00=dgumbel_p1(y,t0,ymn=vv0[1],slope=vv0[2],sigma=vv0[3])
				Fp1=dgumbel_p1(y,t0,ymn=vvp[1],slope=vvp[2],sigma=vvp[3])
				f2[i,i,]=(Fp1-2*F00+Fm1)/(dd[i]*dd[i])
			} else if(i<j) {
				vvmm=vv
				vvmp=vv
				vvpm=vv
				vvpp=vv
				vvmm[i]=vv[i]-dd[i];vvmm[j]=vv[j]-dd[j]
				vvmp[i]=vv[i]-dd[i];vvmp[j]=vv[j]+dd[j]
				vvpm[i]=vv[i]+dd[i];vvpm[j]=vv[j]-dd[j]
				vvpp[i]=vv[i]+dd[i];vvpp[j]=vv[j]+dd[j]
				Fm1m1=dgumbel_p1(y,t0,ymn=vvmm[1],slope=vvmm[2],sigma=vvmm[3])
				Fm1p1=dgumbel_p1(y,t0,ymn=vvmp[1],slope=vvmp[2],sigma=vvmp[3])
				Fp1m1=dgumbel_p1(y,t0,ymn=vvpm[1],slope=vvpm[2],sigma=vvpm[3])
				Fp1p1=dgumbel_p1(y,t0,ymn=vvpp[1],slope=vvpp[2],sigma=vvpp[3])
				f2[i,j,]=(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j])
				f2[j,i,]=f2[i,j,]
			}
		}
	}
	return(f2)
}
#' DMGS equation 2.1, p2 term
#' @inheritParams	manf
gumbel_p1_p2f=function(y,t0,v1,d1,v2,d2,v3,fd3){
	d3=fd3*v3
# new method
	dd=c(d1,d2,d3)
	vv=c(v1,v2,v3)
	p2=array(0,c(3,3,length(y)))
	for (i in 1:3){
		for (j in 1:3){
			if(i==j){
				vvm=vv
				vv0=vv
				vvp=vv
				vvm[i]=vv[i]-dd[i]
				vvp[i]=vv[i]+dd[i]
				Fm1=dgumbel_p1(y,t0,ymn=vvm[1],slope=vvm[2],sigma=vvm[3])
				F00=dgumbel_p1(y,t0,ymn=vv0[1],slope=vv0[2],sigma=vv0[3])
				Fp1=dgumbel_p1(y,t0,ymn=vvp[1],slope=vvp[2],sigma=vvp[3])
				p2[i,i,]=(Fp1-2*F00+Fm1)/(dd[i]*dd[i])
			} else if(i<j) {
				vvmm=vv
				vvmp=vv
				vvpm=vv
				vvpp=vv
				vvmm[i]=vv[i]-dd[i];vvmm[j]=vv[j]-dd[j]
				vvmp[i]=vv[i]-dd[i];vvmp[j]=vv[j]+dd[j]
				vvpm[i]=vv[i]+dd[i];vvpm[j]=vv[j]-dd[j]
				vvpp[i]=vv[i]+dd[i];vvpp[j]=vv[j]+dd[j]
				Fm1m1=dgumbel_p1(y,t0,ymn=vvmm[1],slope=vvmm[2],sigma=vvmm[3])
				Fm1p1=dgumbel_p1(y,t0,ymn=vvmp[1],slope=vvmp[2],sigma=vvmp[3])
				Fp1m1=dgumbel_p1(y,t0,ymn=vvpm[1],slope=vvpm[2],sigma=vvpm[3])
				Fp1p1=dgumbel_p1(y,t0,ymn=vvpp[1],slope=vvpp[2],sigma=vvpp[3])
				p2[i,j,]=(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j])
				p2[j,i,]=p2[i,j,]
			}
		}
	}
	return(p2)
}
#' DMGS equation 3.3, mu2 term
#' @inheritParams	manf
gumbel_p1_mu2f=function(alpha,t0,v1,d1,v2,d2,v3,fd3){
	q00=qgumbel_p1((1-alpha),t0,ymn=v1,slope=v2,sigma=v3)
	d3=fd3*v3
# new method
	dd=c(d1,d2,d3)
	vv=c(v1,v2,v3)
	mu2=array(0,c(3,3,length(alpha)))
	for (i in 1:3){
		for (j in 1:3){
			if(i==j){
				vvm=vv
				vv0=vv
				vvp=vv
				vvm[i]=vv[i]-dd[i]
				vvp[i]=vv[i]+dd[i]
				Fm1=pgumbel_p1(q00,t0,ymn=vvm[1],slope=vvm[2],sigma=vvm[3])
				F00=pgumbel_p1(q00,t0,ymn=vv0[1],slope=vv0[2],sigma=vv0[3])
				Fp1=pgumbel_p1(q00,t0,ymn=vvp[1],slope=vvp[2],sigma=vvp[3])
				mu2[i,i,]=-(Fp1-2*F00+Fm1)/(dd[i]*dd[i])
			} else if(i<j) {
				vvmm=vv
				vvmp=vv
				vvpm=vv
				vvpp=vv
				vvmm[i]=vv[i]-dd[i];vvmm[j]=vv[j]-dd[j]
				vvmp[i]=vv[i]-dd[i];vvmp[j]=vv[j]+dd[j]
				vvpm[i]=vv[i]+dd[i];vvpm[j]=vv[j]-dd[j]
				vvpp[i]=vv[i]+dd[i];vvpp[j]=vv[j]+dd[j]
				Fm1m1=pgumbel_p1(q00,t0,ymn=vvmm[1],slope=vvmm[2],sigma=vvmm[3])
				Fm1p1=pgumbel_p1(q00,t0,ymn=vvmp[1],slope=vvmp[2],sigma=vvmp[3])
				Fp1m1=pgumbel_p1(q00,t0,ymn=vvpm[1],slope=vvpm[2],sigma=vvpm[3])
				Fp1p1=pgumbel_p1(q00,t0,ymn=vvpp[1],slope=vvpp[2],sigma=vvpp[3])
				mu2[i,j,]=-(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j])
				mu2[j,i,]=mu2[i,j,]
			}
		}
	}
	return(mu2)
}
#' Gumbel distribution: RHP mean
#' @inheritParams	manf
gumbel_p1_means=function(means,t0,ml_params,lddi,lddd,lambdad_rhp,nx,dim){

	if(means){
# intro
		eulerconstant=0.57721566490153286060651209008240243104215933593992
		v1=ml_params[1]
		v2=ml_params[2]
		v3=ml_params[3]

# ml mean
		eulerconstant=0.57721566490153286060651209008240243104215933593992
		mu_hat=v1+v2*t0
		ml_mean=mu_hat+v3*eulerconstant

# rhp mean
		meand1=array(0,c(3,1))
		meand1[1,1]=1
		meand1[2,1]=t0
		meand1[3,1]=eulerconstant
		meand2=array(0,c(3,3,1)) #but all zero for gumbel
		dmean=dmgs(lddi,lddd,meand1,lambdad_rhp,meand2,dim=3)
		rh_mean=ml_mean+dmean/nx
	}else{
		ml_mean="means not selected"
		rh_mean="means not selected"
	}

	list(ml_mean=ml_mean,rh_mean=rh_mean)

}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inheritParams	manf
gumbel_p1_logscores=function(logscores,x,t,d1,d2,fd3,aderivs=TRUE){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){

			x1=x[-i]
			t1=t[-i]

			dd=dgumbel_p1sub(x1,t1,x[i],t[i],d1,d2,fd3,aderivs)

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
dgumbel_p1sub=function(x,t,y,t0,d1,d2,fd3,aderivs=TRUE){

		nx=length(x)

		lm=lm(x~t)
		v1start=lm$coefficients[1]
		v2start=lm$coefficients[2]
		xhat=v1start+v2start*t
		v3start=sqrt((sum((x-xhat)^2))/(nx-1))
		opt1=optim(c(v1start,v2start,v3start),gumbel_p1_loglik,x=x,t=t,
			control=list(fnscale=-1))
		v1hat=opt1$par[1]
		v2hat=opt1$par[2]
		v3hat=opt1$par[3]
		ml_params=c(v1hat,v2hat,v3hat)

# ml
		muhat=v1hat+v2hat*t0
		ml_pdf=dgumbel(y,mu=muhat,sigma=v3hat)
		ml_cdf=pgumbel(y,mu=muhat,sigma=v3hat)

# rhp
		if(aderivs) ldd=gumbel_p1_ldda(x,t,v1hat,v2hat,v3hat)
		if(!aderivs)ldd=gumbel_p1_ldd(x,t,v1hat,d1,v2hat,d2,v3hat,fd3)
		lddi=solve(ldd)

		if(aderivs) lddd=gumbel_p1_lddda(x,t,v1hat,v2hat,v3hat)
		if(!aderivs)lddd=gumbel_p1_lddd(x,t,v1hat,d1,v2hat,d2,v3hat,fd3)

		if(aderivs) f1=gumbel_p1_f1fa(y,t0,v1hat,v2hat,v3hat)
		if(!aderivs)f1=gumbel_p1_f1f(y,t0,v1hat,d1,v2hat,d2,v3hat,fd3)

		if(aderivs) f2=gumbel_p1_f2fa(y,t0,v1hat,v2hat,v3hat)
		if(!aderivs)f2=gumbel_p1_f2f(y,t0,v1hat,d1,v2hat,d2,v3hat,fd3)

		if(aderivs) p1=gumbel_p1_p1fa(y,t0,v1hat,v2hat,v3hat)
		if(!aderivs)p1=gumbel_p1_p1f(y,t0,v1hat,d1,v2hat,d2,v3hat,fd3)

		if(aderivs) p2=gumbel_p1_p2fa(y,t0,v1hat,v2hat,v3hat)
		if(!aderivs)p2=gumbel_p1_p2f(y,t0,v1hat,d1,v2hat,d2,v3hat,fd3)

		lambdad_rhp=c(0,0,-1/v3hat)
		df=dmgs(lddi,lddd,f1,lambdad_rhp,f2,dim=3)
		dp=dmgs(lddi,lddd,p1,lambdad_rhp,p2,dim=3)
		rh_pdf=pmax(ml_pdf+df/nx,0)
		rh_cdf=pmin(pmax(ml_cdf+dp/nx,0),1)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					rh_pdf=rh_pdf,
					ml_cdf=ml_cdf,
					rh_cdf=rh_cdf)
}
