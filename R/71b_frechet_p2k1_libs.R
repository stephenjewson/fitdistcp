#' Waic
#' @inheritParams manf
frechet_p2k1_waic=function(waicscores,x,t,v1hat,d1,v2hat,d2,v3hat,fd3,kloc,
	lddi,lddd,lambdad,aderivs){
		if(waicscores){
			if(aderivs){
				f1f=frechet_p2k1_f1fa(x,t,v1hat,v2hat,v3hat,kloc=kloc)
				f2f=frechet_p2k1_f2fa(x,t,v1hat,v2hat,v3hat,kloc=kloc)
			} else {
				f1f=frechet_p2k1_f1f(x,t,v1hat,d1,v2hat,d2,v3hat,fd3,kloc=kloc)
				f2f=frechet_p2k1_f2f(x,t,v1hat,d1,v2hat,d2,v3hat,fd3,kloc=kloc)
			}
			fhatx=dfrechet_p2k1(x,t,ymn=v1hat,slope=v2hat,lambda=v3hat,log=FALSE,kloc=kloc)
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
frechet_p2k1_predictordata=function(predictordata,x,t,t0,params,kloc){
	if(predictordata){
#
# calculate the probabilities of the data using the fited model
#
		a=params[1]
		b=params[2]
		s=params[3]
		sg=exp(a+b*t)
		px=pfrechet(x,mu=kloc,sigma=sg,lambda=s)
#
# calculate the quantiles for those probabilities at t0
#
		sg0=exp(a+b*t0)
		qx=qfrechet(px,mu=kloc,sigma=sg0,lambda=s)
	} else{
		predictedparameter="predictordata not selected"
		adjustedx="predictordata not selected"
	}

	list(predictedparameter=sg,adjustedx=qx)
}
#' Logf for RUST
#' @inheritParams manf
frechet_p2k1_logf=function(params,x,t,kloc){
#	a=params[1]
#	b=params[2]
#	s=params[3]
#	sg=exp(a+b*t)
#	if(s>0){
#		logf=sum(dfrechet(x,mu=kloc,sigma=sg,lambda=s,log=TRUE))-log(s)
#	}else{
#		logf=-Inf
#	}
	a=params[1]
	b=params[2]
	s=pmax(params[3],.Machine$double.eps)
	sg=exp(a+b*t)
	logf=sum(dfrechet(x,mu=kloc,sigma=sg,lambda=s,log=TRUE))-log(s)
	return(logf)
}
#'  observed log-likelihood function
#' @inheritParams	manf
frechet_p2k1_loglik=function(vv,x,t,kloc){
	mu=vv[1]+vv[2]*t
	loglik=sum(dfrechet(x,mu=kloc,sigma=exp(mu),lambda=max(vv[3],.Machine$double.eps),log=TRUE))
	return(loglik)
}
#' Frechet_k1-with-p2 quantile function
#' @inheritParams	manf
qfrechet_p2k1=function(p,t0,ymn,slope,lambda,kloc){

	mu=(ymn+slope*t0)
	return(qfrechet(p,mu=kloc,sigma=exp(mu),lambda=lambda))

}
#' Frechet_k1-with-p2 density function
#' @inheritParams	manf
dfrechet_p2k1=function(x,t0,ymn,slope,lambda,log=FALSE,kloc){

	mu=(ymn+slope*t0)
	return(dfrechet(x,mu=kloc,sigma=exp(mu),lambda=lambda,log=log))

}
#' Frechet_k1-with-p2 distribution function
#' @inheritParams	manf
pfrechet_p2k1=function(x,t0,ymn,slope,lambda,kloc){

	mu=(ymn+slope*t0)
	return(pfrechet(x,mu=kloc,sigma=exp(mu),lambda=lambda))

}
#' One component of the second derivative of the normalized log-likelihood
#' @inheritParams manf
frechet_p2k1_lmn=function(x,t,v1,d1,v2,d2,v3,fd3,kloc,mm,nn){
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
			lmn[i]=sum(dfrechet_p2k1(x,t,ymn=vvd[1],slope=vvd[2],lambda=vvd[3],kloc=kloc,log=TRUE))/nx
		}
		dld=(lmn[1]-lmn[2]-lmn[3]+lmn[4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:3){
			for (j in 1:3){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[i]=sum(dfrechet_p2k1(x,t,ymn=vvd[1],slope=vvd[2],lambda=vvd[3],kloc=kloc,log=TRUE))/nx
		}
		dld=(lmn[1]-2*lmn[2]+lmn[3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' Second derivative matrix of the normalized log-likelihood
#' @inheritParams	manf
frechet_p2k1_ldd=function(x,t,v1,d1,v2,d2,v3,fd3,kloc){
	ldd=matrix(0,3,3)
	for (i in 1:3){
		for (j in i:3){
			ldd[i,j]=frechet_p2k1_lmn(x,t,v1,d1,v2,d2,v3,fd3,kloc,i,j)
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
frechet_p2k1_lmnp=function(x,t,v1,d1,v2,d2,v3,fd3,kloc,mm,nn,rr){
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
			lmn[i]=sum(dfrechet_p2k1(x,t,ymn=vvd[1],slope=vvd[2],lambda=vvd[3],kloc=kloc,log=TRUE))/nx
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
			lmn[i]=sum(dfrechet_p2k1(x,t,ymn=vvd[1],slope=vvd[2],lambda=vvd[3],kloc=kloc,log=TRUE))/nx
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
			lmn[i]=sum(dfrechet_p2k1(x,t,ymn=vvd[1],slope=vvd[2],lambda=vvd[3],kloc=kloc,log=TRUE))/nx
		}
		dld1=(lmn[3]-2*lmn[2]+lmn[1])/(dd[m2]*dd[m2])
		dld2=(lmn[6]-2*lmn[5]+lmn[4])/(dd[m2]*dd[m2])
		dld=(dld2-dld1)/(2*dd[n2])
	}
	return(dld)
}
#' Third derivative tensor of the normalized log-likelihood
#' @inheritParams	manf
frechet_p2k1_lddd=function(x,t,v1,d1,v2,d2,v3,fd3,kloc){
	lddd=array(0,c(3,3,3))
	for (i in 1:3){
		for (j in i:3){
			for (k in j:3){
				lddd[i,j,k]=frechet_p2k1_lmnp(x,t,v1,d1,v2,d2,v3,fd3,kloc,i,j,k)
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
frechet_p2k1_f1f=function(y,t0,v1,d1,v2,d2,v3,fd3,kloc){
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
	F1m1=dfrechet_p2k1(y,t0,ymn=v1m1,slope=v200,lambda=v3,kloc=kloc)
	F1p1=dfrechet_p2k1(y,t0,ymn=v1p1,slope=v200,lambda=v3,kloc=kloc)
# v2 derivatives
	F2m1=dfrechet_p2k1(y,t0,ymn=v100,slope=v2m1,lambda=v3,kloc=kloc)
	F2p1=dfrechet_p2k1(y,t0,ymn=v100,slope=v2p1,lambda=v3,kloc=kloc)
# v3 derivatives
	F3m1=dfrechet_p2k1(y,t0,ymn=v100,slope=v200,lambda=v3m1,kloc=kloc)
	F3p1=dfrechet_p2k1(y,t0,ymn=v100,slope=v200,lambda=v3p1,kloc=kloc)
	f1=matrix(0,3,length(y))
	f1[1,]=(F1p1-F1m1)/(2*d1)
	f1[2,]=(F2p1-F2m1)/(2*d2)
	f1[3,]=(F3p1-F3m1)/(2*d3)
	return(f1)
}
#' DMGS equation 2.1, p1 term
#' @inheritParams	manf
frechet_p2k1_p1f=function(y,t0,v1,d1,v2,d2,v3,fd3,kloc){
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
	F1m1=pfrechet_p2k1(y,t0,ymn=v1m1,slope=v200,lambda=v3,kloc=kloc)
	F1p1=pfrechet_p2k1(y,t0,ymn=v1p1,slope=v200,lambda=v3,kloc=kloc)
# v2 derivatives
	F2m1=pfrechet_p2k1(y,t0,ymn=v100,slope=v2m1,lambda=v3,kloc=kloc)
	F2p1=pfrechet_p2k1(y,t0,ymn=v100,slope=v2p1,lambda=v3,kloc=kloc)
# v3 derivatives
	F3m1=pfrechet_p2k1(y,t0,ymn=v100,slope=v200,lambda=v3m1,kloc=kloc)
	F3p1=pfrechet_p2k1(y,t0,ymn=v100,slope=v200,lambda=v3p1,kloc=kloc)
	p1=matrix(0,3,length(y))
	p1[1,]=(F1p1-F1m1)/(2*d1)
	p1[2,]=(F2p1-F2m1)/(2*d2)
	p1[3,]=(F3p1-F3m1)/(2*d3)
	return(p1)
}
#' DMGS equation 3.3, mu1 term
#' @inheritParams	manf
frechet_p2k1_mu1f=function(alpha,t0,v1,d1,v2,d2,v3,fd3,kloc){
	q00=qfrechet_p2k1((1-alpha),t0,ymn=v1,slope=v2,lambda=v3,kloc=kloc)
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
	F1m1=pfrechet_p2k1(q00,t0,ymn=v1m1,slope=v200,lambda=v3,kloc=kloc)
	F1p1=pfrechet_p2k1(q00,t0,ymn=v1p1,slope=v200,lambda=v3,kloc=kloc)
# v2 derivatives
	F2m1=pfrechet_p2k1(q00,t0,ymn=v100,slope=v2m1,lambda=v3,kloc=kloc)
	F2p1=pfrechet_p2k1(q00,t0,ymn=v100,slope=v2p1,lambda=v3,kloc=kloc)
# v3 derivatives
	F3m1=pfrechet_p2k1(q00,t0,ymn=v100,slope=v200,lambda=v3m1,kloc=kloc)
	F3p1=pfrechet_p2k1(q00,t0,ymn=v100,slope=v200,lambda=v3p1,kloc=kloc)
	mu1=matrix(0,3,length(alpha))
	mu1[1,]=-(F1p1-F1m1)/(2*d1)
	mu1[2,]=-(F2p1-F2m1)/(2*d2)
	mu1[3,]=-(F3p1-F3m1)/(2*d3)
	return(mu1)
}
#' DMGS equation 2.1, f2 term
#' @inheritParams	manf
frechet_p2k1_f2f=function(y,t0,v1,d1,v2,d2,v3,fd3,kloc){
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
				Fm1=dfrechet_p2k1(y,t0,ymn=vvm[1],slope=vvm[2],lambda=vvm[3],kloc=kloc)
				F00=dfrechet_p2k1(y,t0,ymn=vv0[1],slope=vv0[2],lambda=vv0[3],kloc=kloc)
				Fp1=dfrechet_p2k1(y,t0,ymn=vvp[1],slope=vvp[2],lambda=vvp[3],kloc=kloc)
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
				Fm1m1=dfrechet_p2k1(y,t0,ymn=vvmm[1],slope=vvmm[2],lambda=vvmm[3],kloc=kloc)
				Fm1p1=dfrechet_p2k1(y,t0,ymn=vvmp[1],slope=vvmp[2],lambda=vvmp[3],kloc=kloc)
				Fp1m1=dfrechet_p2k1(y,t0,ymn=vvpm[1],slope=vvpm[2],lambda=vvpm[3],kloc=kloc)
				Fp1p1=dfrechet_p2k1(y,t0,ymn=vvpp[1],slope=vvpp[2],lambda=vvpp[3],kloc=kloc)
				f2[i,j,]=(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j])
				f2[j,i,]=f2[i,j,]
			}
		}
	}
	return(f2)
}
#' DMGS equation 2.1, p2 term
#' @inheritParams	manf
frechet_p2k1_p2f=function(y,t0,v1,d1,v2,d2,v3,fd3,kloc){
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
				Fm1=pfrechet_p2k1(y,t0,ymn=vvm[1],slope=vvm[2],lambda=vvm[3],kloc=kloc)
				F00=pfrechet_p2k1(y,t0,ymn=vv0[1],slope=vv0[2],lambda=vv0[3],kloc=kloc)
				Fp1=pfrechet_p2k1(y,t0,ymn=vvp[1],slope=vvp[2],lambda=vvp[3],kloc=kloc)
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
				Fm1m1=pfrechet_p2k1(y,t0,ymn=vvmm[1],slope=vvmm[2],lambda=vvmm[3],kloc=kloc)
				Fm1p1=pfrechet_p2k1(y,t0,ymn=vvmp[1],slope=vvmp[2],lambda=vvmp[3],kloc=kloc)
				Fp1m1=pfrechet_p2k1(y,t0,ymn=vvpm[1],slope=vvpm[2],lambda=vvpm[3],kloc=kloc)
				Fp1p1=pfrechet_p2k1(y,t0,ymn=vvpp[1],slope=vvpp[2],lambda=vvpp[3],kloc=kloc)
				p2[i,j,]=(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j])
				p2[j,i,]=p2[i,j,]
			}
		}
	}
	return(p2)
}
#' DMGS equation 3.3, mu2 term
#' @inheritParams	manf
frechet_p2k1_mu2f=function(alpha,t0,v1,d1,v2,d2,v3,fd3,kloc){
	q00=qfrechet_p2k1((1-alpha),t0,ymn=v1,slope=v2,lambda=v3,kloc=kloc)
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
				Fm1=pfrechet_p2k1(q00,t0,ymn=vvm[1],slope=vvm[2],lambda=vvm[3],kloc=kloc)
				F00=pfrechet_p2k1(q00,t0,ymn=vv0[1],slope=vv0[2],lambda=vv0[3],kloc=kloc)
				Fp1=pfrechet_p2k1(q00,t0,ymn=vvp[1],slope=vvp[2],lambda=vvp[3],kloc=kloc)
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
				Fm1m1=pfrechet_p2k1(q00,t0,ymn=vvmm[1],slope=vvmm[2],lambda=vvmm[3],kloc=kloc)
				Fm1p1=pfrechet_p2k1(q00,t0,ymn=vvmp[1],slope=vvmp[2],lambda=vvmp[3],kloc=kloc)
				Fp1m1=pfrechet_p2k1(q00,t0,ymn=vvpm[1],slope=vvpm[2],lambda=vvpm[3],kloc=kloc)
				Fp1p1=pfrechet_p2k1(q00,t0,ymn=vvpp[1],slope=vvpp[2],lambda=vvpp[3],kloc=kloc)
				mu2[i,j,]=-(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j])
				mu2[j,i,]=mu2[i,j,]
			}
		}
	}
	return(mu2)
}
#' frechet_k1 distribution: RHP mean
#' @inheritParams	manf
frechet_p2k1_means=function(means,t0,ml_params,lddi,lddd,lambdad_rhp,nx,dim,kloc){

	if(means){
# intro
		v1=ml_params[1]
		v2=ml_params[2]
		v3=ml_params[3]

# ml mean
		mu_hat=v1+v2*t0
		if(mu_hat>1){
			ml_mean=mu_hat*gamma(1-(1/v3))
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
frechet_p2k1_logscores=function(logscores,x,t,d1,d2,fd3,kloc,aderivs=TRUE){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){

			x1=x[-i]
			t1=t[-i]

			dd=dfrechet_p2k1sub(x1,t1,x[i],t[i],d1,d2,fd3,kloc,aderivs)

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
dfrechet_p2k1sub=function(x,t,y,t0,d1,d2,fd3,kloc,aderivs=TRUE){

		nx=length(x)

		lm=lm(x~t)
		v1start=lm$coefficients[1]
		v2start=lm$coefficients[2]
		xhat=v1start+v2start*t
		v3start=sqrt((sum((x-xhat)^2))/(nx-1))
		v1start=kloc
		v2start=0
		v3start=2
		opt1=optim(c(v1start,v2start,v3start),frechet_p2k1_loglik,x=x,t=t,
			kloc=kloc,control=list(fnscale=-1))
		v1hat=opt1$par[1]
		v2hat=opt1$par[2]
		v3hat=opt1$par[3]
		ml_params=c(v1hat,v2hat,v3hat)

# ml
		muhat=v1hat+v2hat*t0
		ml_pdf=dfrechet(y,mu=kloc,sigma=exp(muhat),lambda=v3hat)
		ml_cdf=pfrechet(y,mu=kloc,sigma=exp(muhat),lambda=v3hat)

# rhp
		if(aderivs) ldd=frechet_p2k1_ldda(x,t,v1hat,v2hat,v3hat,kloc=kloc)
		if(!aderivs)ldd=frechet_p2k1_ldd(x,t,v1hat,d1,v2hat,d2,v3hat,fd3,kloc=kloc)
		lddi=solve(ldd)

		if(aderivs) lddd=frechet_p2k1_lddda(x,t,v1hat,v2hat,v3hat,kloc=kloc)
		if(!aderivs)lddd=frechet_p2k1_lddd(x,t,v1hat,d1,v2hat,d2,v3hat,fd3,kloc=kloc)

		if(aderivs) f1=frechet_p2k1_f1fa(y,t0,v1hat,v2hat,v3hat,kloc=kloc)
		if(!aderivs)f1=frechet_p2k1_f1f(y,t0,v1hat,d1,v2hat,d2,v3hat,fd3,kloc=kloc)

		if(aderivs) f2=frechet_p2k1_f2fa(y,t0,v1hat,v2hat,v3hat,kloc=kloc)
		if(!aderivs)f2=frechet_p2k1_f2f(y,t0,v1hat,d1,v2hat,d2,v3hat,fd3,kloc=kloc)

		if(aderivs) p1=frechet_p2k1_p1fa(y,t0,v1hat,v2hat,v3hat,kloc=kloc)
		if(!aderivs)p1=frechet_p2k1_p1f(y,t0,v1hat,d1,v2hat,d2,v3hat,fd3,kloc=kloc)

		if(aderivs) p2=frechet_p2k1_p2fa(y,t0,v1hat,v2hat,v3hat,kloc=kloc)
		if(!aderivs)p2=frechet_p2k1_p2f(y,t0,v1hat,d1,v2hat,d2,v3hat,fd3,kloc=kloc)

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
