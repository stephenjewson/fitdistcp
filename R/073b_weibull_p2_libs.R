#' Waic
#' @inherit manwaic return
#' @inheritParams manf
weibull_p2_waic=function(waicscores,x,t,v1hat,fd1,v2hat,d2,v3hat,d3,
	lddi,lddd,lambdad,aderivs=TRUE){
		if(waicscores){
			if(aderivs){
				f1f=weibull_p2_f1fa(x,t,v1hat,v2hat,v3hat)
				f2f=weibull_p2_f2fa(x,t,v1hat,v2hat,v3hat)
			} else {
				f1f=weibull_p2_f1f(x,t,v1hat,fd1,v2hat,d2,v3hat,d3)
				f2f=weibull_p2_f2f(x,t,v1hat,fd1,v2hat,d2,v3hat,d3)
			}
			fhatx=dweibull_p2(x,t,shape=v1hat,ymn=v2hat,slope=v3hat,log=FALSE)
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
#' @inherit manpredictor return
#' @inheritParams manf
weibull_p2_predictordata=function(predictordata,x,t,t0,params){
	if(predictordata){
#
# calculate the probabilities of the data using the fited model
#
		sh=min(20,params[1])
		a=params[2]
		b=params[3]
		sc=exp(a+b*t)
		px=pweibull(x,shape=sh,scale=sc)
#
# calculate the quantiles for those probabilities at t0
#
		sc0=exp(a+b*t0)
		qx=qweibull(px,shape=sh,scale=sc0)
	} else{
		predictedparameter="predictordata not selected"
		adjustedx="predictordata not selected"
	}

	list(predictedparameter=sc,adjustedx=qx)
}
#' Logf for RUST
#' @inherit manlogf return
#' @inheritParams manf
weibull_p2_logf=function(params,x,t){
#	sh=min(50,params[1]) #high values of shape give NaNs in dweibull
#	a=params[2]
#	b=params[3]
#	sc=exp(a+b*t)
#	if((min(sh)>sqrt(.Machine$double.eps))&(min(sc)>sqrt(.Machine$double.eps))){
## if sc or s get too small, dweibull crashes, it seems
#	logf=sum(dweibull(x,shape=sh,scale=sc,log=TRUE))-log(sh)
#	}else{
#		logf=-Inf
#	}
	sh=pmax(min(20,params[1]),sqrt(.Machine$double.eps)) #high values of shape give NaNs in dweibull
	a=params[2]
	b=params[3]
	sc=pmax(exp(a+b*t),sqrt(.Machine$double.eps))
	logf=sum(dweibull(x,shape=sh,scale=sc,log=TRUE))-log(sh)
	if(is.na(logf)){
		message("dweibull is giving NaNs again...let's have a look why.")
		message("logf=",logf)
		message("a=",a)
		message("b=",b)
		message("sh=",sh)
		message("sc=",sc)
		for(i in 1:length(x)){
			message("x,d=",x[i],dweibull(x[i],shape=sh,scale=sc,log=TRUE))
		}
	}
	return(logf)
}
#'  observed log-likelihood function
#' @inherit manloglik return
#' @inheritParams	manf
weibull_p2_loglik=function(vv,x,t){
	sh=pmax(min(20,vv[1]),.Machine$double.eps)
	sc=pmax(exp(vv[2]+vv[3]*t),.Machine$double.eps)
	loglik=sum(dweibull(x,shape=sh,scale=sc,log=TRUE))
	if(is.na(loglik))message("\n B:",loglik,sh,sc)
	return(loglik)
}
#' Weibull-with-p1 quantile function
#' @inherit manvector return
#' @inheritParams	manf
qweibull_p2=function(p,t0,shape,ymn,slope){

	sc=exp(ymn+slope*t0)
	return(qweibull(p,shape=shape,scale=sc))

}
#' Weibull-with-p1 density function
#' @inherit manvector return
#' @inheritParams	manf
dweibull_p2=function(x,t0,shape,ymn,slope,log=FALSE){

	shape=pmax(min(20,shape),.Machine$double.eps)
	sc=pmax(exp(ymn+slope*t0),.Machine$double.eps)
	return(dweibull(x,shape=shape,scale=sc,log=log))

}
#' Weibull-with-p1 distribution function
#' @inherit manvector return
#' @inheritParams	manf
pweibull_p2=function(x,t0,shape,ymn,slope){

	sc=exp(ymn+slope*t0)
	return(pweibull(x,shape=shape,scale=sc))

}
#' One component of the second derivative of the normalized log-likelihood
#' @inherit manlnn return
#' @inheritParams manf
weibull_p2_lmn=function(x,t,v1,fd1,v2,d2,v3,d3,mm,nn){
	d1=fd1*v1
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
			lmn[i]=sum(dweibull_p2(x,t,shape=vvd[1],ymn=vvd[2],slope=vvd[3],log=TRUE))/nx
		}
		dld=(lmn[1]-lmn[2]-lmn[3]+lmn[4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:3){
			for (j in 1:3){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[i]=sum(dweibull_p2(x,t,shape=vvd[1],ymn=vvd[2],slope=vvd[3],log=TRUE))/nx
		}
		dld=(lmn[1]-2*lmn[2]+lmn[3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' Second derivative matrix of the normalized log-likelihood
#' @inherit manldd return
#' @inheritParams	manf
weibull_p2_ldd=function(x,t,v1,fd1,v2,d2,v3,d3){
	ldd=matrix(0,3,3)
	for (i in 1:3){
		for (j in i:3){
			ldd[i,j]=weibull_p2_lmn(x,t,v1,fd1,v2,d2,v3,d3,i,j)
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
#' @inherit manlnnn return
#' @inheritParams manf
weibull_p2_lmnp=function(x,t,v1,fd1,v2,d2,v3,d3,mm,nn,rr){
	d1=fd1*v1
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
			lmn[i]=sum(dweibull_p2(x,t,shape=vvd[1],ymn=vvd[2],slope=vvd[3],log=TRUE))/nx
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
			lmn[i]=sum(dweibull_p2(x,t,shape=vvd[1],ymn=vvd[2],slope=vvd[3],log=TRUE))/nx
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
			lmn[i]=sum(dweibull_p2(x,t,shape=vvd[1],ymn=vvd[2],slope=vvd[3],log=TRUE))/nx
		}
		dld1=(lmn[3]-2*lmn[2]+lmn[1])/(dd[m2]*dd[m2])
		dld2=(lmn[6]-2*lmn[5]+lmn[4])/(dd[m2]*dd[m2])
		dld=(dld2-dld1)/(2*dd[n2])
	}
	return(dld)
}
#' Third derivative tensor of the normalized log-likelihood
#' @inherit manlddd return
#' @inheritParams	manf
weibull_p2_lddd=function(x,t,v1,fd1,v2,d2,v3,d3){
	lddd=array(0,c(3,3,3))
	for (i in 1:3){
		for (j in i:3){
			for (k in j:3){
				lddd[i,j,k]=weibull_p2_lmnp(x,t,v1,fd1,v2,d2,v3,d3,i,j,k)
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
#' @inherit man1f return
#' @inheritParams	manf
weibull_p2_f1f=function(y,t0,v1,fd1,v2,d2,v3,d3){
	d1=fd1*v1
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
	F1m1=dweibull_p2(y,t0,shape=v1m1,ymn=v200,slope=v300)
	F1p1=dweibull_p2(y,t0,shape=v1p1,ymn=v200,slope=v300)
# v2 derivatives
	F2m1=dweibull_p2(y,t0,shape=v1,ymn=v2m1,slope=v300)
	F2p1=dweibull_p2(y,t0,shape=v1,ymn=v2p1,slope=v300)
# v3 derivatives
	F3m1=dweibull_p2(y,t0,shape=v1,ymn=v200,slope=v3m1)
	F3p1=dweibull_p2(y,t0,shape=v1,ymn=v200,slope=v3p1)
	f1=matrix(0,3,length(y))
	f1[1,]=(F1p1-F1m1)/(2*d1)
	f1[2,]=(F2p1-F2m1)/(2*d2)
	f1[3,]=(F3p1-F3m1)/(2*d3)
	return(f1)
}
#' DMGS equation 2.1, p1 term
#' @inherit man1f return
#' @inheritParams	manf
weibull_p2_p1f=function(y,t0,v1,fd1,v2,d2,v3,d3){
	d1=fd1*v1
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
	F1m1=pweibull_p2(y,t0,shape=v1m1,ymn=v200,slope=v300)
	F1p1=pweibull_p2(y,t0,shape=v1p1,ymn=v200,slope=v300)
# v2 derivatives
	F2m1=pweibull_p2(y,t0,shape=v1,ymn=v2m1,slope=v300)
	F2p1=pweibull_p2(y,t0,shape=v1,ymn=v2p1,slope=v300)
# v3 derivatives
	F3m1=pweibull_p2(y,t0,shape=v1,ymn=v200,slope=v3m1)
	F3p1=pweibull_p2(y,t0,shape=v1,ymn=v200,slope=v3p1)
	p1=matrix(0,3,length(y))
	p1[1,]=(F1p1-F1m1)/(2*d1)
	p1[2,]=(F2p1-F2m1)/(2*d2)
	p1[3,]=(F3p1-F3m1)/(2*d3)
	return(p1)
}
#' DMGS equation 3.3, mu1 term
#' @inherit man1f return
#' @inheritParams	manf
weibull_p2_mu1f=function(alpha,t0,v1,fd1,v2,d2,v3,d3){
	q00=qweibull_p2((1-alpha),t0,shape=v1,ymn=v2,slope=v3)
	d1=fd1*v1
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
	F1m1=pweibull_p2(q00,t0,shape=v1m1,ymn=v200,slope=v300)
	F1p1=pweibull_p2(q00,t0,shape=v1p1,ymn=v200,slope=v300)
# v2 derivatives
	F2m1=pweibull_p2(q00,t0,shape=v1,ymn=v2m1,slope=v300)
	F2p1=pweibull_p2(q00,t0,shape=v1,ymn=v2p1,slope=v300)
# v3 derivatives
	F3m1=pweibull_p2(q00,t0,shape=v1,ymn=v200,slope=v3m1)
	F3p1=pweibull_p2(q00,t0,shape=v1,ymn=v200,slope=v3p1)
	mu1=matrix(0,3,length(alpha))
	mu1[1,]=-(F1p1-F1m1)/(2*d1)
	mu1[2,]=-(F2p1-F2m1)/(2*d2)
	mu1[3,]=-(F3p1-F3m1)/(2*d3)
	return(mu1)
}
#' DMGS equation 2.1, f2 term
#' @inherit man2f return
#' @inheritParams	manf
weibull_p2_f2f=function(y,t0,v1,fd1,v2,d2,v3,d3){
	d1=fd1*v1
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
				Fm1=dweibull_p2(y,t0,shape=vvm[1],ymn=vvm[2],slope=vvm[3])
				F00=dweibull_p2(y,t0,shape=vv0[1],ymn=vv0[2],slope=vv0[3])
				Fp1=dweibull_p2(y,t0,shape=vvp[1],ymn=vvp[2],slope=vvp[3])
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
				Fm1m1=dweibull_p2(y,t0,shape=vvmm[1],ymn=vvmm[2],slope=vvmm[3])
				Fm1p1=dweibull_p2(y,t0,shape=vvmp[1],ymn=vvmp[2],slope=vvmp[3])
				Fp1m1=dweibull_p2(y,t0,shape=vvpm[1],ymn=vvpm[2],slope=vvpm[3])
				Fp1p1=dweibull_p2(y,t0,shape=vvpp[1],ymn=vvpp[2],slope=vvpp[3])
				f2[i,j,]=(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j])
				f2[j,i,]=f2[i,j,]
			}
		}
	}
	return(f2)
}
#' DMGS equation 2.1, p2 term
#' @inherit man2f return
#' @inheritParams	manf
weibull_p2_p2f=function(y,t0,v1,fd1,v2,d2,v3,d3){
	d1=fd1*v1
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
				Fm1=pweibull_p2(y,t0,shape=vvm[1],ymn=vvm[2],slope=vvm[3])
				F00=pweibull_p2(y,t0,shape=vv0[1],ymn=vv0[2],slope=vv0[3])
				Fp1=pweibull_p2(y,t0,shape=vvp[1],ymn=vvp[2],slope=vvp[3])
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
				Fm1m1=pweibull_p2(y,t0,shape=vvmm[1],ymn=vvmm[2],slope=vvmm[3])
				Fm1p1=pweibull_p2(y,t0,shape=vvmp[1],ymn=vvmp[2],slope=vvmp[3])
				Fp1m1=pweibull_p2(y,t0,shape=vvpm[1],ymn=vvpm[2],slope=vvpm[3])
				Fp1p1=pweibull_p2(y,t0,shape=vvpp[1],ymn=vvpp[2],slope=vvpp[3])
				p2[i,j,]=(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j])
				p2[j,i,]=p2[i,j,]
			}
		}
	}
	return(p2)
}
#' DMGS equation 3.3, mu2 term
#' @inherit man2f return
#' @inheritParams	manf
weibull_p2_mu2f=function(alpha,t0,v1,fd1,v2,d2,v3,d3){
	q00=qweibull_p2((1-alpha),t0,shape=v1,ymn=v2,slope=v3)
	d1=fd1*v1
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
				Fm1=pweibull_p2(q00,t0,shape=vvm[1],ymn=vvm[2],slope=vvm[3])
				F00=pweibull_p2(q00,t0,shape=vv0[1],ymn=vv0[2],slope=vv0[3])
				Fp1=pweibull_p2(q00,t0,shape=vvp[1],ymn=vvp[2],slope=vvp[3])
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
				Fm1m1=pweibull_p2(q00,t0,shape=vvmm[1],ymn=vvmm[2],slope=vvmm[3])
				Fm1p1=pweibull_p2(q00,t0,shape=vvmp[1],ymn=vvmp[2],slope=vvmp[3])
				Fp1m1=pweibull_p2(q00,t0,shape=vvpm[1],ymn=vvpm[2],slope=vvpm[3])
				Fp1p1=pweibull_p2(q00,t0,shape=vvpp[1],ymn=vvpp[2],slope=vvpp[3])
				mu2[i,j,]=-(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j])
				mu2[j,i,]=mu2[i,j,]
			}
		}
	}
	return(mu2)
}
#' weibull distribution: RHP mean
#' @inherit manmeans return
#' @inheritParams	manf
weibull_p2_means=function(means,t0,ml_params,lddi,lddd,lambdad_rhp,nx,dim){

	if(means){
# intro
		v1=ml_params[1]
		v2=ml_params[2]
		v3=ml_params[3]

# ml mean
		mu_hat=v2+v3*t0
		arg=1+(1/v1)
		ml_mean=exp(mu_hat)*gamma(arg)

# rhp mean
# starts with derivatives of the log of the mean, to make the algebra easier
#
		lmeand2=array(0,c(3,1))
		lmeand2[1,1]=1
		lmeand2[2,1]=t0
		lmeand2[3,1]=(1/(v1*v1))*digamma(arg)
		meand2=array(0,c(3,1))
		meand2[1,1]=ml_mean*lmeand2[1,1]
		meand2[2,1]=ml_mean*lmeand2[2,1]
		meand2[3,1]=ml_mean*lmeand2[3,1]
#
		lmeand3=array(0,c(3,3,1))
		lmeand3[1,1,1]=0 #alpha-alpha
		lmeand3[1,2,1]=0 #alpha-beta
		lmeand3[1,3,1]=0 #alpha-k
		lmeand3[2,1,1]=0 #beta-alpha
		lmeand3[2,2,1]=0 #beta-beta
		lmeand3[2,2,1]=0 #beta-k
		lmeand3[3,1,1]=0 #alpha-k
		lmeand3[3,2,1]=0 #beta-k
		lmeand3[3,3,1]=2*digamma(arg)/(v1*v1*v1)+trigamma(arg)/(v1*v1*v1*v1)
#
		meand3=array(0,c(3,3,1))
		meand3[1,1,1]=lmeand2[1,1]/ml_mean								#alpha-alpha
		meand3[1,2,1]=lmeand2[1,1]*lmeand2[2,1]/ml_mean	#alpha-beta
		meand3[1,3,1]=lmeand2[1,1]*lmeand2[3,1]/ml_mean	#alpha-k
		meand3[2,1,1]=meand3[1,2,1] 										#beta-alpha
		meand3[2,2,1]=lmeand2[2,1]*lmeand2[2,1]/ml_mean	#beta-beta
		meand3[2,3,1]=lmeand2[2,1]*lmeand2[3,1]/ml_mean	#beta-k
		meand3[3,1,1]=meand3[1,3,1] 										#alpha-k
		meand3[3,2,1]=ml_mean*lmeand3[3,3,1]+lmeand2[3,1]*lmeand2[3,1]/ml_mean #beta-k
#
		dmean=dmgs(lddi,lddd,meand2,lambdad_rhp,meand3,dim=3)
		rh_mean=ml_mean+dmean/nx
	}else{
		ml_mean="means not selected"
		rh_mean="means not selected"
	}


	list(ml_mean=ml_mean,rh_mean=rh_mean)

}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inherit manlogscores return
#' @inheritParams	manf
weibull_p2_logscores=function(logscores,x,t,fd1,d2,d3,aderivs){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){

			x1=x[-i]
			t1=t[-i]

			dd=dweibull_p2sub(x1,t1,x[i],t[i],fd1,d2,d3,aderivs)

			ml_pdf=dd$ml_pdf
			ml_oos_logscore=ml_oos_logscore+log(ml_pdf)

			rh_pdf=dd$rh_pdf
			rh_oos_logscore=rh_oos_logscore+log(rh_pdf)
		}
	}else{
		ml_oos_logscore="extras not selected"
		rh_oos_logscore="extras not selected"
	}
	list(ml_oos_logscore=ml_oos_logscore,rh_oos_logscore=rh_oos_logscore)
}
#' Densities from MLE and RHP
#' @inherit mandsub return
#' @inheritParams	manf
dweibull_p2sub=function(x,t,y,t0,fd1,d2,d3,aderivs=TRUE){

		nx=length(x)

		opt1=optim(c(1,0,0),weibull_p2_loglik,x=x,t=t,
			control=list(fnscale=-1))
		v1hat=opt1$par[1]
		v2hat=opt1$par[2]
		v3hat=opt1$par[3]
		ml_params=c(v1hat,v2hat,v3hat)

# ml
		sc=exp(v2hat+v3hat*t0)
		ml_pdf=dweibull(y,shape=v1hat,scale=sc)
		ml_cdf=pweibull(y,shape=v1hat,scale=sc)

# rhp
		if(aderivs) ldd=weibull_p2_ldda(x,t,v1hat,v2hat,v3hat)
		if(!aderivs)ldd=weibull_p2_ldd(x,t,v1hat,fd1,v2hat,d2,v3hat,d3)
		lddi=solve(ldd)

		if(aderivs) lddd=weibull_p2_lddda(x,t,v1hat,v2hat,v3hat)
		if(!aderivs)lddd=weibull_p2_lddd(x,t,v1hat,fd1,v2hat,d2,v3hat,d3)

		if(aderivs) f1=weibull_p2_f1fa(y,t0,v1hat,v2hat,v3hat)
		if(!aderivs)f1=weibull_p2_f1f(y,t0,v1hat,fd1,v2hat,d2,v3hat,d3)

		if(aderivs) f2=weibull_p2_f2fa(y,t0,v1hat,v2hat,v3hat)
		if(!aderivs)f2=weibull_p2_f2f(y,t0,v1hat,fd1,v2hat,d2,v3hat,d3)

		if(aderivs) p1=weibull_p2_p1fa(y,t0,v1hat,v2hat,v3hat)
		if(!aderivs)p1=weibull_p2_p1f(y,t0,v1hat,fd1,v2hat,d2,v3hat,d3)

		if(aderivs) p2=weibull_p2_p2fa(y,t0,v1hat,v2hat,v3hat)
		if(!aderivs)p2=weibull_p2_p2f(y,t0,v1hat,fd1,v2hat,d2,v3hat,d3)

		lambdad_rhp=c(-1/v1hat,0,0)
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
