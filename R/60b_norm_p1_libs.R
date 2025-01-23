#' Linear regression formula, densities
#' @inheritParams manf
pnorm_p1_formula=function(y,ta,ta0,nx,muhat0,v3hat){
	top=ta0*ta0
	sx=sum(ta*ta)

# convert maxlik sg to sum^2/(n-2), which is used in the standard formula
# noting that maxlik is equivalent to using 1/(n-1)
	sg1=v3hat*sqrt((nx-1)/(nx-2))
	sg2=sg1*sqrt((1+1/nx+top/sx))

	yd=(y-muhat0)/sg2
	rh_cdf=(pt(yd,df=nx-2))
	return(rh_cdf)
}
#' Linear regression formula, densities
#' @inheritParams manf
dnorm_p1_formula=function(y,ta,ta0,nx,muhat0,v3hat){
	top=ta0*ta0
	sx=sum(ta*ta)

# convert maxlik sg to sum^2/(n-2), which is used in the standard formula
	sg1=v3hat*sqrt((nx-1)/(nx-2))
	sg2=sg1*sqrt((1+1/nx+top/sx))

	yd=(y-muhat0)/sg2
	rh_pdf=(dt(yd,df=nx-2))/sg2
	return(rh_pdf)
}
#' Linear regression formula, quantiles
#' @inheritParams manf
qnorm_p1_formula=function(alpha,ta,ta0,nx,muhat0,v3hat){
	top=ta0*ta0
	sx=sum(ta*ta)

# convert maxlik sg to sum^2/(n-2), which is used in the standard formula
	sg1=v3hat*sqrt((nx-1)/(nx-2))
	sg2=sg1*sqrt((1+1/nx+top/sx))

	temp=qt((1-alpha),df=nx-2)
	rh_quantiles=muhat0+temp*sg2
	return(rh_quantiles)
}
#' Waic
#' @inheritParams manf
norm_p1_waic=function(waicscores,x,t,v1hat,d1,v2hat,d2,v3hat,fd3,aderivs=TRUE){
	if(waicscores){
		f1f=norm_p1_f1f(x,t,v1hat,d1,v2hat,d2,v3hat,fd3)
		f2f=norm_p1_f2f(x,t,v1hat,d1,v2hat,d2,v3hat,fd3)
		if(aderivs) ldd=norm_p1_ldda(x,t,v1hat,v2hat,v3hat)
		if(!aderivs)ldd=norm_p1_ldd(x,t,v1hat,d1,v2hat,d2,v3hat,fd3)
		lddi=solve(ldd)

		if(aderivs) lddd=norm_p1_lddda(x,t,v1hat,v2hat,v3hat)
		if(!aderivs)lddd=norm_p1_lddd(x,t,v1hat,d1,v2hat,d2,v3hat,fd3)
		fhatx=dnorm_p1(x,t,ymn=v1hat,slope=v2hat,sigma=v3hat,log=FALSE)
		lambdad_rhp=c(0,0,-1/v3hat)
		waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad_rhp,f2f,dim=3)
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
norm_p1_predictordata=function(x,t,t0,params){
#
# calculate the probabilities of the data using the fited model
#
# note that t may be centred
# -but the params take that into account
# -so the mu is not affected by centering
		a=params[1]
		b=params[2]
		s=params[3]
		mu=a+b*t
		px=pnorm(x,mean=mu,sd=s)
#
# calculate the quantiles for those probabilities at t0
#
		mu0=a+b*t0
		qx=qnorm(px,mean=mu0,sd=s)

	list(predictedparameter=mu,adjustedx=qx)
}
#' Logf for RUST
#' @inheritParams manf
norm_p1_logf=function(params,x,t){
#	a=params[1]
#	b=params[2]
#	s=params[3]
#	mu=a+b*t
#	if(s>0){
#		logf=sum(dnorm(x,mean=mu,sd=s,log=TRUE))
#	}else{
#		logf=-Inf
#	}
	a=params[1]
	b=params[2]
	s=pmax(params[3],.Machine$double.eps)
	mu=a+b*t
	logf=sum(dnorm(x,mean=mu,sd=s,log=TRUE))
	return(logf)
}
#' Normal-with-p1  observed log-likelihood function
#' @inheritParams manf
norm_p1_loglik=function(vv,x,t){
	n=length(x)
	mean=vv[1]+vv[2]*t #so mean is a vector, just like x
	loglik=sum(dnorm(x,mean=mean,sd=max(vv[3],0),log=TRUE))
	return(loglik)
}
#' Maximum likelihood estimator
#' @inheritParams manf
norm_p1_mlparams=function(x,t){
	mlparams=matrix(0,3)
	reg=lm(x~t)
	mlparams[1]=reg$coefficients[1]
	mlparams[2]=reg$coefficients[2]
	mlparams[3]=sd(reg$residuals)
	return(mlparams)
}
#' Normal-with-p1 quantile function
#' @inheritParams manf
qnorm_p1=function(p,t0,ymn,slope,sigma){

	return(qnorm(p,mean=(ymn+slope*t0),sd=sigma))

}
#' Normal-with-p1 density function
#' @inheritParams manf
dnorm_p1=function(x,t0,ymn,slope,sigma,log=FALSE){

	return(dnorm(x,mean=(ymn+slope*t0),sd=sigma,log=log))

}
#' Normal-with-p1 distribution function
#' @inheritParams manf
pnorm_p1=function(x,t0,ymn,slope,sigma){

	return(pnorm(x,mean=(ymn+slope*t0),sd=sigma))

}
#' One component of the second derivative of the normalized log-likelihood
#' @inheritParams manf
norm_p1_lmn=function(x,t,v1,d1,v2,d2,v3,fd3,mm,nn){
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
			lmn[i]=sum(dnorm_p1(x,t,ymn=vvd[1],slope=vvd[2],sigma=vvd[3],log=TRUE))/nx
		}
		dld=(lmn[1]-lmn[2]-lmn[3]+lmn[4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:3){
			for (j in 1:3){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[i]=sum(dnorm_p1(x,t,ymn=vvd[1],slope=vvd[2],sigma=vvd[3],log=TRUE))/nx
		}
		dld=(lmn[1]-2*lmn[2]+lmn[3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' Second derivative matrix of the normalized log-likelihood
#' @inheritParams manf
norm_p1_ldd=function(x,t,v1,d1,v2,d2,v3,fd3){
	nx=length(x)
	ldd=matrix(0,3,3)
	for (i in 1:3){
		for (j in i:3){
			ldd[i,j]=norm_p1_lmn(x,t,v1,d1,v2,d2,v3,fd3,i,j)
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
norm_p1_lmnp=function(x,t,v1,d1,v2,d2,v3,fd3,mm,nn,rr){
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
			lmn[i]=sum(dnorm_p1(x,t,ymn=vvd[1],slope=vvd[2],sigma=vvd[3],log=TRUE))/nx
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
			lmn[i]=sum(dnorm_p1(x,t,ymn=vvd[1],slope=vvd[2],sigma=vvd[3],log=TRUE))/nx
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
			lmn[i]=sum(dnorm_p1(x,t,ymn=vvd[1],slope=vvd[2],sigma=vvd[3],log=TRUE))/nx
		}
		dld1=(lmn[3]-2*lmn[2]+lmn[1])/(dd[m2]*dd[m2])
		dld2=(lmn[6]-2*lmn[5]+lmn[4])/(dd[m2]*dd[m2])
		dld=(dld2-dld1)/(2*dd[n2])
	}
	return(dld)
}
#' Third derivative tensor of the normalized log-likelihood
#' @inheritParams manf
norm_p1_lddd=function(x,t,v1,d1,v2,d2,v3,fd3){
# calculate the unique values
	lddd=array(0,c(3,3,3))
	for (i in 1:3){
		for (j in i:3){
			for (k in j:3){
				lddd[i,j,k]=norm_p1_lmnp(x,t,v1,d1,v2,d2,v3,fd3,i,j,k)
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
#' @inheritParams manf
norm_p1_f1f=function(y,t0,v1,d1,v2,d2,v3,fd3){
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
	F1m1=dnorm_p1(y,t0,ymn=v1m1,slope=v200,sigma=v3)
	F1p1=dnorm_p1(y,t0,ymn=v1p1,slope=v200,sigma=v3)
# v2 derivatives
	F2m1=dnorm_p1(y,t0,ymn=v100,slope=v2m1,sigma=v3)
	F2p1=dnorm_p1(y,t0,ymn=v100,slope=v2p1,sigma=v3)
# v3 derivatives
	F3m1=dnorm_p1(y,t0,ymn=v100,slope=v200,sigma=v3m1)
	F3p1=dnorm_p1(y,t0,ymn=v100,slope=v200,sigma=v3p1)
	f1=matrix(0,3,length(y))
	f1[1,]=(F1p1-F1m1)/(2*d1)
	f1[2,]=(F2p1-F2m1)/(2*d2)
	f1[3,]=(F3p1-F3m1)/(2*d3)
	return(f1)
}
#' DMGS equation 2.1, f2 term
#' @inheritParams manf
norm_p1_f2f=function(y,t0,v1,d1,v2,d2,v3,fd3){
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
				Fm1=dnorm_p1(y,t0,ymn=vvm[1],slope=vvm[2],sigma=vvm[3])
				F00=dnorm_p1(y,t0,ymn=vv0[1],slope=vv0[2],sigma=vv0[3])
				Fp1=dnorm_p1(y,t0,ymn=vvp[1],slope=vvp[2],sigma=vvp[3])
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
				Fm1m1=dnorm_p1(y,t0,ymn=vvmm[1],slope=vvmm[2],sigma=vvmm[3])
				Fm1p1=dnorm_p1(y,t0,ymn=vvmp[1],slope=vvmp[2],sigma=vvmp[3])
				Fp1m1=dnorm_p1(y,t0,ymn=vvpm[1],slope=vvpm[2],sigma=vvpm[3])
				Fp1p1=dnorm_p1(y,t0,ymn=vvpp[1],slope=vvpp[2],sigma=vvpp[3])
				f2[i,j,]=(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j])
				f2[j,i,]=f2[i,j,]
			}
		}
	}
	return(f2)
}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inheritParams manf
norm_p1_logscores=function(logscores,x,t,aderivs=TRUE){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){
			x1=x[-i]
			t1=t[-i]

			dd=dnorm_p1sub(x1,t1,x[i],t[i],aderivs)
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
#' @inheritParams manf
dnorm_p1sub=function(x,t,y,t0,aderivs=TRUE){

		nx=length(x)

		meant=mean(t)
		ta=t-meant
		ta0=t0-meant
# we have to centre this here again because this is used in a cross-validation loop
# note that t0 should never be centred

		ml_params=norm_p1_mlparams(x,t)
		v1hat=ml_params[1]
		v2hat=ml_params[2]
		v3hat=ml_params[3]
		muhat0=v1hat+v2hat*t0
# ml
		ml_pdf=dnorm(y,mean=muhat0,sd=v3hat)
		ml_cdf=pnorm(y,mean=muhat0,sd=v3hat)
# rhp

		rh_pdf=dnorm_p1_formula(y,ta,ta0,nx,muhat0,v3hat)
		rh_pdf=pmax(rh_pdf,0)

		rh_cdf=pnorm_p1_formula(y,ta,ta0,nx,muhat0,v3hat)
		rh_cdf=pmin(pmax(rh_cdf,0),1)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					rh_pdf=rh_pdf,
					ml_cdf=ml_cdf,
					rh_cdf=rh_cdf)
}

