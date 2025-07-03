#' Waic
#' @inherit manwaic return
#' @inheritParams manf
norm_waic=function(waicscores,x,v1hat,d1,v2hat,fd2,aderivs){
	if(waicscores){

		if(aderivs) f1f=norm_f1fa(x,v1hat,v2hat)
		if(!aderivs)f1f=norm_f1f(x,v1hat,d1,v2hat,fd2)

		if(aderivs) f2f=norm_f2fa(x,v1hat,v2hat)
		if(!aderivs)f2f=norm_f2f(x,v1hat,d1,v2hat,fd2)

		if(aderivs) ldd=norm_ldda(x,v1hat,v2hat)
		if(!aderivs)ldd=norm_ldd(x,v1hat,d1,v2hat,fd2)
		lddi=solve(ldd)
		if(aderivs)	lddd=norm_lddda(x,v1hat,v2hat)
		if(!aderivs)lddd=norm_lddd(x,v1hat,d1,v2hat,fd2)
		fhatx=dnorm(x,mean=v1hat,sd=v2hat)
		lambdad_rhp=c(0,-1/v2hat)
		waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad_rhp,f2f,dim=2)
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
norm_logf=function(params,x){
	m=params[1]
	s=pmax(params[2],.Machine$double.eps)
	logf=sum(dnorm(x,mean=m,sd=s,log=TRUE))-log(s)
	return(logf)
}
#' Maximum likelihood estimator
#' @inherit manloglik return
#' @inheritParams manf
norm_ml_params=function(x){
	mlparams=matrix(0,2)
	nx=length(x)
	mlparams[1]=mean(x)
	mlparams[2]=sqrt((nx-1)/nx)*sd(x)
	return(mlparams)
}
#' Method of moments estimator
#' @return Vector
#' @inheritParams manf
norm_unbiasedv_params=function(x){
	params=matrix(0,2)
	nx=length(x)
	params[1]=mean(x)
	params[2]=sd(x)
	return(params)
}
#' One component of the second derivative of the normalized log-likelihood
#' @inherit manlnn return
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
#' Second derivative matrix of the normalized log-likelihood
#' @inherit manldd return
#' @inheritParams manf
norm_ldd=function(x,v1,d1,v2,fd2){
	ldd=matrix(0,2,2)
	for (i in 1:2){
		for (j in i:2){
			ldd[i,j]=norm_lmn(x,v1,d1,v2,fd2,i,j)
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
#' @inherit manlnnn return
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
#' Third derivative tensor of the normalized log-likelihood
#' @inherit manlddd return
#' @inheritParams manf
norm_lddd=function(x,v1,d1,v2,fd2){
	lddd=array(0,c(2,2,2))
# calculate the unique values
	for (i in 1:2){
		for (j in i:2){
			for (k in j:2){
				lddd[i,j,k]=norm_lmnp(x,v1,d1,v2,fd2,i,j,k)
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
#' DMGS equation 3.3, f1 term
#' @inherit man1f return
#' @inheritParams manf
norm_f1f=function(y,v1,d1,v2,fd2){
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
	F1m1=dnorm(y,mean=v1m1,sd=v200)
	F1p1=dnorm(y,mean=v1p1,sd=v200)
# v2 derivatives
	F2m1=dnorm(y,mean=v100,sd=v2m1)
	F2p1=dnorm(y,mean=v100,sd=v2p1)
	f1=matrix(0,2,length(y))
	f1[1,]=(F1p1-F1m1)/(2*d1)
	f1[2,]=(F2p1-F2m1)/(2*d2)
	return(f1)
}
#' DMGS equation 3.3, f2 term
#' @inherit man2f return
#' @inheritParams manf
norm_f2f=function(y,v1,d1,v2,fd2){
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
	F1m2=dnorm(y,mean=v1m2,sd=v200)
	F1m1=dnorm(y,mean=v1m1,sd=v200)
	F100=dnorm(y,mean=v100,sd=v200)
	F1p1=dnorm(y,mean=v1p1,sd=v200)
	F1p2=dnorm(y,mean=v1p2,sd=v200)
# v2 derivative
	F2m2=dnorm(y,mean=v100,sd=v2m2)
	F2m1=dnorm(y,mean=v100,sd=v2m1)
	F200=dnorm(y,mean=v100,sd=v200)
	F2p1=dnorm(y,mean=v100,sd=v2p1)
	F2p2=dnorm(y,mean=v100,sd=v2p2)
# cross derivative
	Fcm1m1=dnorm(y,mean=v1m1,sd=v2m1)
	Fcm1p1=dnorm(y,mean=v1m1,sd=v2p1)
	Fcp1m1=dnorm(y,mean=v1p1,sd=v2m1)
	Fcp1p1=dnorm(y,mean=v1p1,sd=v2p1)
	f2=array(0,c(2,2,length(y)))
	f2[1,1,]=(F1p1-2*F100+F1m1)/(d1*d1)
	f2[2,2,]=(F2p1-2*F200+F2m1)/(d2*d2)
	f2[1,2,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2)
	# copy
	f2[2,1,]=f2[1,2,]
	return(f2)
}
#' One component of the second derivative of the expected log-likelihood
#' @inherit manlnn return
#' @inheritParams manf
norm_gmn=function(alpha,v1,d1,v2,fd2,mm,nn){
	nx=length(alpha)
	d2=fd2*v2
  x=qnorm((1-alpha),mean=v1,sd=v2)
	net3=matrix(0,3,2)
	net4=matrix(0,4,2)
	lmn=matrix(0,nx,4)
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
			lmn[,i]=dnorm(x,mean=vvd[1],sd=vvd[2],log=TRUE)
		}
		dld=(lmn[,1]-lmn[,2]-lmn[,3]+lmn[,4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:3){
			for (j in 1:2){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[,i]=dnorm(x,mean=vvd[1],sd=vvd[2],log=TRUE)
		}
		dld=(lmn[,1]-2*lmn[,2]+lmn[,3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' Second derivative matrix of the expected per-observation log-likelihood
#' @inherit manldd return
#' @inheritParams manf
norm_gg=function(nx,v1,d1,v2,fd2){
	expinfmat=matrix(0,2,2)
	for (i in 1:2){
		for (j in i:2){
			expinfmat[i,j]=-quad(norm_gmn,xa=0,xb=1,v1=v1,d1=d1,v2=v2,fd2=fd2,mm=i,nn=j)
		}
	}
	for (i in 2:2){
		for (j in 1:(i-1)){
			expinfmat[i,j]=expinfmat[j,i]
		}
	}
 return(expinfmat)
}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inherit manlogscores return
#' @inheritParams manf
norm_logscores=function(logscores,x){

	if(logscores){

		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){
			x1=x[-i]

			dd=dnormsub(x1,x[i])
			ml_params=dd$ml_params

			ml_pdf=dd$ml_pdf
			ml_oos_logscore=ml_oos_logscore+log(ml_pdf)

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
dnormsub=function(x,y,aderivs=TRUE){

		nx=length(x)

# ml
		ml_params=norm_ml_params(x)
		ml_pdf=dnorm(y,mean=ml_params[1],sd=ml_params[2])
		ml_cdf=pnorm(y,mean=ml_params[1],sd=ml_params[2])

# rhp pdf
# -first, convert sigma from maxlik to unbiased
		sgu=ml_params[2]*sqrt(nx/(nx-1))
# then, convert sigma to predictive sigma
		sg1=sgu*sqrt((nx+1)/nx)

		yd=(y-ml_params[1])/sg1
		rh_pdf=(dt(yd,df=nx-1))/sg1
		rh_pdf=pmax(rh_pdf,0)

# rhp cdf
		rh_cdf=pt(yd,df=nx-1)
		rh_cdf=pmin(pmax(rh_cdf,0),1)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					rh_pdf=rh_pdf,
					ml_cdf=ml_cdf,
					rh_cdf=rh_cdf)
}
