#' Waic for RUST
#' @inherit manwaic return
#' @inheritParams manf
lnorm_waic=function(waicscores,x,v1hat,d1,v2hat,fd2,aderivs){
	if(waicscores){

		if(aderivs) f1f=lnorm_f1fa(x,v1hat,v2hat)
		if(!aderivs)f1f=lnorm_f1f(x,v1hat,d1,v2hat,fd2)

		if(aderivs) f2f=lnorm_f2fa(x,v1hat,v2hat)
		if(!aderivs)f2f=lnorm_f2f(x,v1hat,d1,v2hat,fd2)

		if(aderivs)	ldd=lnorm_ldda(x,v1hat,v2hat)
		if(!aderivs)ldd=lnorm_ldd(x,v1hat,d1,v2hat,fd2)
		lddi=solve(ldd)
		if(aderivs)	lddd=lnorm_lddda(x,v1hat,v2hat)
		if(!aderivs)lddd=lnorm_lddd(x,v1hat,d1,v2hat,fd2)
		fhatx=dlnorm(x,meanlog=v1hat,sdlog=v2hat)
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
lnorm_logf=function(params,x){
	m=params[1]
	s=pmax(params[2],.Machine$double.eps)
	logf=sum(dlnorm(x,meanlog=m,sdlog=s,log=TRUE))-log(s)
	return(logf)
}
#' One component of the second derivative of the normalized log-likelihood
#' @inherit manlnn return
#' @inheritParams manf
lnorm_lmn=function(x,v1,d1,v2,fd2,mm,nn){
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
			lmn[i]=sum(dlnorm(x,meanlog=vvd[1],sdlog=vvd[2],log=TRUE))/nx
		}
		dld=(lmn[1]-lmn[2]-lmn[3]+lmn[4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:3){
			for (j in 1:2){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[i]=sum(dlnorm(x,meanlog=vvd[1],sdlog=vvd[2],log=TRUE))/nx
		}
		dld=(lmn[1]-2*lmn[2]+lmn[3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' Second derivative matrix of the lnormalized log-likelihood
#' @inherit manldd return
#' @inheritParams manf
lnorm_ldd=function(x,v1,d1,v2,fd2){
	ldd=matrix(0,2,2)
	for (i in 1:2){
		for (j in i:2){
			ldd[i,j]=lnorm_lmn(x,v1,d1,v2,fd2,i,j)
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
lnorm_lmnp=function(x,v1,d1,v2,fd2,mm,nn,rr){
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
			lmn[i]=sum(dlnorm(x,meanlog=vvd[1],sdlog=vvd[2],log=TRUE))/nx
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
			lmn[i]=sum(dlnorm(x,meanlog=vvd[1],sdlog=vvd[2],log=TRUE))/nx
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
			lmn[i]=sum(dlnorm(x,meanlog=vvd[1],sdlog=vvd[2],log=TRUE))/nx
		}
		dld1=(lmn[3]-2*lmn[2]+lmn[1])/(dd[m2]*dd[m2])
		dld2=(lmn[6]-2*lmn[5]+lmn[4])/(dd[m2]*dd[m2])
		dld=(dld2-dld1)/(2*dd[n2])
	}
	return(dld)
}
#' Third derivative tensor of the lnormalized log-likelihood
#' @inherit manlddd return
#' @inheritParams manf
lnorm_lddd=function(x,v1,d1,v2,fd2){
	lddd=array(0,c(2,2,2))
	for (i in 1:2){
		for (j in i:2){
			for (k in j:2){
				lddd[i,j,k]=lnorm_lmnp(x,v1,d1,v2,fd2,i,j,k)
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
lnorm_f1f=function(y,v1,d1,v2,fd2){
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
	F1m1=dlnorm(y,meanlog=v1m1,sdlog=v200)
	F1p1=dlnorm(y,meanlog=v1p1,sdlog=v200)
# v2 derivatives
	F2m1=dlnorm(y,meanlog=v100,sdlog=v2m1)
	F2p1=dlnorm(y,meanlog=v100,sdlog=v2p1)
	f1=matrix(0,2,length(y))
	f1[1,]=(F1p1-F1m1)/(2*d1)
	f1[2,]=(F2p1-F2m1)/(2*d2)
	return(f1)
}
#' DMGS equation 3.3, f2 term
#' @inherit man2f return
#' @inheritParams manf
lnorm_f2f=function(y,v1,d1,v2,fd2){
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
	F1m2=dlnorm(y,meanlog=v1m2,sdlog=v200)
	F1m1=dlnorm(y,meanlog=v1m1,sdlog=v200)
	F100=dlnorm(y,meanlog=v100,sdlog=v200)
	F1p1=dlnorm(y,meanlog=v1p1,sdlog=v200)
	F1p2=dlnorm(y,meanlog=v1p2,sdlog=v200)
# v2 derivative
	F2m2=dlnorm(y,meanlog=v100,sdlog=v2m2)
	F2m1=dlnorm(y,meanlog=v100,sdlog=v2m1)
	F200=dlnorm(y,meanlog=v100,sdlog=v200)
	F2p1=dlnorm(y,meanlog=v100,sdlog=v2p1)
	F2p2=dlnorm(y,meanlog=v100,sdlog=v2p2)
# cross derivative
	Fcm1m1=dlnorm(y,meanlog=v1m1,sdlog=v2m1)
	Fcm1p1=dlnorm(y,meanlog=v1m1,sdlog=v2p1)
	Fcp1m1=dlnorm(y,meanlog=v1p1,sdlog=v2m1)
	Fcp1p1=dlnorm(y,meanlog=v1p1,sdlog=v2p1)
	f2=array(0,c(2,2,length(y)))
	f2[1,1,]=(F1p1-2*F100+F1m1)/(d1*d1)
	f2[2,2,]=(F2p1-2*F200+F2m1)/(d2*d2)
	f2[1,2,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2)
	# copy
	f2[2,1,]=f2[1,2,]
	return(f2)
}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inherit manlogscores return
#' @inheritParams manf
lnorm_logscores=function(logscores,x){

	if(logscores){
		y=log(x)
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){

			x1=x[-i]

			dd=dlnormsub(x1,x[i])

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
dlnormsub=function(x,y,aderivs=TRUE){

# great potential for confusion here
# x and y are both lognormal
# logx and logy are both normal

		nx=length(x)
		logx=log(x)
		logy=log(y)
# ml
		ml_params=norm_ml_params(logx) #this really should be norm not lnorm
		ml_pdf=dlnorm(y,meanlog=ml_params[1],sdlog=ml_params[2])
		ml_cdf=plnorm(y,meanlog=ml_params[1],sdlog=ml_params[2])

# rhp pdf
# -first, convert sigma from maxlik to unbiased
		sgu=ml_params[2]*sqrt(nx/(nx-1))
# then, convert sigma to predictive sigma
		sg1=sgu*sqrt((nx+1)/nx)

		logyd=(logy-ml_params[1])/sg1
		rh_pdf=(dt(logyd,df=nx-1))/(y*sg1) #note the extra exp term here, to convert to loglnormal predictive density
		rh_pdf=pmax(rh_pdf,0)

# rhp cdf
		rh_cdf=pt(logyd,df=nx-1)
		rh_cdf=pmin(pmax(rh_cdf,0),1)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					rh_pdf=rh_pdf,
					ml_cdf=ml_cdf,
					rh_cdf=rh_cdf)
}
