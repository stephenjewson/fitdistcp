#' Waic
#' @inherit manwaic return
#' @inheritParams manf
logis_waic=function(waicscores,x,v1hat,d1,v2hat,fd2,lddi,lddd,lambdad,aderivs){
		if(waicscores){

			if(aderivs) f1f=logis_f1fa(x,v1hat,v2hat)
			if(!aderivs)f1f=logis_f1f(x,v1hat,d1,v2hat,fd2)

			if(aderivs) f2f=logis_f2fa(x,v1hat,v2hat)
			if(!aderivs)f2f=logis_f2f(x,v1hat,d1,v2hat,fd2)

			fhatx=dlogis(x,location=v1hat,scale=v2hat)
			waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad,f2f,dim=2)
			waic1=waic$waic1
			waic2=waic$waic2
		}else{
			waic1="extras not selected"
			waic2="extras not selected"
		}
		list(waic1=waic1,waic2=waic2)
}
#' Logf for RUST
#' @inherit manlogf return
#' @inheritParams manf
logis_logf=function(params,x){
#	l=params[1]
#	s=params[2]
#	if(s>0){
#		logf=sum(dlogis(x,location=l,scale=s,log=TRUE))-log(s)
#	}else{
#		logf=-Inf
#	}
	l=params[1]
	s=pmax(params[2],.Machine$double.eps)
	logf=sum(dlogis(x,location=l,scale=s,log=TRUE))-log(s)
	return(logf)
}
#'  log-likelihood function
#' @inherit manloglik return
#' @inheritParams manf
logis_loglik=function(vv,x){
	n=length(x)
	loglik=sum(dlogis(x,location=vv[1],scale=max(vv[2],.Machine$double.eps),log=TRUE))
	return(loglik)
}
#' One component of the second derivative of the normalized log-likelihood
#' @inherit manlnn return
#' @inheritParams manf
logis_lmn=function(x,v1,d1,v2,fd2,mm,nn){
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
			lmn[i]=sum(dlogis(x,location=vvd[1],scale=vvd[2],log=TRUE))/nx
		}
		dld=(lmn[1]-lmn[2]-lmn[3]+lmn[4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:3){
			for (j in 1:2){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[i]=sum(dlogis(x,location=vvd[1],scale=vvd[2],log=TRUE))/nx
		}
		dld=(lmn[1]-2*lmn[2]+lmn[3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' Second derivative matrix of the normalized log-likelihood
#' @inherit manldd return
#' @inheritParams manf
logis_ldd=function(x,v1,d1,v2,fd2){
	nx=length(x)
	ldd=matrix(0,2,2)
	for (i in 1:2){
		for (j in i:2){
			ldd[i,j]=logis_lmn(x,v1,d1,v2,fd2,i,j)
		}
	}
	for (i in 2:1){
		for (j in 1:(i-1)){
			ldd[i,j]=ldd[j,i]
		}
	}
	return(ldd)
}
#' One component of the third derivative of the normalized log-likelihood
#' @inherit manlnnn return
#' @inheritParams manf
logis_lmnp=function(x,v1,d1,v2,fd2,mm,nn,rr){
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
			lmn[i]=sum(dlogis(x,location=vvd[1],scale=vvd[2],log=TRUE))/nx
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
			lmn[i]=sum(dlogis(x,location=vvd[1],scale=vvd[2],log=TRUE))/nx
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
			lmn[i]=sum(dlogis(x,location=vvd[1],scale=vvd[2],log=TRUE))/nx
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
logis_lddd=function(x,v1,d1,v2,fd2){
# calculate the unique values
	lddd=array(0,c(2,2,2))
	for (i in 1:2){
		for (j in i:2){
			for (k in j:2){
				lddd[i,j,k]=logis_lmnp(x,v1,d1,v2,fd2,i,j,k)
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
logis_f1f=function(y,v1,d1,v2,fd2){
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
	F1m1=dlogis(y,location=v1m1,scale=v200)
	F1p1=dlogis(y,location=v1p1,scale=v200)
# v2 derivatives
	F2m1=dlogis(y,location=v100,scale=v2m1)
	F2p1=dlogis(y,location=v100,scale=v2p1)
	f1=matrix(0,2,length(y))
	f1[1,]=(F1p1-F1m1)/(2*d1)
	f1[2,]=(F2p1-F2m1)/(2*d2)
	return(f1)
}
#' DMGS equation 3.3, p1 term
#' @inherit man1f return
#' @inheritParams manf
logis_p1f=function(y,v1,d1,v2,fd2){
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
	F1m1=plogis(y,location=v1m1,scale=v200)
	F1p1=plogis(y,location=v1p1,scale=v200)
# v2 derivatives
	F2m1=plogis(y,location=v100,scale=v2m1)
	F2p1=plogis(y,location=v100,scale=v2p1)
	p1=matrix(0,2,length(y))
	p1[1,]=(F1p1-F1m1)/(2*d1)
	p1[2,]=(F2p1-F2m1)/(2*d2)
	return(p1)
}
#' DMGS equation 3.3, mu1 term
#' @inherit man1f return
#' @inheritParams manf
logis_mu1f=function(alpha,v1,d1,v2,fd2){
	q00=qlogis((1-alpha),location=v1,scale=v2)
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
	F1m1=plogis(q00,location=v1m1,scale=v200)
	F1p1=plogis(q00,location=v1p1,scale=v200)
# v2 derivatives
	F2m1=plogis(q00,location=v100,scale=v2m1)
	F2p1=plogis(q00,location=v100,scale=v2p1)
	mu1=matrix(0,2,length(alpha))
	mu1[1,]=-(F1p1-F1m1)/(2*d1)
	mu1[2,]=-(F2p1-F2m1)/(2*d2)
	return(mu1)
}
#' DMGS equation 3.3, f2 term
#' @inherit man2f return
#' @inheritParams manf
logis_f2f=function(y,v1,d1,v2,fd2){
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
	F1m2=dlogis(y,location=v1m2,scale=v200)
	F1m1=dlogis(y,location=v1m1,scale=v200)
	F100=dlogis(y,location=v100,scale=v200)
	F1p1=dlogis(y,location=v1p1,scale=v200)
	F1p2=dlogis(y,location=v1p2,scale=v200)
# v2 derivative
	F2m2=dlogis(y,location=v100,scale=v2m2)
	F2m1=dlogis(y,location=v100,scale=v2m1)
	F200=dlogis(y,location=v100,scale=v200)
	F2p1=dlogis(y,location=v100,scale=v2p1)
	F2p2=dlogis(y,location=v100,scale=v2p2)
# cross derivative
	Fcm1m1=dlogis(y,location=v1m1,scale=v2m1)
	Fcm1p1=dlogis(y,location=v1m1,scale=v2p1)
	Fcp1m1=dlogis(y,location=v1p1,scale=v2m1)
	Fcp1p1=dlogis(y,location=v1p1,scale=v2p1)
	f2=array(0,c(2,2,length(y)))
	f2[1,1,]=(F1p1-2*F100+F1m1)/(d1*d1)
	f2[2,2,]=(F2p1-2*F200+F2m1)/(d2*d2)
	f2[1,2,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2)
	# copy
	f2[2,1,]=f2[1,2,]
	return(f2)
}
#' DMGS equation 3.3, p2 term
#' @inherit man2f return
#' @inheritParams manf
logis_p2f=function(y,v1,d1,v2,fd2){
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
	F1m2=plogis(y,location=v1m2,scale=v200)
	F1m1=plogis(y,location=v1m1,scale=v200)
	F100=plogis(y,location=v100,scale=v200)
	F1p1=plogis(y,location=v1p1,scale=v200)
	F1p2=plogis(y,location=v1p2,scale=v200)
# v2 derivative
	F2m2=plogis(y,location=v100,scale=v2m2)
	F2m1=plogis(y,location=v100,scale=v2m1)
	F200=plogis(y,location=v100,scale=v200)
	F2p1=plogis(y,location=v100,scale=v2p1)
	F2p2=plogis(y,location=v100,scale=v2p2)
# cross derivative
	Fcm1m1=plogis(y,location=v1m1,scale=v2m1)
	Fcm1p1=plogis(y,location=v1m1,scale=v2p1)
	Fcp1m1=plogis(y,location=v1p1,scale=v2m1)
	Fcp1p1=plogis(y,location=v1p1,scale=v2p1)
	p2=array(0,c(2,2,length(y)))
	p2[1,1,]=(F1p1-2*F100+F1m1)/(d1*d1)
	p2[2,2,]=(F2p1-2*F200+F2m1)/(d2*d2)
	p2[1,2,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2)
	# copy
	p2[2,1,]=p2[1,2,]
	return(p2)
}
#' DMGS equation 3.3, mu2 term
#' @inherit man2f return
#' @inheritParams manf
logis_mu2f=function(alpha,v1,d1,v2,fd2){
	q00=qlogis((1-alpha),location=v1,scale=v2)
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
	F1m2=plogis(q00,location=v1m2,scale=v200)
	F1m1=plogis(q00,location=v1m1,scale=v200)
	F100=plogis(q00,location=v100,scale=v200)
	F1p1=plogis(q00,location=v1p1,scale=v200)
	F1p2=plogis(q00,location=v1p2,scale=v200)
# v2 derivative
	F2m2=plogis(q00,location=v100,scale=v2m2)
	F2m1=plogis(q00,location=v100,scale=v2m1)
	F200=plogis(q00,location=v100,scale=v200)
	F2p1=plogis(q00,location=v100,scale=v2p1)
	F2p2=plogis(q00,location=v100,scale=v2p2)
# cross derivative
	Fcm1m1=plogis(q00,location=v1m1,scale=v2m1)
	Fcm1p1=plogis(q00,location=v1m1,scale=v2p1)
	Fcp1m1=plogis(q00,location=v1p1,scale=v2m1)
	Fcp1p1=plogis(q00,location=v1p1,scale=v2p1)
	mu2[1,1,]=-(F1p1-2*F100+F1m1)/(d1*d1)
	mu2[2,2,]=-(F2p1-2*F200+F2m1)/(d2*d2)
	mu2[1,2,]=-(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2)
	# copy
	mu2[2,1,]=mu2[1,2,]
	return(mu2)
}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inherit manlogscores return
#' @inheritParams manf
logis_logscores=function(logscores,x,d1=0.01,fd2=0.01,aderivs=TRUE){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){
			x1=x[-i]

			dd=dlogis2sub(x1,x[i],d1,fd2,aderivs)

			ml_params=dd$ml_params

			ml_pdf=dd$ml_pdf
			ml_oos_logscore=ml_oos_logscore+log(ml_pdf)

			rh_pdf=dd$rh_pdf
			rh_oos_logscore=rh_oos_logscore+log(max(rh_pdf,.Machine$double.eps))
		}
	}else{
		ml_oos_logscore="extras not selected"
		rh_oos_logscore="extras not selected"
	}
	list(ml_oos_logscore=ml_oos_logscore,rh_oos_logscore=rh_oos_logscore)
}
#' Densities from MLE and RHP
#' @inherit mandsub return
#' @inheritParams manf
dlogis2sub=function(x,y,d1=0.01,fd2=0.01,aderivs=TRUE){

		nx=length(x)

		v1start=mean(x)
		v2start=sd(x)
		opt=optim(c(v1start,v2start),logis_loglik,x=x,control=list(fnscale=-1))
		v1hat=opt$par[1]
		v2hat=opt$par[2]
		ml_params=c(v1hat,v2hat)

# mle
		ml_pdf=dlogis(y,location=v1hat,scale=v2hat)
		ml_cdf=plogis(y,location=v1hat,scale=v2hat)

# rhp
		if(aderivs)	ldd=logis_ldda(x,v1hat,v2hat)
		if(!aderivs)ldd=logis_ldd(x,v1hat,d1,v2hat,fd2)
		lddi=solve(ldd)
		if(aderivs)lddd=logis_lddda(x,v1hat,v2hat)
		if(!aderivs)lddd=logis_lddd(x,v1hat,d1,v2hat,fd2)

		if(aderivs) f1=logis_f1fa(y,v1hat,v2hat)
		if(!aderivs)f1=logis_f1f(y,v1hat,d1,v2hat,fd2)

		if(aderivs) f2=logis_f2fa(y,v1hat,v2hat)
		if(!aderivs)f2=logis_f2f(y,v1hat,d1,v2hat,fd2)

		if(aderivs) p1=logis_p1fa(y,v1hat,v2hat)
		if(!aderivs)p1=logis_p1f(y,v1hat,d1,v2hat,fd2)

		if(aderivs) p2=logis_p2fa(y,v1hat,v2hat)
		if(!aderivs)p2=logis_p2f(y,v1hat,d1,v2hat,fd2)

		lambdad_rhp=c(0,-1/v2hat)
		df1=dmgs(lddi,lddd,f1,lambdad_rhp,f2,dim=2)
		dp1=dmgs(lddi,lddd,p1,lambdad_rhp,p2,dim=2)
		rh_pdf=pmax(ml_pdf+df1/nx,0)
		rh_cdf=pmin(pmax(ml_cdf+dp1/nx,0),1)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					rh_pdf=rh_pdf,
					ml_cdf=ml_cdf,
					rh_cdf=rh_cdf)
}
