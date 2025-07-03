#' Waic for RUST
#' @inherit manwaic return
#' @inheritParams manf
gnorm_waic=function(waicscores,x,v1hat,d1,v2hat,fd2,kbeta,lddi,lddd,lambdad,aderivs){
		if(waicscores){

			if(aderivs) f1f=gnorm_k3_f1fa(x,v1hat,v2hat,kbeta)
			if(!aderivs)f1f=gnorm_k3_f1f(x,v1hat,d1,v2hat,fd2,kbeta)

			if(aderivs) f2f=gnorm_k3_f2fa(x,v1hat,v2hat,kbeta)
			if(!aderivs)f2f=gnorm_k3_f2f(x,v1hat,d1,v2hat,fd2,kbeta)

			fhatx=dgnorm(x,mu=v1hat,alpha=v2hat,beta=kbeta)
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
gnorm_k3_logf=function(params,x,kbeta){
	m=params[1]
	s=pmax(params[2],.Machine$double.eps)
	logf=sum(dgnorm(x,mu=m,alpha=s,beta=kbeta,log=TRUE))-log(s)
	return(logf)
}
#'  log-likelihood function
#' @inherit manloglik return
#' @inheritParams manf
gnorm_k3_loglik=function(vv,x,kbeta){
	n=length(x)
	loglik=sum(dgnorm(x,mu=vv[1],alpha=max(vv[2],.Machine$double.eps),beta=kbeta,log=TRUE))
	return(loglik)
}
#' One component of the second derivative of the normalized log-likelihood
#' @inherit manlnn return
#' @inheritParams manf
gnorm_k3_lmn=function(x,v1,d1,v2,fd2,kbeta,mm,nn){
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
			lmn[i]=sum(dgnorm(x,mu=vvd[1],alpha=vvd[2],beta=kbeta,log=TRUE))/nx
		}
		dld=(lmn[1]-lmn[2]-lmn[3]+lmn[4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:3){
			for (j in 1:2){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[i]=sum(dgnorm(x,mu=vvd[1],alpha=vvd[2],beta=kbeta,log=TRUE))/nx
		}
		dld=(lmn[1]-2*lmn[2]+lmn[3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' Second derivative matrix of the normalized log-likelihood
#' @inherit manldd return
#' @inheritParams manf
gnorm_k3_ldd=function(x,v1,d1,v2,fd2,kbeta){
	nx=length(x)
	ldd=matrix(0,2,2)
	for (i in 1:2){
		for (j in i:2){
			ldd[i,j]=gnorm_k3_lmn(x,v1,d1,v2,fd2,kbeta,i,j)
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
gnorm_lmnp=function(x,v1,d1,v2,fd2,kbeta,mm,nn,rr){
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
			lmn[i]=sum(dgnorm(x,mu=vvd[1],alpha=vvd[2],beta=kbeta,log=TRUE))/nx
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
			lmn[i]=sum(dgnorm(x,mu=vvd[1],alpha=vvd[2],beta=kbeta,log=TRUE))/nx
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
			lmn[i]=sum(dgnorm(x,mu=vvd[1],alpha=vvd[2],beta=kbeta,log=TRUE))/nx
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
gnorm_k3_lddd=function(x,v1,d1,v2,fd2,kbeta){
	lddd=array(0,c(2,2,2))
	for (i in 1:2){
		for (j in i:2){
			for (k in j:2){
				lddd[i,j,k]=gnorm_lmnp(x,v1,d1,v2,fd2,kbeta,i,j,k)
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
gnorm_k3_f1f=function(y,v1,d1,v2,fd2,kbeta){
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
	F1m1=dgnorm(y,mu=v1m1,alpha=v200,beta=kbeta)
	F1p1=dgnorm(y,mu=v1p1,alpha=v200,beta=kbeta)
# v2 derivatives
	F2m1=dgnorm(y,mu=v100,alpha=v2m1,beta=kbeta)
	F2p1=dgnorm(y,mu=v100,alpha=v2p1,beta=kbeta)
	f1=matrix(0,2,length(y))
	f1[1,]=(F1p1-F1m1)/(2*d1)
	f1[2,]=(F2p1-F2m1)/(2*d2)
	return(f1)
}
#' DMGS equation 3.3, p1 term
#' @inherit man1f return
#' @inheritParams manf
gnorm_k3_p1f=function(y,v1,d1,v2,fd2,kbeta){
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
	F1m1=pgnorm(y,mu=v1m1,alpha=v200,beta=kbeta)
	F1p1=pgnorm(y,mu=v1p1,alpha=v200,beta=kbeta)
# v2 derivatives
	F2m1=pgnorm(y,mu=v100,alpha=v2m1,beta=kbeta)
	F2p1=pgnorm(y,mu=v100,alpha=v2p1,beta=kbeta)
	p1=matrix(0,2,length(y))
	p1[1,]=(F1p1-F1m1)/(2*d1)
	p1[2,]=(F2p1-F2m1)/(2*d2)
	return(p1)
}
#' DMGS equation 3.3, mu1 term
#' @inherit man1f return
#' @inheritParams manf
gnorm_k3_mu1f=function(alpha,v1,d1,v2,fd2,kbeta){
	q00=qgnorm((1-alpha),mu=v1,alpha=v2,beta=kbeta)
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
	F1m1=pgnorm(q00,mu=v1m1,alpha=v200,beta=kbeta)
	F1p1=pgnorm(q00,mu=v1p1,alpha=v200,beta=kbeta)
# v2 derivatives
	F2m1=pgnorm(q00,mu=v100,alpha=v2m1,beta=kbeta)
	F2p1=pgnorm(q00,mu=v100,alpha=v2p1,beta=kbeta)
	mu1=matrix(0,2,length(alpha))
	mu1[1,]=-(F1p1-F1m1)/(2*d1)
	mu1[2,]=-(F2p1-F2m1)/(2*d2)
	return(mu1)
}
#' DMGS equation 3.3, f2 term
#' @inherit man2f return
#' @inheritParams manf
gnorm_k3_f2f=function(y,v1,d1,v2,fd2,kbeta){
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
	F1m2=dgnorm(y,mu=v1m2,alpha=v200,beta=kbeta)
	F1m1=dgnorm(y,mu=v1m1,alpha=v200,beta=kbeta)
	F100=dgnorm(y,mu=v100,alpha=v200,beta=kbeta)
	F1p1=dgnorm(y,mu=v1p1,alpha=v200,beta=kbeta)
	F1p2=dgnorm(y,mu=v1p2,alpha=v200,beta=kbeta)
# v2 derivative
	F2m2=dgnorm(y,mu=v100,alpha=v2m2,beta=kbeta)
	F2m1=dgnorm(y,mu=v100,alpha=v2m1,beta=kbeta)
	F200=dgnorm(y,mu=v100,alpha=v200,beta=kbeta)
	F2p1=dgnorm(y,mu=v100,alpha=v2p1,beta=kbeta)
	F2p2=dgnorm(y,mu=v100,alpha=v2p2,beta=kbeta)
# cross derivative
	Fcm1m1=dgnorm(y,mu=v1m1,alpha=v2m1,beta=kbeta)
	Fcm1p1=dgnorm(y,mu=v1m1,alpha=v2p1,beta=kbeta)
	Fcp1m1=dgnorm(y,mu=v1p1,alpha=v2m1,beta=kbeta)
	Fcp1p1=dgnorm(y,mu=v1p1,alpha=v2p1,beta=kbeta)
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
gnorm_k3_p2f=function(y,v1,d1,v2,fd2,kbeta){
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
	F1m2=pgnorm(y,mu=v1m2,alpha=v200,beta=kbeta)
	F1m1=pgnorm(y,mu=v1m1,alpha=v200,beta=kbeta)
	F100=pgnorm(y,mu=v100,alpha=v200,beta=kbeta)
	F1p1=pgnorm(y,mu=v1p1,alpha=v200,beta=kbeta)
	F1p2=pgnorm(y,mu=v1p2,alpha=v200,beta=kbeta)
# v2 derivative
	F2m2=pgnorm(y,mu=v100,alpha=v2m2,beta=kbeta)
	F2m1=pgnorm(y,mu=v100,alpha=v2m1,beta=kbeta)
	F200=pgnorm(y,mu=v100,alpha=v200,beta=kbeta)
	F2p1=pgnorm(y,mu=v100,alpha=v2p1,beta=kbeta)
	F2p2=pgnorm(y,mu=v100,alpha=v2p2,beta=kbeta)
# cross derivative
	Fcm1m1=pgnorm(y,mu=v1m1,alpha=v2m1,beta=kbeta)
	Fcm1p1=pgnorm(y,mu=v1m1,alpha=v2p1,beta=kbeta)
	Fcp1m1=pgnorm(y,mu=v1p1,alpha=v2m1,beta=kbeta)
	Fcp1p1=pgnorm(y,mu=v1p1,alpha=v2p1,beta=kbeta)
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
gnorm_k3_mu2f=function(alpha,v1,d1,v2,fd2,kbeta){
	q00=qgnorm((1-alpha),mu=v1,alpha=v2,beta=kbeta)
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
	F1m2=pgnorm(q00,mu=v1m2,alpha=v200,beta=kbeta)
	F1m1=pgnorm(q00,mu=v1m1,alpha=v200,beta=kbeta)
	F100=pgnorm(q00,mu=v100,alpha=v200,beta=kbeta)
	F1p1=pgnorm(q00,mu=v1p1,alpha=v200,beta=kbeta)
	F1p2=pgnorm(q00,mu=v1p2,alpha=v200,beta=kbeta)
# v2 derivative
	F2m2=pgnorm(q00,mu=v100,alpha=v2m2,beta=kbeta)
	F2m1=pgnorm(q00,mu=v100,alpha=v2m1,beta=kbeta)
	F200=pgnorm(q00,mu=v100,alpha=v200,beta=kbeta)
	F2p1=pgnorm(q00,mu=v100,alpha=v2p1,beta=kbeta)
	F2p2=pgnorm(q00,mu=v100,alpha=v2p2,beta=kbeta)
# cross derivative
	Fcm1m1=pgnorm(q00,mu=v1m1,alpha=v2m1,beta=kbeta)
	Fcm1p1=pgnorm(q00,mu=v1m1,alpha=v2p1,beta=kbeta)
	Fcp1m1=pgnorm(q00,mu=v1p1,alpha=v2m1,beta=kbeta)
	Fcp1p1=pgnorm(q00,mu=v1p1,alpha=v2p1,beta=kbeta)
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
gnorm_k3_logscores=function(logscores,x,d1=0.01,fd2=0.01,kbeta,aderivs){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){
			x1=x[-i]
			dd=dgnorm_k3sub(x1,x[i],d1,fd2,kbeta,aderivs)

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
dgnorm_k3sub=function(x,y,d1=0.01,fd2=0.01,kbeta,aderivs=TRUE){

		nx=length(x)

		v1start=mean(x)
		v2start=sd(x)
		opt=optim(c(v1start,v2start),gnorm_k3_loglik,x=x,kbeta=kbeta,control=list(fnscale=-1))
		v1hat=opt$par[1]
		v2hat=opt$par[2]
		ml_params=c(v1hat,v2hat)

# mle
		ml_pdf=dgnorm(y,mu=v1hat,alpha=v2hat,beta=kbeta)
		ml_cdf=pgnorm(y,mu=v1hat,alpha=v2hat,beta=kbeta)

# rhp
		if(aderivs)	ldd=gnorm_k3_ldda(x,v1hat,v2hat,kbeta)
		if(!aderivs)ldd=gnorm_k3_ldd(x,v1hat,d1,v2hat,fd2,kbeta)
		lddi=solve(ldd)
		if(aderivs)	lddd=gnorm_k3_lddda(x,v1hat,v2hat,kbeta)
		if(!aderivs)lddd=gnorm_k3_lddd(x,v1hat,d1,v2hat,fd2,kbeta)

		if(aderivs) f1=gnorm_k3_f1fa(y,v1hat,v2hat,kbeta)
		if(!aderivs)f1=gnorm_k3_f1f(y,v1hat,d1,v2hat,fd2,kbeta)

		if(aderivs) f2=gnorm_k3_f2fa(y,v1hat,v2hat,kbeta)
		if(!aderivs)f2=gnorm_k3_f2f(y,v1hat,d1,v2hat,fd2,kbeta)

# no aderiv for the cdf in gnorm
		p1=gnorm_k3_p1f(y,v1hat,d1,v2hat,fd2,kbeta)
		p2=gnorm_k3_p2f(y,v1hat,d1,v2hat,fd2,kbeta)
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



