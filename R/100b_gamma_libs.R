# I use v1-=shape,v2=scale
#' Waic
#' @inherit manwaic return
#' @inheritParams manf
gamma_waic=function(waicscores,x,v1hat,fd1,v2hat,fd2,lddi,lddd,lambdad,aderivs){
		if(waicscores){

			if(aderivs) f1f=gamma_f1fa(x,v1hat,v2hat)
			if(!aderivs)f1f=gamma_f1f(x,v1hat,fd1,v2hat,fd2)

			if(aderivs) f2f=gamma_f2fa(x,v1hat,v2hat)
			if(!aderivs)f2f=gamma_f2f(x,v1hat,fd1,v2hat,fd2)

			fhatx=dgamma(x,shape=v1hat,scale=v2hat)
			waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad,f2f,dim=2)
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
gamma_logf=function(params,x){
	sh=pmax(params[1],.Machine$double.eps)
	sc=pmax(params[2],.Machine$double.eps)
	logf=sum(dgamma(x,shape=sh,scale=sc,log=TRUE))-log(sh)-log(sc)
	return(logf)
}#'  log-likelihood function
#' @inherit manloglik return
#' @inheritParams manf
gamma_loglik=function(vv,x){
	n=length(x)
	loglik=sum(dgamma(x,shape=max(vv[1],.Machine$double.eps),scale=max(vv[2],.Machine$double.eps),log=TRUE))
	return(loglik)
}
#' One component of the second derivative of the normalized log-likelihood
#' @inherit manlnn return
#' @inheritParams manf
gamma_lmn=function(x,v1,fd1,v2,fd2,mm,nn){
	d1=fd1*v1
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
			lmn[i]=sum(dgamma(x,shape=vvd[1],scale=vvd[2],log=TRUE))/nx
		}
		dld=(lmn[1]-lmn[2]-lmn[3]+lmn[4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:3){
			for (j in 1:2){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[i]=sum(dgamma(x,shape=vvd[1],scale=vvd[2],log=TRUE))/nx
		}
		dld=(lmn[1]-2*lmn[2]+lmn[3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' Second derivative matrix of the normalized log-likelihood
#' @inherit manldd return
#' @inheritParams manf
gamma_ldd=function(x,v1,fd1,v2,fd2){
	ldd=matrix(0,2,2)
	for (i in 1:2){
		for (j in i:2){
			ldd[i,j]=gamma_lmn(x,v1,fd1,v2,fd2,i,j)
		}
	}
	for (i in 2:1){
		for (j in 1:(i-1)){
			ldd[i,j]=ldd[j,i]
		}
	}
	return(ldd)
}
#' One component of the second derivative of the expected log-likelihood
#' @inherit manlnn return
#' @inheritParams manf
gamma_gmn=function(alpha,v1,fd1,v2,fd2,mm,nn){
	nx=length(alpha)
	d1=fd1*v1
	d2=fd2*v2
  x=qgamma((1-alpha),shape=v1,scale=v2)
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
			lmn[,i]=dgamma(x,shape=vvd[1],scale=vvd[2],log=TRUE)
		}
		dld=(lmn[,1]-lmn[,2]-lmn[,3]+lmn[,4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:3){
			for (j in 1:2){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[,i]=dgamma(x,shape=vvd[1],scale=vvd[2],log=TRUE)
		}
		dld=(lmn[,1]-2*lmn[,2]+lmn[,3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' Second derivative matrix of the expected log-likelihood
#' @inherit manldd return
#' @inheritParams manf
gamma_gg=function(v1,fd1,v2,fd2){
	expinfmat=matrix(0,2,2)
	for (i in 1:2){
		for (j in i:2){
			expinfmat[i,j]=-quad(gamma_gmn,xa=0,xb=1,v1=v1,fd1=fd1,v2=v2,fd2=fd2,mm=i,nn=j)
		}
	}
	for (i in 2:2){
		for (j in 1:(i-1)){
			expinfmat[i,j]=expinfmat[j,i]
		}
	}
 return(expinfmat)
}
#' One component of the second derivative of the normalized log-likelihood
#' @inherit manlnnn return
#' @inheritParams manf
gamma_lmnp=function(x,v1,fd1,v2,fd2,mm,nn,rr){
	d1=fd1*v1
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
			lmn[i]=sum(dgamma(x,shape=vvd[1],scale=vvd[2],log=TRUE))/nx
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
			lmn[i]=sum(dgamma(x,shape=vvd[1],scale=vvd[2],log=TRUE))/nx
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
			lmn[i]=sum(dgamma(x,shape=vvd[1],scale=vvd[2],log=TRUE))/nx
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
gamma_lddd=function(x,v1,fd1,v2,fd2){
# calculate the unique values
	lddd=array(0,c(2,2,2))
	for (i in 1:2){
		for (j in i:2){
			for (k in j:2){
				lddd[i,j,k]=gamma_lmnp(x,v1,fd1,v2,fd2,i,j,k)
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
gamma_f1f=function(y,v1,fd1,v2,fd2){
	d1=fd1*v1
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
	F1m1=dgamma(y,shape=v1m1,scale=v200)
	F1p1=dgamma(y,shape=v1p1,scale=v200)
# v2 derivatives
	F2m1=dgamma(y,shape=v100,scale=v2m1)
	F2p1=dgamma(y,shape=v100,scale=v2p1)
	f1=matrix(0,2,length(y))
	f1[1,]=(F1p1-F1m1)/(2*d1)
	f1[2,]=(F2p1-F2m1)/(2*d2)
	return(f1)
}
#' DMGS equation 3.3, p1 term
#' @inherit man1f return
#' @inheritParams manf
gamma_p1f=function(y,v1,fd1,v2,fd2){
	d1=fd1*v1
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
	F1m1=pgamma(y,shape=v1m1,scale=v200)
	F1p1=pgamma(y,shape=v1p1,scale=v200)
# v2 derivatives
	F2m1=pgamma(y,shape=v100,scale=v2m1)
	F2p1=pgamma(y,shape=v100,scale=v2p1)
	p1=matrix(0,2,length(y))
	p1[1,]=(F1p1-F1m1)/(2*d1)
	p1[2,]=(F2p1-F2m1)/(2*d2)
	return(p1)
}
#' DMGS equation 3.3, mu1 term
#' @inherit man1f return
#' @inheritParams manf
gamma_mu1f=function(alpha,v1,fd1,v2,fd2){
	q00=qgamma((1-alpha),shape=v1,scale=v2)
	d1=fd1*v1
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
	F1m1=pgamma(q00,shape=v1m1,scale=v200)
	F1p1=pgamma(q00,shape=v1p1,scale=v200)
# v2 derivatives
	F2m1=pgamma(q00,shape=v100,scale=v2m1)
	F2p1=pgamma(q00,shape=v100,scale=v2p1)
	mu1=matrix(0,2,length(alpha))
	mu1[1,]=-(F1p1-F1m1)/(2*d1)
	mu1[2,]=-(F2p1-F2m1)/(2*d2)
	return(mu1)
}
#' DMGS equation 3.3, f2 term
#' @inherit man2f return
#' @inheritParams manf
gamma_f2f=function(y,v1,fd1,v2,fd2){
	d1=fd1*v1
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
	F1m2=dgamma(y,shape=v1m2,scale=v200)
	F1m1=dgamma(y,shape=v1m1,scale=v200)
	F100=dgamma(y,shape=v100,scale=v200)
	F1p1=dgamma(y,shape=v1p1,scale=v200)
	F1p2=dgamma(y,shape=v1p2,scale=v200)
# v2 derivative
	F2m2=dgamma(y,shape=v100,scale=v2m2)
	F2m1=dgamma(y,shape=v100,scale=v2m1)
	F200=dgamma(y,shape=v100,scale=v200)
	F2p1=dgamma(y,shape=v100,scale=v2p1)
	F2p2=dgamma(y,shape=v100,scale=v2p2)
# cross derivative
	Fcm1m1=dgamma(y,shape=v1m1,scale=v2m1)
	Fcm1p1=dgamma(y,shape=v1m1,scale=v2p1)
	Fcp1m1=dgamma(y,shape=v1p1,scale=v2m1)
	Fcp1p1=dgamma(y,shape=v1p1,scale=v2p1)
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
gamma_p2f=function(y,v1,fd1,v2,fd2){
	d1=fd1*v1
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
	F1m2=pgamma(y,shape=v1m2,scale=v200)
	F1m1=pgamma(y,shape=v1m1,scale=v200)
	F100=pgamma(y,shape=v100,scale=v200)
	F1p1=pgamma(y,shape=v1p1,scale=v200)
	F1p2=pgamma(y,shape=v1p2,scale=v200)
# v2 derivative
	F2m2=pgamma(y,shape=v100,scale=v2m2)
	F2m1=pgamma(y,shape=v100,scale=v2m1)
	F200=pgamma(y,shape=v100,scale=v200)
	F2p1=pgamma(y,shape=v100,scale=v2p1)
	F2p2=pgamma(y,shape=v100,scale=v2p2)
# cross derivative
	Fcm1m1=pgamma(y,shape=v1m1,scale=v2m1)
	Fcm1p1=pgamma(y,shape=v1m1,scale=v2p1)
	Fcp1m1=pgamma(y,shape=v1p1,scale=v2m1)
	Fcp1p1=pgamma(y,shape=v1p1,scale=v2p1)
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
gamma_mu2f=function(alpha,v1,fd1,v2,fd2){
	q00=qgamma((1-alpha),shape=v1,scale=v2)
	d1=fd1*v1
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
	F1m2=pgamma(q00,shape=v1m2,scale=v200)
	F1m1=pgamma(q00,shape=v1m1,scale=v200)
	F100=pgamma(q00,shape=v100,scale=v200)
	F1p1=pgamma(q00,shape=v1p1,scale=v200)
	F1p2=pgamma(q00,shape=v1p2,scale=v200)
# v2 derivative
	F2m2=pgamma(q00,shape=v100,scale=v2m2)
	F2m1=pgamma(q00,shape=v100,scale=v2m1)
	F200=pgamma(q00,shape=v100,scale=v200)
	F2p1=pgamma(q00,shape=v100,scale=v2p1)
	F2p2=pgamma(q00,shape=v100,scale=v2p2)
# cross derivative
	Fcm1m1=pgamma(q00,shape=v1m1,scale=v2m1)
	Fcm1p1=pgamma(q00,shape=v1m1,scale=v2p1)
	Fcp1m1=pgamma(q00,shape=v1p1,scale=v2m1)
	Fcp1p1=pgamma(q00,shape=v1p1,scale=v2p1)
	mu2[1,1,]=-(F1p1-2*F100+F1m1)/(d1*d1)
	mu2[2,2,]=-(F2p1-2*F200+F2m1)/(d2*d2)
	mu2[1,2,]=-(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2)
	# copy
	mu2[2,1,]=mu2[1,2,]
	return(mu2)
}
#' MLE and RHP predictive means
#' @inherit manmeans return
#' @inheritParams manf
gamma_means=function(means,ml_params,lddi,lddd,lambdad_cp,nx,dim=2){
# v1 is shape
# v2 is scale
	if(means){
		v1=ml_params[1]
		v2=ml_params[2]

# ml mean
		ml_mean=v1*v2

# cp mean

		meand1=array(0,c(2,1))
		meand1[1,1]=v2
		meand1[2,1]=v1

		meand2=array(0,c(2,2,1))
		meand2[1,1,1]=0
		meand2[1,2,1]=1
		meand2[2,1,1]=1
		meand2[2,2,1]=0

		dmean=dmgs(lddi,lddd,meand1,lambdad_cp,meand2,dim=2)
		cp_mean=ml_mean+dmean/nx
	}else{
		ml_mean="means not selected"
		cp_mean="means not selected"
	}


	list(ml_mean=ml_mean,cp_mean=cp_mean)

}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inherit manlogscores return
#' @inheritParams manf
gamma_logscores=function(logscores,x,fd1=0.01,fd2=0.01,aderivs=TRUE){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		cp_oos_logscore=0
		for (i in 1:nx){
			x1=x[-i]

			dd=dgammasub(x1,x[i],fd1,fd2,aderivs)

			ml_params=dd$ml_params

			ml_pdf=dd$ml_pdf
			ml_oos_logscore=ml_oos_logscore+log(ml_pdf)

			cp_pdf=dd$cp_pdf
			cp_oos_logscore=cp_oos_logscore+log(max(cp_pdf,.Machine$double.eps))
		}
	}else{
		ml_oos_logscore="logscores not selected"
		cp_oos_logscore="logscores not selected"
	}
	list(ml_oos_logscore=ml_oos_logscore,cp_oos_logscore=cp_oos_logscore)
}
#' Densities from MLE and RHP
#' @inherit mandsub return
#' @inheritParams manf
dgammasub=function(x,y,fd1=0.01,fd2=0.01,aderivs=TRUE){

		nx=length(x)

		opt=optim(c(1,1),gamma_loglik,x=x,control=list(fnscale=-1))
		v1hat=opt$par[1]
		v2hat=opt$par[2]
		ml_params=c(v1hat,v2hat)

# mle
		ml_pdf=dgamma(y,shape=v1hat,scale=v2hat)
		ml_cdf=pgamma(y,shape=v1hat,scale=v2hat)

# cp
		if(aderivs) ldd=gamma_ldda(x,v1hat,v2hat)
		if(!aderivs)ldd=gamma_ldd(x,v1hat,fd1,v2hat,fd2)
		lddi=solve(ldd)

		if(aderivs) lddd=gamma_lddda(x,v1hat,v2hat)
		if(!aderivs)lddd=gamma_lddd(x,v1hat,fd1,v2hat,fd2)

		if(aderivs) f1=gamma_f1fa(y,v1hat,v2hat)
		if(!aderivs)f1=gamma_f1f(y,v1hat,fd1,v2hat,fd2)

		if(aderivs) f2=gamma_f2fa(y,v1hat,v2hat)
		if(!aderivs)f2=gamma_f2f(y,v1hat,fd1,v2hat,fd2)

		p1=gamma_p1f(y,v1hat,fd1,v2hat,fd2)
		p2=gamma_p2f(y,v1hat,fd1,v2hat,fd2)
		lambdad_cp=c(-1/v1hat,-1/v2hat)
		df1=dmgs(lddi,lddd,f1,lambdad_cp,f2,dim=2)
		dp1=dmgs(lddi,lddd,p1,lambdad_cp,p2,dim=2)
		cp_pdf=pmax(ml_pdf+df1/nx,0)
		cp_cdf=pmin(pmax(ml_cdf+dp1/nx,0),1)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					cp_pdf=cp_pdf,
					ml_cdf=ml_cdf,
					cp_cdf=cp_cdf)
}
