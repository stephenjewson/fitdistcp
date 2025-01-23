#' rgpd for gpd_k1 but with maxlik xi within bounds
#' @inheritParams manf
rgpd_k1_minmax=function(nx,kloc,sigma,xi,minxi=-0.45,maxxi=0.45){
	xihat=-999
	while((xihat<minxi)||(xihat>maxxi)){ #0.46 also works...0.47 doesn't
		xx=extraDistr::rgpd(nx,mu=kloc,sigma=sigma,xi=xi)
		ics=gpd_k1_setics(xx,c(0,0))
		opt1=optim(ics,gpd_k1_loglik,x=xx,kloc=kloc,control=list(fnscale=-1))
		xihat=opt1$par[2]
	}
	return(xx)
}
#' Waic
#' @inheritParams manf
gpd_k1_waic=function(waicscores,x,v1hat,fd1,v2hat,d2,kloc,lddi,lddd,lambdad,
	aderivs){
		if(waicscores){

			if(aderivs) f1f=gpd_k1_f1fa(x,v1hat,v2hat,kloc)
			if(!aderivs)f1f=gpd_k1_f1f(x,v1hat,fd1,v2hat,d2,kloc)

			if(aderivs) f2f=gpd_k1_f2fa(x,v1hat,v2hat,kloc)
			if(!aderivs)f2f=gpd_k1_f2f(x,v1hat,fd1,v2hat,d2,kloc)

			fhatx=extraDistr::dgpd(x,mu=kloc,sigma=v1hat,xi=v2hat)
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
#' @inheritParams manf
gpd_k1_logf=function(params,x,kloc){
#	sc=params[1]
#	sh=params[2]
#	if(sc>0){
#		logf=sum(dgpd(x,mu=kloc,sigma=sc,xi=sh,log=TRUE))-log(sc)
#	}else{
#		logf=-Inf
#	}
	sc=pmax(params[1],.Machine$double.eps)
	sh=params[2]
	logf=sum(dgpd(x,mu=kloc,sigma=sc,xi=sh,log=TRUE))-log(sc)
	return(logf)
}#' Set initial conditions
#' @inheritParams manf
gpd_k1_setics=function(x,ics){
	if((ics[1]==0)&&(ics[2]==0)){
		ics[1]=sd(x)
		ics[2]=0
	}
	return(ics)
}
#'  log-likelihood function
#' @inheritParams manf
gpd_k1_loglik=function(vv,x,kloc){
	n=length(x)
	loglik=sum(extraDistr::dgpd(x,mu=kloc,sigma=max(vv[1],.Machine$double.eps),xi=vv[2],log=TRUE))
	return(loglik)
}
#' Check MLE
#' @inheritParams manf
gpd_k1_checkmle=function(ml_params,kloc,minxi,maxxi){
	v1hat=ml_params[1]
	v2hat=ml_params[2]
	if(is.na(v1hat))stop()
	if(is.na(v2hat))stop()
#
# min xi
#
##	minxi=0
	if(v2hat<minxi){cat("\n***v2hat=",v2hat,"=> execution halted because maxlik shape parameter <",minxi,"***\n");stop()}
#
# max xi
#
	if(v2hat>maxxi){cat("\n***v2hat=",v2hat,"=> execution halted because maxlik shape parameter >",maxxi,"***\n");stop()}
# This max value is ad-hoc
# If it's lowered to 1, then the ppm results for xi=0.6 go wrong, which I understand.
# If it's increased to 100, then in about 1 in a billion cases, for nx=25,
# the xi-hat value is very large and the code crashes because lddi can't be calculated.
# I suspect there is more to understand about that latter part, but for now
# this is a compromise.
}
#' One component of the second derivative of the normalized log-likelihood
#' @inheritParams manf
gpd_k13_l11=function(x,v1,fd1,v2,kloc){
  nx=length(x)
	d1=fd1*v1
	v1m1=v1-1*d1
	v100=v1
	v1p1=v1+1*d1
	lm1=sum(log(extraDistr::dgpd(x,mu=kloc,sigma=v1m1,xi=v2)))/nx
	l00=sum(log(extraDistr::dgpd(x,mu=kloc,sigma=v100,xi=v2)))/nx
	lp1=sum(log(extraDistr::dgpd(x,mu=kloc,sigma=v1p1,xi=v2)))/nx
	dld22=(lp1-2*l00+lm1)/(d1*d1)
	return(dld22)
}
#' One component of the second derivative of the normalized log-likelihood
#' @inheritParams manf
gpd_k1_lmn=function(x,v1,fd1,v2,d2,kloc,mm,nn){
	d1=fd1*v1
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
			lmn[i]=sum(dgpd(x,mu=kloc,sigma=vvd[1],xi=vvd[2],log=TRUE))/nx
		}
		dld=(lmn[1]-lmn[2]-lmn[3]+lmn[4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:3){
			for (j in 1:2){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[i]=sum(dgpd(x,mu=kloc,sigma=vvd[1],xi=vvd[2],log=TRUE))/nx
		}
		dld=(lmn[1]-2*lmn[2]+lmn[3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' Second derivative matrix of the normalized log-likelihood, with fixed shape
#' @inheritParams manf
gpd_k13_ldd=function(x,v1,fd1,v2,kloc){
	ldd=matrix(0,1,1)
	ldd[1,1]=gpd_k13_l11(x,v1,fd1,v2,kloc)
	return(ldd)
}
#' Second derivative matrix of the normalized log-likelihood
#' @inheritParams manf
gpd_k1_ldd=function(x,v1,fd1,v2,d2,kloc){
	ldd=matrix(0,2,2)
	for (i in 1:2){
		for (j in i:2){
			ldd[i,j]=gpd_k1_lmn(x,v1,fd1,v2,d2,kloc,i,j)
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
#' @inheritParams manf
gpd_k13_l111=function(x,v1,fd1,v2,kloc){
  nx=length(x)
	d1=fd1*v1
	v1m2=v1-2*d1
	v1m1=v1-1*d1
	v1p1=v1+1*d1
	v1p2=v1+2*d1
	lm2=sum(extraDistr::dgpd(x,mu=kloc,sigma=v1m2,xi=v2,log=TRUE))/nx
	lm1=sum(extraDistr::dgpd(x,mu=kloc,sigma=v1m1,xi=v2,log=TRUE))/nx
	lp1=sum(extraDistr::dgpd(x,mu=kloc,sigma=v1p1,xi=v2,log=TRUE))/nx
	lp2=sum(extraDistr::dgpd(x,mu=kloc,sigma=v1p2,xi=v2,log=TRUE))/nx
	dld111=(lp2-2*lp1+2*lm1-lm2)/(2*d1*d1*d1)
	return(dld111)
}
#' Third derivative tensor of the normalized log-likelihood
#' @inheritParams manf
gpd_k13_lddd=function(x,v1,fd1,v2,kloc){
	lddd=array(0,c(1,1,1))
	lddd[1,1,1]=gpd_k13_l111(x,v1,fd1,v2,kloc)
	return(lddd)
}
#' One component of the second derivative of the normalized log-likelihood
#' @inheritParams manf
gpd_k1_lmnp=function(x,v1,fd1,v2,d2,kloc,mm,nn,rr){
	d1=fd1*v1
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
			lmn[i]=sum(dgpd(x,mu=kloc,sigma=vvd[1],xi=vvd[2],log=TRUE))/nx
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
			lmn[i]=sum(dgpd(x,mu=kloc,sigma=vvd[1],xi=vvd[2],log=TRUE))/nx
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
			lmn[i]=sum(dgpd(x,mu=kloc,sigma=vvd[1],xi=vvd[2],log=TRUE))/nx
		}
		dld1=(lmn[3]-2*lmn[2]+lmn[1])/(dd[m2]*dd[m2])
		dld2=(lmn[6]-2*lmn[5]+lmn[4])/(dd[m2]*dd[m2])
		dld=(dld2-dld1)/(2*dd[n2])
	}
	return(dld)
}
#' Third derivative tensor of the normalized log-likelihood
#' @inheritParams manf
gpd_k1_lddd=function(x,v1,fd1,v2,d2,kloc){
	lddd=array(0,c(2,2,2))
	for (i in 1:2){
		for (j in i:2){
			for (k in j:2){
				lddd[i,j,k]=gpd_k1_lmnp(x,v1,fd1,v2,d2,kloc,i,j,k)
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
#' @inheritParams manf
gpd_k13_f1f=function(y,v1,fd1,v2,kloc){
  d1=fd1*v1
	ymax=ifelse(v2<0,(d1-v1)/v2,Inf)
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v1 derivatives
	F1m1=extraDistr::dgpd(y,mu=kloc,sigma=v1m1,xi=v2)
	F1p1=extraDistr::dgpd(y,mu=kloc,sigma=v1p1,xi=v2)
	f1=matrix(0,1,length(y))
	f1[1,]=ifelse(y<ymax,(F1p1-F1m1)/(2*d1),0)
	return(f1)
}
#' DMGS equation 3.3, f1 term
#' @inheritParams manf
gpd_k1_f1f=function(y,v1,fd1,v2,d2,kloc){
  d1=fd1*v1
	ymax=ifelse(v2<0,0-(v1-2*d1)/(v2-2*d2),Inf)
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v2 stuff
	v2m1=v2-1*d2
	v200=v2+0*d2
	v2p1=v2+1*d2
# v1 derivatives
	F1m1=extraDistr::dgpd(y,mu=kloc,sigma=v1m1,xi=v200)
	F1p1=extraDistr::dgpd(y,mu=kloc,sigma=v1p1,xi=v200)
# v2 derivatives
	F2m1=extraDistr::dgpd(y,mu=kloc,sigma=v100,xi=v2m1)
	F2p1=extraDistr::dgpd(y,mu=kloc,sigma=v100,xi=v2p1)
	f1=matrix(0,2,length(y))
	f1[1,]=ifelse(y<ymax,(F1p1-F1m1)/(2*d1),0)
	f1[2,]=ifelse(y<ymax,(F2p1-F2m1)/(2*d2),0)
	return(f1)
}
#' DMGS equation 3.3, p1 term
#' @inheritParams manf
gpd_k13_p1f=function(y,v1,fd1,v2,kloc){
  d1=fd1*v1
	ymax=ifelse(v2<0,0-(v1-d1)/v2,Inf)
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v1 derivatives
	F1m1=extraDistr::pgpd(y,mu=kloc,sigma=v1m1,xi=v2)
	F1p1=extraDistr::pgpd(y,mu=kloc,sigma=v1p1,xi=v2)
	p1=matrix(0,1,length(y))
	p1[1,]=ifelse(y<ymax,(F1p1-F1m1)/(2*d1),0)
	return(p1)
}
#' DMGS equation 3.3, p1 term
#' @inheritParams manf
gpd_k1_p1f=function(y,v1,fd1,v2,d2,kloc){
  d1=fd1*v1
	ymax=ifelse(v2<0,0-(v1-2*d1)/(v2-2*d2),Inf)
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v2 stuff
	v2m1=v2-1*d2
	v200=v2+0*d2
	v2p1=v2+1*d2
# v1 derivatives
	F1m1=extraDistr::pgpd(y,mu=kloc,sigma=v1m1,xi=v200)
	F1p1=extraDistr::pgpd(y,mu=kloc,sigma=v1p1,xi=v200)
# v2 derivatives
	F2m1=extraDistr::pgpd(y,mu=kloc,sigma=v100,xi=v2m1)
	F2p1=extraDistr::pgpd(y,mu=kloc,sigma=v100,xi=v2p1)
	p1=matrix(0,2,length(y))
	p1[1,]=ifelse(y<ymax,(F1p1-F1m1)/(2*d1),0)
	p1[2,]=ifelse(y<ymax,(F2p1-F2m1)/(2*d2),0)
	return(p1)
}
#' DMGS equation 3.3, mu1 term
#' @inheritParams manf
gpd_k13_mu1f=function(alpha,v1,fd1,v2,kloc){
	q00=extraDistr::qgpd((1-alpha),mu=kloc,sigma=v1,xi=v2)
  d1=fd1*v1
	qmax=ifelse(v2<0,0-(v1-d1)/v2,Inf)
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v1 derivatives
	F1m1=extraDistr::pgpd(q00,mu=kloc,sigma=v1m1,xi=v2)
	F1p1=extraDistr::pgpd(q00,mu=kloc,sigma=v1p1,xi=v2)
	mu1=matrix(0,1,length(alpha))
	mu1[1,]=-(F1p1-F1m1)/(2*d1)
	return(mu1)
}
#' DMGS equation 3.3, mu1 term
#' @inheritParams manf
gpd_k1_mu1f=function(alpha,v1,fd1,v2,d2,kloc){
	q00=extraDistr::qgpd((1-alpha),mu=kloc,sigma=v1,xi=v2)
  d1=fd1*v1
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v2 stuff
	v2m1=v2-1*d2
	v200=v2+0*d2
	v2p1=v2+1*d2
# v1 derivatives
	F1m1=extraDistr::pgpd(q00,mu=kloc,sigma=v1m1,xi=v200)
	F1p1=extraDistr::pgpd(q00,mu=kloc,sigma=v1p1,xi=v200)
# v2 derivatives
	F2m1=extraDistr::pgpd(q00,mu=kloc,sigma=v100,xi=v2m1)
	F2p1=extraDistr::pgpd(q00,mu=kloc,sigma=v100,xi=v2p1)
	mu1=matrix(0,2,length(alpha))
	mu1[1,]=-(F1p1-F1m1)/(2*d1)
	mu1[2,]=-(F2p1-F2m1)/(2*d2)
	return(mu1)
}
#' DMGS equation 3.3, f2 term
#' @inheritParams manf
gpd_k13_f2f=function(y,v1,fd1,v2,kloc){
  d1=fd1*v1
	ymax=ifelse(v2<0,0-(v1-d1)/v2,Inf)
# v1 stuff
	v1m2=v1-2*d1
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
	v1p2=v1+2*d1
	f2=array(0,c(1,1,length(y)))
# v1
	F1m2=extraDistr::dgpd(y,mu=kloc,sigma=v1m2,xi=v2)
	F1m1=extraDistr::dgpd(y,mu=kloc,sigma=v1m1,xi=v2)
	F100=extraDistr::dgpd(y,mu=kloc,sigma=v100,xi=v2)
	F1p1=extraDistr::dgpd(y,mu=kloc,sigma=v1p1,xi=v2)
	F1p2=extraDistr::dgpd(y,mu=kloc,sigma=v1p2,xi=v2)
	f2[1,1,]=ifelse(y<ymax,(F1p1-2*F100+F1m1)/(d1*d1),0)
	return(f2)
}
#' DMGS equation 3.3, f2 term
#' @inheritParams manf
gpd_k1_f2f=function(y,v1,fd1,v2,d2,kloc){
  d1=fd1*v1
	ymax=ifelse(v2<0,0-(v1-2*d1)/(v2-2*d2),Inf)
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
	F1m2=extraDistr::dgpd(y,mu=kloc,sigma=v1m2,xi=v2)
	F1m1=extraDistr::dgpd(y,mu=kloc,sigma=v1m1,xi=v2)
	F100=extraDistr::dgpd(y,mu=kloc,sigma=v100,xi=v2)
	F1p1=extraDistr::dgpd(y,mu=kloc,sigma=v1p1,xi=v2)
	F1p2=extraDistr::dgpd(y,mu=kloc,sigma=v1p2,xi=v2)
	f2[1,1,]=ifelse(y<ymax,(F1p1-2*F100+F1m1)/(d1*d1),0)
# v2 derivative
	F2m2=extraDistr::dgpd(y,mu=kloc,sigma=v1,xi=v2m2)
	F2m1=extraDistr::dgpd(y,mu=kloc,sigma=v1,xi=v2m1)
	F200=extraDistr::dgpd(y,mu=kloc,sigma=v1,xi=v200)
	F2p1=extraDistr::dgpd(y,mu=kloc,sigma=v1,xi=v2p1)
	F2p2=extraDistr::dgpd(y,mu=kloc,sigma=v1,xi=v2p2)
	f2[2,2,]=ifelse(y<ymax,(F2p1-2*F200+F2m1)/(d2*d2),0)
# cross derivative12
	Fcm1m1=extraDistr::dgpd(y,mu=kloc,sigma=v1m1,xi=v2m1)
	Fcm1p1=extraDistr::dgpd(y,mu=kloc,sigma=v1m1,xi=v2p1)
	Fcp1m1=extraDistr::dgpd(y,mu=kloc,sigma=v1p1,xi=v2m1)
	Fcp1p1=extraDistr::dgpd(y,mu=kloc,sigma=v1p1,xi=v2p1)
	f2[1,2,]=ifelse(y<ymax,(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2),0)
	f2[2,1,]=f2[1,2,]
	return(f2)
}
#' DMGS equation 3.3, p2 term
#' @inheritParams manf
gpd_k13_p2f=function(y,v1,fd1,v2,kloc){
  d1=fd1*v1
	ymax=ifelse(v2<0,0-(v1-d1)/v2,Inf)
# v1 stuff
	v1m2=v1-2*d1
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
	v1p2=v1+2*d1
	p2=array(0,c(1,1,length(y)))
# v1
	F1m2=extraDistr::pgpd(y,mu=kloc,sigma=v1m2,xi=v2)
	F1m1=extraDistr::pgpd(y,mu=kloc,sigma=v1m1,xi=v2)
	F100=extraDistr::pgpd(y,mu=kloc,sigma=v100,xi=v2)
	F1p1=extraDistr::pgpd(y,mu=kloc,sigma=v1p1,xi=v2)
	F1p2=extraDistr::pgpd(y,mu=kloc,sigma=v1p2,xi=v2)
	p2[1,1,]=ifelse(y<ymax,(F1p1-2*F100+F1m1)/(d1*d1),0)
	return(p2)
}
#' DMGS equation 3.3, p2 term
#' @inheritParams manf
gpd_k1_p2f=function(y,v1,fd1,v2,d2,kloc){
  d1=fd1*v1
	ymax=ifelse(v2<0,0-(v1-2*d1)/(v2-2*d2),Inf)
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
	F1m2=extraDistr::pgpd(y,mu=kloc,sigma=v1m2,xi=v2)
	F1m1=extraDistr::pgpd(y,mu=kloc,sigma=v1m1,xi=v2)
	F100=extraDistr::pgpd(y,mu=kloc,sigma=v100,xi=v2)
	F1p1=extraDistr::pgpd(y,mu=kloc,sigma=v1p1,xi=v2)
	F1p2=extraDistr::pgpd(y,mu=kloc,sigma=v1p2,xi=v2)
	p2[1,1,]=ifelse(y<ymax,(F1p1-2*F100+F1m1)/(d1*d1),0)
# v2 derivative
	F2m2=extraDistr::pgpd(y,mu=kloc,sigma=v1,xi=v2m2)
	F2m1=extraDistr::pgpd(y,mu=kloc,sigma=v1,xi=v2m1)
	F200=extraDistr::pgpd(y,mu=kloc,sigma=v1,xi=v200)
	F2p1=extraDistr::pgpd(y,mu=kloc,sigma=v1,xi=v2p1)
	F2p2=extraDistr::pgpd(y,mu=kloc,sigma=v1,xi=v2p2)
	p2[2,2,]=ifelse(y<ymax,(F2p1-2*F200+F2m1)/(d2*d2),0)
# cross derivative12
	Fcm1m1=extraDistr::pgpd(y,mu=kloc,sigma=v1m1,xi=v2m1)
	Fcm1p1=extraDistr::pgpd(y,mu=kloc,sigma=v1m1,xi=v2p1)
	Fcp1m1=extraDistr::pgpd(y,mu=kloc,sigma=v1p1,xi=v2m1)
	Fcp1p1=extraDistr::pgpd(y,mu=kloc,sigma=v1p1,xi=v2p1)
	p2[1,2,]=ifelse(y<ymax,(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2),0)
	p2[2,1,]=p2[1,2,]
	return(p2)
}
#' DMGS equation 3.3, mu2 term
#' @inheritParams manf
gpd_k13_mu2f=function(alpha,v1,fd1,v2,kloc){
	q00=extraDistr::qgpd((1-alpha),mu=kloc,sigma=v1,xi=v2)
  d1=fd1*v1
# v1 stuff
	v1m2=v1-2*d1
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
	v1p2=v1+2*d1
	mu2=array(0,c(1,1,length(alpha)))
# v1
	F1m2=extraDistr::pgpd(q00,mu=kloc,sigma=v1m2,xi=v2)
	F1m1=extraDistr::pgpd(q00,mu=kloc,sigma=v1m1,xi=v2)
	F100=extraDistr::pgpd(q00,mu=kloc,sigma=v100,xi=v2)
	F1p1=extraDistr::pgpd(q00,mu=kloc,sigma=v1p1,xi=v2)
	F1p2=extraDistr::pgpd(q00,mu=kloc,sigma=v1p2,xi=v2)
	mu2[1,1,]=-(F1p1-2*F100+F1m1)/(d1*d1)
	return(mu2)
}
#' DMGS equation 3.3, mu2 term
#' @inheritParams manf
gpd_k1_mu2f=function(alpha,v1,fd1,v2,d2,kloc){
	q00=extraDistr::qgpd((1-alpha),mu=kloc,sigma=v1,xi=v2)
  d1=fd1*v1
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
	F1m2=extraDistr::pgpd(q00,mu=kloc,sigma=v1m2,xi=v2)
	F1m1=extraDistr::pgpd(q00,mu=kloc,sigma=v1m1,xi=v2)
	F100=extraDistr::pgpd(q00,mu=kloc,sigma=v100,xi=v2)
	F1p1=extraDistr::pgpd(q00,mu=kloc,sigma=v1p1,xi=v2)
	F1p2=extraDistr::pgpd(q00,mu=kloc,sigma=v1p2,xi=v2)
	mu2[1,1,]=-(F1p1-2*F100+F1m1)/(d1*d1)
# v2 derivative
	F2m2=extraDistr::pgpd(q00,mu=kloc,sigma=v1,xi=v2m2)
	F2m1=extraDistr::pgpd(q00,mu=kloc,sigma=v1,xi=v2m1)
	F200=extraDistr::pgpd(q00,mu=kloc,sigma=v1,xi=v200)
	F2p1=extraDistr::pgpd(q00,mu=kloc,sigma=v1,xi=v2p1)
	F2p2=extraDistr::pgpd(q00,mu=kloc,sigma=v1,xi=v2p2)
	mu2[2,2,]=-(F2p1-2*F200+F2m1)/(d2*d2)
# cross derivative12
	Fcm1m1=extraDistr::pgpd(q00,mu=kloc,sigma=v1m1,xi=v2m1)
	Fcm1p1=extraDistr::pgpd(q00,mu=kloc,sigma=v1m1,xi=v2p1)
	Fcp1m1=extraDistr::pgpd(q00,mu=kloc,sigma=v1p1,xi=v2m1)
	Fcp1p1=extraDistr::pgpd(q00,mu=kloc,sigma=v1p1,xi=v2p1)
	mu2[1,2,]=-(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2)
	mu2[2,1,]=mu2[1,2,]
	return(mu2)
}
#' Derivative of expected information matrix, based on MEV routine gpd.infomat
#' @inheritParams manf
gpd_k1_ggd_mev=function(v1,fd1,v2,d2,kloc){
	ggd=array(0,c(2,2,2))
  d1=fd1*v1

	ggm1=gpd.infomat(c(v1-d1,v2),dat=c(1),method=c("exp"))
	ggp1=gpd.infomat(c(v1+d1,v2),dat=c(1),method=c("exp"))
	ggd[1,,]=(ggp1-ggm1)/(2*d1)

	ggm2=gpd.infomat(c(v1,v2-d2),dat=c(1),method=c("exp"))
	ggp2=gpd.infomat(c(v1,v2+d2),dat=c(1),method=c("exp"))
	ggd[2,,]=(ggp2-ggm2)/(2*d2)

  return(ggd)
}
#' Analytical Expressions for Predictive Means
#' RHP mean based on the expectation of DMGS equation 2.1
#' @inheritParams manf
gpd_k1_means=function(means,ml_params,lddi,lddi_k2,lddd,lddd_k2,
									lambdad_flat,lambdad_rh_mle,lambdad_rh_flat,lambdad_jp,
									nx,dim=2,kloc=0){

	if(means){
# intro
		eulerconstant=0.57721566490153286060651209008240243104215933593992
		mu=kloc
		sigma=ml_params[1]
		xi=ml_params[2]

# set up derivative arrays
		meand1=array(0,c(2,1))
		meand2=array(0,c(2,2,1))
		meand1_k2=array(0,c(1,1))
		meand2_k2=array(0,c(1,1,1))

# mle
		onemxi=1-xi
		onemxisq=onemxi*onemxi
		onemxicu=onemxi*onemxi*onemxi
		ml_mean=kloc+sigma/onemxi
# calculate first derivative array for bayesian xi!=0 cases
		meand1[1,1]=1/onemxi
		meand1[2,1]=sigma/onemxisq
# and copy for the 2d case
		meand1_k2[1,1]=meand1[1,1]
# calculate second derivative array (only has 1 non-zero term!)
		meand2[1,1,1]=0
		meand2[1,2,1]=1/onemxisq
		meand2[2,2,1]=-2*sigma/onemxicu
		meand2_k2[1,1,1]=0

# I've realized now that when I integrate over xi, the mean is Inf

#	flat_mean			=ml_mean+dmgs(lddi,lddd,meand1,lambdad_flat,meand2,dim=2)/nx
		flat_mean			=Inf

		rh_ml_mean		=ml_mean+dmgs(lddi_k2,lddd_k2,meand1_k2,lambdad_rh_mle,meand2_k2,dim=1)/nx

#	rh_flat_mean	=ml_mean+dmgs(lddi,lddd,meand1,lambdad_rh_flat,meand2,dim=2)/nx
		rh_flat_mean	=Inf

#	jp_mean				=ml_mean+dmgs(lddi,lddd,meand1,lambdad_jp,meand2,dim=2)/nx
		jp_mean				=Inf

	}else{
		flat_mean="means not selected"
		ml_mean="means not selected"
		rh_ml_mean="means not selected"
		rh_flat_mean="means not selected"
		jp_mean="means not selected"
	}

# return
	list(ml_mean=ml_mean,flat_mean=flat_mean,rh_ml_mean=rh_ml_mean,rh_flat_mean=rh_flat_mean,jp_mean=jp_mean)
}
#' Densities for 5 predictions
#' @inheritParams manf
dgpdsub=function(x,y,ics,fd1=0.01,d2=0.01,kloc=0,dlogpi=0,
	minxi,maxxi,extramodels=FALSE,aderivs=TRUE){

		nx=length(x)

		ics=gpd_k1_setics(x,ics)
		opt=optim(ics,gpd_k1_loglik,x=x,kloc=kloc,control=list(fnscale=-1))
		v1hat=opt$par[1]
		v2hat=min(maxxi,max(minxi,opt$par[2])) #just reset in this case
#		v2hat=opt$par[2]
		ml_params=c(v1hat,v2hat)
#		cat("ml_params=",ml_params,"\n")

		y=fixgpdrange(y,kloc,v1hat,v2hat)


# mle
		ml_pdf=extraDistr::dgpd(y,mu=kloc,sigma=v1hat,xi=v2hat)
		ml_cdf=extraDistr::pgpd(y,mu=kloc,sigma=v1hat,xi=v2hat)

# return
		list(
					ml_params=ml_params,
					ml_pdf=ml_pdf,
					ml_cdf=ml_cdf)
}

