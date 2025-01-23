#' rgev for gev_p1 but with maxlik xi within bounds
#' @inheritParams manf
rgev_p1_minmax=function(nx,mu,sigma,xi,tt,minxi=-0.45,maxxi=0.45,centering=TRUE){
	xihat=-999
  if(centering)tt=tt-mean(tt)
	while((xihat<minxi)||(xihat>maxxi)){ #0.46 also works...0.47 doesn't
		xx=extraDistr::rgev(nx,mu=mu,sigma=sigma,xi=xi)
		ics=gev_p1_setics(xx,tt,c(0,0,0,0))
		opt1=optim(ics,gev_p1_loglik,x=xx,t=tt,control=list(fnscale=-1))
		xihat=opt1$par[4]
	}
#	cat("1: xihat=",xihat,":")
#	ics=gev_p1_setics(xx,tt,c(0,0,0,0))
#	opt1=optim(ics,gev_p1_loglik,x=xx,t=tt,control=list(fnscale=-1))
#	xihat=opt1$par[4]
#	cat(xihat,":")
	return(xx)
}
#' Waic
#' @inheritParams manf
gev_p1_waic=function(waicscores,x,t,v1hat,d1,v2hat,d2,v3hat,fd3,v4hat,d4,
	lddi,lddd,lambdad,aderivs){
		if(waicscores){
			if(aderivs){
				f1f=gev_p1_f1fa(x,t,v1hat,v2hat,v3hat,v4hat)
				f2f=gev_p1_f2fa(x,t,v1hat,v2hat,v3hat,v4hat)
			} else {
				f1f=gev_p1_f1f(x,t,v1hat,d1,v2hat,d2,v3hat,fd3,v4hat,d4)
				f2f=gev_p1_f2f(x,t,v1hat,d1,v2hat,d2,v3hat,fd3,v4hat,d4)
			}
			fhatx=dgev_p1(x,t,ymn=v1hat,slope=v2hat,sigma=v3hat,xi=v4hat,log=FALSE)
			waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad,f2f,dim=4)
			waic1=waic$waic1
			waic2=waic$waic2
		}else{
			waic1="waicscores not selected"
			waic2="waicscores not selected"
		}
		list(waic1=waic1,waic2=waic2)
}
#' Predicted Parameter and Generalized Residuals
#' @inheritParams manf
gev_p1_predictordata=function(predictordata,x,t,t0,params){
	if(predictordata){
#
# calculate the probabilities of the data using the fited model
#
		a=params[1]
		b=params[2]
		sc=params[3]
		sh=params[4]
		mu=a+b*t
		px=extraDistr::pgev(x,mu=mu,sigma=sc,xi=sh)
#
# calculate the quantiles for those probabilities at t0
#
		mu0=a+b*t0
		qx=extraDistr::qgev(px,mu=mu0,sigma=sc,xi=sh)
	} else{
		predictedparameter="predictordata not selected"
		adjustedx="predictordata not selected"
	}

	list(predictedparameter=mu,adjustedx=qx)
}#' Logf for RUST
#' @inheritParams manf
gev_p1_logf=function(params,x,t){
#	a=params[1]
#	b=params[2]
#	sc=params[3]
#	sh=params[4]
#	mu=a+b*t
#	if(sc>0){
#		logf=sum(extraDistr::dgev(x,mu=mu,sigma=sc,xi=sh,log=TRUE))-log(sc)
#	}else{
#		logf=-Inf
#	}
	a=params[1]
	b=params[2]
	sc=pmax(params[3],sqrt(.Machine$double.eps))
	sh=params[4]
	mu=a+b*t
	logf=sum(extraDistr::dgev(x,mu=mu,sigma=sc,xi=sh,log=TRUE))-log(sc)
	return(logf)
}#' Set initial conditions
#' @inheritParams manf
gev_p1_setics=function(x,t,ics){
	nx=length(x)
	if((ics[1]==0)&&(ics[2]==0)&&(ics[3]==0)&&(ics[4]==0)){
		lm=lm(x~t)
		ics[1]=lm$coefficients[1]
		ics[2]=lm$coefficients[2]
		xhat=ics[1]+ics[2]*t
		ics[3]=sqrt((sum((x-xhat)^2))/nx)
		ics[4]=0
	}
	return(ics)
}
#'  observed log-likelihood function
#' @inheritParams manf
gev_p1_loglik=function(vv,x,t){
	n=length(x)
	mu=vv[1]+vv[2]*t #so mean is a vector, just like x
	loglik=sum(extraDistr::dgev(x,mu=mu,sigma=max(vv[3],.Machine$double.eps),xi=vv[4],log=TRUE))
	return(loglik)
}
#' Check MLE
#' @inheritParams manf
gev_p1_checkmle=function(ml_params,minxi,maxxi){
	v1hat=ml_params[1]
	v2hat=ml_params[2]
	v3hat=ml_params[3]
	v4hat=ml_params[4]
	if(is.na(v1hat))stop()
	if(is.na(v2hat))stop()
	if(is.na(v3hat))stop()
	if(is.na(v4hat))stop()
	if(v4hat<minxi){cat("\n***v4hat=",v4hat,"=> execution halted because maxlik shape parameter <",minxi,"***\n");stop()}
	if(v4hat>maxxi){cat("\n***v4hat=",v4hat,"=> execution halted because maxlik shape parameter >",maxxi,"***\n");stop()}
#	cat("2: v4hat=",v4hat,"\n")
}
#' GEVD-with-p1: Quantile function
#' @inheritParams manf
qgev_p1=function(p,t0,ymn,slope,sigma,xi){

	return(extraDistr::qgev(p,mu=(ymn+slope*t0),sigma=sigma,xi=xi))

}
#' GEVD-with-p1: Density function
#' @inheritParams manf
dgev_p1=function(x,t0,ymn,slope,sigma,xi,log=FALSE){

	mu=(ymn+slope*t0)
	return(extraDistr::dgev(x,mu=mu,sigma=sigma,xi=xi,log=log))

}
#' GEVD-with-p1: Distribution function
#' @inheritParams manf
pgev_p1=function(y,t0,ymn,slope,sigma,xi){

	return(extraDistr::pgev(y,mu=(ymn+slope*t0),sigma=sigma,xi=xi))

}
#' One component of the second derivative of the normalized log-likelihood
#' @inheritParams manf
gev_p1_lmn=function(x,t,v1,d1,v2,d2,v3,fd3,v4,d4,mm,nn){
	d3=fd3*v3
	net3=matrix(0,3,4)
	net4=matrix(0,4,4)
	lmn=matrix(0,4)
	dd=c(d1,d2,d3,d4)
	vv=c(v1,v2,v3,v4)
	vvd=matrix(0,4)
	nx=length(x)
# different
	if(mm!=nn){
		net4[,mm]=c(-1,-1,1,1)
		net4[,nn]=c(-1,1,-1,1)
		for (i in 1:4){
			for (j in 1:4){
				vvd[j]=vv[j]+net4[i,j]*dd[j]
			}
			lmn[i]=sum(dgev_p1(x,t,ymn=vvd[1],slope=vvd[2],sigma=vvd[3],xi=vvd[4],log=TRUE))/nx
		}
		dld=(lmn[1]-lmn[2]-lmn[3]+lmn[4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:3){
			for (j in 1:4){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[i]=sum(dgev_p1(x,t,ymn=vvd[1],slope=vvd[2],sigma=vvd[3],xi=vvd[4],log=TRUE))/nx
		}
		dld=(lmn[1]-2*lmn[2]+lmn[3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' Second derivative matrix of the normalized log-likelihood, with fixed shape parameter
#' @inheritParams manf
gev_p1k3_ldd=function(x,t,v1,d1,v2,d2,v3,fd3,v4){
	d4=999
	nx=length(x)
	ldd=matrix(0,3,3)
# unique terms
	for (i in 1:3){
		for (j in i:3){
			ldd[i,j]=gev_p1_lmn(x,t,v1,d1,v2,d2,v3,fd3,v4,d4,i,j)
		}
	}
# copies
	for (i in 3:2){
		for (j in 1:(i-1)){
			ldd[i,j]=ldd[j,i]
		}
	}
#	cat("new ldd:",ldd,"\n")
	return(ldd)
}
#' Second derivative matrix of the normalized log-likelihood
#' @inheritParams manf
gev_p1_ldd=function(x,t,v1,d1,v2,d2,v3,fd3,v4,d4){
	nx=length(x)
	ldd=matrix(0,4,4)
	for (i in 1:4){
		for (j in i:4){
			ldd[i,j]=gev_p1_lmn(x,t,v1,d1,v2,d2,v3,fd3,v4,d4,i,j)
		}
	}
	for (i in 4:2){
		for (j in 1:(i-1)){
			ldd[i,j]=ldd[j,i]
		}
	}
#	cat("new:",ldd,"\n")
	return(ldd)
}
#' One component of the second derivative of the normalized log-likelihood
#' @inheritParams manf
gev_p1_lmnp=function(x,t,v1,d1,v2,d2,v3,fd3,v4,d4,mm,nn,rr){
	d3=fd3*v3
	net4=matrix(0,4,4)
	net6=matrix(0,6,4)
	net8=matrix(0,8,4)
	lmn=matrix(0,8)
	dd=c(d1,d2,d3,d4)
	vv=c(v1,v2,v3,v4)
	vvd=matrix(0,4)
	nx=length(x)
# all diff
	if ((mm!=nn)&(nn!=rr)&(rr!=mm)){
		net8[,mm]=c(-1,1,-1,1,-1,1,-1,1)
		net8[,nn]=c(-1,-1,1,1,-1,-1,1,1)
		net8[,rr]=c(-1,-1,-1,-1,1,1,1,1)
		for (i in 1:8){
			for (j in 1:4){
				vvd[j]=vv[j]+net8[i,j]*dd[j]
			}
			lmn[i]=sum(dgev_p1(x,t,ymn=vvd[1],slope=vvd[2],sigma=vvd[3],xi=vvd[4],log=TRUE))/nx
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
			for (j in 1:4){
				vvd[j]=vv[j]+net4[i,j]*dd[j]
			}
			lmn[i]=sum(dgev_p1(x,t,ymn=vvd[1],slope=vvd[2],sigma=vvd[3],xi=vvd[4],log=TRUE))/nx
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
			for (j in 1:4){
				vvd[j]=vv[j]+net6[i,j]*dd[j]
			}
			lmn[i]=sum(dgev_p1(x,t,ymn=vvd[1],slope=vvd[2],sigma=vvd[3],xi=vvd[4],log=TRUE))/nx
		}
		dld1=(lmn[3]-2*lmn[2]+lmn[1])/(dd[m2]*dd[m2])
		dld2=(lmn[6]-2*lmn[5]+lmn[4])/(dd[m2]*dd[m2])
		dld=(dld2-dld1)/(2*dd[n2])
	}
	return(dld)
}
#' Third derivative tensor of the normalized log-likelihood, with fixed shape parameter
#' @inheritParams manf
gev_p1k3_lddd=function(x,t,v1,d1,v2,d2,v3,fd3,v4){
	lddd=array(0,c(3,3,3))
	d4=999
# calculate the unique values
	for (i in 1:3){
		for (j in i:3){
			for (k in j:3){
				lddd[i,j,k]=gev_p1_lmnp(x,t,v1,d1,v2,d2,v3,fd3,v4,d4,i,j,k)
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
#	cat("sum(ldddnew)=",sum(lddd),"\n")
	return(lddd)
}
#' Third derivative tensor of the normalized log-likelihood, with fixed shape parameter
#' @inheritParams manf
gev_p1_lddd=function(x,t,v1,d1,v2,d2,v3,fd3,v4,d4){
	lddd=array(0,c(4,4,4))
# calculate the unique values
	for (i in 1:4){
		for (j in i:4){
			for (k in j:4){
				lddd[i,j,k]=gev_p1_lmnp(x,t,v1,d1,v2,d2,v3,fd3,v4,d4,i,j,k)
			}
		}
	}
# steves dumb algorithm for filling in the non-unique values
	for (i in 1:4){
		for (j in 1:4){
			for (k in 1:4){
				a=c(i,j,k)
				b=sort(a)
				lddd[a[1],a[2],a[3]]=lddd[b[1],b[2],b[3]]
			}
		}
	}
#	cat("sum(ldddnew)=",sum(lddd),"\n")
	return(lddd)
}
#' DMGS equation 2.1, f1 term, fixed shape parameter
#' @inheritParams manf
gev_p1k3_f1f=function(y,t0,v1,d1,v2,d2,v3,fd3,v4){
	d3=fd3*v3
	ymax=ifelse(v4<0,v1+v2*t0-2*(d1+d2*t0)-(v3-2*d3)/(v4),Inf)
# new method
	dd=c(d1,d2,d3)
	vv=c(v1,v2,v3)
	f1=matrix(0,3,length(y))
	for (i in 1:3){
		vvm=vv
		vvp=vv
		vvm[i]=vv[i]-dd[i]
		vvp[i]=vv[i]+dd[i]
		Fm1=dgev_p1(y,t0,ymn=vvm[1],slope=vvm[2],sigma=vvm[3],xi=v4)
		Fp1=dgev_p1(y,t0,ymn=vvp[1],slope=vvp[2],sigma=vvp[3],xi=v4)
		f1[i,]=ifelse(y<ymax,(Fp1-Fm1)/(2*dd[i]),0)
	}
#	cat("new f1=",sum(f1),"\n")
	return(f1)
}
#' DMGS equation 2.1, f1 term
#' @inheritParams manf
gev_p1_f1f=function(y,t0,v1,d1,v2,d2,v3,fd3,v4,d4){
	d3=fd3*v3
	ymax=ifelse(v4<0,v1+v2*t0-2*(d1+d2*t0)-(v3-2*d3)/(v4-2*d4),Inf)
# new method
	dd=c(d1,d2,d3,d4)
	vv=c(v1,v2,v3,v4)
	f1=matrix(0,4,length(y))
	for (i in 1:4){
		vvm=vv
		vvp=vv
		vvm[i]=vv[i]-dd[i]
		vvp[i]=vv[i]+dd[i]
		Fm1=dgev_p1(y,t0,ymn=vvm[1],slope=vvm[2],sigma=vvm[3],xi=vvm[4])
		Fp1=dgev_p1(y,t0,ymn=vvp[1],slope=vvp[2],sigma=vvp[3],xi=vvp[4])
		f1[i,]=ifelse(y<ymax,(Fp1-Fm1)/(2*dd[i]),0)
	}
#	cat("new f1=",sum(f1),"\n")
	return(f1)
}
#' GEVD-with-p1: DMGS equation 3.3 mu1 term, fixed shape parameter
#' @inheritParams manf
gev_p1k3_mu1f=function(alpha,t0,v1,d1,v2,d2,v3,fd3,v4){
	d4=999
	q00=qgev_p1((1-alpha),t0,ymn=v1,slope=v2,sigma=v3,xi=v4)
	d3=fd3*v3
# new method
	dd=c(d1,d2,d3)
	vv=c(v1,v2,v3)
	mu1=matrix(0,3,length(alpha))
	for (i in 1:3){
		vvm=vv
		vvp=vv
		vvm[i]=vv[i]-dd[i]
		vvp[i]=vv[i]+dd[i]
		Fm1=pgev_p1(q00,t0,ymn=vvm[1],slope=vvm[2],sigma=vvm[3],xi=v4)
		Fp1=pgev_p1(q00,t0,ymn=vvp[1],slope=vvp[2],sigma=vvp[3],xi=v4)
		mu1[i,]=-(Fp1-Fm1)/(2*dd[i])
	}
#	cat("new mu1=",sum(mu1),"\n")
	return(mu1)
}
#' GEVD-with-p1: DMGS equation 3.3 mu1 term
#' @inheritParams manf
gev_p1_mu1f=function(alpha,t0,v1,d1,v2,d2,v3,fd3,v4,d4){
	q00=qgev_p1((1-alpha),t0,ymn=v1,slope=v2,sigma=v3,xi=v4)
	d3=fd3*v3
# new method
	dd=c(d1,d2,d3,d4)
	vv=c(v1,v2,v3,v4)
	mu1=matrix(0,4,length(alpha))
	for (i in 1:4){
		vvm=vv
		vvp=vv
		vvm[i]=vv[i]-dd[i]
		vvp[i]=vv[i]+dd[i]
		Fm1=pgev_p1(q00,t0,ymn=vvm[1],slope=vvm[2],sigma=vvm[3],xi=vvm[4])
		Fp1=pgev_p1(q00,t0,ymn=vvp[1],slope=vvp[2],sigma=vvp[3],xi=vvp[4])
		mu1[i,]=-(Fp1-Fm1)/(2*dd[i])
	}
#	cat("new mu1=",sum(mu1),"\n")
	return(mu1)
}
#' GEVD-with-p1: DMGS equation 3.3 mu2 term, fixed shape parameter
#' @inheritParams manf
gev_p1k3_mu2f=function(alpha,t0,v1,d1,v2,d2,v3,fd3,v4){
	d4=999
	q00=qgev_p1((1-alpha),t0,ymn=v1,slope=v2,sigma=v3,xi=v4)
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
				Fm1=pgev_p1(q00,t0,ymn=vvm[1],slope=vvm[2],sigma=vvm[3],xi=v4)
				F00=pgev_p1(q00,t0,ymn=vv0[1],slope=vv0[2],sigma=vv0[3],xi=v4)
				Fp1=pgev_p1(q00,t0,ymn=vvp[1],slope=vvp[2],sigma=vvp[3],xi=v4)
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
				Fm1m1=pgev_p1(q00,t0,ymn=vvmm[1],slope=vvmm[2],sigma=vvmm[3],xi=v4)
				Fm1p1=pgev_p1(q00,t0,ymn=vvmp[1],slope=vvmp[2],sigma=vvmp[3],xi=v4)
				Fp1m1=pgev_p1(q00,t0,ymn=vvpm[1],slope=vvpm[2],sigma=vvpm[3],xi=v4)
				Fp1p1=pgev_p1(q00,t0,ymn=vvpp[1],slope=vvpp[2],sigma=vvpp[3],xi=v4)
				mu2[i,j,]=-(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j])
				mu2[j,i,]=mu2[i,j,]
			}
		}
	}
#	cat("new mu2=",sum(mu2),"\n")
	return(mu2)
}
#' GEVD-with-p1: DMGS equation 1.2 f2 term
#' @inheritParams manf
gev_p1k3_f2f=function(y,t0,v1,d1,v2,d2,v3,fd3,v4){
	d3=fd3*v3
	ymax=ifelse(v4<0,v1+v2*t0-2*(d1+d2*t0)-(v3-2*d3)/(v4),Inf)
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
				Fm1=dgev_p1(y,t0,ymn=vvm[1],slope=vvm[2],sigma=vvm[3],xi=v4)
				F00=dgev_p1(y,t0,ymn=vv0[1],slope=vv0[2],sigma=vv0[3],xi=v4)
				Fp1=dgev_p1(y,t0,ymn=vvp[1],slope=vvp[2],sigma=vvp[3],xi=v4)
				f2[i,i,]=ifelse(y<ymax,(Fp1-2*F00+Fm1)/(dd[i]*dd[i]),0)
			} else if(i<j) {
				vvmm=vv
				vvmp=vv
				vvpm=vv
				vvpp=vv
				vvmm[i]=vv[i]-dd[i];vvmm[j]=vv[j]-dd[j]
				vvmp[i]=vv[i]-dd[i];vvmp[j]=vv[j]+dd[j]
				vvpm[i]=vv[i]+dd[i];vvpm[j]=vv[j]-dd[j]
				vvpp[i]=vv[i]+dd[i];vvpp[j]=vv[j]+dd[j]
				Fm1m1=dgev_p1(y,t0,ymn=vvmm[1],slope=vvmm[2],sigma=vvmm[3],xi=v4)
				Fm1p1=dgev_p1(y,t0,ymn=vvmp[1],slope=vvmp[2],sigma=vvmp[3],xi=v4)
				Fp1m1=dgev_p1(y,t0,ymn=vvpm[1],slope=vvpm[2],sigma=vvpm[3],xi=v4)
				Fp1p1=dgev_p1(y,t0,ymn=vvpp[1],slope=vvpp[2],sigma=vvpp[3],xi=v4)
				f2[i,j,]=ifelse(y<ymax,(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j]),0)
				f2[j,i,]=f2[i,j,]
			}
		}
	}
#	cat("new f2=",sum(f2),"\n")
	return(f2)
}
#' GEVD-with-p1: DMGS equation 1.2 f2 term
#' @inheritParams manf
gev_p1_f2f=function(y,t0,v1,d1,v2,d2,v3,fd3,v4,d4){
	d3=fd3*v3
	ymax=ifelse(v4<0,v1+v2*t0-2*(d1+d2*t0)-(v3-2*d3)/(v4-2*d4),Inf)
# new method
	dd=c(d1,d2,d3,d4)
	vv=c(v1,v2,v3,v4)
	f2=array(0,c(4,4,length(y)))
	for (i in 1:4){
		for (j in 1:4){
			if(i==j){
				vvm=vv
				vv0=vv
				vvp=vv
				vvm[i]=vv[i]-dd[i]
				vvp[i]=vv[i]+dd[i]
				Fm1=dgev_p1(y,t0,ymn=vvm[1],slope=vvm[2],sigma=vvm[3],xi=vvm[4])
				F00=dgev_p1(y,t0,ymn=vv0[1],slope=vv0[2],sigma=vv0[3],xi=vv0[4])
				Fp1=dgev_p1(y,t0,ymn=vvp[1],slope=vvp[2],sigma=vvp[3],xi=vvp[4])
				f2[i,i,]=ifelse(y<ymax,(Fp1-2*F00+Fm1)/(dd[i]*dd[i]),0)
			} else if(i<j) {
				vvmm=vv
				vvmp=vv
				vvpm=vv
				vvpp=vv
				vvmm[i]=vv[i]-dd[i];vvmm[j]=vv[j]-dd[j]
				vvmp[i]=vv[i]-dd[i];vvmp[j]=vv[j]+dd[j]
				vvpm[i]=vv[i]+dd[i];vvpm[j]=vv[j]-dd[j]
				vvpp[i]=vv[i]+dd[i];vvpp[j]=vv[j]+dd[j]
				Fm1m1=dgev_p1(y,t0,ymn=vvmm[1],slope=vvmm[2],sigma=vvmm[3],xi=vvmm[4])
				Fm1p1=dgev_p1(y,t0,ymn=vvmp[1],slope=vvmp[2],sigma=vvmp[3],xi=vvmp[4])
				Fp1m1=dgev_p1(y,t0,ymn=vvpm[1],slope=vvpm[2],sigma=vvpm[3],xi=vvpm[4])
				Fp1p1=dgev_p1(y,t0,ymn=vvpp[1],slope=vvpp[2],sigma=vvpp[3],xi=vvpp[4])
				f2[i,j,]=ifelse(y<ymax,(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j]),0)
				f2[j,i,]=f2[i,j,]
			}
		}
	}
#	cat("new f2=",sum(f2),"\n")
	return(f2)
}
#' GEVD-with-p1: DMGS equation 3.3 mu2 term
#' @inheritParams manf
gev_p1_mu2f=function(alpha,t0,v1,d1,v2,d2,v3,fd3,v4,d4){
	q00=qgev_p1((1-alpha),t0,ymn=v1,slope=v2,sigma=v3,xi=v4)
	d3=fd3*v3
# new method
	dd=c(d1,d2,d3,d4)
	vv=c(v1,v2,v3,v4)
	mu2=array(0,c(4,4,length(alpha)))
	for (i in 1:4){
		for (j in 1:4){
			if(i==j){
				vvm=vv
				vv0=vv
				vvp=vv
				vvm[i]=vv[i]-dd[i]
				vvp[i]=vv[i]+dd[i]
				Fm1=pgev_p1(q00,t0,ymn=vvm[1],slope=vvm[2],sigma=vvm[3],xi=vvm[4])
				F00=pgev_p1(q00,t0,ymn=vv0[1],slope=vv0[2],sigma=vv0[3],xi=vv0[4])
				Fp1=pgev_p1(q00,t0,ymn=vvp[1],slope=vvp[2],sigma=vvp[3],xi=vvp[4])
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
				Fm1m1=pgev_p1(q00,t0,ymn=vvmm[1],slope=vvmm[2],sigma=vvmm[3],xi=vvmm[4])
				Fm1p1=pgev_p1(q00,t0,ymn=vvmp[1],slope=vvmp[2],sigma=vvmp[3],xi=vvmp[4])
				Fp1m1=pgev_p1(q00,t0,ymn=vvpm[1],slope=vvpm[2],sigma=vvpm[3],xi=vvpm[4])
				Fp1p1=pgev_p1(q00,t0,ymn=vvpp[1],slope=vvpp[2],sigma=vvpp[3],xi=vvpp[4])
				mu2[i,j,]=-(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j])
				mu2[j,i,]=mu2[i,j,]
			}
		}
	}
#	cat("new mu2=",sum(mu2),"\n")
	return(mu2)
}
#' Derivative of information matrix, based on ldd
#' @inheritParams manf
gev_p1_ggd=function(x,t,v1,d1,v2,d2,v3,fd3,v4,d4){
	ggd=array(0,c(4,4,4))
  d3=fd3*v3

	ggm1=gev_p1_ldd(x,t,v1-d1,d1,v2,d2,v3,fd3,v4,d4)
	ggp1=gev_p1_ldd(x,t,v1+d1,d1,v2,d2,v3,fd3,v4,d4)
	ggd[1,,]=(ggp1-ggm1)/(2*d1)

	ggm2=gev_p1_ldd(x,t,v1,d1,v2-d2,d2,v3,fd3,v4,d4)
	ggp2=gev_p1_ldd(x,t,v1,d1,v2+d2,d2,v3,fd3,v4,d4)
	ggd[2,,]=(ggp2-ggm2)/(2*d2)

	ggm3=gev_p1_ldd(x,t,v1,d1,v2,d2,v3-d3,fd3,v4,d4)
	ggp3=gev_p1_ldd(x,t,v1,d1,v2,d2,v3+d3,fd3,v4,d4)
	ggd[3,,]=(ggp3-ggm3)/(2*d3)

	ggm4=gev_p1_ldd(x,t,v1,d1,v2,d2,v3,fd3,v4-d4,d4)
	ggp4=gev_p1_ldd(x,t,v1,d1,v2,d2,v3,fd3,v4+d4,d4)
	ggd[4,,]=(ggp4-ggm4)/(2*d4)

  return(ggd)
}
#' Analytical expressions for Predictive Means
#' RHP mean based on the expectation of DMGS equation 2.1
#' @inheritParams manf
gev_p1_means=function(means,t0,ml_params,lddi,lddi_k4,lddd,lddd_k4,
									lambdad_flat,lambdad_rh_mle,lambdad_rh_flat,lambdad_jp,
									nx,dim=4){

	if(means){
# intro
		eulerconstant=0.57721566490153286060651209008240243104215933593992
		ymn=ml_params[1]
		slope=ml_params[2]
		sigma=ml_params[3]
		xi=ml_params[4]

# set up derivative arrays
		meand1=array(0,c(4,1))
		meand2=array(0,c(4,4,1)) #but all zero for gumbel
		meand1_k4=array(0,c(4,1))
		meand2_k4=array(0,c(3,3,1)) #but all zero for gumbel

		if(ml_params[4]==0){
# xi=0 case

# mle
			ml_mean=ymn+slope*t0+sigma*eulerconstant
# calculate first derivative array for bayesian xi=0 cases
			meand1[1,1]=1
			meand1[2,1]=t0
			meand1[3,1]=eulerconstant
			meand1[4,1]=0
# and copy for the 3d case
			meand1_k4[1,1]=meand1[1,1]
			meand1_k4[2,1]=meand1[2,1]
			meand1_k4[3,1]=meand1[3,1]
# meand2 is all zero as initialized

		} else{
# non-gumbel case

			g0=gamma(1-xi)
			g1=g0*digamma(1-xi)
			g2=(trigamma(1-xi)*g0*g0+g1*g1)/g0
# mle
			ml_mean=ymn+slope*t0+sigma*(g0-1)/xi
# calculate first derivative array for bayesian xi!=0 cases
			meand1[1,1]=1
			meand1[2,1]=t0
			meand1[3,1]=(g0-1)/xi
			meand1[4,1]=(1-g0-xi*g1)*sigma/(xi*xi)
# and copy for the 2d case
			meand1_k4[1,1]=meand1[1,1]
			meand1_k4[2,1]=meand1[2,1]
			meand1_k4[3,1]=meand1[3,1]
# calculate second derivative array (only has 1 non-zero term!)
			meand2[3,4,1]=(1-g0-xi*g1)/(xi*xi)
			meand2[4,3,1]=meand2[3,4,1]
			meand2[4,4,1]=sigma*(-2+2*g0+2*xi*g1+xi*xi*g2)/(xi*xi*xi)
		}
# I've realized now that when I integrate over xi, the mean in Inf

#	flat_mean			=ml_mean+dmgs(lddi,lddd,meand1,lambdad_flat,meand2,dim=4)/nx
		flat_mean			=Inf

		rh_ml_mean		=ml_mean+dmgs(lddi_k4,lddd_k4,meand1_k4,lambdad_rh_mle,meand2_k4,dim=3)/nx

#	rh_flat_mean	=ml_mean+dmgs(lddi,lddd,meand1,lambdad_rh_flat,meand2,dim=4)/nx
		rh_flat_mean	=Inf

#	jp_mean				=ml_mean+dmgs(lddi,lddd,meand1,lambdad_jp,meand2,dim=4)/nx
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
dgev_p1sub=function(x,t,y,t0,ics,d1=0.01,d2=0.01,fd3=0.01,d4=0.01,
	minxi,maxxi,extramodels=FALSE,aderivs=TRUE){

		nx=length(x)

		ics=gev_p1_setics(x,t,ics)
		opt=optim(ics,gev_p1_loglik,x=x,t=t,control=list(fnscale=-1))
		v1hat=opt$par[1]
		v2hat=opt$par[2]
		v3hat=opt$par[3]
		v4hat=min(maxxi,max(minxi,opt$par[4]))
		ml_params=c(v1hat,v2hat,v3hat,v4hat)
#		cat("ml_params=",ml_params,"\n")

# now that I've dropped dmgs, not sure I need this anymore
#		muhat0=v1hat+v2hat*t0
#		y=fixgevrange(y,muhat0,v3hat,v4hat)

# mle
		ml_pdf=dgev_p1(y,t0,ymn=v1hat,slope=v2hat,sigma=v3hat,xi=v4hat)
		ml_cdf=pgev_p1(y,t0,ymn=v1hat,slope=v2hat,sigma=v3hat,xi=v4hat)

# return
		list(
					ml_params=ml_params,
					ml_pdf=ml_pdf,
					ml_cdf=ml_cdf)
}

