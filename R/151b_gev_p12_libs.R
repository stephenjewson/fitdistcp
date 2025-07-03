#
# NOTE THAT PARAMETER 4 HAS A DEFAULT VALUE OF ZERO AND RANGE FROM -INF TO INF
# BECAUSE ITS INSIDE AN EXPONENTIAL
#
#' rgev for gev_p12 but with maxlik xi within bounds
#' @return Vector
#' @inheritParams manf
rgev_p12_minmax=function(nx,mu,sigma,xi,t1,t2,minxi=-0.45,maxxi=0.45,centering=TRUE){
	xihat=-9999
  if(centering){
  	t1=t1-mean(t1)
  	t2=t2-mean(t2)
  }
	while((xihat<minxi)||(xihat>maxxi)){ #0.46 also works...0.47 doesn't
		xx=extraDistr::rgev(nx,mu=mu,sigma=sigma,xi=xi)
		ics=gev_p12_setics(xx,t1,t2,c(0,0,0,0,0))
		opt1=optim(ics,gev_p12_loglik,x=xx,t1=t1,t2=t2,control=list(fnscale=-1))
		xihat=opt1$par[5]
	}
	return(xx)
}
#' Waic
#' @inherit manwaic return
#' @inheritParams manf
gev_p12_waic=function(waicscores,x,t1,t2,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4,v5hat,d5,
	lddi,lddd,lambdad,aderivs){
		if(waicscores){
			if(aderivs){
				f1f=gev_p12_f1fa(x,t1,t2,v1hat,v2hat,v3hat,v4hat,v5hat)
				f2f=gev_p12_f2fa(x,t1,t2,v1hat,v2hat,v3hat,v4hat,v5hat)
			} else {
				f1f=gev_p12_f1f(x,t1,t2,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4,v5hat,d5)
				f2f=gev_p12_f2f(x,t1,t2,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4,v5hat,d5)
			}
			fhatx=dgev_p12(x,t1,t2,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat,xi=v5hat,log=FALSE)
			waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad,f2f,dim=5)
			waic1=waic$waic1
			waic2=waic$waic2
		}else{
			waic1="waicscores not selected"
			waic2="waicscores not selected"
		}
		list(waic1=waic1,waic2=waic2)
}
#' Predicted Parameter and Generalized Residuals
#' @inherit manpredictor return
#' @inheritParams manf
gev_p12_predictordata=function(predictordata,x,t1,t2,t01,t02,params){
	if(predictordata){
#
# calculate the probabilities of the data using the fited model
#
# t is always a matrix
# t0 is always a vector
		a=params[1]
		b=params[2]
		sc1=params[3]
		sc2=params[4]
		sh=params[5]
		mu=a+b*t1
		sigma=exp(sc1+sc2*t2)
		px=extraDistr::pgev(x,mu=mu,sigma=sigma,xi=sh)
#
# calculate the quantiles for those probabilities at t01,t02
#
		mu0=a+b*t01
		sigma0=exp(sc1+sc2*t02)
		qx=extraDistr::qgev(px,mu=mu0,sigma=sigma0,xi=sh)
	} else{
		predictedparameter="predictordata not selected"
		adjustedx="predictordata not selected"
	}


	list(predictedparameter=mu,adjustedx=qx)
}#' Logf for RUST
#' @inherit manlogf return
#' @inheritParams manf
gev_p12_logf=function(params,x,t1,t2){
	a=params[1]
	b=params[2]
	sc1=params[3]
	sc2=params[4]
	sh=params[5]
#	if(is.vector(t)){
#		mu=a+b*t1
#		sigma=exp(sc1+sc2*t2)
#	} else {
#		mu=a+b*t1
#		sigma=exp(sc1+sc2*t2)
#	}
	mu=a+b*t1
	sigma=exp(sc1+sc2*t2)
	logf=sum(extraDistr::dgev(x,mu=mu,sigma=sigma,xi=sh,log=TRUE))
	return(logf)
}#' Set initial conditions
#' @return Vector
#' @inheritParams manf
gev_p12_setics=function(x,t1,t2,ics){
# t is always a matrix
	nx=length(x)
	if((ics[1]==0)&&(ics[2]==0)&&(ics[3]==0)&&(ics[4]==0)&&(ics[5]==0)){
		lm=lm(x~t1)
		ics[1]=lm$coefficients[1]
		ics[2]=lm$coefficients[2]
		xhat=ics[1]+ics[2]*t1
		ics[3]=0
		ics[4]=0 #zero because it's inside an exponential
		ics[5]=0
	}
	return(ics)
}
#'  observed log-likelihood function
#' @inherit manloglik return
#' @inheritParams manf
gev_p12_loglik=function(vv,x,t1,t2){
# t is always a matrix
	n=length(x)
	mu=vv[1]+vv[2]*t1 #so mean is a vector, just like x
	sigma=exp(vv[3]+vv[4]*t2)
	loglik=sum(extraDistr::dgev(x,mu=mu,sigma=sigma,xi=vv[5],log=TRUE))
	return(loglik)
}
#' Check MLE
#' @return No return value (just a message to the screen).
#' @inheritParams manf
gev_p12_checkmle=function(ml_params,minxi,maxxi){
	v1hat=ml_params[1]
	v2hat=ml_params[2]
	v3hat=ml_params[3]
	v4hat=ml_params[4]
	v5hat=ml_params[5]
	if(is.na(v1hat))stop()
	if(is.na(v2hat))stop()
	if(is.na(v3hat))stop()
	if(is.na(v4hat))stop()
	if(is.na(v5hat))stop()
	if(v5hat<minxi){warning("\n***v5hat=",v5hat,"=> execution halted because maxlik shape parameter <",minxi,"***");stop()}
	if(v5hat>maxxi){warning("\n***v5hat=",v5hat,"=> execution halted because maxlik shape parameter >",maxxi,"***");stop()}
}
#' GEVD-with-p1: Quantile function
#' @inherit manvector return
#' @inheritParams manf
qgev_p12=function(p,t1,t2,ymn,slope,sigma1,sigma2,xi){
# t is sometimes a vector, sometimes a matrix

#	if(is.vector(t)){
#		mu=ymn+slope*t1
#		sigma=exp(sigma1+sigma2*t2)
#	} else {
#		mu=ymn+slope*t1
#		sigma=exp(sigma1+sigma2*t2)
#	}
	mu=ymn+slope*t1
	sigma=exp(sigma1+sigma2*t2)

	return(extraDistr::qgev(p,mu=mu,sigma=sigma,xi=xi))

}
#' GEVD-with-p1: Density function
#' @inherit manvector return
#' @inheritParams manf
dgev_p12=function(x,t1,t2,ymn,slope,sigma1,sigma2,xi,log=FALSE){
# t is sometimes a vector, sometimes a matrix

#	if(is.vector(t)){
#		mu=ymn+slope*t1
#		sigma=exp(sigma1+sigma2*t2)
#	} else {
#		mu=ymn+slope*t1
#		sigma=exp(sigma1+sigma2*t2)
#	}
		mu=ymn+slope*t1
		sigma=exp(sigma1+sigma2*t2)

	return(extraDistr::dgev(x,mu=mu,sigma=sigma,xi=xi,log=log))

}
#' GEVD-with-p1: Distribution function
#' @inherit manvector return
#' @inheritParams manf
pgev_p12=function(y,t1,t2,ymn,slope,sigma1,sigma2,xi){
# t is sometimes a vector, sometimes a matrix

#	if(is.vector(t)){
#		mu=ymn+slope*t1
#		sigma=exp(sigma1+sigma2*t2)
#	} else {
#		mu=ymn+slope*t1
#		sigma=exp(sigma1+sigma2*t2)
#	}
		mu=ymn+slope*t1
		sigma=exp(sigma1+sigma2*t2)

	return(extraDistr::pgev(y,mu=mu,sigma=sigma,xi=xi))

}
#' One component of the second derivative of the normalized log-likelihood
#' @inherit manlnn return
#' @inheritParams manf
gev_p12_lmn=function(x,t1,t2,v1,d1,v2,d2,v3,d3,v4,d4,v5,d5,mm,nn){

	nparams=5
	nlmn3=3
	nlmn4=4
	net3=matrix(0,nlmn3,nparams)
	net4=matrix(0,nlmn4,nparams)
	lmn=matrix(0,nlmn4)
	dd=c(d1,d2,d3,d4,d5)
	vv=c(v1,v2,v3,v4,v5)
	vvd=matrix(0,nparams)
	nx=length(x)
# different
	if(mm!=nn){
		net4[,mm]=c(-1,-1,1,1)
		net4[,nn]=c(-1,1,-1,1)
		for (i in 1:nlmn4){
			for (j in 1:nparams){
				vvd[j]=vv[j]+net4[i,j]*dd[j]
			}
			lmn[i]=sum(dgev_p12(x,t1,t2,ymn=vvd[1],slope=vvd[2],sigma1=vvd[3],sigma2=vvd[4],xi=vvd[5],log=TRUE))/nx
		}
		dld=(lmn[1]-lmn[2]-lmn[3]+lmn[4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:(nlmn4-1)){
			for (j in 1:nparams){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[i]=sum(dgev_p12(x,t1,t2,ymn=vvd[1],slope=vvd[2],sigma1=vvd[3],sigma2=vvd[4],xi=vvd[5],log=TRUE))/nx
		}
		dld=(lmn[1]-2*lmn[2]+lmn[3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' Second derivative matrix of the normalized log-likelihood, with fixed shape parameter
#' @inherit manldd return
#' @inheritParams manf
gev_p12k3_ldd=function(x,t1,t2,v1,d1,v2,d2,v3,d3,v4,d4,v5){
	d5=999
	nparams=4
	nx=length(x)
	ldd=matrix(0,nparams,nparams)
# unique terms
	for (i in 1:nparams){
		for (j in i:nparams){
			ldd[i,j]=gev_p12_lmn(x,t1,t2,v1,d1,v2,d2,v3,d3,v4,d4,v5,d5,i,j)
		}
	}
# copies
	for (i in nparams:2){
		for (j in 1:(i-1)){
			ldd[i,j]=ldd[j,i]
		}
	}
	return(ldd)
}
#' Second derivative matrix of the normalized log-likelihood
#' @inherit manldd return
#' @inheritParams manf
gev_p12_ldd=function(x,t1,t2,v1,d1,v2,d2,v3,d3,v4,d4,v5,d5){
	nx=length(x)
	nparams=5
	ldd=matrix(0,nparams,nparams)
	for (i in 1:nparams){
		for (j in i:nparams){
			ldd[i,j]=gev_p12_lmn(x,t1,t2,v1,d1,v2,d2,v3,d3,v4,d4,v5,d5,i,j)
		}
	}
	for (i in nparams:2){
		for (j in 1:(i-1)){
			ldd[i,j]=ldd[j,i]
		}
	}
	return(ldd)
}
#' One component of the second derivative of the normalized log-likelihood
#' @inherit manlnnn return
#' @inheritParams manf
gev_p12_lmnp=function(x,t1,t2,v1,d1,v2,d2,v3,d3,v4,d4,v5,d5,mm,nn,rr){

	nparams=5
	nlmn4=4
	nlmn6=6
	nlmn8=8
	net4=matrix(0,nlmn4,nparams)
	net6=matrix(0,nlmn6,nparams)
	net8=matrix(0,nlmn8,nparams)
	lmn=matrix(0,nlmn8)
	dd=c(d1,d2,d3,d4,d5)
	vv=c(v1,v2,v3,v4,v5)
	vvd=matrix(0,nparams)
	nx=length(x)
# all diff
	if ((mm!=nn)&(nn!=rr)&(rr!=mm)){
		net8[,mm]=c(-1,1,-1,1,-1,1,-1,1)
		net8[,nn]=c(-1,-1,1,1,-1,-1,1,1)
		net8[,rr]=c(-1,-1,-1,-1,1,1,1,1)
		for (i in 1:nlmn8){
			for (j in 1:nparams){
				vvd[j]=vv[j]+net8[i,j]*dd[j]
			}
			lmn[i]=sum(dgev_p12(x,t1,t2,ymn=vvd[1],slope=vvd[2],sigma1=vvd[3],sigma2=vvd[4],xi=vvd[5],log=TRUE))/nx
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
		for (i in 1:nlmn4){
			for (j in 1:nparams){
				vvd[j]=vv[j]+net4[i,j]*dd[j]
			}
			lmn[i]=sum(dgev_p12(x,t1,t2,ymn=vvd[1],slope=vvd[2],sigma1=vvd[3],sigma2=vvd[4],xi=vvd[5],log=TRUE))/nx
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
		for (i in 1:nlmn6){
			for (j in 1:nparams){
				vvd[j]=vv[j]+net6[i,j]*dd[j]
			}
			lmn[i]=sum(dgev_p12(x,t1,t2,ymn=vvd[1],slope=vvd[2],sigma1=vvd[3],sigma2=vvd[4],xi=vvd[5],log=TRUE))/nx
		}
		dld1=(lmn[3]-2*lmn[2]+lmn[1])/(dd[m2]*dd[m2])
		dld2=(lmn[6]-2*lmn[5]+lmn[4])/(dd[m2]*dd[m2])
		dld=(dld2-dld1)/(2*dd[n2])
	}
	return(dld)
}
#' Third derivative tensor of the normalized log-likelihood, with fixed shape parameter
#' @inherit manlddd return
#' @inheritParams manf
gev_p12k3_lddd=function(x,t1,t2,v1,d1,v2,d2,v3,d3,v4,d4,v5){ #no d5 deliberately
	nparams=4
	lddd=array(0,c(nparams,nparams,nparams))
	d5=999
# calculate the unique values
	for (i in 1:nparams){
		for (j in i:nparams){
			for (k in j:nparams){
				lddd[i,j,k]=gev_p12_lmnp(x,t1,t2,v1,d1,v2,d2,v3,d3,v4,d4,v5,d5,i,j,k)
			}
		}
	}
# steves dumb algorithm for filling in the non-unique values
	for (i in 1:nparams){
		for (j in 1:nparams){
			for (k in 1:nparams){
				a=c(i,j,k)
				b=sort(a)
				lddd[a[1],a[2],a[3]]=lddd[b[1],b[2],b[3]]
			}
		}
	}
	return(lddd)
}
#' Third derivative tensor of the normalized log-likelihood, with fixed shape parameter
#' @inherit manlddd return
#' @inheritParams manf
gev_p12_lddd=function(x,t1,t2,v1,d1,v2,d2,v3,d3,v4,d4,v5,d5){
	nparams=5
	lddd=array(0,c(nparams,nparams,nparams))
# calculate the unique values
	for (i in 1:nparams){
		for (j in i:nparams){
			for (k in j:nparams){
				lddd[i,j,k]=gev_p12_lmnp(x,t1,t2,v1,d1,v2,d2,v3,d3,v4,d4,v5,d5,i,j,k)
			}
		}
	}
# steves dumb algorithm for filling in the non-unique values
	for (i in 1:nparams){
		for (j in 1:nparams){
			for (k in 1:nparams){
				a=c(i,j,k)
				b=sort(a)
				lddd[a[1],a[2],a[3]]=lddd[b[1],b[2],b[3]]
			}
		}
	}
	return(lddd)
}
#' DMGS equation 2.1, f1 term, fixed shape parameter
#' @inherit man1f return
#' @inheritParams manf
gev_p12k3_f1f=function(y,t01,t02,v1,d1,v2,d2,v3,d3,v4,d4,v5){ #no d5 deliberately

	nparams=4
	ymax=ifelse(v5<0,v1+v2*t01-2*(d1+d2*t01)-(exp(v3+v4*t02)-2*(d3+d4))/(v5),Inf)
# new method
	dd=c(d1,d2,d3,d4)
	vv=c(v1,v2,v3,d4)
	f1=matrix(0,nparams,length(y))
	for (i in 1:nparams){
		vvm=vv
		vvp=vv
		vvm[i]=vv[i]-dd[i]
		vvp[i]=vv[i]+dd[i]
		Fm1=dgev_p12(y,t01,t02,ymn=vvm[1],slope=vvm[2],sigma1=vvm[3],sigma2=vvm[4],xi=v5)
		Fp1=dgev_p12(y,t01,t02,ymn=vvp[1],slope=vvp[2],sigma1=vvp[3],sigma2=vvp[4],xi=v5)
		f1[i,]=ifelse(y<ymax,(Fp1-Fm1)/(2*dd[i]),0)
	}
	return(f1)
}
#' DMGS equation 2.1, f1 term
#' @inherit man1f return
#' @inheritParams manf
gev_p12_f1f=function(y,t01,t02,v1,d1,v2,d2,v3,d3,v4,d4,v5,d5){

	nparams=5
	ymax=ifelse(v5<0,v1+v2*t01-2*(d1+d2*t01)-(exp(v3+v4*t02)-2*(d3+d4*t01))/(v5),Inf)
#	ymax=ifelse(v4<0,v1+v2*t01-2*(d1+d2*t01)-(v3-2*d3)/(v4-2*d4),Inf)
# new method
	dd=c(d1,d2,d3,d4,d5)
	vv=c(v1,v2,v3,v4,v5)
	f1=matrix(0,nparams,length(y))
	for (i in 1:nparams){
		vvm=vv
		vvp=vv
		vvm[i]=vv[i]-dd[i]
		vvp[i]=vv[i]+dd[i]
		Fm1=dgev_p12(y,t01,t02,ymn=vvm[1],slope=vvm[2],sigma1=vvm[3],sigma2=vvm[4],xi=vvm[5])
		Fp1=dgev_p12(y,t01,t02,ymn=vvp[1],slope=vvp[2],sigma1=vvp[3],sigma2=vvp[4],xi=vvp[5])
		f1[i,]=ifelse(y<ymax,(Fp1-Fm1)/(2*dd[i]),0)
	}
	return(f1)
}
#' GEVD-with-p1: DMGS equation 3.3 mu1 term, fixed shape parameter
#' @inherit man1f return
#' @inheritParams manf
gev_p12k3_mu1f=function(alpha,t01,t02,v1,d1,v2,d2,v3,d3,v4,d4,v5){ #no d5 deliberately

	nparams=4
	d4=999
	q00=qgev_p12((1-alpha),t01,t02,ymn=v1,slope=v2,sigma1=v3,sigma2=v4,xi=v5)

# new method
	dd=c(d1,d2,d3,d4)
	vv=c(v1,v2,v3,v4)
	mu1=matrix(0,nparams,length(alpha))
	for (i in 1:nparams){
		vvm=vv
		vvp=vv
		vvm[i]=vv[i]-dd[i]
		vvp[i]=vv[i]+dd[i]
		Fm1=pgev_p12(q00,t01,t02,ymn=vvm[1],slope=vvm[2],sigma1=vvm[3],sigma2=vvm[4],xi=v5)
		Fp1=pgev_p12(q00,t01,t02,ymn=vvp[1],slope=vvp[2],sigma1=vvp[3],sigma2=vvp[4],xi=v5)
		mu1[i,]=-(Fp1-Fm1)/(2*dd[i])
	}
	return(mu1)
}
#' GEVD-with-p1: DMGS equation 3.3 mu1 term
#' @inherit man1f return
#' @inheritParams manf
gev_p12_mu1f=function(alpha,t01,t02,v1,d1,v2,d2,v3,d3,v4,d4,v5,d5){

	nparams=5
	q00=qgev_p12((1-alpha),t01,t02,ymn=v1,slope=v2,sigma1=v3,sigma2=v4,xi=v5)

# new method
	dd=c(d1,d2,d3,d4,d5)
	vv=c(v1,v2,v3,v4,v5)
	mu1=matrix(0,nparams,length(alpha))
	for (i in 1:nparams){
		vvm=vv
		vvp=vv
		vvm[i]=vv[i]-dd[i]
		vvp[i]=vv[i]+dd[i]
		Fm1=pgev_p12(q00,t01,t02,ymn=vvm[1],slope=vvm[2],sigma1=vvm[3],sigma2=vvm[4],xi=vvm[5])
		Fp1=pgev_p12(q00,t01,t02,ymn=vvp[1],slope=vvp[2],sigma1=vvp[3],sigma2=vvp[4],xi=vvp[5])
		mu1[i,]=-(Fp1-Fm1)/(2*dd[i])
	}
	return(mu1)
}
#' GEVD-with-p1: DMGS equation 3.3 mu2 term, fixed shape parameter
#' @inherit man2f return
#' @inheritParams manf
gev_p12k3_mu2f=function(alpha,t01,t02,v1,d1,v2,d2,v3,d3,v4,d4,v5){ #no d5 deliberately

	nparams=4
	d5=999
	q00=qgev_p12((1-alpha),t01,t02,ymn=v1,slope=v2,sigma1=v3,sigma2=v4,xi=v5)

# new method
	dd=c(d1,d2,d3,d4)
	vv=c(v1,v2,v3,v4)
	mu2=array(0,c(nparams,nparams,length(alpha)))
	for (i in 1:nparams){
		for (j in 1:nparams){
			if(i==j){
				vvm=vv
				vv0=vv
				vvp=vv
				vvm[i]=vv[i]-dd[i]
				vvp[i]=vv[i]+dd[i]
				Fm1=pgev_p12(q00,t01,t02,ymn=vvm[1],slope=vvm[2],sigma1=vvm[3],sigma2=vvm[4],xi=v5)
				F00=pgev_p12(q00,t01,t02,ymn=vv0[1],slope=vv0[2],sigma1=vv0[3],sigma2=vv0[4],xi=v5)
				Fp1=pgev_p12(q00,t01,t02,ymn=vvp[1],slope=vvp[2],sigma1=vvp[3],sigma2=vvp[4],xi=v5)
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
				Fm1m1=pgev_p12(q00,t01,t02,ymn=vvmm[1],slope=vvmm[2],sigma1=vvmm[3],sigma2=vvmm[4],xi=v5)
				Fm1p1=pgev_p12(q00,t01,t02,ymn=vvmp[1],slope=vvmp[2],sigma1=vvmp[3],sigma2=vvmp[4],xi=v5)
				Fp1m1=pgev_p12(q00,t01,t02,ymn=vvpm[1],slope=vvpm[2],sigma1=vvpm[3],sigma2=vvpm[4],xi=v5)
				Fp1p1=pgev_p12(q00,t01,t02,ymn=vvpp[1],slope=vvpp[2],sigma1=vvpp[3],sigma2=vvpp[4],xi=v5)
				mu2[i,j,]=-(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j])
				mu2[j,i,]=mu2[i,j,]
			}
		}
	}
	return(mu2)
}
#' GEVD-with-p1: DMGS equation 1.2 f2 term
#' @inherit man2f return
#' @inheritParams manf
gev_p12k3_f2f=function(y,t01,t02,v1,d1,v2,d2,v3,d3,v4,d4,v5){ #no d5 deliberately

	nparams=4
	ymax=ifelse(v5<0,v1+v2*t01-2*(d1+d2*t01)-(exp(v3+v4*t02)-2*(d3+d4*t01))/(v5),Inf)
#	ymax=ifelse(v4<0,v1+v2*t0-2*(d1+d2*t0)-(v3-2*d3)/(v4),Inf)
# new method
	dd=c(d1,d2,d3,d4)
	vv=c(v1,v2,v3,v4)
	f2=array(0,c(nparams,nparams,length(y)))
	for (i in 1:nparams){
		for (j in 1:nparams){
			if(i==j){
				vvm=vv
				vv0=vv
				vvp=vv
				vvm[i]=vv[i]-dd[i]
				vvp[i]=vv[i]+dd[i]
				Fm1=dgev_p12(y,t01,t02,ymn=vvm[1],slope=vvm[2],sigma1=vvm[3],sigma2=vvm[4],xi=v5)
				F00=dgev_p12(y,t01,t02,ymn=vv0[1],slope=vv0[2],sigma1=vv0[3],sigma2=vv0[4],xi=v5)
				Fp1=dgev_p12(y,t01,t02,ymn=vvp[1],slope=vvp[2],sigma1=vvp[3],sigma2=vvp[4],xi=v5)
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
				Fm1m1=dgev_p12(y,t01,t02,ymn=vvmm[1],slope=vvmm[2],sigma1=vvmm[3],sigma2=vvmm[4],xi=v5)
				Fm1p1=dgev_p12(y,t01,t02,ymn=vvmp[1],slope=vvmp[2],sigma1=vvmp[3],sigma2=vvmp[4],xi=v5)
				Fp1m1=dgev_p12(y,t01,t02,ymn=vvpm[1],slope=vvpm[2],sigma1=vvpm[3],sigma2=vvpm[4],xi=v5)
				Fp1p1=dgev_p12(y,t01,t02,ymn=vvpp[1],slope=vvpp[2],sigma1=vvpp[3],sigma2=vvpp[4],xi=v5)
				f2[i,j,]=ifelse(y<ymax,(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j]),0)
				f2[j,i,]=f2[i,j,]
			}
		}
	}
	return(f2)
}
#' GEVD-with-p1: DMGS equation 1.2 f2 term
#' @inherit man2f return
#' @inheritParams manf
gev_p12_f2f=function(y,t01,t02,v1,d1,v2,d2,v3,d3,v4,d4,v5,d5){

	nparams=5
	ymax=ifelse(v5<0,v1+v2*t01-2*(d1+d2*t01)-(exp(v3+v4*t02)-2*(d3+d4*t01))/(v5),Inf)
#	ymax=ifelse(v4<0,v1+v2*t0-2*(d1+d2*t0)-(v3-2*d3)/(v4-2*d4),Inf)
# new method
	dd=c(d1,d2,d3,d4,d5)
	vv=c(v1,v2,v3,v4,v5)
	f2=array(0,c(nparams,nparams,length(y)))
	for (i in 1:nparams){
		for (j in 1:nparams){
			if(i==j){
				vvm=vv
				vv0=vv
				vvp=vv
				vvm[i]=vv[i]-dd[i]
				vvp[i]=vv[i]+dd[i]
				Fm1=dgev_p12(y,t01,t02,ymn=vvm[1],slope=vvm[2],sigma1=vvm[3],sigma2=vvm[4],xi=vvm[5])
				F00=dgev_p12(y,t01,t02,ymn=vv0[1],slope=vv0[2],sigma1=vv0[3],sigma2=vv0[4],xi=vv0[5])
				Fp1=dgev_p12(y,t01,t02,ymn=vvp[1],slope=vvp[2],sigma1=vvp[3],sigma2=vvp[4],xi=vvp[5])
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
				Fm1m1=dgev_p12(y,t01,t02,ymn=vvmm[1],slope=vvmm[2],sigma1=vvmm[3],sigma2=vvmm[4],xi=vvmm[5])
				Fm1p1=dgev_p12(y,t01,t02,ymn=vvmp[1],slope=vvmp[2],sigma1=vvmp[3],sigma2=vvmp[4],xi=vvmp[5])
				Fp1m1=dgev_p12(y,t01,t02,ymn=vvpm[1],slope=vvpm[2],sigma1=vvpm[3],sigma2=vvpm[4],xi=vvpm[5])
				Fp1p1=dgev_p12(y,t01,t02,ymn=vvpp[1],slope=vvpp[2],sigma1=vvpp[3],sigma2=vvpp[4],xi=vvpp[5])
				f2[i,j,]=ifelse(y<ymax,(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j]),0)
				f2[j,i,]=f2[i,j,]
			}
		}
	}
	return(f2)
}
#' GEVD-with-p1: DMGS equation 3.3 mu2 term
#' @inherit man2f return
#' @inheritParams manf
gev_p12_mu2f=function(alpha,t01,t02,v1,d1,v2,d2,v3,d3,v4,d4,v5,d5){
	nparams=5
	q00=qgev_p12((1-alpha),t01,t02,ymn=v1,slope=v2,sigma1=v3,sigma2=v4,xi=v5)

# new method
	dd=c(d1,d2,d3,d4,d5)
	vv=c(v1,v2,v3,v4,v5)
	mu2=array(0,c(nparams,nparams,length(alpha)))
	for (i in 1:nparams){
		for (j in 1:nparams){
			if(i==j){
				vvm=vv
				vv0=vv
				vvp=vv
				vvm[i]=vv[i]-dd[i]
				vvp[i]=vv[i]+dd[i]
				Fm1=pgev_p12(q00,t01,t02,ymn=vvm[1],slope=vvm[2],sigma1=vvm[3],sigma2=vvm[4],xi=vvm[5])
				F00=pgev_p12(q00,t01,t02,ymn=vv0[1],slope=vv0[2],sigma1=vv0[3],sigma2=vv0[4],xi=vv0[5])
				Fp1=pgev_p12(q00,t01,t02,ymn=vvp[1],slope=vvp[2],sigma1=vvp[3],sigma2=vvp[4],xi=vvp[5])
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
				Fm1m1=pgev_p12(q00,t01,t02,ymn=vvmm[1],slope=vvmm[2],sigma1=vvmm[3],sigma2=vvmm[4],xi=vvmm[5])
				Fm1p1=pgev_p12(q00,t01,t02,ymn=vvmp[1],slope=vvmp[2],sigma1=vvmp[3],sigma2=vvmp[4],xi=vvmp[5])
				Fp1m1=pgev_p12(q00,t01,t02,ymn=vvpm[1],slope=vvpm[2],sigma1=vvpm[3],sigma2=vvpm[4],xi=vvpm[5])
				Fp1p1=pgev_p12(q00,t01,t02,ymn=vvpp[1],slope=vvpp[2],sigma1=vvpp[3],sigma2=vvpp[4],xi=vvpp[5])
				mu2[i,j,]=-(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j])
				mu2[j,i,]=mu2[i,j,]
			}
		}
	}
	return(mu2)
}
#' Derivative of information matrix, based on ldd
#' @inherit manldd return
#' @inheritParams manf
gev_p12_ggd=function(x,t1,t2,v1,d1,v2,d2,v3,d3,v4,d4,v5,d5){

	nparams=5
	ggd=array(0,c(nparams,nparams,nparams))


	ggm1=gev_p12_ldd(x,t1,t2,v1-d1,d1,v2,d2,v3,d3,v4,d4,v5,d5)
	ggp1=gev_p12_ldd(x,t1,t2,v1+d1,d1,v2,d2,v3,d3,v4,d4,v5,d5)
	ggd[1,,]=(ggp1-ggm1)/(2*d1)

	ggm2=gev_p12_ldd(x,t1,t2,v1,d1,v2-d2,d2,v3,d3,v4,d4,v5,d5)
	ggp2=gev_p12_ldd(x,t1,t2,v1,d1,v2+d2,d2,v3,d3,v4,d4,v5,d5)
	ggd[2,,]=(ggp2-ggm2)/(2*d2)

	ggm3=gev_p12_ldd(x,t1,t2,v1,d1,v2,d2,v3-d3,d3,v4,d4,v5,d5)
	ggp3=gev_p12_ldd(x,t1,t2,v1,d1,v2,d2,v3+d3,d3,v4,d4,v5,d5)
	ggd[3,,]=(ggp3-ggm3)/(2*d3)

	ggm4=gev_p12_ldd(x,t1,t2,v1,d1,v2,d2,v3,d3,v4-d4,d4,v5,d5)
	ggp4=gev_p12_ldd(x,t1,t2,v1,d1,v2,d2,v3,d3,v4+d4,d4,v5,d5)
	ggd[4,,]=(ggp4-ggm4)/(2*d4)

	ggm5=gev_p12_ldd(x,t1,t2,v1,d1,v2,d2,v3,d3,v4,d4,v5-d5,d5)
	ggp5=gev_p12_ldd(x,t1,t2,v1,d1,v2,d2,v3,d3,v4,d4,v5+d5,d5)
	ggd[5,,]=(ggp5-ggm5)/(2*d5)

	  return(ggd)
}
#' Analytical expressions for Predictive Means
#' RHP mean based on the expectation of DMGS equation 2.1
#' @inherit manmeans return
#' @inheritParams manf
gev_p12_means=function(means,t01,t02,ml_params,nx){

	if(means){
# intro
		eulerconstant=0.57721566490153286060651209008240243104215933593992
		ymn=ml_params[1]
		slope=ml_params[2]
		sigma1=ml_params[3]
		sigma2=ml_params[4]
		sigma=exp(sigma1+sigma2*t02)
		xi=ml_params[5]

		if(ml_params[5]==0){
# xi=0 case
			ml_mean=ymn+slope*t01+sigma*eulerconstant
		} else{
# xi!=0 case
			g0=gamma(1-xi)
			g1=g0*digamma(1-xi)
			g2=(trigamma(1-xi)*g0*g0+g1*g1)/g0
			ml_mean=ymn+slope*t01+sigma*(g0-1)/xi
		}
# return
		crhp_mle_mean="haven't worked it out yet"
		pu_mean=Inf
	}else{
		pu_mean="means not selected"
		ml_mean="means not selected"
		crhp_mle_mean="means not selected"
	}
	list(ml_mean=ml_mean,crhp_mle_mean=crhp_mle_mean,pu_mean=pu_mean)
}
#' Densities for 5 predictions
#' @inherit mandsub return
#' @inheritParams manf
dgev_p12sub=function(x,t1,t2,y,t01,t02,ics,d1=0.01,d2=0.01,d3=0.01,d4=0.01,d5=0.01,
	minxi,maxxi,debug,extramodels=FALSE,aderivs=TRUE){

		nx=length(x)

		ics=gev_p12_setics(x,t1,t2,ics)
		opt=optim(ics,gev_p12_loglik,x=x,t1=t1,t2=t2,control=list(fnscale=-1))
		v1hat=opt$par[1]
		v2hat=opt$par[2]
		v3hat=opt$par[3]
		v4hat=opt$par[4]
		v5hat=min(maxxi,max(minxi,opt$par[5]))
		ml_params=c(v1hat,v2hat,v3hat,v4hat,v5hat)

# now that I've dropped dmgs d and p, I don't think I need this any more
#		muhat0=v1hat+v2hat*t01
#		sghat0=exp(v3hat+v4hat*t02)
#		y=fixgevrange(y,muhat0,sghat0,v5hat)

# mle
		ml_pdf=dgev_p12(y,t01,t02,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat,xi=v5hat)
		ml_cdf=pgev_p12(y,t01,t02,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat,xi=v5hat)

# return
		list(ml_params=ml_params,
					ml_pdf=ml_pdf,
					ml_cdf=ml_cdf)
}

