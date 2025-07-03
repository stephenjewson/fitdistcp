#' rgev but with maxlik xi guaranteed within bounds
#' @return Vector
#' @inheritParams manf
rgev_minmax=function(nx,mu,sigma,xi,minxi=-0.45,maxxi=0.45){
	xihat=-9999
	while((xihat<minxi)||(xihat>maxxi)){
		xx=extraDistr::rgev(nx,mu=mu,sigma=sigma,xi=xi)
		ics=gev_setics(xx,c(0,0,0))
		opt1=optim(ics,gev_loglik,x=xx,control=list(fnscale=-1))
		xihat=opt1$par[3]
	}
	return(xx)
}
#' Waic
#' @inherit manwaic return
#' @inheritParams manf
gev_waic=function(waicscores,x,v1hat,d1,v2hat,fd2,v3hat,d3,lddi,lddd,
	lambdad,aderivs){
		if(waicscores){

			if(aderivs) f1f=gev_f1fa(x,v1hat,v2hat,v3hat)
			if(!aderivs)f1f=gev_f1f(x,v1hat,d1,v2hat,fd2,v3hat,d3)

			if(aderivs) f2f=gev_f2fa(x,v1hat,v2hat,v3hat)
			if(!aderivs)f2f=gev_f2f(x,v1hat,d1,v2hat,fd2,v3hat,d3)

			fhatx=extraDistr::dgev(x,mu=v1hat,sigma=v2hat,xi=v3hat)
			waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad,f2f,dim=3)
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
gev_logf=function(params,x){
	mu=params[1]
	sc=pmax(params[2],.Machine$double.eps)
	sh=params[3]
	logf=sum(extraDistr::dgev(x,mu=mu,sigma=sc,xi=sh,log=TRUE))-log(sc)
	return(logf)
}
#' Set initial conditions
#' @return Vector
#' @inheritParams manf
gev_setics=function(x,ics){
	if((ics[1]==0)&&(ics[2]==0)&&(ics[3]==0)){
		ics[1]=mean(x)
		ics[2]=sd(x)
		ics[3]=0
	}
	return(ics)
}
#' PWM parameter estimation
#' @return Vector
#' @inheritParams manf
gev_pwm_params=function(x){
		pw_params=matrix(0,3)
		g=fExtremes::gevFit(x,type="pwm")
		temp=g@fit$par.ests
# reorder to match the mu, sigma, xi ordering that I use
		pw_params[1]=temp[2]
		pw_params[2]=temp[3]
		pw_params[3]=temp[1]
return(pw_params)
}
#'  log-likelihood function
#' @inherit manloglik return
#' @inheritParams manf
gev_loglik=function(vv,x){
	n=length(x)
	loglik=sum(extraDistr::dgev(x,mu=vv[1],sigma=max(vv[2],.Machine$double.eps),xi=vv[3],log=TRUE))
	return(loglik)
}
#' Check MLE
#' @return No return value (just a message to the screen).
#' @inheritParams manf
gev_checkmle=function(ml_params,minxi,maxxi){
	v1hat=ml_params[1]
	v2hat=ml_params[2]
	v3hat=ml_params[3]
	if(is.na(v1hat))stop()
	if(is.na(v2hat))stop()
	if(is.na(v3hat))stop()
	if(v3hat<minxi){warning("\n***v3hat=",v3hat,"=> execution halted because maxlik shape parameter <",minxi,"***\n");stop()}
	if(v3hat>maxxi){warning("\n***v3hat=",v3hat,"=> execution halted because maxlik shape parameter >",maxxi,"***\n");stop()}
}
#' One component of the second derivative of the normalized log-likelihood
#' @inherit manlnn return
#' @inheritParams manf
gev_lmn=function(x,v1,d1,v2,fd2,v3,d3,mm,nn){
	d2=fd2*v2
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
			lmn[i]=sum(extraDistr::dgev(x,mu=vvd[1],sigma=vvd[2],xi=vvd[3],log=TRUE))/nx
		}
		dld=(lmn[1]-lmn[2]-lmn[3]+lmn[4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:3){
			for (j in 1:3){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[i]=sum(extraDistr::dgev(x,mu=vvd[1],sigma=vvd[2],xi=vvd[3],log=TRUE))/nx
		}
		dld=(lmn[1]-2*lmn[2]+lmn[3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' Second derivative matrix of the normalized log-likelihood
#' @inherit manldd return
#' @inheritParams manf
gev_ldd=function(x,v1,d1,v2,fd2,v3,d3){
#
	nx=length(x)
	ldd=matrix(0,3,3)
	for (i in 1:3){
		for (j in i:3){
			ldd[i,j]=gev_lmn(x,v1,d1,v2,fd2,v3,d3,i,j)
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
gev_lmnp=function(x,v1,d1,v2,fd2,v3,d3,mm,nn,rr){
	d2=fd2*v2
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
			lmn[i]=sum(extraDistr::dgev(x,mu=vvd[1],sigma=vvd[2],xi=vvd[3],log=TRUE))/nx
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
			lmn[i]=sum(extraDistr::dgev(x,mu=vvd[1],sigma=vvd[2],xi=vvd[3],log=TRUE))/nx
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
			lmn[i]=sum(extraDistr::dgev(x,mu=vvd[1],sigma=vvd[2],xi=vvd[3],log=TRUE))/nx
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
gev_lddd=function(x,v1,d1,v2,fd2,v3,d3){
	lddd=array(0,c(3,3,3))
# calculate the unique values
	for (i in 1:3){
		for (j in i:3){
			for (k in j:3){
				lddd[i,j,k]=gev_lmnp(x,v1,d1,v2,fd2,v3,d3,i,j,k)
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
#' DMGS equation 3.3, f1 term
#' @inherit man1f return
#' @inheritParams manf
gev_f1f=function(y,v1,d1,v2,fd2,v3,d3){
  d2=fd2*v2
	ymax=ifelse(v3<0,v1-2*d1-(v2-2*d2)/(v3-2*d3),Inf)
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
	F1m1=extraDistr::dgev(y,mu=v1m1,sigma=v200,xi=v3)
	F1p1=extraDistr::dgev(y,mu=v1p1,sigma=v200,xi=v3)
# v2 derivatives
	F2m1=extraDistr::dgev(y,mu=v100,sigma=v2m1,xi=v3)
	F2p1=extraDistr::dgev(y,mu=v100,sigma=v2p1,xi=v3)
# v3 derivatives
	F3m1=extraDistr::dgev(y,mu=v100,sigma=v200,xi=v3m1)
	F3p1=extraDistr::dgev(y,mu=v100,sigma=v200,xi=v3p1)
	f1=matrix(0,3,length(y))
	f1[1,]=ifelse(y<ymax,(F1p1-F1m1)/(2*d1),0)
	f1[2,]=ifelse(y<ymax,(F2p1-F2m1)/(2*d2),0)
	f1[3,]=ifelse(y<ymax,(F3p1-F3m1)/(2*d3),0)
	return(f1)
}
#' DMGS equation 3.3, mu1 term
#' @inherit man1f return
#' @inheritParams manf
gev_mu1f=function(alpha,v1,d1,v2,fd2,v3,d3){
	q00=extraDistr::qgev((1-alpha),mu=v1,sigma=v2,xi=v3)
  d2=fd2*v2
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
	F1m1=extraDistr::pgev(q00,mu=v1m1,sigma=v200,xi=v3)
	F1p1=extraDistr::pgev(q00,mu=v1p1,sigma=v200,xi=v3)
# v2 derivatives
	F2m1=extraDistr::pgev(q00,mu=v100,sigma=v2m1,xi=v3)
	F2p1=extraDistr::pgev(q00,mu=v100,sigma=v2p1,xi=v3)
# v3 derivatives
	F3m1=extraDistr::pgev(q00,mu=v100,sigma=v200,xi=v3m1)
	F3p1=extraDistr::pgev(q00,mu=v100,sigma=v200,xi=v3p1)
	mu1=matrix(0,3,length(alpha))
	mu1[1,]=-(F1p1-F1m1)/(2*d1)
	mu1[2,]=-(F2p1-F2m1)/(2*d2)
	mu1[3,]=-(F3p1-F3m1)/(2*d3)
	return(mu1)
}
#' DMGS equation 3.3, f2 term
#' @inherit man2f return
#' @inheritParams manf
gev_f2f=function(y,v1,d1,v2,fd2,v3,d3){
  d2=fd2*v2
	ymax=ifelse(v3<0,v1-2*d1-(v2-2*d2)/(v3-2*d3),Inf)
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
				Fm1=extraDistr::dgev(y,mu=vvm[1],sigma=vvm[2],xi=vvm[3])
				F00=extraDistr::dgev(y,mu=vv0[1],sigma=vv0[2],xi=vv0[3])
				Fp1=extraDistr::dgev(y,mu=vvp[1],sigma=vvp[2],xi=vvp[3])
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
				Fm1m1=extraDistr::dgev(y,mu=vvmm[1],sigma=vvmm[2],xi=vvmm[3])
				Fm1p1=extraDistr::dgev(y,mu=vvmp[1],sigma=vvmp[2],xi=vvmp[3])
				Fp1m1=extraDistr::dgev(y,mu=vvpm[1],sigma=vvpm[2],xi=vvpm[3])
				Fp1p1=extraDistr::dgev(y,mu=vvpp[1],sigma=vvpp[2],xi=vvpp[3])
				f2[i,j,]=ifelse(y<ymax,(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j]),0)
				f2[j,i,]=f2[i,j,]
			}
		}
	}
	return(f2)
}
#' DMGS equation 3.3, mu2 term
#' @inherit man2f return
#' @inheritParams manf
gev_mu2f=function(alpha,v1,d1,v2,fd2,v3,d3){
	q00=extraDistr::qgev((1-alpha),mu=v1,sigma=v2,xi=v3)
  d2=fd2*v2
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
				Fm1=extraDistr::pgev(q00,mu=vvm[1],sigma=vvm[2],xi=vvm[3])
				F00=extraDistr::pgev(q00,mu=vv0[1],sigma=vv0[2],xi=vv0[3])
				Fp1=extraDistr::pgev(q00,mu=vvp[1],sigma=vvp[2],xi=vvp[3])
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
				Fm1m1=extraDistr::pgev(q00,mu=vvmm[1],sigma=vvmm[2],xi=vvmm[3])
				Fm1p1=extraDistr::pgev(q00,mu=vvmp[1],sigma=vvmp[2],xi=vvmp[3])
				Fp1m1=extraDistr::pgev(q00,mu=vvpm[1],sigma=vvpm[2],xi=vvpm[3])
				Fp1p1=extraDistr::pgev(q00,mu=vvpp[1],sigma=vvpp[2],xi=vvpp[3])
				mu2[i,j,]=-(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j])
				mu2[j,i,]=mu2[i,j,]
			}
		}
	}
	return(mu2)
}
#' Derivative of expected information matrix, based on MEV routine gev.infomat
#' @inherit manldd return
#' @inheritParams manf
gev_ggd_mev=function(v1,d1,v2,fd2,v3,d3){
	ggd=array(0,c(3,3,3))
  d2=fd2*v2

	ggm1=gev.infomat(c(v1-d1,v2,v3),dat=c(1),method=c("exp"))
	ggp1=gev.infomat(c(v1+d1,v2,v3),dat=c(1),method=c("exp"))
	ggd[1,,]=(ggp1-ggm1)/(2*d1)

	ggm2=gev.infomat(c(v1,v2-d2,v3),dat=c(1),method=c("exp"))
	ggp2=gev.infomat(c(v1,v2+d2,v3),dat=c(1),method=c("exp"))
	ggd[2,,]=(ggp2-ggm2)/(2*d2)

	ggm3=gev.infomat(c(v1,v2,v3-d3),dat=c(1),method=c("exp"))
	ggp3=gev.infomat(c(v1,v2,v3+d3),dat=c(1),method=c("exp"))
	ggd[3,,]=(ggp3-ggm3)/(2*d3)

  return(ggd)
}
#' Derivative of inverse expected information matrix, based on MEV routine gev.infomat
#' @inherit manlddd return
#' @inheritParams manf
gev_ggid_mev=function(v1,d1,v2,fd2,v3,d3){
	ggid=array(0,c(3,3,3))
  d2=fd2*v2

	ggm1=gev.infomat(c(v1-d1,v2,v3),dat=c(1),method=c("exp"))
	ggp1=gev.infomat(c(v1+d1,v2,v3),dat=c(1),method=c("exp"))
	ggim1=solve(ggm1)
	ggip1=solve(ggp1)
	ggid[1,,]=(ggip1-ggim1)/(2*d1)

	ggm2=gev.infomat(c(v1,v2-d2,v3),dat=c(1),method=c("exp"))
	ggp2=gev.infomat(c(v1,v2+d2,v3),dat=c(1),method=c("exp"))
	ggim2=solve(ggm2)
	ggip2=solve(ggp2)
	ggid[2,,]=(ggip2-ggim2)/(2*d2)

	ggm3=gev.infomat(c(v1,v2,v3-d3),dat=c(1),method=c("exp"))
	ggp3=gev.infomat(c(v1,v2,v3+d3),dat=c(1),method=c("exp"))
	ggim3=solve(ggm3)
	ggip3=solve(ggp3)
	ggid[3,,]=(ggip3-ggim3)/(2*d3)

  return(ggid)
}
#' Analytical Expressions for Predictive Means
#' RHP mean based on the expectation of DMGS equation 2.1
#' @inherit manmeans return
#' @inheritParams manf
gev_means=function(means,ml_params,lddi,lddi_k3,lddd,lddd_k3,
									lambdad_flat,lambdad_rh_mle,
									lambdad_rh_flat,lambdad_jp,lambdad_custom,
									nx,dim=3){

	if(means){
# intro
		eulerconstant=0.57721566490153286060651209008240243104215933593992
		mu=ml_params[1]
		sigma=ml_params[2]
		xi=ml_params[3]

# set up derivative arrays
		meand1=array(0,c(3,1))
		meand2=array(0,c(3,3,1)) #but all zero for gumbel
		meand1_k3=array(0,c(2,1))
		meand2_k3=array(0,c(2,2,1)) #but all zero for gumbel

		if(ml_params[3]==0){
# xi=0 case

# mle
			ml_mean=mu+sigma*eulerconstant
# calculate first derivative array for bayesian xi=0 cases
			meand1[1,1]=1
			meand1[2,1]=eulerconstant
# and copy for the 2d case
			meand1_k3[1,1]=meand1[1,1]
			meand1_k3[2,1]=meand1[2,1]
# meand2 is all zero as initialized

		} else{
# non-gumbel case

			g0=gamma(1-xi)
			g1=g0*digamma(1-xi)
			g2=(trigamma(1-xi)*g0*g0+g1*g1)/g0
# mle
			ml_mean=mu+sigma*(g0-1)/xi
# calculate first derivative array for bayesian xi!=0 cases
			meand1[1,1]=1
			meand1[2,1]=(g0-1)/xi
			meand1[3,1]=(1-g0-xi*g1)*sigma/(xi*xi)
# and copy for the 2d case
			meand1_k3[1,1]=meand1[1,1]
			meand1_k3[2,1]=meand1[2,1]
# calculate second derivative array (only has 1 non-zero term!)
			meand2[2,3,1]=(1-g0-xi*g1)/(xi*xi)
			meand2[3,2,1]=meand2[2,3,1]
			meand2[3,3,1]=sigma*(-2+2*g0+2*xi*g1+xi*xi*g2)/(xi*xi*xi)
		}
# I've realized now that when we integrate over xi, the mean is Inf

#	flat_mean			=ml_mean+dmgs(lddi,lddd,meand1,lambdad_flat,meand2,dim=3)/nx
	flat_mean			=Inf

	rh_ml_mean		=ml_mean+dmgs(lddi_k3,lddd_k3,meand1_k3,lambdad_rh_mle,meand2_k3,dim=2)/nx

#	rh_flat_mean	=ml_mean+dmgs(lddi,lddd,meand1,lambdad_rh_flat,meand2,dim=3)/nx
	rh_flat_mean	=Inf

#	jp_mean				=ml_mean+dmgs(lddi,lddd,meand1,lambdad_jp,meand2,dim=3)/nx
	jp_mean				=Inf

# custom mean
	custom_mean		=Inf

	}else{
		flat_mean="means not selected"
		ml_mean="means not selected"
		rh_ml_mean="means not selected"
		rh_flat_mean="means not selected"
		jp_mean="means not selected"
		custom_mean="means not selected"
	}

# return
	list(ml_mean=ml_mean,flat_mean=flat_mean,rh_ml_mean=rh_ml_mean,
				rh_flat_mean=rh_flat_mean,jp_mean=jp_mean,custom_mean=custom_mean)
}
#' Densities for 5 predictions
#' @inherit mandsub return
#' @inheritParams manf
dgevsub=function(x,y,ics,d1=0.01,fd2=0.01,d3=0.01,customprior,
	minxi,maxxi,extramodels=FALSE,aderivs=TRUE){

		nx=length(x)

		ics=gev_setics(x,ics)
		opt=optim(ics,gev_loglik,x=x,control=list(fnscale=-1))
		v1hat=opt$par[1]
		v2hat=opt$par[2]
		v3hat=opt$par[3]
#don't want to stop the code here, so just reset the value in this case
#		v3hat=min(maxxi,max(minxi,opt$par[3])) #don't want to stop the code here, so just reset the value in this case
# notes on this:
# -I've commented this out now because I think it's ok with the aderivs
# -but I'm not quite sure

		ml_params=c(v1hat,v2hat,v3hat)

# now that I've dropped dmgs, don't think I need this
#		y=fixgevrange(y,v1hat,v2hat,v3hat)

# mle
		ml_pdf=extraDistr::dgev(y,mu=v1hat,sigma=v2hat,xi=v3hat)
		ml_cdf=extraDistr::pgev(y,mu=v1hat,sigma=v2hat,xi=v3hat)

# return
		list(
					ml_params=ml_params,
					ml_pdf=ml_pdf,
					ml_cdf=ml_cdf)
}

