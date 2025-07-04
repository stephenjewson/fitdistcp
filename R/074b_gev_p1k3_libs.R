#' Waic
#' @inherit manwaic return
#' @inheritParams manf
gev_p1k3_waic=function(waicscores,x,t,v1hat,d1,v2hat,d2,v3hat,fd3,kshape,
	lddi,lddd,lambdad,aderivs){
		if(waicscores){
			if(aderivs){
				f1f=gev_p1k3_f1fa(x,t,v1hat,v2hat,v3hat,kshape=kshape)
				f2f=gev_p1k3_f2fa(x,t,v1hat,v2hat,v3hat,kshape=kshape)
			} else {
				f1f=gev_p1k3_f1f(x,t,v1hat,d1,v2hat,d2,v3hat,fd3,kshape=kshape)
				f2f=gev_p1k3_f2f(x,t,v1hat,d1,v2hat,d2,v3hat,fd3,kshape=kshape)
			}
			fhatx=dgev_p1k3(x,t,ymn=v1hat,slope=v2hat,sigma=v3hat,log=FALSE,kshape=kshape)
			waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad,f2f,dim=3)
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
gev_p1k3_predictordata=function(predictordata,x,t,t0,params,kshape){
	if(predictordata){
#
# calculate the probabilities of the data using the fited model
#
		a=params[1]
		b=params[2]
		s=params[3]
		mu=a+b*t
		px=pgev(x,mu=mu,sigma=s,xi=kshape)
#
# calculate the quantiles for those probabilities at t0
#
		mu0=a+b*t0
		qx=qgev(px,mu=mu0,sigma=s,xi=kshape)
	} else{
		predictedparameter="predictordata not selected"
		adjustedx="predictordata not selected"
	}

	list(predictedparameter=mu,adjustedx=qx)
}
#' Logf for RUST
#' @inherit manlogf return
#' @inheritParams manf
gev_p1k3_logf=function(params,x,t,kshape){
#	a=params[1]
#	b=params[2]
#	s=params[3]
#	mu=a+b*t
#	if(s>0){
#		logf=sum(dgev(x,mu=mu,sigma=s,xi=kshape,log=TRUE))-log(s)
#	}else{
#		logf=-Inf
#	}
	a=params[1]
	b=params[2]
	s=pmax(params[3],.Machine$double.eps)
	mu=a+b*t
	logf=sum(dgev(x,mu=mu,sigma=s,xi=kshape,log=TRUE))-log(s)
	return(logf)
}
#' GEV-with-known-shape-with-p1  observed log-likelihood function
#' @inherit manloglik return
#' @inheritParams	manf
gev_p1k3_loglik=function(vv,x,t,kshape){
	n=length(x)
	mu=vv[1]+vv[2]*t
	loglik=sum(dgev(x,mu=mu,sigma=max(vv[3],.Machine$double.eps),log=TRUE,xi=kshape))
	if(loglik==Inf)loglik=100000000 #just a large negative number
	return(loglik)
}
#' GEV-with-known-shape-with-p1 quantile function
#' @inherit manvector return
#' @inheritParams	manf
qgev_p1k3=function(p,t0,ymn,slope,sigma,kshape){

	return(qgev(p,mu=(ymn+slope*t0),sigma=sigma,xi=kshape))

}
#' GEV-with-known-shape-with-p1 density function
#' @inherit manvector return
#' @inheritParams	manf
dgev_p1k3=function(x,t0,ymn,slope,sigma,log=FALSE,kshape){

	return(dgev(x,mu=(ymn+slope*t0),sigma=sigma,log=log,xi=kshape))

}
#' GEV-with-known-shape-with-p1 distribution function
#' @inherit manvector return
#' @inheritParams	manf
pgev_p1k3=function(x,t0,ymn,slope,sigma,kshape){

	return(pgev(x,mu=(ymn+slope*t0),sigma=sigma,xi=kshape))
}
#' One component of the second derivative of the normalized log-likelihood
#' @inherit manlnn return
#' @inheritParams manf
gev_p1k3_lmn=function(x,t,v1,d1,v2,d2,v3,fd3,kshape,mm,nn){
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
			lmn[i]=sum(dgev_p1k3(x,t,ymn=vvd[1],slope=vvd[2],sigma=vvd[3],kshape=kshape,log=TRUE))/nx
		}
		dld=(lmn[1]-lmn[2]-lmn[3]+lmn[4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:3){
			for (j in 1:3){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[i]=sum(dgev_p1k3(x,t,ymn=vvd[1],slope=vvd[2],sigma=vvd[3],kshape=kshape,log=TRUE))/nx
		}
		dld=(lmn[1]-2*lmn[2]+lmn[3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' Second derivative matrix of the normalized log-likelihood
#' @inherit manldd return
#' @inheritParams	manf
gev_p1k3_ldd=function(x,t,v1,d1,v2,d2,v3,fd3,kshape){
	ldd=matrix(0,3,3)
	for (i in 1:3){
		for (j in i:3){
			ldd[i,j]=gev_p1k3_lmn(x,t,v1,d1,v2,d2,v3,fd3,kshape,i,j)
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
gev_p1k3_lmnp=function(x,t,v1,d1,v2,d2,v3,fd3,kshape,mm,nn,rr){
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
			lmn[i]=sum(dgev_p1k3(x,t,ymn=vvd[1],slope=vvd[2],sigma=vvd[3],kshape=kshape,log=TRUE))/nx
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
			lmn[i]=sum(dgev_p1k3(x,t,ymn=vvd[1],slope=vvd[2],sigma=vvd[3],kshape=kshape,log=TRUE))/nx
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
			lmn[i]=sum(dgev_p1k3(x,t,ymn=vvd[1],slope=vvd[2],sigma=vvd[3],kshape=kshape,log=TRUE))/nx
		}
		dld1=(lmn[3]-2*lmn[2]+lmn[1])/(dd[m2]*dd[m2])
		dld2=(lmn[6]-2*lmn[5]+lmn[4])/(dd[m2]*dd[m2])
		dld=(dld2-dld1)/(2*dd[n2])
	}
	return(dld)
}#' Third derivative tensor of the normalized log-likelihood
#' @inherit manlddd return
#' @inheritParams	manf
gev_p1k3_lddd=function(x,t,v1,d1,v2,d2,v3,fd3,kshape){
	lddd=array(0,c(3,3,3))
	for (i in 1:3){
		for (j in i:3){
			for (k in j:3){
				lddd[i,j,k]=gev_p1k3_lmnp(x,t,v1,d1,v2,d2,v3,fd3,kshape,i,j,k)
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
#' @inherit man1f return
#' @inheritParams	manf
gev_p1k3_f1f=function(y,t0,v1,d1,v2,d2,v3,fd3,kshape){
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
	F1m1=dgev_p1k3(y,t0,ymn=v1m1,slope=v200,sigma=v3,kshape=kshape)
	F1p1=dgev_p1k3(y,t0,ymn=v1p1,slope=v200,sigma=v3,kshape=kshape)
# v2 derivatives
	F2m1=dgev_p1k3(y,t0,ymn=v100,slope=v2m1,sigma=v3,kshape=kshape)
	F2p1=dgev_p1k3(y,t0,ymn=v100,slope=v2p1,sigma=v3,kshape=kshape)
# v3 derivatives
	F3m1=dgev_p1k3(y,t0,ymn=v100,slope=v200,sigma=v3m1,kshape=kshape)
	F3p1=dgev_p1k3(y,t0,ymn=v100,slope=v200,sigma=v3p1,kshape=kshape)
	f1=matrix(0,3,length(y))
	f1[1,]=(F1p1-F1m1)/(2*d1)
	f1[2,]=(F2p1-F2m1)/(2*d2)
	f1[3,]=(F3p1-F3m1)/(2*d3)
	return(f1)
}
#' DMGS equation 3.3, mu1 term
#' @inherit man1f return
#' @inheritParams	manf
gev_p1k3_mu1f=function(alpha,t0,v1,d1,v2,d2,v3,fd3,kshape){
	q00=qgev_p1k3((1-alpha),t0,ymn=v1,slope=v2,sigma=v3,kshape=kshape)
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
	F1m1=pgev_p1k3(q00,t0,ymn=v1m1,slope=v200,sigma=v3,kshape=kshape)
	F1p1=pgev_p1k3(q00,t0,ymn=v1p1,slope=v200,sigma=v3,kshape=kshape)
# v2 derivatives
	F2m1=pgev_p1k3(q00,t0,ymn=v100,slope=v2m1,sigma=v3,kshape=kshape)
	F2p1=pgev_p1k3(q00,t0,ymn=v100,slope=v2p1,sigma=v3,kshape=kshape)
# v3 derivatives
	F3m1=pgev_p1k3(q00,t0,ymn=v100,slope=v200,sigma=v3m1,kshape=kshape)
	F3p1=pgev_p1k3(q00,t0,ymn=v100,slope=v200,sigma=v3p1,kshape=kshape)
	mu1=matrix(0,3,length(alpha))
	mu1[1,]=-(F1p1-F1m1)/(2*d1)
	mu1[2,]=-(F2p1-F2m1)/(2*d2)
	mu1[3,]=-(F3p1-F3m1)/(2*d3)
	return(mu1)
}
#' DMGS equation 2.1, f2 term
#' @inherit man2f return
#' @inheritParams	manf
gev_p1k3_f2f=function(y,t0,v1,d1,v2,d2,v3,fd3,kshape){
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
				Fm1=dgev_p1k3(y,t0,ymn=vvm[1],slope=vvm[2],sigma=vvm[3],kshape=kshape)
				F00=dgev_p1k3(y,t0,ymn=vv0[1],slope=vv0[2],sigma=vv0[3],kshape=kshape)
				Fp1=dgev_p1k3(y,t0,ymn=vvp[1],slope=vvp[2],sigma=vvp[3],kshape=kshape)
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
				Fm1m1=dgev_p1k3(y,t0,ymn=vvmm[1],slope=vvmm[2],sigma=vvmm[3],kshape=kshape)
				Fm1p1=dgev_p1k3(y,t0,ymn=vvmp[1],slope=vvmp[2],sigma=vvmp[3],kshape=kshape)
				Fp1m1=dgev_p1k3(y,t0,ymn=vvpm[1],slope=vvpm[2],sigma=vvpm[3],kshape=kshape)
				Fp1p1=dgev_p1k3(y,t0,ymn=vvpp[1],slope=vvpp[2],sigma=vvpp[3],kshape=kshape)
				f2[i,j,]=(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j])
				f2[j,i,]=f2[i,j,]
			}
		}
	}
	return(f2)
}
#' DMGS equation 3.3, mu2 term
#' @inherit man2f return
#' @inheritParams	manf
gev_p1k3_mu2f=function(alpha,t0,v1,d1,v2,d2,v3,fd3,kshape){
	q00=qgev_p1k3((1-alpha),t0,ymn=v1,slope=v2,sigma=v3,kshape=kshape)
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
				Fm1=pgev_p1k3(q00,t0,ymn=vvm[1],slope=vvm[2],sigma=vvm[3],kshape=kshape)
				F00=pgev_p1k3(q00,t0,ymn=vv0[1],slope=vv0[2],sigma=vv0[3],kshape=kshape)
				Fp1=pgev_p1k3(q00,t0,ymn=vvp[1],slope=vvp[2],sigma=vvp[3],kshape=kshape)
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
				Fm1m1=pgev_p1k3(q00,t0,ymn=vvmm[1],slope=vvmm[2],sigma=vvmm[3],kshape=kshape)
				Fm1p1=pgev_p1k3(q00,t0,ymn=vvmp[1],slope=vvmp[2],sigma=vvmp[3],kshape=kshape)
				Fp1m1=pgev_p1k3(q00,t0,ymn=vvpm[1],slope=vvpm[2],sigma=vvpm[3],kshape=kshape)
				Fp1p1=pgev_p1k3(q00,t0,ymn=vvpp[1],slope=vvpp[2],sigma=vvpp[3],kshape=kshape)
				mu2[i,j,]=-(Fp1p1-Fm1p1-Fp1m1+Fm1m1)/(4*dd[i]*dd[j])
				mu2[j,i,]=mu2[i,j,]
			}
		}
	}
	return(mu2)
}
#' Analytical expressions for Predictive Means
#' RHP mean based on the expectation of DMGS equation 2.1
#' @inherit manmeans return
#' @inheritParams manf
gev_p1k3_means=function(means,t0,ml_params,kshape,nx){

	if(means){
# intro
		eulerconstant=0.57721566490153286060651209008240243104215933593992
		ymn=ml_params[1]
		slope=ml_params[2]
		sigma=ml_params[3]
		xi=kshape

		if(xi==0){
# xi=0 case
			ml_mean=ymn+slope*t0+sigma*eulerconstant
		} else{
# xi!=0 case
			g0=gamma(1-xi)
			g1=g0*digamma(1-xi)
			g2=(trigamma(1-xi)*g0*g0+g1*g1)/g0
			ml_mean=ymn+slope*t0+sigma*(g0-1)/xi
		}
# return
		pu_mean=Inf
	}else{
		pu_mean="means not selected"
		ml_mean="means not selected"
	}
	list(ml_mean=ml_mean,pu_mean=pu_mean)
}
#' Densities from MLE and RHP
#' @inherit mandsub return
#' @inheritParams	manf
dgev_p1k3sub=function(x,t,y,t0,d1,d2,fd3,kshape,aderivs=TRUE){

		nx=length(x)

		lm=lm(x~t)
		v1start=lm$coefficients[1]
		v2start=lm$coefficients[2]
		xhat=v1start+v2start*t
		v3start=100000000
		opt1=optim(c(v1start,v2start,v3start),gev_p1k3_loglik,x=x,t=t,
			kshape=kshape,control=list(fnscale=-1))
		v1hat=opt1$par[1]
		v2hat=opt1$par[2]
		v3hat=opt1$par[3]
		ml_params=c(v1hat,v2hat,v3hat)

# ml
		muhat0=v1hat+v2hat*t0
		y=fixgevrange(y,muhat0,v3hat,kshape)
		ml_pdf=dgev(y,mu=muhat0,sigma=v3hat,xi=kshape)
		ml_cdf=pgev(y,mu=muhat0,sigma=v3hat,xi=kshape)


# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					ml_cdf=ml_cdf)
}
