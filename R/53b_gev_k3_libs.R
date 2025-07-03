#' Waic
#' @inherit manwaic return
#' @inheritParams manf
gev_k3_waic=function(waicscores,x,v1hat,d1,v2hat,fd2,kshape,lddi,lddd,lambdad,aderivs){
		if(waicscores){

			if(aderivs) f1f=gev_k3_f1fa(x,v1hat,v2hat,kshape)
			if(!aderivs)f1f=gev_k3_f1f(x,v1hat,d1,v2hat,fd2,kshape)

			if(aderivs) f2f=gev_k3_f2fa(x,v1hat,v2hat,kshape)
			if(!aderivs)f2f=gev_k3_f2f(x,v1hat,d1,v2hat,fd2,kshape)

			fhatx=extraDistr::dgev(x,mu=v1hat,sigma=v2hat,xi=kshape)
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
gev_k3_logf=function(params,x,kshape){
	l=params[1]
	s=pmax(params[2],.Machine$double.eps)
	logf=sum(extraDistr::dgev(x,mu=l,sigma=s,xi=kshape,log=TRUE))-log(s)
	return(logf)
}
#'  log-likelihood function
#' @inherit manloglik return
#' @inheritParams manf
gev_k3_loglik=function(vv,x,kshape){
	n=length(x)
	loglik=sum(extraDistr::dgev(x,mu=vv[1],sigma=max(vv[2],.Machine$double.eps),xi=kshape,log=TRUE))
	return(loglik)
}
#' One component of the second derivative of the normalized log-likelihood
#' @inherit manlnn return
#' @inheritParams manf
gev_k3_lmn=function(x,v1,d1,v2,fd2,kshape,mm,nn){
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
			lmn[i]=sum(extraDistr::dgev(x,mu=vvd[1],sigma=vvd[2],xi=kshape,log=TRUE))/nx
		}
		dld=(lmn[1]-lmn[2]-lmn[3]+lmn[4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:3){
			for (j in 1:2){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[i]=sum(extraDistr::dgev(x,mu=vvd[1],sigma=vvd[2],xi=kshape,log=TRUE))/nx
		}
		dld=(lmn[1]-2*lmn[2]+lmn[3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' Second derivative matrix of the normalized log-likelihood
#' @inherit manldd return
#' @inheritParams manf
gev_k3_ldd=function(x,v1,d1,v2,fd2,kshape){
	nx=length(x)
	ldd=matrix(0,2,2)
	for (i in 1:2){
		for (j in i:2){
			ldd[i,j]=gev_k3_lmn(x,v1,d1,v2,fd2,kshape,i,j)
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
gev_k3_lmnp=function(x,v1,d1,v2,fd2,kshape,mm,nn,rr){
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
			lmn[i]=sum(extraDistr::dgev(x,mu=vvd[1],sigma=vvd[2],xi=kshape,log=TRUE))/nx
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
			lmn[i]=sum(extraDistr::dgev(x,mu=vvd[1],sigma=vvd[2],xi=kshape,log=TRUE))/nx
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
			lmn[i]=sum(extraDistr::dgev(x,mu=vvd[1],sigma=vvd[2],xi=kshape,log=TRUE))/nx
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
gev_k3_lddd=function(x,v1,d1,v2,fd2,kshape){
# calculate the unique values
	lddd=array(0,c(2,2,2))
	for (i in 1:2){
		for (j in i:2){
			for (k in j:2){
				lddd[i,j,k]=gev_k3_lmnp(x,v1,d1,v2,fd2,kshape,i,j,k)
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
gev_k3_f1f=function(y,v1,d1,v2,fd2,kshape){
	d2=fd2*v2
	v3=kshape
	ymax=ifelse(v3<0,v1-2*d1-(v2-2*d2)/v3,Inf)
# v1 stuff
	v1m1=v1-1*d1
	v100=v1+0*d1
	v1p1=v1+1*d1
# v2 stuff
	v2m1=v2-1*d2
	v200=v2+0*d2
	v2p1=v2+1*d2
# v1 derivatives
	F1m1=extraDistr::dgev(y,mu=v1m1,sigma=v200,xi=kshape)
	F1p1=extraDistr::dgev(y,mu=v1p1,sigma=v200,xi=kshape)
# v2 derivatives
	F2m1=extraDistr::dgev(y,mu=v100,sigma=v2m1,xi=kshape)
	F2p1=extraDistr::dgev(y,mu=v100,sigma=v2p1,xi=kshape)
	f1=matrix(0,2,length(y))
	f1[1,]=ifelse(y<ymax,(F1p1-F1m1)/(2*d1),0)
	f1[2,]=ifelse(y<ymax,(F2p1-F2m1)/(2*d2),0)
	return(f1)
}
#' DMGS equation 3.3, mu1 term
#' @inherit man1f return
#' @inheritParams manf
gev_k3_mu1f=function(alpha,v1,d1,v2,fd2,kshape){
	q00=extraDistr::qgev((1-alpha),mu=v1,sigma=v2,xi=kshape)
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
	F1m1=extraDistr::pgev(q00,mu=v1m1,sigma=v200,xi=kshape)
	F1p1=extraDistr::pgev(q00,mu=v1p1,sigma=v200,xi=kshape)
# v2 derivatives
	F2m1=extraDistr::pgev(q00,mu=v100,sigma=v2m1,xi=kshape)
	F2p1=extraDistr::pgev(q00,mu=v100,sigma=v2p1,xi=kshape)
	mu1=matrix(0,2,length(alpha))
	mu1[1,]=-(F1p1-F1m1)/(2*d1)
	mu1[2,]=-(F2p1-F2m1)/(2*d2)
	return(mu1)
}
#' DMGS equation 3.3, f2 term
#' @inherit man2f return
#' @inheritParams manf
gev_k3_f2f=function(y,v1,d1,v2,fd2,kshape){
	d2=fd2*v2
	v3=kshape
	ymax=ifelse(v3<0,v1-2*d1-(v2-2*d2)/v3,Inf)
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
	F1m2=extraDistr::dgev(y,mu=v1m2,sigma=v200,xi=kshape)
	F1m1=extraDistr::dgev(y,mu=v1m1,sigma=v200,xi=kshape)
	F100=extraDistr::dgev(y,mu=v100,sigma=v200,xi=kshape)
	F1p1=extraDistr::dgev(y,mu=v1p1,sigma=v200,xi=kshape)
	F1p2=extraDistr::dgev(y,mu=v1p2,sigma=v200,xi=kshape)
# v2 derivative
	F2m2=extraDistr::dgev(y,mu=v100,sigma=v2m2,xi=kshape)
	F2m1=extraDistr::dgev(y,mu=v100,sigma=v2m1,xi=kshape)
	F200=extraDistr::dgev(y,mu=v100,sigma=v200,xi=kshape)
	F2p1=extraDistr::dgev(y,mu=v100,sigma=v2p1,xi=kshape)
	F2p2=extraDistr::dgev(y,mu=v100,sigma=v2p2,xi=kshape)
# cross derivative
	Fcm1m1=extraDistr::dgev(y,mu=v1m1,sigma=v2m1,xi=kshape)
	Fcm1p1=extraDistr::dgev(y,mu=v1m1,sigma=v2p1,xi=kshape)
	Fcp1m1=extraDistr::dgev(y,mu=v1p1,sigma=v2m1,xi=kshape)
	Fcp1p1=extraDistr::dgev(y,mu=v1p1,sigma=v2p1,xi=kshape)
	f2=array(0,c(2,2,length(y)))
	f2[1,1,]=ifelse(y<ymax,(F1p1-2*F100+F1m1)/(d1*d1),0)
	f2[2,2,]=ifelse(y<ymax,(F2p1-2*F200+F2m1)/(d2*d2),0)
	f2[1,2,]=ifelse(y<ymax,(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2),0)
	# copy
	f2[2,1,]=f2[1,2,]
	return(f2)
}
#' DMGS equation 3.3, mu2 term
#' @inherit man2f return
#' @inheritParams manf
gev_k3_mu2f=function(alpha,v1,d1,v2,fd2,kshape){
	q00=extraDistr::qgev((1-alpha),mu=v1,sigma=v2,xi=kshape)
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
	F1m2=extraDistr::pgev(q00,mu=v1m2,sigma=v200,xi=kshape)
	F1m1=extraDistr::pgev(q00,mu=v1m1,sigma=v200,xi=kshape)
	F100=extraDistr::pgev(q00,mu=v100,sigma=v200,xi=kshape)
	F1p1=extraDistr::pgev(q00,mu=v1p1,sigma=v200,xi=kshape)
	F1p2=extraDistr::pgev(q00,mu=v1p2,sigma=v200,xi=kshape)
# v2 derivative
	F2m2=extraDistr::pgev(q00,mu=v100,sigma=v2m2,xi=kshape)
	F2m1=extraDistr::pgev(q00,mu=v100,sigma=v2m1,xi=kshape)
	F200=extraDistr::pgev(q00,mu=v100,sigma=v200,xi=kshape)
	F2p1=extraDistr::pgev(q00,mu=v100,sigma=v2p1,xi=kshape)
	F2p2=extraDistr::pgev(q00,mu=v100,sigma=v2p2,xi=kshape)
# cross derivative
	Fcm1m1=extraDistr::pgev(q00,mu=v1m1,sigma=v2m1,xi=kshape)
	Fcm1p1=extraDistr::pgev(q00,mu=v1m1,sigma=v2p1,xi=kshape)
	Fcp1m1=extraDistr::pgev(q00,mu=v1p1,sigma=v2m1,xi=kshape)
	Fcp1p1=extraDistr::pgev(q00,mu=v1p1,sigma=v2p1,xi=kshape)
	mu2[1,1,]=-(F1p1-2*F100+F1m1)/(d1*d1)
	mu2[2,2,]=-(F2p1-2*F200+F2m1)/(d2*d2)
	mu2[1,2,]=-(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2)
	# copy
	mu2[2,1,]=mu2[1,2,]
	return(mu2)
}
#' MLE and RHP means
#' @inherit manmeans return
#' @inheritParams manf
gev_k3_means=function(means,ml_params,lddi,lddd,lambdad_rhp,nx,dim=2,kshape){

	if(means){
# intro
		eulerconstant=0.57721566490153286060651209008240243104215933593992

		if(kshape==0){
# gumbel case
			means=gumbel_means(means,ml_params,lddi,lddd,lambdad_rhp,nx,dim=2)
			ml_mean=means$ml_mean
			rh_mean=means$rh_mean
		} else{
# non-gumbel case

# -ml
			g1=gamma(1-kshape)
			ml_mean=ml_params[1]+ml_params[2]*(g1-1)/kshape
# -rhp
			meand1=array(0,c(2,1))
			meand1[1,1]=1
			meand1[2,1]=(g1-1)/kshape
			meand2=array(0,c(2,2,1)) #and they are all zero it seems
			dmean=dmgs(lddi,lddd,meand1,lambdad_rhp,meand2,dim=2)
			rh_mean=ml_mean+dmean/nx
		}
	}else{
		ml_mean="means not selected"
		rh_mean="means not selected"
	}

	list(ml_mean=ml_mean,rh_mean=rh_mean)
}
#' Densities from MLE and RHP
#' @inherit mandsub return
#' @inheritParams manf
dgev_k3sub=function(x,y,d1=0.01,fd2=0.01,kshape,aderivs=TRUE){

		nx=length(x)

		v1start=mean(x)
		v2start=sd(x)
		opt=optim(c(v1start,v2start),gev_k3_loglik,x=x,kshape=kshape,
			control=list(fnscale=-1))
		v1hat=opt$par[1]
		v2hat=opt$par[2]
		ml_params=c(v1hat,v2hat)

		y=fixgevrange(y,v1hat,v2hat,kshape)

# mle
		ml_pdf=extraDistr::dgev(y,mu=v1hat,sigma=v2hat,xi=kshape)
		ml_cdf=extraDistr::pgev(y,mu=v1hat,sigma=v2hat,xi=kshape)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					ml_cdf=ml_cdf)
}



