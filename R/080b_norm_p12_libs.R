#' Waic
#' @inheritParams manf
norm_p12_waic=function(waicscores,x,t1,t2,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4,lddi,lddd,lambdad){
		if(waicscores){
			f1f=norm_p12_f1f(x,t1,t2,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4)
			f2f=norm_p12_f2f(x,t1,t2,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4)
			fhatx=dnorm_p12(x,t1,t2,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat,log=FALSE)
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
norm_p12_predictordata=function(predictordata,x,t1,t2,t10,t20,params){
	if(predictordata){
#
# calculate the probabilities of the data using the fitted model
#
		a=params[1]
		b=params[2]
		g=params[3]
		d=params[4]
		mu=a+b*t1
		sg=exp(g+d*t2)
		px=pnorm(x,mean=mu,sd=sg)
#
# calculate the quantiles for those probabilities at t0
#
		mu0=a+b*t10
		sg0=exp(g+d*t20)
		qx=qnorm(px,mean=mu0,sd=sg0)
	} else{
		predictedparameter="predictordata not selected"
		adjustedx="predictordata not selected"
	}

	list(predictedparameter=list(mu=mu,sg=sg),adjustedx=qx)
}
#' Logf for RUST
#' @inheritParams manf
norm_p12_logf=function(params,x,t1,t2,nonnegslopesonly=FALSE){
#	a=params[1]
#	b=params[2]
#	g=params[3]
#	d=params[4]
#	mu=a+b*t1
#	sg=exp(g+d*t2)
#	logf=sum(dnorm(x,mean=mu,sd=sg,log=TRUE))
	a=params[1]
	b=params[2]
	g=params[3]
	d=params[4]
	mu=a+b*t1
	if(nonnegslopesonly&&((b<0)||(d<0))){
		logf=-Inf
	} else{
		sg=pmax(exp(g+d*t2),.Machine$double.eps)
		logf=sum(dnorm(x,mean=mu,sd=sg,log=TRUE))
	}


	return(logf)
}
#' Set initial conditions
#' @inheritParams manf
norm_p12_setics=function(x,t1,t2,ics){
	nx=length(x)
	if((ics[1]==0)&&(ics[2]==0)&&(ics[3]==0)&&(ics[4]==0)){
		lm=lm(x~t1)
		ics[1]=lm$coefficients[1]
		ics[2]=lm$coefficients[2]
#		xhat=ics[1]+ics[2]*t1
#		ics[3]=1
		ics[3]=log(sqrt(mean((lm$residuals)**2))) #try this to avoid the small number of maxlik fails
		ics[4]=0
	}
	return(ics)
}
#'  observed log-likelihood function
#' @inheritParams manf
norm_p12_loglik=function(vv,x,t1,t2){
	n=length(x)
	mu=vv[1]+vv[2]*t1 #so mean is a vector, just like x
	sg=exp(vv[3]+vv[4]*t2)
	loglik=sum(dnorm(x,mean=mu,sd=sg,log=TRUE))
	return(loglik)
}
#' Check MLE
#' @inheritParams manf
norm_p12_checkmle=function(ml_params){
	v1hat=ml_params[1]
	v2hat=ml_params[2]
	v3hat=ml_params[3]
	v4hat=ml_params[4]
	if(is.na(v1hat))stop()
	if(is.na(v2hat))stop()
	if(is.na(v3hat))stop()
	if(is.na(v4hat))stop()
}
#' Normal-with-p12: Quantile function
#' @inheritParams manf
qnorm_p12=function(p,t10,t20,ymn,slope,sigma1,sigma2){

	mu=(ymn+slope*t10)
	sg=exp(sigma1+sigma2*t20)
	return(qnorm(p,mean=mu,sd=sg))

}
#' Normal-with-p12: Density function
#' @inheritParams manf
dnorm_p12=function(x,t10,t20,ymn,slope,sigma1,sigma2,log=FALSE){

	mu=(ymn+slope*t10)
	sg=exp(sigma1+sigma2*t20)
	return(dnorm(x,mean=mu,sd=sg,log=log))

}
#' Normal-with-p12: Distribution function
#' @inheritParams manf
pnorm_p12=function(y,t10,t20,ymn,slope,sigma1,sigma2){

	mu=(ymn+slope*t10)
	sg=exp(sigma1+sigma2*t20)
	return(pnorm(y,mean=mu,sd=sg))

}
#' Bootstrap
#' @inheritParams manf
norm_p12_boot=function(x,t1,t2,n){

	sim_vals=matrix(0,n,4)
	ics=c(0,0,0,0)
	base=c(1:length(x))
	for (i in 1:n){
		index=sample(base,replace=TRUE)
		bx=x[index]
		bt1=t1[index]
		bt2=t2[index]
		ics=norm_p12_setics(bx,bt1,bt2,ics)
		opt1=optim(ics,norm_p12_loglik,x=bx,t1=bt1,t2=bt2,control=list(fnscale=-1))
		sim_vals[i,]=opt1$par[]
	}

	return(list(sim_vals=sim_vals))

}
#' One component of the second derivative of the normalized log-likelihood
#' @inheritParams manf
norm_p12_lmn=function(x,t1,t2,v1,d1,v2,d2,v3,d3,v4,d4,mm,nn){
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
			lmn[i]=sum(dnorm_p12(x,t1,t2,ymn=vvd[1],slope=vvd[2],sigma1=vvd[3],sigma2=vvd[4],log=TRUE))/nx
		}
		dld=(lmn[1]-lmn[2]-lmn[3]+lmn[4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:3){
			for (j in 1:4){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[i]=sum(dnorm_p12(x,t1,t2,ymn=vvd[1],slope=vvd[2],sigma1=vvd[3],sigma2=vvd[4],log=TRUE))/nx
		}
		dld=(lmn[1]-2*lmn[2]+lmn[3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' Second derivative matrix of the normalized log-likelihood
#' @inheritParams manf
norm_p12_ldd=function(x,t1,t2,v1,d1,v2,d2,v3,d3,v4,d4){
	ldd=matrix(0,4,4)
	for (i in 1:4){
		for (j in i:4){
			ldd[i,j]=norm_p12_lmn(x,t1,t2,v1,d1,v2,d2,v3,d3,v4,d4,i,j)
		}
	}
	for (i in 4:2){
		for (j in 1:(i-1)){
			ldd[i,j]=ldd[j,i]
		}
	}
	return(ldd)
}
#' One component of the second derivative of the normalized log-likelihood
#' @inheritParams manf
norm_p12_lmnp=function(x,t1,t2,v1,d1,v2,d2,v3,d3,v4,d4,mm,nn,rr){
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
			lmn[i]=sum(dnorm_p12(x,t1,t2,ymn=vvd[1],slope=vvd[2],sigma1=vvd[3],sigma2=vvd[4],log=TRUE))/nx
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
			lmn[i]=sum(dnorm_p12(x,t1,t2,ymn=vvd[1],slope=vvd[2],sigma1=vvd[3],sigma2=vvd[4],log=TRUE))/nx
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
			lmn[i]=sum(dnorm_p12(x,t1,t2,ymn=vvd[1],slope=vvd[2],sigma1=vvd[3],sigma2=vvd[4],log=TRUE))/nx
		}
		dld1=(lmn[3]-2*lmn[2]+lmn[1])/(dd[m2]*dd[m2])
		dld2=(lmn[6]-2*lmn[5]+lmn[4])/(dd[m2]*dd[m2])
		dld=(dld2-dld1)/(2*dd[n2])
	}
	return(dld)
}
#' Third derivative tensor of the normalized log-likelihood, with fixed shape parameter
#' @inheritParams manf
norm_p12_lddd=function(x,t1,t2,v1,d1,v2,d2,v3,d3,v4,d4){
	lddd=array(0,c(4,4,4))
	for (i in 1:4){
		for (j in i:4){
			for (k in j:4){
				lddd[i,j,k]=norm_p12_lmnp(x,t1,t2,v1,d1,v2,d2,v3,d3,v4,d4,i,j,k)
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
	return(lddd)
}
#' DMGS equation 2.1, f1 term
#' @inheritParams manf
norm_p12_f1f=function(y,t10,t20,v1,d1,v2,d2,v3,d3,v4,d4){
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
# v4 stuff
	v4m1=v4-1*d4
	v400=v4+0*d4
	v4p1=v4+1*d4
# v1 derivatives
	F1m1=dnorm_p12(y,t10,t20,ymn=v1m1,slope=v200,sigma1=v300,sigma2=v400)
	F1p1=dnorm_p12(y,t10,t20,ymn=v1p1,slope=v200,sigma1=v300,sigma2=v400)
# v2 derivatives
	F2m1=dnorm_p12(y,t10,t20,ymn=v100,slope=v2m1,sigma1=v300,sigma2=v400)
	F2p1=dnorm_p12(y,t10,t20,ymn=v100,slope=v2p1,sigma1=v300,sigma2=v400)
# v3 derivatives
	F3m1=dnorm_p12(y,t10,t20,ymn=v100,slope=v200,sigma1=v3m1,sigma2=v400)
	F3p1=dnorm_p12(y,t10,t20,ymn=v100,slope=v200,sigma1=v3p1,sigma2=v400)
# v4 derivatives
	F4m1=dnorm_p12(y,t10,t20,ymn=v100,slope=v200,sigma1=v300,sigma2=v4m1)
	F4p1=dnorm_p12(y,t10,t20,ymn=v100,slope=v200,sigma1=v300,sigma2=v4p1)
	f1=matrix(0,4,length(y))
	f1[1,]=(F1p1-F1m1)/(2*d1)
	f1[2,]=(F2p1-F2m1)/(2*d2)
	f1[3,]=(F3p1-F3m1)/(2*d3)
	f1[4,]=(F4p1-F4m1)/(2*d4)
	return(f1)
}
#' DMGS equation 2.1, p1 term
#' @inheritParams manf
norm_p12_p1f=function(y,t10,t20,v1,d1,v2,d2,v3,d3,v4,d4){
	#d3=d3*v3
	#ymax=ifelse(v4<0,v1+v2*t0-2*(d1+d2*t0)-(v3-2*d3)/(v4-2*d4),Inf)
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
# v4 stuff
	v4m1=v4-1*d4
	v400=v4+0*d4
	v4p1=v4+1*d4
# v1 derivatives
	F1m1=pnorm_p12(y,t10,t20,ymn=v1m1,slope=v200,sigma1=v3,sigma2=v4)
	F1p1=pnorm_p12(y,t10,t20,ymn=v1p1,slope=v200,sigma1=v3,sigma2=v4)
# v2 derivatives
	F2m1=pnorm_p12(y,t10,t20,ymn=v100,slope=v2m1,sigma1=v3,sigma2=v4)
	F2p1=pnorm_p12(y,t10,t20,ymn=v100,slope=v2p1,sigma1=v3,sigma2=v4)
# v3 derivatives
	F3m1=pnorm_p12(y,t10,t20,ymn=v100,slope=v200,sigma1=v3m1,sigma2=v4)
	F3p1=pnorm_p12(y,t10,t20,ymn=v100,slope=v200,sigma1=v3p1,sigma2=v4)
# v4 derivatives
	F4m1=pnorm_p12(y,t10,t20,ymn=v100,slope=v200,sigma1=v3,sigma2=v4m1)
	F4p1=pnorm_p12(y,t10,t20,ymn=v100,slope=v200,sigma1=v3,sigma2=v4p1)
	p1=matrix(0,4,length(y))
	p1[1,]=(F1p1-F1m1)/(2*d1)
	p1[2,]=(F2p1-F2m1)/(2*d2)
	p1[3,]=(F3p1-F3m1)/(2*d3)
	p1[4,]=(F4p1-F4m1)/(2*d4)
	return(p1)
}
#' Normal-with-p12: DMGS equation 3.3 mu1 term
#' @inheritParams manf
norm_p12_mu1f=function(alpha,t10,t20,v1,d1,v2,d2,v3,d3,v4,d4){
	q00=qnorm_p12((1-alpha),t10,t20,ymn=v1,slope=v2,sigma1=v3,sigma2=v4)
	#d3=d3*v3
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
# v4 stuff
	v4m1=v4-1*d4
	v400=v4+0*d4
	v4p1=v4+1*d4
# v1 derivatives
	F1m1=pnorm_p12(q00,t10,t20,ymn=v1m1,slope=v200,sigma1=v3,sigma2=v4)
	F1p1=pnorm_p12(q00,t10,t20,ymn=v1p1,slope=v200,sigma1=v3,sigma2=v4)
# v2 derivatives
	F2m1=pnorm_p12(q00,t10,t20,ymn=v100,slope=v2m1,sigma1=v3,sigma2=v4)
	F2p1=pnorm_p12(q00,t10,t20,ymn=v100,slope=v2p1,sigma1=v3,sigma2=v4)
# v3 derivatives
	F3m1=pnorm_p12(q00,t10,t20,ymn=v100,slope=v200,sigma1=v3m1,sigma2=v4)
	F3p1=pnorm_p12(q00,t10,t20,ymn=v100,slope=v200,sigma1=v3p1,sigma2=v4)
# v4 derivatives
	F4m1=pnorm_p12(q00,t10,t20,ymn=v100,slope=v200,sigma1=v3,sigma2=v4m1)
	F4p1=pnorm_p12(q00,t10,t20,ymn=v100,slope=v200,sigma1=v3,sigma2=v4p1)
	mu1=matrix(0,4,length(alpha))
	mu1[1,]=-(F1p1-F1m1)/(2*d1)
	mu1[2,]=-(F2p1-F2m1)/(2*d2)
	mu1[3,]=-(F3p1-F3m1)/(2*d3)
	mu1[4,]=-(F4p1-F4m1)/(2*d4)
	return(mu1)
}
#' Normal-with-p12: DMGS equation 1.2 f2 term
#' @inheritParams manf
norm_p12_f2f=function(y,t10,t20,v1,d1,v2,d2,v3,d3,v4,d4){
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
# v3 stuff
	v3m2=v3-2*d3
	v3m1=v3-1*d3
	v300=v3+0*d3
	v3p1=v3+1*d3
	v3p2=v3+2*d3
# v4 stuff
	v4m2=v4-2*d4
	v4m1=v4-1*d4
	v400=v4+0*d4
	v4p1=v4+1*d4
	v4p2=v4+2*d4
	f2=array(0,c(4,4,length(y)))
# v1
	F1m2=dnorm_p12(y,t10,t20,ymn=v1m2,slope=v2,sigma1=v3,sigma2=v4)
	F1m1=dnorm_p12(y,t10,t20,ymn=v1m1,slope=v2,sigma1=v3,sigma2=v4)
	F100=dnorm_p12(y,t10,t20,ymn=v100,slope=v2,sigma1=v3,sigma2=v4)
	F1p1=dnorm_p12(y,t10,t20,ymn=v1p1,slope=v2,sigma1=v3,sigma2=v4)
	F1p2=dnorm_p12(y,t10,t20,ymn=v1p2,slope=v2,sigma1=v3,sigma2=v4)
	f2[1,1,]=(F1p1-2*F100+F1m1)/(d1*d1)
# v2 derivative
	F2m2=dnorm_p12(y,t10,t20,ymn=v1,slope=v2m2,sigma1=v3,sigma2=v4)
	F2m1=dnorm_p12(y,t10,t20,ymn=v1,slope=v2m1,sigma1=v3,sigma2=v4)
	F200=dnorm_p12(y,t10,t20,ymn=v1,slope=v200,sigma1=v3,sigma2=v4)
	F2p1=dnorm_p12(y,t10,t20,ymn=v1,slope=v2p1,sigma1=v3,sigma2=v4)
	F2p2=dnorm_p12(y,t10,t20,ymn=v1,slope=v2p2,sigma1=v3,sigma2=v4)
	f2[2,2,]=(F2p1-2*F200+F2m1)/(d2*d2)
# v3 derivative
	F3m2=dnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3m2,sigma2=v4)
	F3m1=dnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3m1,sigma2=v4)
	F300=dnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v300,sigma2=v4)
	F3p1=dnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3p1,sigma2=v4)
	F3p2=dnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3p2,sigma2=v4)
	f2[3,3,]=(F3p1-2*F300+F3m1)/(d3*d3)
# v4 derivative
	F4m2=dnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3,sigma2=v4m2)
	F4m1=dnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3,sigma2=v4m1)
	F400=dnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3,sigma2=v400)
	F4p1=dnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3,sigma2=v4p1)
	F4p2=dnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3,sigma2=v4p2)
	f2[4,4,]=(F4p1-2*F400+F4m1)/(d4*d4)
# cross derivative 12
	Fcm1m1=dnorm_p12(y,t10,t20,ymn=v1m1,slope=v2m1,sigma1=v3,sigma2=v4)
	Fcm1p1=dnorm_p12(y,t10,t20,ymn=v1m1,slope=v2p1,sigma1=v3,sigma2=v4)
	Fcp1m1=dnorm_p12(y,t10,t20,ymn=v1p1,slope=v2m1,sigma1=v3,sigma2=v4)
	Fcp1p1=dnorm_p12(y,t10,t20,ymn=v1p1,slope=v2p1,sigma1=v3,sigma2=v4)
	f2[1,2,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2)
	f2[2,1,]=f2[1,2,]
# cross derivative 13
	Fcm1m1=dnorm_p12(y,t10,t20,ymn=v1m1,slope=v2,sigma1=v3m1,sigma2=v4)
	Fcm1p1=dnorm_p12(y,t10,t20,ymn=v1m1,slope=v2,sigma1=v3p1,sigma2=v4)
	Fcp1m1=dnorm_p12(y,t10,t20,ymn=v1p1,slope=v2,sigma1=v3m1,sigma2=v4)
	Fcp1p1=dnorm_p12(y,t10,t20,ymn=v1p1,slope=v2,sigma1=v3p1,sigma2=v4)
	f2[1,3,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d3)
	f2[3,1,]=f2[1,3,]
# cross derivative 14
	Fcm1m1=dnorm_p12(y,t10,t20,ymn=v1m1,slope=v2,sigma1=v3,sigma2=v4m1)
	Fcm1p1=dnorm_p12(y,t10,t20,ymn=v1m1,slope=v2,sigma1=v3,sigma2=v4p1)
	Fcp1m1=dnorm_p12(y,t10,t20,ymn=v1p1,slope=v2,sigma1=v3,sigma2=v4m1)
	Fcp1p1=dnorm_p12(y,t10,t20,ymn=v1p1,slope=v2,sigma1=v3,sigma2=v4p1)
	f2[1,4,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d4)
	f2[4,1,]=f2[1,4,]
# cross derivative 23
	Fcm1m1=dnorm_p12(y,t10,t20,ymn=v1,slope=v2m1,sigma1=v3m1,sigma2=v4)
	Fcm1p1=dnorm_p12(y,t10,t20,ymn=v1,slope=v2p1,sigma1=v3m1,sigma2=v4)
	Fcp1m1=dnorm_p12(y,t10,t20,ymn=v1,slope=v2m1,sigma1=v3p1,sigma2=v4)
	Fcp1p1=dnorm_p12(y,t10,t20,ymn=v1,slope=v2p1,sigma1=v3p1,sigma2=v4)
	f2[2,3,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d2*d3)
	f2[3,2,]=f2[2,3,]
# cross derivative 24
	Fcm1m1=dnorm_p12(y,t10,t20,ymn=v1,slope=v2m1,sigma1=v3,sigma2=v4m1)
	Fcm1p1=dnorm_p12(y,t10,t20,ymn=v1,slope=v2p1,sigma1=v3,sigma2=v4m1)
	Fcp1m1=dnorm_p12(y,t10,t20,ymn=v1,slope=v2m1,sigma1=v3,sigma2=v4p1)
	Fcp1p1=dnorm_p12(y,t10,t20,ymn=v1,slope=v2p1,sigma1=v3,sigma2=v4p1)
	f2[2,4,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d2*d4)
	f2[4,2,]=f2[2,4,]
# cross derivative 34
	Fcm1m1=dnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3m1,sigma2=v4m1)
	Fcm1p1=dnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3p1,sigma2=v4m1)
	Fcp1m1=dnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3m1,sigma2=v4p1)
	Fcp1p1=dnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3p1,sigma2=v4p1)
	f2[3,4,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d3*d4)
	f2[4,3,]=f2[3,4,]
	return(f2)
}
#' Normal-with-p12: DMGS equation 1.2 p2 term
#' @inheritParams manf
norm_p12_p2f=function(y,t10,t20,v1,d1,v2,d2,v3,d3,v4,d4){
	#d3=d3*v3
	##ymax=ifelse(v4<0,v1+v2*t0-2*(d1+d2*t0)-(v3-2*d3)/(v4-2*d4),Inf)
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
# v3 stuff
	v3m2=v3-2*d3
	v3m1=v3-1*d3
	v300=v3+0*d3
	v3p1=v3+1*d3
	v3p2=v3+2*d3
# v4 stuff
	v4m2=v4-2*d4
	v4m1=v4-1*d4
	v400=v4+0*d4
	v4p1=v4+1*d4
	v4p2=v4+2*d4
	p2=array(0,c(4,4,length(y)))
# v1
	F1m2=pnorm_p12(y,t10,t20,ymn=v1m2,slope=v2,sigma1=v3,sigma2=v4)
	F1m1=pnorm_p12(y,t10,t20,ymn=v1m1,slope=v2,sigma1=v3,sigma2=v4)
	F100=pnorm_p12(y,t10,t20,ymn=v100,slope=v2,sigma1=v3,sigma2=v4)
	F1p1=pnorm_p12(y,t10,t20,ymn=v1p1,slope=v2,sigma1=v3,sigma2=v4)
	F1p2=pnorm_p12(y,t10,t20,ymn=v1p2,slope=v2,sigma1=v3,sigma2=v4)
	p2[1,1,]=(F1p1-2*F100+F1m1)/(d1*d1)
# v2 derivative
	F2m2=pnorm_p12(y,t10,t20,ymn=v1,slope=v2m2,sigma1=v3,sigma2=v4)
	F2m1=pnorm_p12(y,t10,t20,ymn=v1,slope=v2m1,sigma1=v3,sigma2=v4)
	F200=pnorm_p12(y,t10,t20,ymn=v1,slope=v200,sigma1=v3,sigma2=v4)
	F2p1=pnorm_p12(y,t10,t20,ymn=v1,slope=v2p1,sigma1=v3,sigma2=v4)
	F2p2=pnorm_p12(y,t10,t20,ymn=v1,slope=v2p2,sigma1=v3,sigma2=v4)
	p2[2,2,]=(F2p1-2*F200+F2m1)/(d2*d2)
# v3 derivative
	F3m2=pnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3m2,sigma2=v4)
	F3m1=pnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3m1,sigma2=v4)
	F300=pnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v300,sigma2=v4)
	F3p1=pnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3p1,sigma2=v4)
	F3p2=pnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3p2,sigma2=v4)
	p2[3,3,]=(F3p1-2*F300+F3m1)/(d3*d3)
# v4 derivative
	F4m2=pnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3,sigma2=v4m2)
	F4m1=pnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3,sigma2=v4m1)
	F400=pnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3,sigma2=v400)
	F4p1=pnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3,sigma2=v4p1)
	F4p2=pnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3,sigma2=v4p2)
	p2[4,4,]=(F4p1-2*F400+F4m1)/(d4*d4)
# cross derivative12
	Fcm1m1=pnorm_p12(y,t10,t20,ymn=v1m1,slope=v2m1,sigma1=v3,sigma2=v4)
	Fcm1p1=pnorm_p12(y,t10,t20,ymn=v1m1,slope=v2p1,sigma1=v3,sigma2=v4)
	Fcp1m1=pnorm_p12(y,t10,t20,ymn=v1p1,slope=v2m1,sigma1=v3,sigma2=v4)
	Fcp1p1=pnorm_p12(y,t10,t20,ymn=v1p1,slope=v2p1,sigma1=v3,sigma2=v4)
	p2[1,2,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2)
	p2[2,1,]=p2[1,2,]
# cross derivative13
	Fcm1m1=pnorm_p12(y,t10,t20,ymn=v1m1,slope=v2,sigma1=v3m1,sigma2=v4)
	Fcm1p1=pnorm_p12(y,t10,t20,ymn=v1m1,slope=v2,sigma1=v3p1,sigma2=v4)
	Fcp1m1=pnorm_p12(y,t10,t20,ymn=v1p1,slope=v2,sigma1=v3m1,sigma2=v4)
	Fcp1p1=pnorm_p12(y,t10,t20,ymn=v1p1,slope=v2,sigma1=v3p1,sigma2=v4)
	p2[1,3,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d3)
	p2[3,1,]=p2[1,3,]
# cross derivative14
	Fcm1m1=pnorm_p12(y,t10,t20,ymn=v1m1,slope=v2,sigma1=v3,sigma2=v4m1)
	Fcm1p1=pnorm_p12(y,t10,t20,ymn=v1m1,slope=v2,sigma1=v3,sigma2=v4p1)
	Fcp1m1=pnorm_p12(y,t10,t20,ymn=v1p1,slope=v2,sigma1=v3,sigma2=v4m1)
	Fcp1p1=pnorm_p12(y,t10,t20,ymn=v1p1,slope=v2,sigma1=v3,sigma2=v4p1)
	p2[1,4,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d4)
	p2[4,1,]=p2[1,4,]
# cross derivative23
	Fcm1m1=pnorm_p12(y,t10,t20,ymn=v1,slope=v2m1,sigma1=v3m1,sigma2=v4)
	Fcm1p1=pnorm_p12(y,t10,t20,ymn=v1,slope=v2p1,sigma1=v3m1,sigma2=v4)
	Fcp1m1=pnorm_p12(y,t10,t20,ymn=v1,slope=v2m1,sigma1=v3p1,sigma2=v4)
	Fcp1p1=pnorm_p12(y,t10,t20,ymn=v1,slope=v2p1,sigma1=v3p1,sigma2=v4)
	p2[2,3,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d2*d3)
	p2[3,2,]=p2[2,3,]
# cross derivative24
	Fcm1m1=pnorm_p12(y,t10,t20,ymn=v1,slope=v2m1,sigma1=v3,sigma2=v4m1)
	Fcm1p1=pnorm_p12(y,t10,t20,ymn=v1,slope=v2p1,sigma1=v3,sigma2=v4m1)
	Fcp1m1=pnorm_p12(y,t10,t20,ymn=v1,slope=v2m1,sigma1=v3,sigma2=v4p1)
	Fcp1p1=pnorm_p12(y,t10,t20,ymn=v1,slope=v2p1,sigma1=v3,sigma2=v4p1)
	p2[2,4,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d2*d4)
	p2[4,2,]=p2[2,4,]
# cross derivative34
	Fcm1m1=pnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3m1,sigma2=v4m1)
	Fcm1p1=pnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3p1,sigma2=v4m1)
	Fcp1m1=pnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3m1,sigma2=v4p1)
	Fcp1p1=pnorm_p12(y,t10,t20,ymn=v1,slope=v2,sigma1=v3p1,sigma2=v4p1)
	p2[3,4,]=(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d3*d4)
	p2[4,3,]=p2[3,4,]
	return(p2)
}
#' Normal-with-p12: DMGS equation 3.3 mu2 term
#' @inheritParams manf
norm_p12_mu2f=function(alpha,t10,t20,v1,d1,v2,d2,v3,d3,v4,d4){
	q00=qnorm_p12((1-alpha),t10,t20,ymn=v1,slope=v2,sigma1=v3,sigma2=v4)
	#d3=d3*v3
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
# v3 stuff
	v3m2=v3-2*d3
	v3m1=v3-1*d3
	v300=v3+0*d3
	v3p1=v3+1*d3
	v3p2=v3+2*d3
# v4 stuff
	v4m2=v4-2*d4
	v4m1=v4-1*d4
	v400=v4+0*d4
	v4p1=v4+1*d4
	v4p2=v4+2*d4
	mu2=array(0,c(4,4,length(alpha)))
# v1
	F1m2=pnorm_p12(q00,t10,t20,ymn=v1m2,slope=v2,sigma1=v3,sigma2=v4)
	F1m1=pnorm_p12(q00,t10,t20,ymn=v1m1,slope=v2,sigma1=v3,sigma2=v4)
	F100=pnorm_p12(q00,t10,t20,ymn=v100,slope=v2,sigma1=v3,sigma2=v4)
	F1p1=pnorm_p12(q00,t10,t20,ymn=v1p1,slope=v2,sigma1=v3,sigma2=v4)
	F1p2=pnorm_p12(q00,t10,t20,ymn=v1p2,slope=v2,sigma1=v3,sigma2=v4)
	mu2[1,1,]=-(F1p1-2*F100+F1m1)/(d1*d1)
# v2 derivative
	F2m2=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2m2,sigma1=v3,sigma2=v4)
	F2m1=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2m1,sigma1=v3,sigma2=v4)
	F200=pnorm_p12(q00,t10,t20,ymn=v1,slope=v200,sigma1=v3,sigma2=v4)
	F2p1=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2p1,sigma1=v3,sigma2=v4)
	F2p2=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2p2,sigma1=v3,sigma2=v4)
	mu2[2,2,]=-(F2p1-2*F200+F2m1)/(d2*d2)
# v3 derivative
	F3m2=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2,sigma1=v3m2,sigma2=v4)
	F3m1=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2,sigma1=v3m1,sigma2=v4)
	F300=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2,sigma1=v300,sigma2=v4)
	F3p1=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2,sigma1=v3p1,sigma2=v4)
	F3p2=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2,sigma1=v3p2,sigma2=v4)
	mu2[3,3,]=-(F3p1-2*F300+F3m1)/(d3*d3)
# v4 derivative
	F4m2=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2,sigma1=v3,sigma2=v4m2)
	F4m1=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2,sigma1=v3,sigma2=v4m1)
	F400=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2,sigma1=v3,sigma2=v400)
	F4p1=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2,sigma1=v3,sigma2=v4p1)
	F4p2=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2,sigma1=v3,sigma2=v4p2)
	mu2[4,4,]=-(F4p1-2*F400+F4m1)/(d4*d4)
# cross derivative12
	Fcm1m1=pnorm_p12(q00,t10,t20,ymn=v1m1,slope=v2m1,sigma1=v3,sigma2=v4)
	Fcm1p1=pnorm_p12(q00,t10,t20,ymn=v1m1,slope=v2p1,sigma1=v3,sigma2=v4)
	Fcp1m1=pnorm_p12(q00,t10,t20,ymn=v1p1,slope=v2m1,sigma1=v3,sigma2=v4)
	Fcp1p1=pnorm_p12(q00,t10,t20,ymn=v1p1,slope=v2p1,sigma1=v3,sigma2=v4)
	mu2[1,2,]=-(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d2)
	mu2[2,1,]=mu2[1,2,]
# cross derivative13
	Fcm1m1=pnorm_p12(q00,t10,t20,ymn=v1m1,slope=v2,sigma1=v3m1,sigma2=v4)
	Fcm1p1=pnorm_p12(q00,t10,t20,ymn=v1m1,slope=v2,sigma1=v3p1,sigma2=v4)
	Fcp1m1=pnorm_p12(q00,t10,t20,ymn=v1p1,slope=v2,sigma1=v3m1,sigma2=v4)
	Fcp1p1=pnorm_p12(q00,t10,t20,ymn=v1p1,slope=v2,sigma1=v3p1,sigma2=v4)
	mu2[1,3,]=-(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d3)
	mu2[3,1,]=mu2[1,3,]
# cross derivative14
	Fcm1m1=pnorm_p12(q00,t10,t20,ymn=v1m1,slope=v2,sigma1=v3,sigma2=v4m1)
	Fcm1p1=pnorm_p12(q00,t10,t20,ymn=v1m1,slope=v2,sigma1=v3,sigma2=v4p1)
	Fcp1m1=pnorm_p12(q00,t10,t20,ymn=v1p1,slope=v2,sigma1=v3,sigma2=v4m1)
	Fcp1p1=pnorm_p12(q00,t10,t20,ymn=v1p1,slope=v2,sigma1=v3,sigma2=v4p1)
	mu2[1,4,]=-(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d1*d4)
	mu2[4,1,]=mu2[1,4,]
# cross derivative23
	Fcm1m1=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2m1,sigma1=v3m1,sigma2=v4)
	Fcm1p1=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2p1,sigma1=v3m1,sigma2=v4)
	Fcp1m1=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2m1,sigma1=v3p1,sigma2=v4)
	Fcp1p1=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2p1,sigma1=v3p1,sigma2=v4)
	mu2[2,3,]=-(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d2*d3)
	mu2[3,2,]=mu2[2,3,]
# cross derivative24
	Fcm1m1=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2m1,sigma1=v3,sigma2=v4m1)
	Fcm1p1=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2p1,sigma1=v3,sigma2=v4m1)
	Fcp1m1=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2m1,sigma1=v3,sigma2=v4p1)
	Fcp1p1=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2p1,sigma1=v3,sigma2=v4p1)
	mu2[2,4,]=-(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d2*d4)
	mu2[4,2,]=mu2[2,4,]
# cross derivative34
	Fcm1m1=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2,sigma1=v3m1,sigma2=v4m1)
	Fcm1p1=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2,sigma1=v3p1,sigma2=v4m1)
	Fcp1m1=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2,sigma1=v3m1,sigma2=v4p1)
	Fcp1p1=pnorm_p12(q00,t10,t20,ymn=v1,slope=v2,sigma1=v3p1,sigma2=v4p1)
	mu2[3,4,]=-(Fcp1p1-Fcm1p1-Fcp1m1+Fcm1m1)/(4*d3*d4)
	mu2[4,3,]=mu2[3,4,]
	return(mu2)
}
#' Log scores for 5 predictions calculated using leave-one-out
#' @inheritParams manf
norm_p12_logscores=function(logscores,x,t1,t2,ics,d1=0.01,d2=0.01,d3=0.01,d4=0.01){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		cp_oos_logscore=0
		for (i in 1:nx){

			x1=x[-i]
			t1m1=t1[-i]
			t2m1=t2[-i]

			dd=dnorm_p12dmgs(x1,t1m1,t2m1,x[i],t1[i],t2[i],ics,d1,d2,d3,d4)

			ml_params=dd$ml_params

			ml_pdf=dd$ml_pdf
			ml_oos_logscore=ml_oos_logscore+log(max(ml_pdf,.Machine$double.eps))

			cp_pdf=dd$cp_pdf
			cp_oos_logscore=cp_oos_logscore+log(max(cp_pdf,.Machine$double.eps))

		}

	}else{
		ml_oos_logscore="logscores not selected"
		cp_oos_logscore="logscores not selected"
	}
	list(	ml_oos_logscore			=ml_oos_logscore,
				cp_oos_logscore			=cp_oos_logscore)
}
#' Densities for 5 predictions
#' @inheritParams manf
dnorm_p12dmgs=function(x,t1,t2,y,t10,t20,ics,d1=0.01,d2=0.01,d3=0.01,d4=0.01,
	aderivs=TRUE){

		debug=TRUE
		debug=FALSE

		nx=length(x)

		ics=norm_p12_setics(x,t1,t2,ics)
		opt=optim(ics,norm_p12_loglik,x=x,t1=t1,t2=t2,control=list(fnscale=-1))
		v1hat=opt$par[1]
		v2hat=opt$par[2]
		v3hat=opt$par[3]
		v4hat=opt$par[4]
		ml_params=c(v1hat,v2hat,v3hat,v4hat)
		if(debug)cat("inside dnorm_p12dmgs:ml_params=",ml_params,"\n")

# mle
		ml_pdf=dnorm_p12(y,t10,t20,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat)
		ml_cdf=pnorm_p12(y,t10,t20,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat)

# bayesian ingredients
		if(aderivs) ldd		=norm_p12_ldda(x,t1,t2,v1hat,v2hat,v3hat,v4hat)
		if(!aderivs)ldd		=norm_p12_ldd		(x,t1,t2,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4)
		lddi=solve(ldd)

		if(debug)cat("det(ldd)=",det(ldd),"\n")
		if(aderivs) lddd	=norm_p12_lddda(x,t1,t2,v1hat,v2hat,v3hat,v4hat)
		if(!aderivs)lddd	=norm_p12_lddd	(x,t1,t2,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4)

		if(aderivs) f1=norm_p12_f1fa(y,t10,t20,v1hat,v2hat,v3hat,v4hat)
		if(!aderivs)f1=norm_p12_f1f(y,t10,t20,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4)

		if(aderivs) f2=norm_p12_f2fa(y,t10,t20,v1hat,v2hat,v3hat,v4hat)
		if(!aderivs)f2=norm_p12_f2f(y,t10,t20,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4)

		if(aderivs) p1=norm_p12_p1fa(y,t10,t20,v1hat,v2hat,v3hat,v4hat)
		if(!aderivs)p1=norm_p12_p1f(y,t10,t20,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4)

		if(aderivs) p2=norm_p12_p2fa(y,t10,t20,v1hat,v2hat,v3hat,v4hat)
		if(!aderivs)p2=norm_p12_p2f(y,t10,t20,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4)

# pu
		lambdad_cp=c(0,0,0,0)
		if(debug)cat("lddi=",lddi,"\n")
		if(debug)cat("lddd=",lddd,"\n")
		if(debug)cat("f1=",f1,"\n")
		if(debug)cat("f2=",f2,"\n")
		df1=dmgs(lddi,lddd,f1,lambdad_cp,f2,dim=4)
		dp1=dmgs(lddi,lddd,p1,lambdad_cp,p2,dim=4)
#		cp_pdf=pmax(ml_pdf+df1/nx,0)
		cp_pdf=pmax(ml_pdf+df1/nx,sqrt(.Machine$double.eps)) #this avoids zeroes in the density
		cp_cdf=pmin(pmax(ml_cdf+dp1/nx,0),1)

# return
		list(	ml_params=ml_params,
					ml_value=opt$val,
					ml_pdf=ml_pdf,
					ml_cdf=ml_cdf,
					cp_pdf=cp_pdf,
					cp_cdf=cp_cdf)
}

