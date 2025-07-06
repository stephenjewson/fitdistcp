#' Waic
#' @inherit manwaic return
#' @inheritParams manf
norm_waic=function(waicscores,x,v1hat,d1,v2hat,fd2){
	if(waicscores){

		f1f=norm_f1fa(x,v1hat,v2hat)
		f2f=norm_f2fa(x,v1hat,v2hat)
		ldd=norm_ldda(x,v1hat,v2hat)
		lddi=solve(ldd)
		lddd=norm_lddda(x,v1hat,v2hat)
		fhatx=dnorm(x,mean=v1hat,sd=v2hat)
		lambdad_rhp=c(0,-1/v2hat)
		waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad_rhp,f2f,dim=2)
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
norm_logf=function(params,x){
	m=params[1]
	s=pmax(params[2],.Machine$double.eps)
	logf=sum(dnorm(x,mean=m,sd=s,log=TRUE))-log(s)
	return(logf)
}
#' Maximum likelihood estimator
#' @inherit manloglik return
#' @inheritParams manf
norm_ml_params=function(x){
	mlparams=matrix(0,2)
	nx=length(x)
	mlparams[1]=mean(x)
	mlparams[2]=sqrt((nx-1)/nx)*sd(x)
	return(mlparams)
}
#' Method of moments estimator
#' @return Vector
#' @inheritParams manf
norm_unbiasedv_params=function(x){
	params=matrix(0,2)
	nx=length(x)
	params[1]=mean(x)
	params[2]=sd(x)
	return(params)
}
#' One component of the second derivative of the expected log-likelihood
#' @inherit manlnn return
#' @inheritParams manf
norm_gmn=function(alpha,v1,d1,v2,fd2,mm,nn){
# saved for future testing of the mpd theory
	nx=length(alpha)
	d2=fd2*v2
  x=qnorm((1-alpha),mean=v1,sd=v2)
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
			lmn[,i]=dnorm(x,mean=vvd[1],sd=vvd[2],log=TRUE)
		}
		dld=(lmn[,1]-lmn[,2]-lmn[,3]+lmn[,4])/(4*dd[mm]*dd[nn])
# same
	} else {
		net3[,mm]=c(-1,0,1)
		for (i in 1:3){
			for (j in 1:2){
				vvd[j]=vv[j]+net3[i,j]*dd[j]
			}
			lmn[,i]=dnorm(x,mean=vvd[1],sd=vvd[2],log=TRUE)
		}
		dld=(lmn[,1]-2*lmn[,2]+lmn[,3])/(dd[mm]*dd[mm])
	}
	return(dld)
}
#' Second derivative matrix of the expected per-observation log-likelihood
#' @inherit manldd return
#' @inheritParams manf
norm_gg=function(nx,v1,d1,v2,fd2){
# saved for future testing of the mpd theory
	expinfmat=matrix(0,2,2)
	for (i in 1:2){
		for (j in i:2){
			expinfmat[i,j]=-quad(norm_gmn,xa=0,xb=1,v1=v1,d1=d1,v2=v2,fd2=fd2,mm=i,nn=j)
		}
	}
	for (i in 2:2){
		for (j in 1:(i-1)){
			expinfmat[i,j]=expinfmat[j,i]
		}
	}
 return(expinfmat)
}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inherit manlogscores return
#' @inheritParams manf
norm_logscores=function(logscores,x){

	if(logscores){

		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){
			x1=x[-i]

			dd=dnormsub(x1,x[i])
			ml_params=dd$ml_params

			ml_pdf=dd$ml_pdf
			ml_oos_logscore=ml_oos_logscore+log(ml_pdf)

			rh_pdf=dd$rh_pdf
			rh_oos_logscore=rh_oos_logscore+log(rh_pdf)
		}
	}else{
		ml_oos_logscore="logscores not selected"
		rh_oos_logscore="logscores not selected"
	}
	list(ml_oos_logscore=ml_oos_logscore,rh_oos_logscore=rh_oos_logscore)
}
#' Densities from MLE and RHP
#' @inherit mandsub return
#' @inheritParams manf
dnormsub=function(x,y){

		nx=length(x)

# ml
		ml_params=norm_ml_params(x)
		ml_pdf=dnorm(y,mean=ml_params[1],sd=ml_params[2])
		ml_cdf=pnorm(y,mean=ml_params[1],sd=ml_params[2])

# rhp pdf
# -first, convert sigma from maxlik to unbiased
		sgu=ml_params[2]*sqrt(nx/(nx-1))
# then, convert sigma to predictive sigma
		sg1=sgu*sqrt((nx+1)/nx)

		yd=(y-ml_params[1])/sg1
		rh_pdf=(dt(yd,df=nx-1))/sg1
		rh_pdf=pmax(rh_pdf,0)

# rhp cdf
		rh_cdf=pt(yd,df=nx-1)
		rh_cdf=pmin(pmax(rh_cdf,0),1)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					rh_pdf=rh_pdf,
					ml_cdf=ml_cdf,
					rh_cdf=rh_cdf)
}
