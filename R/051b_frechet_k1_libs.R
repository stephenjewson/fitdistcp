#
# note that the parameter labelling is confusing
# R uses lambda, mu, sigma
# I use mu, sigma, lambda
# I'm fixing location
# so I label the routine _k1_
# and I say sigma=v1, lambda=v2
#' Waic
#' @inherit manwaic return
#' @inheritParams manf
frechet_k1_waic=function(waicscores,x,v1hat,v2hat,kloc,lddi,lddd,
	lambdad){
		if(waicscores){

			f1f=frechet_k1_f1fa(x,v1hat,v2hat,kloc)
			f2f=frechet_k1_f2fa(x,v1hat,v2hat,kloc)

			fhatx=dfrechet(x,mu=kloc,sigma=v1hat,lambda=v2hat)
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
frechet_k1_logf=function(params,x,kloc){
	s=pmax(params[1],.Machine$double.eps)
	l=pmax(params[2],.Machine$double.eps)
	logf=sum(dfrechet(x,mu=kloc,sigma=s,lambda=l,log=TRUE))-log(s)-log(l)
	return(logf)
}
#'  log-likelihood function
#' @inherit manloglik return
#' @inheritParams manf
frechet_loglik=function(vv,x,kloc){

	loglik=sum(dfrechet(x,mu=kloc,sigma=max(vv[1],.Machine$double.eps),lambda=max(vv[2],.Machine$double.eps),log=TRUE))

	return(loglik)
}
#' MLE and RHP predictive means
#' @inherit manmeans return
#' @inheritParams manf
frechet_means=function(means,ml_params,lddi,lddd,lambdad_rhp,nx,dim=2,kloc){

	if(means){
		v1=ml_params[1]
		v2=ml_params[2]
		f2=1-1/v2
		iv22=1/(v2*v2)
		iv23=1/(v2*v2*v2)
		iv24=1/(v2*v2*v2*v2)

# ml mean
		ml_mean=v1*gamma(f2)+kloc

# rhp mean
		meand1=array(0,c(2,1))
		meand1[1,1]=gamma(f2)
		meand1[1,1]=v1*gamma(f2)*digamma(f2)*iv22

		meand2=array(0,c(2,2,1))
		meand2[1,1,1]=0
		meand2[1,2,1]=gamma(f2)*digamma(f2)*iv22
		meand2[2,1,1]=meand2[1,2,1]
		meand2[2,2,1]=v1*(gamma(f2)*trigamma(f2)*iv24-2*gamma(f2)*digamma(f2)*iv23)
		dmean=dmgs(lddi,lddd,meand1,lambdad_rhp,meand2,dim=2)
		rh_mean=ml_mean+dmean/nx

	}else{
		ml_mean="means not selected"
		rh_mean="means not selected"
	}

	list(ml_mean=ml_mean,rh_mean=rh_mean)

}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inherit manlogscores return
#' @inheritParams manf
frechet_logscores=function(logscores,x,kloc){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){
			x1=x[-i]
			dd=dfrechetsub(x1,x[i],kloc)

			ml_params=dd$ml_params

			ml_pdf=dd$ml_pdf
			ml_oos_logscore=ml_oos_logscore+log(ml_pdf)

			rh_pdf=dd$rh_pdf
			rh_oos_logscore=rh_oos_logscore+log(max(rh_pdf,.Machine$double.eps))
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
dfrechetsub=function(x,y,kloc){

		nx=length(x)

		opt=optim(c(1,1),frechet_loglik,x=x,kloc=kloc,control=list(fnscale=-1))
		v1hat=opt$par[1]
		v2hat=opt$par[2]
		ml_params=c(v1hat,v2hat)

# mle
		ml_pdf=dfrechet(y,mu=kloc,sigma=v1hat,lambda=v2hat)
		ml_cdf=pfrechet(y,mu=kloc,sigma=v1hat,lambda=v2hat)

# rhp
		ldd=frechet_k1_ldda(x,v1hat,v2hat,kloc)
		lddi=solve(ldd)
		lddd=frechet_k1_lddda(x,v1hat,v2hat,kloc)

		f1=frechet_k1_f1fa(y,v1hat,v2hat,kloc)
		f2=frechet_k1_f2fa(y,v1hat,v2hat,kloc)

		p1=frechet_k1_p1fa(y,v1hat,v2hat,kloc)
		p2=frechet_k1_p2fa(y,v1hat,v2hat,kloc)

		lambdad_rhp=c(-1/v1hat,-1/v2hat)
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



