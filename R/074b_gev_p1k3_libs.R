#' Waic
#' @inherit manwaic return
#' @inheritParams manf
gev_p1k3_waic=function(waicscores,x,t,v1hat,v2hat,v3hat,kshape,
	lddi,lddd,lambdad){
		if(waicscores){
			f1f=gev_p1k3_f1fw(x,t,v1hat,v2hat,v3hat,kshape=kshape)
			f2f=gev_p1k3_f2fw(x,t,v1hat,v2hat,v3hat,kshape=kshape)
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
dgev_p1k3sub=function(x,t,y,t0,kshape){

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
