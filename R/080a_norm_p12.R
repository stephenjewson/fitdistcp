#' Normal Distribution with Predictors on Mean and SD, with Parameter Uncertainty
#'
#' @inherit man description author references seealso return
#' @inheritParams man
#'
#' @inheritSection man Optional Return Values
#' @inheritSection man Details (homogeneous models)
# #' @inheritSection man Details (non-homogeneous models)
# #' @inheritSection man Details (analytic integration)
#' @inheritSection man Details (DMGS integration)
#'
#' @example man/examples/example_080_norm_p12.R
#'
#' @name norm_p12_cp
NULL
#' @rdname norm_p12_cp
#' @inheritParams man
#' @export
#'
# Centering the input data for x and t1 is often critical
#
qnorm_p12_cp=function(x,t1,t2,t01=NA,t02=NA,n01=NA,n02=NA,p=seq(0.1,0.9,0.1),ics=c(0,0,0,0),
	means=FALSE,waicscores=FALSE,logscores=FALSE,
	dmgs=TRUE,rust=FALSE,nrust=100000,
	extramodels=FALSE,predictordata=TRUE,centering=TRUE,
	debug=FALSE){

	stopifnot(	is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1,
							is.finite(t1),!is.na(t1),
							is.finite(t2),!is.na(t2),
							length(ics)==4)
#
# 1 intro
#
	alpha=1-p
	nx=length(x)
	nalpha=length(alpha)
	t01=maket0(t0=t01,n0=n01,t=t1)
	t02=maket0(t0=t02,n0=n02,t=t2)
	if(debug)message("t01,t02=",t01,t02)
#
# 2 centering
#
	if(centering){
		meant1=mean(t1)
		meant2=mean(t2)
		t1=t1-meant1
		t2=t2-meant2
		t01=t01-meant1
		t02=t02-meant2
	}
#
# 3 ml param estimate
#
	if(debug)message("  maxlik")
	ics=norm_p12_setics(x,t1,t2,ics)
	opt1=optim(ics,norm_p12_loglik,x=x,t1=t1,t2=t2,control=list(fnscale=-1))
	v1hat=opt1$par[1]
	v2hat=opt1$par[2]
	v3hat=opt1$par[3]
	v4hat=opt1$par[4]
	ml_params=c(v1hat,v2hat,v3hat,v4hat)
	muhat=ml_params[1]+ml_params[2]*t1
	muhat0=makemuhat0(t01,n01,t1,ml_params)
	sghat=exp(ml_params[3]+ml_params[4]*t2)
	residuals=x-muhat
	norm_p12_checkmle(ml_params)
	if(debug)message("  ml_params=",ml_params)
#
# 4 predictordata
#
	prd=norm_p12_predictordata(predictordata,x,t1,t2,t01,t02,ml_params)
	predictedparameter=prd$predictedparameter
	adjustedx=prd$adjustedx
#
# 5 aic
#
	if(debug)message("3")
	ml_value=opt1$val
	maic=make_maic(ml_value,nparams=4)
#
# 6 calc ml quantiles and density
#
	if(debug)message("4")
	ml_quantiles=qnorm_p12((1-alpha),t01,t02,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat)
	fhat=dnorm_p12(ml_quantiles,t01,t02,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat,log=FALSE)
#
# 7 ldd (two versions)
#
#
# dmgs
#
	standard_errors="dmgs not selected"
	rh_quantiles="dmgs not selected"
	ru_quantiles="dmgs not selected"
	waic1="dmgs not selected"
	waic2="dmgs not selected"
	ml_oos_logscore="dmgs not selected"
	rh_oos_logscore="dmgs not selected"
	cp_oos_logscore="dmgs not selected"
	ml_mean="dmgs not selected"
	rh_mean="dmgs not selected"
	cp_mean="dmgs not selected"
	cp_method="dmgs not selected"
	if(dmgs){
#
# 8 calculate ldd
#
		ldd=norm_p12_ldda(x,t1,t2,v1hat,v2hat,v3hat,v4hat)
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
#
# 9 calculate lddd
#
		if(debug)message("  lddd")
		lddd=norm_p12_lddda(x,t1,t2,v1hat,v2hat,v3hat,v4hat)
#
# 10 mu1
#
		if(debug)message("calculate mu1")
		mu1=norm_p12_mu1fa(alpha,t01,t02,v1hat,v2hat,v3hat,v4hat)
#
# 11 mu2
#
		if(debug)message("calculate mu2")
		mu2=norm_p12_mu2fa(alpha,t01,t02,v1hat,v2hat,v3hat,v4hat)
#
# 12 model 4: pu
#
		lambdad_rhp=c(0,0,0,0)
		dq=dmgs(lddi,lddd,mu1,lambdad_rhp,mu2,dim=4)
		rh_quantiles=ml_quantiles+dq/(nx*fhat)
#
# 13 model 6: Laplace's method
#
		if(extramodels|means){
			lambdad_lp=c(0,0,0,0)
			lddd_lp=array(0,c(4,4,4))
			dq=dmgs(lddi,lddd_lp,mu1,lambdad_lp,mu2,dim=4)
			lp_quantiles=ml_quantiles+dq/(nx*fhat)
		} else {
			lp_quantiles="extramodels not selected"
		}
#
# 14 means (might as well always calculate)
#
		ml_mean=muhat0
		rh_mean=ml_mean
#
# 15 waicscores
#
		waic=norm_p12_waic(waicscores,x,t1,t2,v1hat,v2hat,v3hat,v4hat,lddi,lddd,lambdad_rhp)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 16 logscores
#
		logscores=norm_p12_logscores(logscores,x,t1,t2,ics)
		ml_oos_logscore				=logscores$ml_oos_logscore
		cp_oos_logscore		=logscores$cp_oos_logscore
#
# 17 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rnorm_p12_cp(nrust,x,t1=t1,t2=t2,t01=t01,t02=t02,n01=NA,n02=NA,rust=TRUE,mlcp=FALSE)
			ru_quantiles=makeq(rustsim$ru_deviates,p)
		}
	} #end of if(dmgs)
#
# 18 decentering
#
	if(centering){
		ml_params[1]=ml_params[1]-ml_params[2]*meant1
		ml_params[3]=ml_params[3]-ml_params[4]*meant2
		if(predictordata){
			predictedparameter$mu=predictedparameter$mu-ml_params[2]*meant1
			predictedparameter$sg=predictedparameter$sg*exp(-ml_params[4]*meant2)
		}
	}

	list(	ml_params=ml_params,
				ml_value=ml_value,
				predictedparameter=predictedparameter,
				adjustedx=adjustedx,
#				expinfmat=expinfmat,
#				expinfmati=expinfmati,
				standard_errors=standard_errors,
				ml_quantiles=ml_quantiles,
				cp_quantiles=rh_quantiles,
				ru_quantiles=ru_quantiles,
				lp_quantiles=lp_quantiles,
				p1=muhat,
				residuals=residuals,
				maic=maic,
				waic1=waic1,
				waic2=waic2,
				ml_oos_logscore=ml_oos_logscore,
				cp_oos_logscore=cp_oos_logscore,
				ml_mean=ml_mean,
				cp_mean=rh_mean,
				cp_method=rhp_dmgs_cpmethod())

}
#' @rdname norm_p12_cp
#' @inheritParams man
#' @export
rnorm_p12_cp=function(n,x,t1,t2,n01=NA,n02=NA,t01=NA,t02=NA,ics=c(0,0,0,0),
	rust=FALSE,mlcp=TRUE,
	debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),
#						is.finite(t1),!is.na(t1),
#						is.finite(t2),!is.na(t2),
#						length(t1)==length(x),length(t1)==length(t2),
#						length(ics)==4)
	stopifnot(is.finite(x),!is.na(x),
						is.finite(t1),!is.na(t1),
						is.finite(t2),!is.na(t2),
						length(ics)==4)

	if(debug)message("t01,n01=",t01,n01)
	if(debug)message("t02,n02=",t02,n02)
	t01=maket0(t0=t01,n0=n01,t=t1)
	t02=maket0(t0=t02,n0=n02,t=t2)

#
# 2 centering
#
	meant1=mean(t1)
	meant2=mean(t2)
	t1=t1-meant1
	t2=t2-meant2
	t01=t01-meant1
	t02=t02-meant2

	ml_params="mlcp not selected"
	ml_deviates="mlcp not selected"
	cp_deviates="mlcp not selected"
	ru_deviates="rust not selected"

	if(mlcp){
		q=qnorm_p12_cp(x,t1=t1,t2=t2,t01=t01,t02=t02,n01=NA,n02=NA,
			p=runif(n),ics=ics)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}


	if(rust){
		th=tnorm_p12_cp("rust",n,x,t1,t2)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			mu=th[i,1]+t01*th[i,2]
			sd=exp(th[i,3]+t02*th[i,4])
			ru_deviates[i]=rnorm(1,mean=mu,sd=sd)
		}
	}

#
# decentering
#
	if(mlcp){
		ml_params[1]=ml_params[1]-ml_params[2]*meant1
		ml_params[3]=ml_params[3]-ml_params[4]*meant2
	}

	op=list(ml_params=ml_params,
			 ml_deviates=ml_deviates,
			 cp_deviates=cp_deviates,
			 ru_deviates=ru_deviates,
			 cp_method=rhp_dmgs_cpmethod())

	return(op)

}
#' @rdname norm_p12_cp
#' @inheritParams man
#' @export
dnorm_p12_cp=function(x,t1,t2,t01=NA,t02=NA,n01=NA,n02=NA,
	y=x,ics=c(0,0,0,0),
	rust=FALSE,nrust=1000,
	boot=FALSE,nboot=10,
	centering=TRUE,rnonnegslopesonly=FALSE,debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
						is.finite(t1),is.finite(t2),
						!is.na(t1),!is.na(t2),
						length(ics)==4)

	if(debug)message("inside dnorm_p12_cp")

	t01=maket0(t0=t01,n0=n01,t=t1)
	t02=maket0(t0=t02,n0=n02,t=t2)

#
# centering
#
	if(centering){
		meant1=mean(t1)
		meant2=mean(t2)
		t1=t1-meant1
		t2=t2-meant2
		t01=t01-meant1
		t02=t02-meant2
	}

	dd=dnorm_p12dmgs(x=x,t1=t1,t2=t2,y=y,t01=t01,t02=t02,ics=ics)
	ru_pdf="rust not selected"
	bs_pdf="boot not selected"
	ml_params=dd$ml_params
	ml_value=dd$ml_value

	if(rust){
		th=tnorm_p12_cp("rust",nrust,x,t1,t2,nonnegslopesonly=rnonnegslopesonly)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+t01*th[ir,2]
			sg=exp(th[ir,3]+t02*th[ir,4])
			ru_pdf=ru_pdf+dnorm(y,mean=mu,sd=sg)
		}
		ru_pdf=ru_pdf/nrust
	}

	if(boot){
		th=tnorm_p12_cp("boot",nboot,x,t1,t2,nonnegslopesonly=rnonnegslopesonly)$theta_samples
		bs_pdf=numeric(length(y))
		for (ir in 1:nboot){
			mu=th[ir,1]+t01*th[ir,2]
			sg=exp(th[ir,3]+t02*th[ir,4])
			bs_pdf=bs_pdf+dnorm(y,mean=mu,sd=sg)
		}
		bs_pdf=bs_pdf/nboot
	}

#
# decentering
#
	if(centering){
		ml_params[1]=ml_params[1]-ml_params[2]*meant1
		ml_params[3]=ml_params[3]-ml_params[4]*meant2
	}

	op=list(ml_params=ml_params,
					ml_value=ml_value,
					ml_pdf=dd$ml_pdf,
					cp_pdf=dd$cp_pdf,
					ml_cdf=dd$ml_cdf,
					cp_cdf=dd$cp_cdf,
					ru_pdf=ru_pdf,
					bs_pdf=bs_pdf,
					cp_method=rhp_dmgs_cpmethod())
	return(op)

}
#' @rdname norm_p12_cp
#' @inheritParams man
#' @export
pnorm_p12_cp=function(x,t1,t2,t01=NA,t02=NA,n01=NA,n02=NA,
	y=x,ics=c(0,0,0,0),
	rust=FALSE,nrust=1000,
	boot=FALSE,nboot=10,
	centering=TRUE,rnonnegslopesonly=FALSE,debug=FALSE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
						is.finite(t1),is.finite(t2),
						!is.na(t1),!is.na(t2),
						length(ics)==4)

	t01=maket0(t0=t01,n0=n01,t=t1)
	t02=maket0(t0=t02,n0=n02,t=t2)

#
# centering
#
	if(centering){
		meant1=mean(t1)
		meant2=mean(t2)
		t1=t1-meant1
		t2=t2-meant2
		t01=t01-meant1
		t02=t02-meant2
	}

	dd=dnorm_p12dmgs(x=x,t1=t1,t2=t2,y=y,t01=t01,t02=t02,ics=ics)
	ru_cdf="rust not selected"
	bs_cdf="boot not selected"
	ml_params=dd$ml_params

	if(rust){
		th=tnorm_p12_cp("rust",nrust,x,t1,t2,nonnegslopesonly=rnonnegslopesonly)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+t01*th[ir,2]
			sg=exp(th[ir,3]+t02*th[ir,4])
			ru_cdf=ru_cdf+pnorm(y,mean=mu,sd=sg)
		}
		ru_cdf=ru_cdf/nrust
	}

	if(boot){
		th=tnorm_p12_cp("boot",nboot,x,t1,t2,nonnegslopesonly=rnonnegslopesonly)$theta_samples
		bs_cdf=numeric(length(y))
		for (ir in 1:nboot){
			mu=th[ir,1]+t01*th[ir,2]
			sg=exp(th[ir,3]+t02*th[ir,4])
			bs_cdf=bs_cdf+pnorm(y,mean=mu,sd=sg)
		}
		bs_cdf=bs_cdf/nboot
	}

#
# decentering
#
	if(centering){
		ml_params[1]=ml_params[1]-ml_params[2]*meant1
		ml_params[3]=ml_params[3]-ml_params[4]*meant2
	}

	op=list(ml_params=ml_params,
					ml_cdf=dd$ml_cdf,
					cp_cdf=dd$cp_cdf,
					ru_cdf=ru_cdf,
					bs_cdf=bs_cdf,
					cp_method=rhp_dmgs_cpmethod())
	return(op)
}
#' @rdname norm_p12_cp
#' @inheritParams man
#' @export
tnorm_p12_cp=function(method,n,x,t1,t2,nonnegslopesonly=FALSE,ics=c(0,0,0,0),
	debug=FALSE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),
#						is.finite(t1),!is.na(t1),
#						is.finite(t2),!is.na(t2),
#						length(t1)==length(x),length(t1)==length(t2),
#						length(ics)==4)
	stopifnot(is.finite(x),!is.na(x),
						is.finite(t1),!is.na(t1),
						is.finite(t2),!is.na(t2),
						length(ics)==4)

#
# centering
#
		meant1=mean(t1)
		meant2=mean(t2)
		t1=t1-meant1
		t2=t2-meant2

	if(method=="rust"){
		th=ru(norm_p12_logf,x=x,t1=t1,t2=t2,nonnegslopesonly=nonnegslopesonly,n=n,d=4,init=c(0,0,0,0))
	} else if (method=="boot"){
		th=norm_p12_boot(x=x,t1=t1,t2=t2,n=n)
	} else{
		message("tnorm_p12 method not valid so stopping.")
		stop()
	}
	theta_samples=th$sim_vals
#
# decentering
#
	theta_samples[,1]=theta_samples[,1]-theta_samples[,2]*meant1
	theta_samples[,3]=theta_samples[,3]-theta_samples[,4]*meant2

	list(theta_samples=theta_samples)

}
