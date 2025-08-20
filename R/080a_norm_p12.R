#' Normal Distribution with Predictors on Mean and SD, with Parameter Uncertainty
#'
#' @inherit man description author references seealso
#' @inheritParams man
#'
#' @inheritSection man Default Return Values
#' @inheritSection man Optional Return Values
#' @inheritSection man Optional Return Values (EVD models only)
#' @inheritSection man Optional Return Values (non-RHP models only)
# #' @inheritSection man Details (homogeneous models)
#' @inheritSection man Details (non-homogeneous models)
# #' @inheritSection man Details (analytic integration)
#' @inheritSection man Details (DMGS integration)
#'
#' @example examples/example_80_norm_p12_cp.R
#'
#' @name Normal_p12_cp
NULL
#' @rdname Normal_p12_cp
#' @inheritParams man
#' @export
#'
# Centering the input data for x and t1 is often critical
#
qnorm_p12_cp=function(x,t1,t2,t10=NA,t20=NA,n10=NA,n20=NA,p=seq(0.1,0.9,0.1),ics=c(0,0,0,0),
	d1=0.01,d2=0.01,d3=0.01,d4=0.01,fdalpha=0.01,
	means=FALSE,waicscores=FALSE,logscores=FALSE,
	dmgs=TRUE,rust=FALSE,nrust=100000,
	extramodels=FALSE,predictordata=TRUE,centering=TRUE,
	debug=FALSE,aderivs=TRUE){

	stopifnot(	is.finite(x),!is.na(x),is.finite(p),!is.na(p),p>0,p<1,
							is.finite(t1),!is.na(t1),
							is.finite(t2),!is.na(t2),
							length(ics)==4)
#
# 1 intro
#
	alpha=1-p
#	cat("alpha1=",alpha,"\n")
	nx=length(x)
	nalpha=length(alpha)
	t10=maket0(t0=t10,n0=n10,t=t1)
	t20=maket0(t0=t20,n0=n20,t=t2)
	if(debug)cat("t10,t20=",t10,t20,"\n")
#
# 2 centering
#
	if(centering){
		meant1=mean(t1)
		meant2=mean(t2)
		t1=t1-meant1
		t2=t2-meant2
		t10=t10-meant1
		t20=t20-meant2
	}
#
# 3 ml param estimate
#
	if(debug)cat("  maxlik\n")
	ics=norm_p12_setics(x,t1,t2,ics)
	opt1=optim(ics,norm_p12_loglik,x=x,t1=t1,t2=t2,control=list(fnscale=-1))
	v1hat=opt1$par[1]
	v2hat=opt1$par[2]
	v3hat=opt1$par[3]
	v4hat=opt1$par[4]
	ml_params=c(v1hat,v2hat,v3hat,v4hat)
	muhat=ml_params[1]+ml_params[2]*t1
	muhat0=makemuhat0(t10,n10,t1,ml_params)
	sghat=exp(ml_params[3]+ml_params[4]*t2)
	residuals=x-muhat
	norm_p12_checkmle(ml_params)
	if(debug)cat("  ml_params=",ml_params,"\n")
#
# 4 predictordata
#
	prd=norm_p12_predictordata(predictordata,x,t1,t2,t10,t20,ml_params)
	predictedparameter=prd$predictedparameter
	adjustedx=prd$adjustedx
#
# 5 aic
#
	if(debug)cat("3\n")
	ml_value=opt1$val
	maic=make_maic(ml_value,nparams=4)
#
# 6 calc ml quantiles and density
#
	if(debug)cat("4\n")
	ml_quantiles=qnorm_p12((1-alpha),t10,t20,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat)
	fhat=dnorm_p12(ml_quantiles,t10,t20,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat,log=FALSE)
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
		if(aderivs) ldd=norm_p12_ldda(x,t1,t2,v1hat,v2hat,v3hat,v4hat)
		if(!aderivs)ldd=norm_p12_ldd(x,t1,t2,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4)
		lddi=solve(ldd)
		standard_errors=make_se(nx,lddi)
#
# 9 calculate lddd
#
		if(debug)cat("  lddd\n")
		if(aderivs) lddd=norm_p12_lddda(x,t1,t2,v1hat,v2hat,v3hat,v4hat)
		if(!aderivs)lddd=norm_p12_lddd(x,t1,t2,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4)
#
# 10 mu1
#
		if(debug)cat("calculate mu1\n")
		if(aderivs) mu1=norm_p12_mu1fa(alpha,t10,t20,v1hat,v2hat,v3hat,v4hat)
		if(!aderivs)mu1=norm_p12_mu1f(alpha,t10,t20,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4)
#
# 11 mu2
#
		if(debug)cat("calculate mu2\n")
		if(aderivs) mu2=norm_p12_mu2fa(alpha,t10,t20,v1hat,v2hat,v3hat,v4hat)
		if(!aderivs)mu2=norm_p12_mu2f(alpha,t10,t20,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4)
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
		waic=norm_p12_waic(waicscores,x,t1,t2,v1hat,d1,v2hat,d2,v3hat,d3,v4hat,d4,lddi,lddd,lambdad_rhp)
		waic1=waic$waic1
		waic2=waic$waic2
#
# 16 logscores
#
		logscores=norm_p12_logscores(logscores,x,t1,t2,ics,d1,d2,d3,d4)
		ml_oos_logscore				=logscores$ml_oos_logscore
		cp_oos_logscore		=logscores$cp_oos_logscore
#
# 17 rust
#
		ru_quantiles="rust not selected"
		if(rust){
			rustsim=rnorm_p12_cp(nrust,x,t1=t1,t2=t2,t10=t10,t20=t20,n10=NA,n20=NA,rust=TRUE,mlpu=FALSE)
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
#' @rdname Normal_p12_cp
#' @inheritParams man
#' @export
rnorm_p12_cp=function(n,x,t1,t2,n10=NA,n20=NA,t10=NA,t20=NA,ics=c(0,0,0,0),
	d1=0.01,d2=0.01,d3=0.01,d4=0.01,rust=FALSE,mlpu=TRUE,
	debug=FALSE,aderivs=TRUE){

#	stopifnot(is.finite(n),!is.na(n),is.finite(x),!is.na(x),
#						is.finite(t1),!is.na(t1),
#						is.finite(t2),!is.na(t2),
#						length(t1)==length(x),length(t1)==length(t2),
#						length(ics)==4)
	stopifnot(is.finite(x),!is.na(x),
						is.finite(t1),!is.na(t1),
						is.finite(t2),!is.na(t2),
						length(ics)==4)

	if(debug)cat("t10,n10=",t10,n10,"\n")
	if(debug)cat("t20,n20=",t20,n20,"\n")
	t10=maket0(t0=t10,n0=n10,t=t1)
	t20=maket0(t0=t20,n0=n20,t=t2)

#
# 2 centering
#
	meant1=mean(t1)
	meant2=mean(t2)
	t1=t1-meant1
	t2=t2-meant2
	t10=t10-meant1
	t20=t20-meant2

	ml_params="mlpu not selected"
	ml_deviates="mlpu not selected"
	cp_deviates="mlpu not selected"
	ru_deviates="rust not selected"

	if(mlpu){
		q=qnorm_p12_cp(x,t1=t1,t2=t2,t10=t10,t20=t20,n10=NA,n20=NA,
			p=runif(n),ics=ics,d1=d1,d2=d2,d3=d3,d4=d4,aderivs=aderivs)
		ml_params=q$ml_params
		ml_deviates=q$ml_quantiles
		cp_deviates=q$cp_quantiles
	}


	if(rust){
		th=tnorm_p12_cp("rust",n,x,t1,t2)$theta_samples
		ru_deviates=numeric(0)
		for (i in 1:n){
			mu=th[i,1]+t10*th[i,2]
			sd=exp(th[i,3]+t20*th[i,4])
			ru_deviates[i]=rnorm(1,mean=mu,sd=sd)
		}
	}

#
# decentering
#
	if(mlpu){
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
#' @rdname Normal_p12_cp
#' @inheritParams man
#' @export
dnorm_p12_cp=function(x,t1,t2,t10=NA,t20=NA,n10=NA,n20=NA,
	y=x,ics=c(0,0,0,0),d1=0.01,d2=0.01,d3=0.01,d4=0.01,
	rust=FALSE,nrust=1000,
	boot=FALSE,nboot=10,
	centering=TRUE,rnonnegslopesonly=FALSE,debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
						is.finite(t1),is.finite(t2),
						!is.na(t1),!is.na(t2),
						length(ics)==4)

	if(debug)cat("inside dnorm_p12_cp\n")

	t10=maket0(t0=t10,n0=n10,t=t1)
	t20=maket0(t0=t20,n0=n20,t=t2)

#
# centering
#
	if(centering){
		meant1=mean(t1)
		meant2=mean(t2)
		t1=t1-meant1
		t2=t2-meant2
		t10=t10-meant1
		t20=t20-meant2
	}

	dd=dnorm_p12dmgs(x=x,t1=t1,t2=t2,y=y,t10=t10,t20=t20,ics=ics,d1,d2,d3,d4,
		aderivs=aderivs)
	ru_pdf="rust not selected"
	bs_pdf="boot not selected"
	ml_params=dd$ml_params
	ml_value=dd$ml_value

	if(rust){
		th=tnorm_p12_cp("rust",nrust,x,t1,t2,nonnegslopesonly=rnonnegslopesonly)$theta_samples
		ru_pdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+t10*th[ir,2]
			sg=exp(th[ir,3]+t20*th[ir,4])
			ru_pdf=ru_pdf+dnorm(y,mean=mu,sd=sg)
		}
		ru_pdf=ru_pdf/nrust
	}

	if(boot){
		th=tnorm_p12_cp("boot",nrust,x,t1,t2,nonnegslopesonly=rnonnegslopesonly)$theta_samples
		bs_pdf=numeric(length(y))
		for (ir in 1:nboot){
			mu=th[ir,1]+t10*th[ir,2]
			sg=exp(th[ir,3]+t20*th[ir,4])
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
#' @rdname Normal_p12_cp
#' @inheritParams man
#' @export
pnorm_p12_cp=function(x,t1,t2,t10=NA,t20=NA,n10=NA,n20=NA,
	y=x,ics=c(0,0,0,0),d1=0.01,d2=0.01,d3=0.01,d4=0.01,
	rust=FALSE,nrust=1000,
	boot=FALSE,nboot=10,
	centering=TRUE,rnonnegslopesonly=FALSE,debug=FALSE,aderivs=TRUE){

	stopifnot(is.finite(x),!is.na(x),is.finite(y),!is.na(y),
						is.finite(t1),is.finite(t2),
						!is.na(t1),!is.na(t2),
						length(ics)==4)

	t10=maket0(t0=t10,n0=n10,t=t1)
	t20=maket0(t0=t20,n0=n20,t=t2)

#
# centering
#
	if(centering){
		meant1=mean(t1)
		meant2=mean(t2)
		t1=t1-meant1
		t2=t2-meant2
		t10=t10-meant1
		t20=t20-meant2
	}

	dd=dnorm_p12dmgs(x=x,t1=t1,t2=t2,y=y,t10=t10,t20=t20,ics=ics,d1,d2,d3,d4,
		aderivs=aderivs)
	ru_cdf="rust not selected"
	bs_cdf="boot not selected"
	ml_params=dd$ml_params

	if(rust){
		th=tnorm_p12_cp("rust",nrust,x,t1,t2,nonnegslopesonly=rnonnegslopesonly)$theta_samples
		ru_cdf=numeric(length(y))
		for (ir in 1:nrust){
			mu=th[ir,1]+t10*th[ir,2]
			sg=exp(th[ir,3]+t20*th[ir,4])
			ru_cdf=ru_cdf+pnorm(y,mean=mu,sd=sg)
		}
		ru_cdf=ru_cdf/nrust
	}

	if(boot){
		th=tnorm_p12_cp("boot",nrust,x,t1,t2,nonnegslopesonly=rnonnegslopesonly)$theta_samples
		bs_cdf=numeric(length(y))
		for (ir in 1:nboot){
			mu=th[ir,1]+t10*th[ir,2]
			sg=exp(th[ir,3]+t20*th[ir,4])
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
#' @rdname Normal_p12_cp
#' @inheritParams man
#' @export
tnorm_p12_cp=function(method,n,x,t1,t2,nonnegslopesonly=FALSE,ics=c(0,0,0,0),
	d1=0.01,d2=0.01,d3=0.01,d4=0.01,debug=FALSE){

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
		th=bnorm_p12(x=x,t1=t1,t2=t2,n=n)
	} else{
		message("tnorm_p12 method not valid so stopping.\n")
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
