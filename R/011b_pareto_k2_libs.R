# extraDistr is very confusing:
# a = shape parameter (they call scale), that I'm varying here
# b = scale parameter (they call location), that I'm calling kscale
#
#' Waic
#' @inherit manwaic return
#' @inheritParams manf
pareto_k2_waic=function(waicscores,x,v1hat,kscale){
	if(waicscores){
		f1f=pareto_k2_f1fa(x,v1hat,kscale=kscale)
		f2f=pareto_k2_f2fa(x,v1hat,kscale=kscale)
		ldd=pareto_k2_ldda(x,v1hat,kscale=kscale)
		lddi=solve(ldd)
		lddd=pareto_k2_lddda(x,v1hat,kscale=kscale)
		fhatx=extraDistr::dpareto(x,a=v1hat,b=kscale)
		lambdad_rhp=-1/v1hat
		waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad_rhp,f2f,dim=1)
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
pareto_k2_logf=function(params,x,kscale){
	sh=pmax(params[1],.Machine$double.eps)
	logf=sum(extraDistr::dpareto(x,a=sh,b=kscale,log=TRUE))-log(sh)
	return(logf)
}
#' Maximum likelihood estimator
#' @inherit manloglik return
#' @inheritParams manf
pareto_k2_ml_params=function(x,kscale){
	nx=length(x)
#	lnx=log(x)
	lnx=log(x/kscale)
	slnx=sum(lnx)
	mlparams=nx/slnx
	return(mlparams)
}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inherit manlogscores return
#' @inheritParams manf
pareto_k2_logscores=function(logscores,x,kscale){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){
			x1=x[-i]
			dd=dpareto_k2_sub(x1,x[i],kscale)
			ml_params1=dd$ml_params
# ml
			ml_pdf=dd$ml_pdf
			ml_oos_logscore=ml_oos_logscore+log(ml_pdf)
# rhp
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
dpareto_k2_sub=function(x,y,kscale){

		nx=length(x)

# ml
		ml_params=pareto_k2_ml_params(x,kscale)
		ml_pdf=extraDistr::dpareto(y,a=ml_params,b=kscale)
		ml_cdf=extraDistr::ppareto(y,a=ml_params,b=kscale)

# rhp
		slnx=sum(log(x/kscale))
		bigx=slnx
		top=nx*(bigx**nx)
		bigz=log(y/kscale)+slnx
		bot=y*((bigz**(nx+1)))
		rh_pdf=top/bot
		rh_pdf=pmax(rh_pdf,0)

		rh_cdf=1-(bigx**nx)/(bigz**nx)
		rh_cdf=pmin(pmax(rh_cdf,0),1)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					rh_pdf=rh_pdf,
					ml_cdf=ml_cdf,
					rh_cdf=rh_cdf)
}

