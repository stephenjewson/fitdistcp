#' Waicscores
#' @inherit manwaic return
#' @inheritParams manf
exp_waic=function(waicscores,x,v1hat){
	if(waicscores){
		f1f=exp_f1fa(x,v1hat)
		f2f=exp_f2fa(x,v1hat)

		ldd=exp_ldda(x,v1hat)
		lddi=solve(ldd)
		lddd=exp_lddda(x,v1hat)
		fhatx=dexp(x,rate=v1hat)
		lambdad_rhp=-1/v1hat
		waic=make_waic(x,fhatx,lddi,lddd,f1f,lambdad_rhp,f2f,dim=1)
		waic1=waic$waic1
		waic2=waic$waic2
	}else{
		waic1="waicscores not selected"
		waic2="waicscores not selected"
	}
	return(list(waic1=waic1,waic2=waic2))
}
#' Logf for RUST
#' @inherit manlogf return
#' @inheritParams manf
exp_logf=function(params,x){
	l=pmax(params[1],.Machine$double.eps)
	logf=sum(dexp(x,rate=l,log=TRUE))-log(l)
	return(logf)
}
#' Log scores for MLE and RHP predictions calculated using leave-one-out
#' @inherit manlogscores return
#' @inheritParams manf
exp_logscores=function(logscores,x){

	if(logscores){
# I could put the logs inside dexpsub, but I'd have toa actually calculate the log for the rhp case

		nx=length(x)
		ml_oos_logscore=0
		rh_oos_logscore=0
		for (i in 1:nx){
			x1=x[-i]
			dd=dexpsub(x1,x[i])
			ml_params1=(nx-1)/sum(x1)
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
dexpsub=function(x,y){

		nx=length(x)

# ml
		ml_params=nx/sum(x)
		ml_pdf=dexp(y,rate=ml_params)
		ml_cdf=pexp(y,rate=ml_params)

# rhp pdf
		sx=sum(x)
		top=sx**nx
		bot1=(sx+y)**(nx+1)
		rh_pdf=nx*top/bot1
		rh_pdf=pmax(rh_pdf,0)

# rhp cdf
		bot2=(sx+y)**nx
		rh_cdf=1-top/bot2
		rh_cdf=pmin(pmax(rh_cdf,0),1)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					rh_pdf=rh_pdf,
					ml_cdf=ml_cdf,
					rh_cdf=rh_cdf)
}
