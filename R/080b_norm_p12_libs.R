#' Norm_p12 example data
#' @returns A list containing data to run an example
#' @param  iseed The random seed
norm_p12_exampledata=function(iseed){
	set.seed(iseed)
	nx=100
	t1=c(1:nx)/nx
	t2=rnorm(nx)
	n01=nx
	n02=nx
	t01=t1[nx]
	t02=t2[nx]
	params=c(0,1,0,1)
	x=rnorm(nx,mean=params[1]+params[2]*t1,sd=exp(params[3]+params[4]*t2))
	y=sort(params[1]+params[2]*t1) #shouldn't be sorted because the x are ordered
	return(list(nx=nx,x=x,t1=t1,t2=t2,n01=n01,n02=n02,t01=t01,t02=t02,params=params,y=y))
}
#' Waic
#' @inherit manwaic return
#' @inheritParams manf
norm_p12_waic=function(waicscores,x,t1,t2,v1hat,v2hat,v3hat,v4hat,lddi,lddd,lambdad){
		if(waicscores){
 			f1f=norm_p12_f1fw(x,t1,t2,v1hat,v2hat,v3hat,v4hat)
			f2f=norm_p12_f2fw(x,t1,t2,v1hat,v2hat,v3hat,v4hat)
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
#' @inherit manpredictor return
#' @inheritParams manf
norm_p12_predictordata=function(predictordata,x,t1,t2,t01,t02,params){
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
		mu0=a+b*t01
		sg0=exp(g+d*t02)
		qx=qnorm(px,mean=mu0,sd=sg0)
	} else{
		predictedparameter="predictordata not selected"
		adjustedx="predictordata not selected"
	}

	list(predictedparameter=list(mu=mu,sg=sg),adjustedx=qx)
}
#' Logf for RUST
#' @inherit manlogf return
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
#' @return Vector
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
#' @inherit manloglik return
#' @inheritParams manf
norm_p12_loglik=function(vv,x,t1,t2){
	n=length(x)
	mu=vv[1]+vv[2]*t1 #so mean is a vector, just like x
	sg=exp(vv[3]+vv[4]*t2)
	loglik=sum(dnorm(x,mean=mu,sd=sg,log=TRUE))
	return(loglik)
}
#' Check MLE
#' @inherit mancheckmle return
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
#' @returns Vector
#' @inheritParams manf
qnorm_p12=function(p,t01,t02,ymn,slope,sigma1,sigma2){

	mu=(ymn+slope*t01)
	sg=exp(sigma1+sigma2*t02)
	return(qnorm(p,mean=mu,sd=sg))

}
#' Normal-with-p12: Density function
#' @returns Vector
#' @inheritParams manf
dnorm_p12=function(x,t01,t02,ymn,slope,sigma1,sigma2,log=FALSE){

	mu=(ymn+slope*t01)
	sg=exp(sigma1+sigma2*t02)
	return(dnorm(x,mean=mu,sd=sg,log=log))

}
#' Normal-with-p12: Distribution function
#' @returns Vector
#' @inheritParams manf
pnorm_p12=function(y,t01,t02,ymn,slope,sigma1,sigma2){

	mu=(ymn+slope*t01)
	sg=exp(sigma1+sigma2*t02)
	return(pnorm(y,mean=mu,sd=sg))

}
#' Bootstrap
#' @inherit manboot return
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
#' Log scores for 5 predictions calculated using leave-one-out
#' @return Two scalars
#' @inheritParams manf
norm_p12_logscores=function(logscores,x,t1,t2,ics){

	if(logscores){
		nx=length(x)
		ml_oos_logscore=0
		cp_oos_logscore=0
		for (i in 1:nx){

			x1=x[-i]
			t1m1=t1[-i]
			t2m1=t2[-i]

			dd=dnorm_p12dmgs(x1,t1m1,t2m1,x[i],t1[i],t2[i],ics)

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
#' @inherit mandsub return
#' @inheritParams manf
dnorm_p12dmgs=function(x,t1,t2,y,t01,t02,ics){

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
		if(debug)message("inside dnorm_p12dmgs:ml_params=",ml_params)

# mle
		ml_pdf=dnorm_p12(y,t01,t02,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat)
		ml_cdf=pnorm_p12(y,t01,t02,ymn=v1hat,slope=v2hat,sigma1=v3hat,sigma2=v4hat)

# bayesian ingredients
		ldd		=norm_p12_ldda(x,t1,t2,v1hat,v2hat,v3hat,v4hat)
		lddi=solve(ldd)

		if(debug)message("det(ldd)=",det(ldd))
		lddd	=norm_p12_lddda(x,t1,t2,v1hat,v2hat,v3hat,v4hat)

		f1=norm_p12_f1fa(y,t01,t02,v1hat,v2hat,v3hat,v4hat)
		f2=norm_p12_f2fa(y,t01,t02,v1hat,v2hat,v3hat,v4hat)

		norm_p12_f2fa(y,t01,t02,v1hat,v2hat,v3hat,v4hat)

		p1=norm_p12_p1fa(y,t01,t02,v1hat,v2hat,v3hat,v4hat)
		p2=norm_p12_p2fa(y,t01,t02,v1hat,v2hat,v3hat,v4hat)

# pu
		lambdad_cp=c(0,0,0,0)
		if(debug)message("lddi=",lddi)
		if(debug)message("lddd=",lddd)
		if(debug)message("f1=",f1)
		if(debug)message("f2=",f2)
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

