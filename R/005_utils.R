#' If a variable is a vector, convert it to a matrix
#' @inheritParams manf
#'
#' @return Vector
ifvectorthenmatrix=function(t){

		if(is.vector(t))t=matrix(t)

	return(t)
}
#' Find the number of predictors in the predictor
#' @inheritParams manf
#'
#' @return Vector
findnt=function(t){

		if(is.vector(t))nt=1

		if(is.matrix(t)){
			nt=dim(t)[2]
		}

		if(nt>3){
			message("Max. number of columns in the predictor is 3.\n")
			stop()
		}

	return(nt)
}
#' Calculate the location parameter when there are predictors (single time point)
#' @inheritParams manf
#'
#' @return Vector
makebetat0=function(nt,params,t0){
		add=0
		for (i in 1:nt){
			add=add+params[i+1]*t0[i]
		}
		betat0=add
	return(betat0)
}
#' Calculate the location parameter when there are predictors (multiple time points)
#' @inheritParams manf
#'
#' @return Vector
makebetatm=function(nt,params,t){
		add=0
		for (i in 1:nt){
			add=add+params[i+1]*t[,i]
		}
		betat=add
	return(betat)
}
#' Deal with situations in which the user wants d or p outside the GEV range
#' @inheritParams manf
#'
#' @return Vector
fixgevrange=function(y,v1,v2,v3){
#
# I've added minxi now, so that even in the case where I adjust xi slightly
# in the 'movexiawayfromzero' routine, to avoid being too near zero when
# calculating analytical gradients, we still adjust the y and we don't
# get NA errors from the analytical gradient routines, from logp1
	minxi=10^-7
	if(v3<0){
#		ximax=v1-v2/v3
		ximax=v1-v2/(v3-minxi)
		for (i in 1:length(y)){
			if(y[i]>ximax)y[i]=ximax
		}
	}
	if(v3>0){
#		ximin=v1-v2/v3
		ximin=v1-v2/(v3+minxi)
		for (i in 1:length(y)){
			if(y[i]<ximin)y[i]=ximin
		}
	}
	return(y)
}
#' Deal with situations in which the user wants d or p outside the GPD range
#' @inheritParams manf
#'
#' @return Vector
fixgpdrange=function(y,v1,v2,v3){
	if(v3<0){
		ximin=v1
		ximax=v1-v2/v3
		for (i in 1:length(y)){
			if(y[i]<ximin)y[i]=ximin
			if(y[i]>ximax)y[i]=ximax
		}
	} else{
		for (i in 1:length(y)){
			if(y[i]<v1)y[i]=v1
		}
	}
	return(y)
}
#' Extract the results from derivatives and put them into f2
#' @param temp1 					output from derivative calculations
#' @param nx 							number of x values
#' @param dim 						number of parameters
#'
#' @return 3d array
deriv_copyfdd=function(temp1,nx,dim){
	f2=array(0,c(dim,dim,nx))
	for (i in 1:dim){
		for (j in 1:dim){
			f2[i,j,]=temp1[(i-1)*dim+j,]
		}
	}
	return(f2)
}
#' Extract the results from derivatives and put them into ldd
#' @param temp1 					output from derivative calculations
#' @param nx 							number of x values
#' @param dim 						number of parameters
#'
#' @return Matrix
deriv_copyldd=function(temp1,nx,dim){
	temp2=apply(temp1,1,sum)/nx
	ldd=matrix(0,dim,dim)
	for (i in 1:dim){
		for (j in 1:dim){
			ldd[i,j]=temp2[(i-1)*dim+j]
		}
	}
	return(ldd)
}
#' Extract the results from derivatives and put them into lddd
#' @param temp1 					output from derivative calculations
#' @param nx 							number of x values
#' @param dim 						number of parameters
#'
#' @return 3d array
deriv_copylddd=function(temp1,nx,dim){
	temp2=apply(temp1,1,sum)/nx
	lddd=array(0,c(dim,dim,dim))
	for (i in 1:dim){
		for (j in 1:dim){
			for (k in 1:dim){
				lddd[i,j,k]=temp2[(i-1)*dim*dim+(j-1)*dim+k]
			}
		}
	}
	return(lddd)
}
#' Extract the results from derivatives and put them into ldd
#' @param temp1 					output from derivative calculations
#' @param nx 							number of x values
#' @param dim 						number of parameters
#'
#' @return 3d array
deriv_copyld2=function(temp1,nx,dim){
	ld2=array(0,c(dim,dim,nx))
	for (i in 1:dim){
		for (j in 1:dim){
			ld2[i,j,]=temp1[(i-1)*dim+j,]
		}
	}
	return(ld2)
}
#' Move xi away from zero a bit
#' @param xi 					xi
#'
#' @return Scalar
movexiawayfromzero=function(xi){
	minxi=10^-7
	if(abs(xi)<minxi){
		if(xi>=0){
# the problem with this idea could be that increasing xi, when xi is positive,
# increases the minimum of the gev, and some y values might then lie below the minimum
# but I've made a change to the 'fixgevrange' routine to deal with that
# and that seems to fix the problem now
			xi=minxi
		} else {
			xi=-minxi
		}
	}
	return(xi)
}
#' Determine t0
#' @inheritParams manf
#'
#' @return Scalar
maket0=function(t0,n0,t){

# if t0 is specified, does nothing
# if t0 isn't specified, calculates t0 from n0

	if( (is.na(t0))  && (is.na(n0))  ){message("Either t0 or n0 must be specified")}
	if( (!is.na(t0)) && (!is.na(n0)) ){message("Only one of t0 or n0 must be specified")}

	if(is.na(t0))t0=t[n0]

	return(t0)
}
#' Make muhat0
#' @param t0					the value of the predictor vector at which to make the prediction (if n0 not specified)
#' @param n0					the position in the predictor vector at which to make the prediction (positive integer less than or equal to the length of \eqn{x}) (if t0 not specified)
#' @param t 					predictor
#' @param mle_params	MLE params
#'
#' @return Scalar
makemuhat0=function(t0,n0,t,mle_params){

	muhat=mle_params[1]+mle_params[2]*t

	if(is.na(t0)){
		muhat0=muhat[n0]
	} else {
		muhat0=mle_params[1]+mle_params[2]*t0
	}

	return(muhat0)
}
#' Make ta0
#' @param t0					the value of the predictor vector at which to make the prediction (if n0 not specified)
#' @param n0					the position in the predictor vector at which to make the prediction (positive integer less than or equal to the length of \eqn{x}) (if t0 not specified)
#' @param t 					predictor
#'
#' @return Scalar
maketresid0=function(t0,n0,t){

	tresid=t-mean(t)

	if(is.na(t0)){
		tresid0=tresid[n0]
	} else {
		tresid0=t0-mean(t)
	}

	return(tresid0)
}
#' Make Standard Errors from lddi
#' @param nx					length of training data
#' @param lddi				the inverse log-likelihood matrix
#'
#' @return Vector
make_se=function(nx,lddi){
	nd=dim(lddi)[1]
	standard_errors=matrix(0,nd)
	for (i in 1:nd){
		if(lddi[i,i]>0){
			standard_errors[i]="square root not possible"
		} else{
			standard_errors[i]=sqrt(-lddi[i,i]/nx)
		}
	}
	return(standard_errors)

}
#' Make WAIC
#' @param 	x						the training data
#' @param 	fhatx				density of x at the maximum likelihood parameters
#' @param 	lddi				inverse of the second derivative log-likelihood matrix
#' @param 	lddd				the third derivative log-likelihood tensor
#' @param 	f1f					the f1 term from DMGS equation 2.1
#' @param 	lambdad	the slope of the log prior
#' @param 	f2f					the f2 term from DMGS equation 2.1
#' @param 	dim					number of free parameters
#'
#' @return Two scalars
make_waic=function(x,fhatx,lddi,lddd,f1f,lambdad,f2f,dim){

# waic seems to be a bust because it doesn't make sense for the mle models...it doesn't penalize them at all for parameters
# the variance version always penalizes twice as much as the difference version...is there a factor of 2 wrong somewhere?

# x is implicit in f1f and f2f

	nx=length(x)
	lfhatx=log(fhatx)
	fhatxsq=fhatx*fhatx
	lfhatxof=lfhatx/fhatx
	omlfhatx=1-lfhatx

	outerf1f=array(0,c(dim,dim,nx))
	for (i in 1:nx){
		for (j in 1:dim){
			for (k in 1:dim){
				outerf1f[j,k,i]=outer(f1f[j,i],f1f[k,i])
			}
		}
	}

# first term is the llpd
# which is just the in-sample log-likelihood
# for mle, should match the log-likelihood at the max
# calc the probabilities, take the log, add them together
	t1=array(0,c(dim,nx))
	t2=array(0,c(dim,dim,nx))
	for (i in 1:nx){
		t1[,i]	=f1f[,i]
		t2[,,i]	=f2f[,,i]
	}
	dq=dmgs(lddi,lddd,t1,lambdad,t2,dim=dim)/nx
	dqq=pmax(fhatx+dq,.Machine$double.eps)
	logf1_rhp=log(dqq)
	sumlogf1_rhp=sum(logf1_rhp)
	lppd_rhp=sumlogf1_rhp
#	logf1_mle=log(fhatx)
#	sumlogf1_mle=sum(logf1_mle)
#	lppd_mle=sumlogf1_mle #no point in calculating for mle because the adjustment is zero

# second term is the penalty term
# this version of it is a difference
# seems to be the same for maxlik
	t1=array(0,c(dim,nx))
	t2=array(0,c(dim,dim,nx))
	for (i in 1:nx){
		t1[,i]	=f1f[,i]/fhatx[i]
		t2[,,i]	=f2f[,,i]/fhatx[i]-outerf1f[,,i]/fhatxsq[i]
	}
	logf2_rhp=lfhatx+dmgs(lddi,lddd,t1,lambdad,t2,dim=dim)/nx
	sumlogf2_rhp=sum(logf2_rhp)
	pwaic1_rhp=2*(sumlogf1_rhp-sumlogf2_rhp)
#	logf2_mle=log(fhatx)
#	sumlogf2_mle=sum(logf2_mle)
#	pwaic1_mle=2*(sumlogf1_mle-sumlogf2_mle)

# alternative version uses the variance thing for the second term

# this version of it is a variance
# it penalizes more...is there really a factor of 2 needed
	t2inner=array(0,c(dim,dim,nx))
	for (i in 1:nx){
		t1[,i]=2*lfhatxof[i]*f1f[,i]
		t2inner[,,i]=f2f[,,i]*lfhatx[i]+outerf1f[,,i]*omlfhatx[i]/fhatx[i]
		t2[,,i]=2*t2inner[,,i]/fhatx[i]
	}
	dq=dmgs(lddi,lddd,t1,lambdad,t2,dim=dim)
	sumvlogf_rhp=sum(lfhatx*lfhatx+dq/nx-logf2_rhp*logf2_rhp)
	pwaic2_rhp=sumvlogf_rhp
#	sumvlogf_mle=sum(lfhatx*lfhatx-logf2_mle*logf2_mle)
#	pwaic2_mle=sumvlogf_mle
#	waic_mle=matrix(0,5)

	waic1=numeric(3)
	waic2=numeric(3)
	waic1[1]=lppd_rhp
	waic1[2]=-pwaic1_rhp
	waic1[3]=waic1[1]+waic1[2]
	waic2[1]=lppd_rhp
	waic2[2]=-pwaic2_rhp
	waic2[3]=waic2[1]+waic2[2]

	list(waic1=waic1,waic2=waic2)
}
#' Make WAIC
#' @param 	x						the training data
#' @param 	fhatx				density of x at the maximum likelihood parameters
#' @param 	lddi				inverse of the second derivative log-likelihood matrix
#' @param 	lddd				the third derivative log-likelihood tensor
#' @param 	f1f					the f1 term from DMGS equation 2.1
#' @param 	lambdad	the slope of the log prior
#' @param 	f2f					the f2 term from DMGS equation 2.1
#' @param 	dim					number of free parameters
#'
#' @return Two scalars
make_cwaic=function(x,fhatx,lddi,lddd,f1f,lambdad,f2f,dim){

# waic seems to be a bust because it doesn't make sense for the mle models...it doesn't penalize them at all for parameters
# the variance version always penalizes twice as much as the difference version...is there a factor of 2 wrong somewhere?

# x is implicit in f1f and f2f

	nx=length(x)
	lfhatx=log(fhatx)
	fhatxsq=fhatx*fhatx
	lfhatxof=lfhatx/fhatx
	omlfhatx=1-lfhatx

	outerf1f=array(0,c(dim,dim,nx))
	for (i in 1:nx){
		for (j in 1:dim){
			for (k in 1:dim){
				outerf1f[j,k,i]=outer(f1f[j,i],f1f[k,i])
			}
		}
	}

# first term is the llpd
# which is just the in-sample log-likelihood
# for mle, should match the log-likelihood at the max
# calc the probabilities, take the log, add them together
	t1=array(0,c(dim,nx))
	t2=array(0,c(dim,dim,nx))
	for (i in 1:nx){
		t1[,i]	=f1f[,i]
		t2[,,i]	=f2f[,,i]
	}
	dq=dmgs(lddi,lddd,t1,lambdad,t2,dim=dim)/nx
	dqq=pmax(fhatx+dq,.Machine$double.eps)
	logf1_rhp=log(dqq)
	sumlogf1_rhp=sum(logf1_rhp)
	lppd_rhp=sumlogf1_rhp

# second term is the penalty term
# this version of it is a difference
# seems to be the same for maxlik
	t1=array(0,c(dim,nx))
	t2=array(0,c(dim,dim,nx))
	for (i in 1:nx){
		t1[,i]	=f1f[,i]/fhatx[i]
		t2[,,i]	=f2f[,,i]/fhatx[i]-outerf1f[,,i]/fhatxsq[i]
	}
	logf2_rhp=lfhatx+dmgs(lddi,lddd,t1,lambdad,t2,dim=dim)/nx
	sumlogf2_rhp=sum(logf2_rhp)
	pwaic1_rhp=2*(sumlogf1_rhp-sumlogf2_rhp)
#	logf2_mle=log(fhatx)
#	sumlogf2_mle=sum(logf2_mle)
#	pwaic1_mle=2*(sumlogf1_mle-sumlogf2_mle)

# alternative version uses the variance thing for the second term

# this version of it is a variance
# it penalizes more...is there really a factor of 2 needed
	t2inner=array(0,c(dim,dim,nx))
	for (i in 1:nx){
		t1[,i]=2*lfhatxof[i]*f1f[,i]
		t2inner[,,i]=f2f[,,i]*lfhatx[i]+outerf1f[,,i]*omlfhatx[i]/fhatx[i]
		t2[,,i]=2*t2inner[,,i]/fhatx[i]
	}
	dq=dmgs(lddi,lddd,t1,lambdad,t2,dim=dim)
	sumvlogf_rhp=sum(lfhatx*lfhatx+dq/nx-logf2_rhp*logf2_rhp)
	pwaic2_rhp=sumvlogf_rhp
#	sumvlogf_mle=sum(lfhatx*lfhatx-logf2_mle*logf2_mle)
#	pwaic2_mle=sumvlogf_mle
#	waic_mle=matrix(0,5)

	waic1=numeric(3)
	waic2=numeric(3)
	waic1[1]=lppd_rhp
	waic1[2]=-pwaic1_rhp
	waic1[3]=waic1[1]+waic1[2]
	waic2[1]=lppd_rhp
	waic2[2]=-pwaic2_rhp
	waic2[3]=waic2[1]+waic2[2]

	list(waic1=waic1,waic2=waic2)
}
#' Calculate MAIC
#' @param ml_value	maximum of the likelihood
#' @param nparams		number of parameters
#'
#' @return Vector of 3 values
#' Returns the two compoments of MAIC, and their sum
make_maic=function(ml_value,nparams){
	maic=numeric(3)
	maic[1]=ml_value
	maic[2]=-nparams
	maic[3]=maic[1]+maic[2]
	return(maic)
}

#' Generates a comment about the method
#'
#' @return String
rust_pumethod=function(){

  method="The cp results are based posterior simulation using
  				ratio of uniforms sampling (RUST),
  				using the predictive matching right Haar prior."
	method=strwrap(method,width=10000,simplify=TRUE)

	return(method)
}
#' Generates a comment about the method
#' @return String
analytic_cpmethod=function(){

  method="The cp results are based on an analytic solution of
					the Bayesian prediction integral,
					using the predictive matching right Haar prior."
	method=strwrap(method,width=10000,simplify=TRUE)

	return(method)
}
#' Generates a comment about the method
#' @return String
rhp_dmgs_cpmethod=function(){

  method="The cp results are based on the DMGS approximation of
					the Bayesian prediction integral,
					using the predictive matching right Haar prior."
	method=strwrap(method,width=10000,simplify=TRUE)

	return(method)
}
#' Generates a comment about the method
#' @return String
crhpflat_dmgs_cpmethod=function(){

  method="The cp results are based on the DMGS approximation of
					the Bayesian prediction integral,
					using the CRHP/flat prior."
	method=strwrap(method,width=10000,simplify=TRUE)

	return(method)
}
#' Generates a comment about the method
#' @return String
adhoc_dmgs_cpmethod=function(){

  method="The cp results are based on the DMGS approximation of
					the Bayesian prediction integral,
					using an ad-hoc prior."
	method=strwrap(method,width=10000,simplify=TRUE)

	return(method)
}
#' Calculates quantiles from simulations by inverting the Hazen CDF
#' @inheritParams manf
#' @return Vector
makeq=function(yy,pp){
	nyy=length(yy)
	npp=length(pp)
	syy=sort(yy)
	lambda=nyy*pp+0.5
	qq=numeric(npp)
	aa=floor(lambda)
	eps=lambda-aa
	for (ii in 1:npp){
		if(lambda[ii]<=1){
			qq[ii]=syy[1]
		}else if(lambda[ii]>=nyy){
			qq[ii]=syy[nyy]
		}else{
			qq[ii]=(1-eps[ii])*syy[aa[ii]]+eps[ii]*syy[aa[ii]+1]
		}
	}
	return(qq)
}
#' Message to explain why GEV and GPD \code{d***} and \code{p***} routines
#' don't return DMGS pdfs and cdfs
#' @inheritParams manf
#' @return String
nopdfcdfmsg=function(yy,pp){

	msg="For the pdf and cdf for the GEV and GPD, considered as a function of the rv, DMGS doesn't return anything."
	msg=paste(msg,"That's because DMGS can't work outside the range limits.")
	msg=paste(msg,"Instead, you could use:")
	msg=paste(msg,"(a) q***, for the inverse CDF, or")
	msg=paste(msg,"(b) q***, with pdf=TRUE, which returns the pdf as a function of probability, or")
	msg=paste(msg,"(c) d*** or p***, with rust=TRUE, which uses rust to generate the pdf and cdf.")

	return(msg)
}

