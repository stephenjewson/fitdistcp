#' Set default params for the chosen model
#' @param model			which distribution to test. Possibles values are
#'	"\code{exp}",
#'	"\code{pareto_k2}",
#'	"\code{halfnorm}",
#'	"\code{unif}",
#'	"\code{norm}",
#'	"\code{norm_dmgs}",
#'	"\code{gnorm_k3}",
#'	"\code{lnorm}",
#'	"\code{lnorm_dmgs}",
#'	"\code{logis}",
#'	"\code{lst_k3}",
#'	"\code{cauchy}",
#'	"\code{gumbel}",
#'	"\code{frechet_k1}",
#'	"\code{weibull}",
#'	"\code{gev_k3}",
#'	"\code{exp_p1}",
#'	"\code{pareto_p1k2}",
#'	"\code{norm_p1}",
#'	"\code{lnorm_p1}",
#'	"\code{logis_p1}",
#'	"\code{lst_p1k3}",
#'	"\code{cauchy_p1}",
#'	"\code{gumbel_p1}",
#'	"\code{frechet_p2k1}",
#'	"\code{weibull_p2}",
#'	"\code{gev_p1k3}",
#'	"\code{norm_p12}",
#'	"\code{lst_p12k3}",
#'	"\code{gamma}",
#'	"\code{invgamma}",
#'	"\code{invgauss}",
#'	"\code{gev}",
#'	"\code{gpd_k1}",
#'	"\code{gev_p1}".
#'	"\code{gev_p12}".
#'	"\code{gev_p123}".
#' @param params		values for the parameters for the specified distribution
#'
#' @return Vector
reltest_params=function(model="exp",params){

	if(is.na(params[1])){
		if(model=="exp")							{params=c(1)}
		if(model=="pareto_k1")				{params=c(1)}
		if(model=="halfnorm")					{params=c(1)}
		if(model=="norm")							{params=c(0,1)}
		if(model=="norm_dmgs")				{params=c(0,1)}
		if(model=="gnorm_k3")					{params=c(0,1,4)}
    if(model=="lnorm")						{params=c(0,1)}
    if(model=="lnorm_dmgs")				{params=c(0,1)}
    if(model=="logis")						{params=c(0,1)}
    if(model=="lst_k3")						{params=c(0,1,10)}
    if(model=="cauchy")						{params=c(0,1)}
    if(model=="gumbel")						{params=c(0,1)}
    if(model=="frechet_k1")				{params=c(1,2)} #lambda =1 has infinite mean
    if(model=="weibull")					{params=c(1,1)}
    if(model=="gev_k3")						{params=c(0,1,0.1)}
    if(model=="exp_p1")						{params=c(0,1)}
    if(model=="pareto_p1k3")			{params=c(0,1)}
    if(model=="norm_p1")					{params=c(0,1,1)}
    if(model=="lnorm_p1")					{params=c(0,1,1)}
    if(model=="logis_p1")					{params=c(0,1,2)}
    if(model=="lst_p1k3") 				{params=c(0,1,1,10)}
    if(model=="cauchy_p1") 				{params=c(0,1,2)}
    if(model=="gumbel_p1")				{params=c(0,1,1)}
    if(model=="frechet_p2k1")			{params=c(0,0.001,2)}
    if(model=="weibull_p1")				{params=c(0,0,1)}
    if(model=="weibull_p2")				{params=c(2,0,0)}
    if(model=="gev_p1k3")					{params=c(0,1,1,0.1)}
    if(model=="cweibull_p2")			{params=c(2,0,0)}
    if(model=="norm_p12")					{params=c(0,1,0,1)} 	# params 2 and 4 can't be zero
    if(model=="lst_p12k5")				{params=c(0,1,0,1,50)}# params 2 and 4 can't be zero...but rust gives errors
    if(model=="weibull_p12")			{params=c(0,1,0,1)}		# params 2 and 4 can't be zero
    if(model=="gamma")						{params=c(1,1)}
    if(model=="invgamma")					{params=c(1,1)}
    if(model=="invgauss")					{params=c(2,2)}
    if(model=="gev")							{params=c(0,1,0.1)}
    if(model=="gpd_k1")						{params=c(1,0.1)}
    if(model=="gev_p1")						{params=c(0,1,1,0.1)}
    if(model=="gev_p12")					{params=c(0,0.0001,0.01,0.01,0.1)}
    if(model=="gev_p123")					{params=c(0,1,0,0.01,0.1,0)}
	}

	return(params)
}
#' Random training data from one model
#' @param model			which distribution to test. Possibles values are
#'	"\code{exp}",
#'	"\code{pareto_k2}",
#'	"\code{halfnorm}",
#'	"\code{unif}",
#'	"\code{norm}",
#'	"\code{norm_dmgs}",
#'	"\code{gnorm_k3}",
#'	"\code{lnorm}",
#'	"\code{lnorm_dmgs}",
#'	"\code{logis}",
#'	"\code{lst_k3}",
#'	"\code{cauchy}",
#'	"\code{gumbel}",
#'	"\code{frechet_k1}",
#'	"\code{weibull}",
#'	"\code{gev_k3}",
#'	"\code{exp_p1}",
#'	"\code{pareto_p1k2}",
#'	"\code{norm_p1}",
#'	"\code{lnorm_p1}",
#'	"\code{logis_p1}",
#'	"\code{lst_p1k3}",
#'	"\code{cauchy_p1}",
#'	"\code{gumbel_p1}",
#'	"\code{frechet_p2k1}",
#'	"\code{weibull_p2}",
#'	"\code{gev_p1k3}",
#'	"\code{norm_p12}",
#'	"\code{lst_p12k3}",
#'	"\code{gamma}",
#'	"\code{invgamma}",
#'	"\code{invgauss}",
#'	"\code{gev}",
#'	"\code{gpd_k1}",
#'	"\code{gev_p1}".
#'	"\code{gev_p12}".
#'	"\code{gev_p123}".
#' @param nx				the length of the training data to use.
#' @param tt				predictor vector
#' @param tt1				predictor vector 1
#' @param tt2				predictor vector 2
#' @param tt3				predictor vector 2
#' @param params		values for the parameters for the specified distribution
#' @param minxi			minimum value for EVT shape parameter
#' @param maxxi			maximum value for EVT shape parameter
#'
#' @return Vector
reltest_simulate=function(model="exp",nx=20,tt,tt1,tt2,tt3,params,minxi=-10,maxxi=-10){

	if(model=="exp")						xx=rexp			(nx,rate=params[1])
	if(model=="pareto_k2")			xx=rpareto	(nx,a=params[2],b=params[1])
	if(model=="halfnorm")				xx=rhalfnorm(nx,theta=params[1])
	if(model=="unif")						xx=runif		(nx,min=params[1],max=params[2])
	if(model=="norm")						xx=rnorm		(nx,mean=params[1],sd=params[2])
	if(model=="norm_dmgs")			xx=rnorm		(nx,mean=params[1],sd=params[2])
	if(model=="gnorm_k3")				xx=rgnorm		(nx,mu=params[1],alpha=params[2],beta=4)
	if(model=="lnorm")					xx=rlnorm		(nx,meanlog=params[1],sdlog=params[2])
	if(model=="lnorm_dmgs")			xx=rlnorm		(nx,meanlog=params[1],sdlog=params[2])
	if(model=="logis")					xx=rlogis		(nx,location=params[1],scale=params[2])
	if(model=="lst_k3")					xx=rlst			(nx,mu=params[1],sigma=params[2],df=params[3])
	if(model=="cauchy")					xx=rcauchy	(nx,location=params[1],scale=params[2])
	if(model=="gumbel")					xx=rgumbel	(nx,mu=params[1],sigma=params[2])
	if(model=="frechet_k1")			xx=rfrechet	(nx,lambda=params[3],sigma=params[2],mu=params[1]) #param ordering confusion
	if(model=="weibull")				xx=rweibull	(nx,shape=params[1],scale=params[2])
	if(model=="gev_k3")					xx=rgev			(nx,mu=params[1],sigma=params[2],xi=params[3])
	if(model=="gamma")					xx=rgamma		(nx,shape=params[1],scale=params[2])
	if(model=="invgamma")				xx=rinvgamma(nx,shape=params[1],scale=params[2])
	if(model=="invgauss")				xx=rinvgauss(nx,mean=params[1],shape=params[2])
  if(model=="gev")						xx=rgev_minmax(nx,mu=params[1],sigma=params[2],xi=params[3],minxi=minxi,maxxi=maxxi)
 # note that the way this routine works, params[1] really is kloc!! Don't think that's wrong.
	if(model=="gpd_k1")					xx=rgpd_k1_minmax(nx,kloc=params[1],sigma=params[2],xi=params[3],minxi=minxi,maxxi=maxxi)

#
# slightly differently for p1 models
#
	if(model=="exp_p1")				xx=rexp			(nx,rate=1/exp(params[1]+params[2]*tt))
	if(model=="pareto_p1k2")	xx=rpareto	(nx,a=1/exp(params[2]+params[3]*tt),b=params[1])
	if(model=="norm_p1")			xx=rnorm		(nx,mean=params[1]+params[2]*tt,sd=params[3])
	if(model=="lnorm_p1")			xx=rlnorm		(nx,meanlog=params[1]+params[2]*tt,sdlog=params[3])
	if(model=="logis_p1")			xx=rlogis		(nx,location=params[1]+params[2]*tt,scale=params[3])
	if(model=="lst_p1k3")			xx=rlst			(nx,mu=params[1]+params[2]*tt,sigma=params[3],df=params[4])
	if(model=="cauchy_p1")		xx=rcauchy	(nx,location=params[1]+params[2]*tt,scale=params[3])
	if(model=="gumbel_p1")		xx=rgumbel	(nx,mu=params[1]+params[2]*tt,sigma=params[3])
	if(model=="frechet_p2k1")	xx=rfrechet	(nx,mu=params[1],sigma=exp(params[2]+params[3]*tt),lambda=params[4])
	if(model=="weibull_p2")		xx=rweibull	(nx,shape=params[1],scale=exp(params[2]+params[3]*tt))
# Note that params[4] really is kxi...don't think that's a mistake
	if(model=="gev_p1k3")			xx=rgev_p1_minmax(nx,mu=params[1]+params[2]*tt,sigma=params[3],xi=params[4],tt,minxi=minxi,maxxi=maxxi)
#	if(model=="norm_p12")			xx=rnorm		(nx,mean=params[1]+params[2]*tt1,sd   =exp(params[3]+params[4]*tt2))
#	if(model=="lst_p12k3")		xx=rlst			(nx,mu  =params[1]+params[2]*tt1,sigma=exp(params[3]+params[4]*tt2),df=params[5])
	if(model=="gev_p1")				xx=rgev_p1_minmax(nx,mu=params[1]+params[2]*tt1,sigma=params[3],xi=params[4],tt1,minxi=minxi,maxxi=maxxi)

	if(model=="gev_p12"){
		xx=rgev_p12_minmax(nx,mu=params[1]+params[2]*tt1,sigma=exp(params[3]+params[4]*tt2),xi=params[5],tt1,tt2,minxi=minxi,maxxi=maxxi)
	}

	if(model=="gev_p123")			xx=rgev_p123_minmax(nx,mu=params[1]+params[2]*tt1,sigma=exp(params[3]+params[4]*tt2),xi=params[5]+params[6]*tt3,tt1,tt2,tt3,minxi=minxi,maxxi=maxxi)

	return(xx)
}
#' Make prediction from one model
#' @param model			which distribution to test. Possibles values are
#'	"\code{exp}",
#'	"\code{pareto_k2}",
#'	"\code{halfnorm}",
#'	"\code{unif}",
#'	"\code{norm}",
#'	"\code{norm_dmgs}",
#'	"\code{gnorm_k3}",
#'	"\code{lnorm}",
#'	"\code{lnorm_dmgs}",
#'	"\code{logis}",
#'	"\code{lst_k3}",
#'	"\code{cauchy}",
#'	"\code{gumbel}",
#'	"\code{frechet_k1}",
#'	"\code{weibull}",
#'	"\code{gev_k3}",
#'	"\code{exp_p1}",
#'	"\code{pareto_p1k2}",
#'	"\code{norm_p1}",
#'	"\code{lnorm_p1}",
#'	"\code{logis_p1}",
#'	"\code{lst_p1k3}",
#'	"\code{cauchy_p1}",
#'	"\code{gumbel_p1}",
#'	"\code{frechet_p2k1}",
#'	"\code{weibull_p2}",
#'	"\code{exp_p1k4}",
#'	"\code{norm_p12}",
#'	"\code{lst_p12k3}",
#'	"\code{gamma}",
#'	"\code{invgamma}",
#'	"\code{invgauss}",
#'	"\code{gev}",
#'	"\code{gpd_k1}",
#'	"\code{gev_p1}".
#'	"\code{gev_p12}".
#'	"\code{gev_p123}".
#' @param xx				training data
#' @param tt				predictor vector
#' @param tt1				predictor vector 1
#' @param tt2				predictor vector 2
#' @param tt3				predictor vector 3
#' @param n0				index for predictor vector
#' @param n10				index for predictor vector 1
#' @param n20				index for predictor vector 2
#' @param n30				index for predictor vector 2
#' @param pp				probabilites at which to make quantile predictions
#' @param params		model parameters
#' @param dmgs			flag for whether to run dmgs calculations or not
#' @param	debug			flag for turning debug messages on
#' @param aderivs		a logical for whether to use analytic derivatives (instead of numerical)
#' @param unbiasedv a logical for whether to use the unbiased variance instead of maxlik (for the normal)
#' @param pwm				a logical for whether to use PWM instead of maxlik (for the GEV)
#' @param minxi			minimum value for EVT shape parameter
#' @param maxxi			maximum value for EVT shape parameter
#'
#' @return Two vectors
reltest_predict=function(model,xx,tt,tt1,tt2,tt3,n0,n10,n20,n30,pp,params,dmgs=TRUE,
	debug=FALSE,aderivs=TRUE,unbiasedv=FALSE,pwm=FALSE,minxi=-10,maxxi=10){

	nalpha=length(pp)
	nmethods=2
#	revert2ml=FALSE
#
# make predictions
#
	if(model=="exp")						pred0=qexp_cp(xx,p=pp,debug=debug,aderivs=aderivs)
	if(model=="pareto_k2")			pred0=qpareto_k2_cp(xx,p=pp,kscale=params[1],debug=debug,aderivs=aderivs)
	if(model=="halfnorm")				pred0=qhalfnorm_cp(xx,p=pp,dmgs=dmgs,debug=debug,aderivs=aderivs)

	if(model=="unif")						pred0=qunif_cp(xx,p=pp,debug=debug,aderivs=aderivs)
	if(model=="norm")						pred0=qnorm_cp(xx,p=pp,debug=debug,aderivs=aderivs,unbiasedv=unbiasedv)
	if(model=="norm_dmgs")			pred0=qnorm_dmgs_cp(xx,p=pp,dmgs=dmgs,debug=debug,aderivs=aderivs)
	if(model=="lnorm")					pred0=qlnorm_cp(xx,p=pp,debug=debug,aderivs=aderivs)
	if(model=="lnorm_dmgs")			pred0=qlnorm_dmgs_cp(xx,p=pp,dmgs=dmgs,debug=debug,aderivs=aderivs)
	if(model=="gnorm_k3")				pred0=qgnorm_k3_cp(xx,p=pp,kbeta=4,debug=debug,aderivs=aderivs)

	if(model=="logis")					pred0=qlogis_cp(xx,p=pp,dmgs=dmgs,debug=debug,aderivs=aderivs)
	if(model=="lst_k3")					pred0=qlst_k3_cp(xx,p=pp,kdf=params[3],dmgs=dmgs,debug=debug,aderivs=aderivs)
	if(model=="cauchy")					pred0=qcauchy_cp(xx,p=pp,dmgs=dmgs,debug=debug,aderivs=aderivs)

	if(model=="gumbel")					pred0=qgumbel_cp(xx,p=pp,dmgs=dmgs,debug=debug,aderivs=aderivs)
	if(model=="frechet_k1")			pred0=qfrechet_k1_cp(xx,p=pp,kloc=params[1],
																dmgs=dmgs,debug=debug,aderivs=aderivs)
	if(model=="weibull")				pred0=qweibull_cp(xx,p=pp,dmgs=dmgs,debug=debug,aderivs=aderivs)
	if(model=="gev_k3")					pred0=qgev_k3_cp(xx,p=pp,kshape=params[3],
																dmgs=dmgs,debug=debug)

	if(model=="exp_p1")					pred0=qexp_p1_cp(xx,tt,n0=n0,p=pp,dmgs=dmgs,aderivs=aderivs)
	if(model=="pareto_p1k2")		pred0=qpareto_p1k2_cp(xx,tt,n0=n0,p=pp,
																kscale=params[1],dmgs=dmgs,debug=debug,aderivs=aderivs)

	if(model=="norm_p1")				pred0=qnorm_p1_cp(xx,tt,n0=n0,p=pp,debug=debug,aderivs=aderivs)
	if(model=="lnorm_p1")				pred0=qlnorm_p1_cp(xx,tt,n0=n0,p=pp,debug=debug,aderivs=aderivs)

	if(model=="logis_p1")				pred0=qlogis_p1_cp(xx,tt,n0=n0,p=pp,dmgs=dmgs,debug=debug,aderivs=aderivs)
	if(model=="lst_p1k3")				pred0=qlst_p1k3_cp(xx,tt,n0=n0,p=pp,kdf=params[4],
																dmgs=dmgs,debug=debug,aderivs=aderivs)
	if(model=="cauchy_p1")			pred0=qcauchy_p1_cp(xx,tt,n0=n0,p=pp,dmgs=dmgs,debug=debug,aderivs=aderivs)

	if(model=="gumbel_p1")			pred0=qgumbel_p1_cp(xx,tt,n0=n0,p=pp,dmgs=dmgs,debug=debug,aderivs=aderivs)
	if(model=="frechet_p2k1")		pred0=qfrechet_p2k1_cp(xx,tt,n0=n0,p=pp,
																kloc=params[1],dmgs=dmgs,debug=debug,aderivs=aderivs)
	if(model=="weibull_p2")			pred0=qweibull_p2_cp(xx,tt,n0=n0,p=pp,dmgs=dmgs,debug=debug,aderivs=aderivs)

	if(model=="gev_p1k3")				pred0=qgev_p1k3_cp(xx,p=pp,kshape=params[4],
																dmgs=dmgs,debug=debug)
#	if(model=="norm_p12")				pred0=qnorm_p12_cp(xx,tt1,tt2,n10=n0,n20=n0,p=pp,
#																dmgs=dmgs,debug=debug,aderivs=aderivs)
#	if(model=="lst_p12k3")			pred0=qlst_p12k3_cp(xx,tt1,tt2,n10=n0,n20=n0,p=pp,
#																dmgs=dmgs,kdf=params[5],debug=debug,aderivs=aderivs)

	if(model=="gamma")					pred0=qgamma_cp(xx,p=pp,dmgs=dmgs,debug=debug,aderivs=aderivs)
	if(model=="invgamma")				pred0=qinvgamma_cp(xx,p=pp,dmgs=dmgs,debug=debug,aderivs=aderivs)
	if(model=="invgauss")				pred0=qinvgauss_cp(xx,p=pp,dmgs=dmgs,debug=debug,aderivs=aderivs)
	if(model=="gev")						pred0=qgev_cp(xx,pp,dmgs=dmgs,debug=debug,pwm=pwm,minxi=minxi,maxxi=maxxi)
	if(model=="gpd_k1")					pred0=qgpd_k1_cp(xx,pp,kloc=params[1],dmgs=dmgs,debug=debug,minxi=minxi,maxxi=maxxi)

	if(model=="gev_p1")					pred0=qgev_p1_cp(x=xx,t=tt1,n0=n10,p=pp,dmgs=dmgs,debug=debug,minxi=minxi,maxxi=maxxi)

	nx=length(tt)
	if(model=="gev_p12")				pred0=qgev_p12_cp(x=xx,t1=tt1,t2=tt2,n01=n0,n02=n0,p=pp,dmgs=dmgs,debug=debug,minxi=minxi,maxxi=maxxi)

	if(model=="gev_p123")				pred0=qgev_p123_cp(x=xx,t1=tt1,t2=tt2,t3=tt3,n01=n0,n02=n0,n03=n0,p=pp,dmgs=dmgs,debug=debug,minxi=minxi,maxxi=maxxi)

	pred=matrix(0,nmethods,nalpha)

	pred[1,]=pred0$ml_quantiles

# two special cases
# -uv instead of ml for norm
# -pwm instead of ml for gev
	if(unbiasedv)	pred[1,]=pred0$uv_quantiles
	if(pwm)				pred[1,]=pred0$pw_quantiles

	if(dmgs)pred[2,]=pred0$cp_quantiles

# some more revert2ml stuff that I don't need any more
#### do the revert2ml thing if the flag says so
###	if(dmgs&&revert2ml){
###		pred[2,]=pred0$ml_quantiles
###	} else {
###		pred[2,]=pred0$cp_quantiles
###	}

# for evt, add the max
	ml_max="Only available for some EVT models"
	if((model=="gev")||(model=="gpd_k1")||(model=="gev_p1")||(model=="gev_p12")){
		ml_max=pred0$ml_max
	}

	return(list(pred=pred,ml_max=ml_max))

}
#' Calculate EP from one model
#' @param model			which distribution to test. Possibles values are
#'	"\code{exp}",
#'	"\code{pareto_k2}",
#'	"\code{halfnorm}",
#'	"\code{unif}",
#'	"\code{norm}",
#'	"\code{norm_dmgs}",
#'	"\code{gnorm_k3}",
#'	"\code{lnorm}",
#'	"\code{lnorm_dmgs}",
#'	"\code{logis}",
#'	"\code{lst_k3}",
#'	"\code{cauchy}",
#'	"\code{gumbel}",
#'	"\code{frechet_k1}",
#'	"\code{weibull}",
#'	"\code{gev_k3}",
#'	"\code{exp_p1}",
#'	"\code{pareto_p1k2}",
#'	"\code{norm_p1}",
#'	"\code{lnorm_p1}",
#'	"\code{logis_p1}",
#'	"\code{lst_p1k3}",
#'	"\code{cauchy_p1}",
#'	"\code{gumbel_p1}",
#'	"\code{frechet_p2k1}",
#'	"\code{weibull_p2}",
#'	"\code{gev_p1k3}",
#'	"\code{norm_p12}",
#'	"\code{lst_p12k3}",
#'	"\code{gamma}",
#'	"\code{invgamma}",
#'	"\code{invgauss}",
#'	"\code{gev}",
#'	"\code{gpd_k1}",
#'	"\code{gev_p1}".
#'	"\code{gev_p12}".
#'	"\code{gev_p123}".
#' @param pred1			quantile predictions
#' @param	tt0				value of the predictor
#' @param	tt10			value of predictor 1
#' @param	tt20			value of predictor 2
#' @param	tt30			value of predictor 3
#' @param	params		the model parameters
#'
#' @return Vector
reltest_makeep=function(model,pred1,tt0,tt10,tt20,tt30,params){

	if(model=="exp")						ep11=1-pexp(pred1,rate=params[1])
	if(model=="pareto_k2")			ep11=1-ppareto(pred1,a=params[2],b=params[1])
	if(model=="halfnorm")				ep11=1-phalfnorm(pred1,theta=params[1])

	if(model=="unif")						ep11=1-punif(pred1,min=params[1],max=params[2])
	if(model=="norm")						ep11=1-pnorm(pred1,mean=params[1],sd=params[2])
	if(model=="norm_dmgs")			ep11=1-pnorm(pred1,mean=params[1],sd=params[2])
	if(model=="gnorm_k3")				ep11=1-pgnorm(pred1,mu=params[1],alpha=params[2],beta=4)
	if(model=="lnorm")					ep11=1-plnorm(pred1,meanlog=params[1],sdlog=params[2])
	if(model=="lnorm_dmgs")			ep11=1-plnorm(pred1,meanlog=params[1],sdlog=params[2])

	if(model=="logis")					ep11=1-plogis(pred1,location=params[1],scale=params[2])
	if(model=="lst_k3")					ep11=1-plst(pred1,mu=params[1],sigma=params[2],df=params[3])
	if(model=="cauchy")					ep11=1-pcauchy(pred1,location=params[1],scale=params[2])

	if(model=="gumbel")					ep11=1-pgumbel(pred1,mu=params[1],sigma=params[2])
	if(model=="frechet_k1")			ep11=1-pfrechet(pred1,lambda=params[3],sigma=params[2],mu=params[1])
	if(model=="weibull")				ep11=1-pweibull(pred1,shape=params[1],scale=params[2])
	if(model=="gev_k3")					ep11=1-extraDistr::pgev(pred1,mu=params[1],sigma=params[2],xi=params[3])

	if(model=="exp_p1")					ep11=1-pexp(pred1,rate=1/exp(params[1]+params[2]*tt0))
	if(model=="pareto_p1k2")		ep11=1-ppareto(pred1,a=1/exp(params[2]+params[3]*tt0),b=params[1])

	if(model=="norm_p1")				ep11=1-pnorm(pred1,mean=params[1]+params[2]*tt0,sd=params[3])
	if(model=="lnorm_p1")				ep11=1-plnorm(pred1,meanlog=params[1]+params[2]*tt0,sdlog=params[3])

	if(model=="logis_p1")				ep11=1-plogis(pred1,location=params[1]+params[2]*tt0,scale=params[3])
	if(model=="lst_p1k3")				ep11=1-plst(pred1,mu=params[1]+params[2]*tt0,sigma=params[3],df=params[4])
	if(model=="cauchy_p1")			ep11=1-pcauchy(pred1,location=params[1]+params[2]*tt0,scale=params[3])

	if(model=="gumbel_p1")			ep11=1-pgumbel(pred1,mu=params[1]+params[2]*tt0,sigma=params[3])
	if(model=="frechet_p2k1")		ep11=1-pfrechet(pred1,mu=params[1],sigma=exp(params[2]+params[3]*tt0),lambda=params[4])
	if(model=="weibull_p2")			ep11=1-pweibull(pred1,shape=params[1],scale=exp(params[2]+params[3]*tt0))
	if(model=="gev_p1k3")				ep11=1-extraDistr::pgev(pred1,mu=params[1]+params[2]*tt0,sigma=params[3],xi=params[4])

#	if(model=="norm_p12")				ep11=1-pnorm(pred1,mean=params[1]+params[2]*tt10,sd=exp(params[3]+params[4]*tt20))
#	if(model=="lst_p12k3")			ep11=1-plst(pred1,mu=params[1]+params[2]*tt10,sigma=exp(params[3]+params[4]*tt20),df=params[5])

	if(model=="gamma")					ep11=1-pgamma(pred1,shape=params[1],scale=params[2])
	if(model=="invgamma")				ep11=1-pinvgamma(pred1,shape=params[1],scale=params[2])
	if(model=="invgauss")				ep11=1-pinvgauss(pred1,mean=params[1],shape=params[2])
	if(model=="gev")						ep11=1-extraDistr::pgev(pred1,
		mu=params[1],sigma=params[2],xi=params[3])
	if(model=="gpd_k1")					ep11=1-extraDistr::pgpd(pred1,mu=params[1],sigma=params[2],xi=params[3])

	if(model=="gev_p1")					ep11=1-extraDistr::pgev(pred1,mu=params[1]+params[2]*tt0,sigma=params[3],xi=params[4])
	if(model=="gev_p12")				ep11=1-extraDistr::pgev(pred1,mu=params[1]+params[2]*tt10,sigma=exp(params[3]+params[4]*tt20),xi=params[5])
	if(model=="gev_p123")				ep11=1-extraDistr::pgev(pred1,mu=params[1]+params[2]*tt10,sigma=exp(params[3]+params[4]*tt20),xi=params[5]+params[6]*tt30)

	return(ep11)
}
#' Calculate MaxEP from one model
#' @param model			which distribution to test. Possibles values are
#'	"\code{gev}",
#'	"\code{gpd_k1}",
#'	"\code{gev_p1}".
#'	"\code{gev_p12}".
#'	"\code{gev_p123}".
#' @param ml_max		predicted max value
#' @param	tt0				value of the predictor
#' @param	tt10			value of predictor 1
#' @param	tt20			value of predictor 2
#' @param	tt30			value of predictor 3
#' @param	params		the model parameters
#'
#' @return Vector
reltest_makemaxep=function(model,ml_max,tt0,tt10,tt20,tt30,params){

	if(model=="gev")						maxep=1-extraDistr::pgev(ml_max,mu=params[1],sigma=params[2],xi=params[3])
	if(model=="gpd_k1")					maxep=1-extraDistr::pgpd(ml_max,mu=params[1],sigma=params[2],xi=params[3])

	if(model=="gev_p1")					maxep=1-extraDistr::pgev(ml_max,mu=params[1]+params[2]*tt0,sigma=params[3],xi=params[4])
	if(model=="gev_p12")				maxep=1-extraDistr::pgev(ml_max,mu=params[1]+params[2]*tt10,sigma=exp(params[3]+params[4]*tt20),xi=params[5])
	if(model=="gev_p123")				maxep=1-extraDistr::pgev(ml_max,mu=params[1]+params[2]*tt10,sigma=exp(params[3]+params[4]*tt20),xi=params[5]+params[6]*tt30)

	return(maxep)
}

