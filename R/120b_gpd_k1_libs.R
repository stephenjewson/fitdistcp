#' rgpd for gpd_k1 but with maxlik xi within bounds
#' @return Vector
#' @inheritParams manf
rgpd_k1_minmax=function(nx,kloc,sigma,xi,minxi=-0.45,maxxi=0.45){
	xihat=-9999
	while((xihat<minxi)||(xihat>maxxi)){ #0.46 also works...0.47 doesn't
		xx=extraDistr::rgpd(nx,mu=kloc,sigma=sigma,xi=xi)
		ics=gpd_k1_setics(xx,c(0,0))
		opt1=optim(ics,gpd_k1_loglik,x=xx,kloc=kloc,control=list(fnscale=-1))
		xihat=opt1$par[2]
	}
	return(xx)
}
#' Waic
#' @inherit manwaic return
#' @inheritParams manf
gpd_k1_waic=function(waicscores,x,v1hat,v2hat,kloc,lddi,lddd,lambdad){
		if(waicscores){

			f1f=gpd_k1_f1fa(x,v1hat,v2hat,kloc)
			f2f=gpd_k1_f2fa(x,v1hat,v2hat,kloc)

			fhatx=extraDistr::dgpd(x,mu=kloc,sigma=v1hat,xi=v2hat)
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
gpd_k1_logf=function(params,x,kloc){
	sc=pmax(params[1],.Machine$double.eps)
	sh=params[2]
	logf=sum(dgpd(x,mu=kloc,sigma=sc,xi=sh,log=TRUE))-log(sc)
	return(logf)
}
#' Set initial conditions
#' @return Vector
#' @inheritParams manf
gpd_k1_setics=function(x,ics){
	if((ics[1]==0)&&(ics[2]==0)){
		ics[1]=sd(x)
		ics[2]=0
	}
	return(ics)
}
#'  log-likelihood function
#' @inherit manloglik return
#' @inheritParams manf
gpd_k1_loglik=function(vv,x,kloc){
	n=length(x)
	loglik=sum(extraDistr::dgpd(x,mu=kloc,sigma=max(vv[1],.Machine$double.eps),xi=vv[2],log=TRUE))
	return(loglik)
}
#' Check MLE
#' @inherit mancheckmle return
#' @inheritParams manf
gpd_k1_checkmle=function(ml_params,kloc,minxi=-1,maxxi=2){
# currently not used, because instead I revert2ml
	v1hat=ml_params[1]
	v2hat=ml_params[2]
	if(is.na(v1hat))stop()
	if(is.na(v2hat))stop()
#
# min xi
#
##	minxi=0
	if(v2hat<minxi){warning("\n***v2hat=",v2hat,"=> execution halted because maxlik shape parameter <",minxi,"***\n");stop()}
#
# max xi
#
	if(v2hat>maxxi){warning("\n***v2hat=",v2hat,"=> execution halted because maxlik shape parameter >",maxxi,"***\n");stop()}
# This max value is ad-hoc
# If it's lowered to 1, then the ppm results for xi=0.6 go wrong, which I understand.
# If it's increased to 100, then in about 1 in a billion cases, for nx=25,
# the xi-hat value is very large and the code crashes because lddi can't be calculated.
# I suspect there is more to understand about that latter part, but for now
# this is a compromise.
}
#' Analytical Expressions for Predictive Means
#' RHP mean based on the expectation of DMGS equation 2.1
#' @inherit manmeans return
#' @inheritParams manf
gpd_k1_means=function(means,ml_params,lddi,lddd,
									lambdad_rh_flat,lambdad_jp,
									nx,dim=2,kloc=0){

	if(means){
# intro
		eulerconstant=0.57721566490153286060651209008240243104215933593992
		mu=kloc
		sigma=ml_params[1]
		xi=ml_params[2]

# set up derivative arrays
		meand1=array(0,c(2,1))
		meand2=array(0,c(2,2,1))

# mle
		onemxi=1-xi
		onemxisq=onemxi*onemxi
		onemxicu=onemxi*onemxi*onemxi
		ml_mean=kloc+sigma/onemxi
# calculate first derivative array for bayesian xi!=0 cases
		meand1[1,1]=1/onemxi
		meand1[2,1]=sigma/onemxisq
# calculate second derivative array (only has 1 non-zero term!)
		meand2[1,1,1]=0
		meand2[1,2,1]=1/onemxisq
		meand2[2,2,1]=-2*sigma/onemxicu

# I've realized now that when I integrate over xi, the mean is Inf

#	rh_flat_mean	=ml_mean+dmgs(lddi,lddd,meand1,lambdad_rh_flat,meand2,dim=2)/nx
		rh_flat_mean	=Inf

	}else{
		ml_mean="means not selected"
		rh_flat_mean="means not selected"
	}

# return
	list(ml_mean=ml_mean,rh_flat_mean=rh_flat_mean)
}
#' Densities for 5 predictions
#' @inherit mandsub return
#' @inheritParams manf
dgpdsub=function(x,y,ics,kloc=0,dlogpi=0,
	minxi,maxxi,extramodels=FALSE){

		nx=length(x)

		ics=gpd_k1_setics(x,ics)
		opt=optim(ics,gpd_k1_loglik,x=x,kloc=kloc,control=list(fnscale=-1))
		v1hat=opt$par[1]
		v2hat=min(maxxi,max(minxi,opt$par[2])) #just reset in this case
#		v2hat=opt$par[2]
		ml_params=c(v1hat,v2hat)

		y=fixgpdrange(y,kloc,v1hat,v2hat)


# mle
		ml_pdf=extraDistr::dgpd(y,mu=kloc,sigma=v1hat,xi=v2hat)
		ml_cdf=extraDistr::pgpd(y,mu=kloc,sigma=v1hat,xi=v2hat)

# return
		list(
					ml_params=ml_params,
					ml_pdf=ml_pdf,
					ml_cdf=ml_cdf)
}

