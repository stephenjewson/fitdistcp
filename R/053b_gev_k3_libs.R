#' Waic
#' @inherit manwaic return
#' @inheritParams manf
gev_k3_waic=function(waicscores,x,v1hat,v2hat,kshape,lddi,lddd,lambdad){
		if(waicscores){

			f1f=gev_k3_f1fa(x,v1hat,v2hat,kshape)

			f2f=gev_k3_f2fa(x,v1hat,v2hat,kshape)

			fhatx=extraDistr::dgev(x,mu=v1hat,sigma=v2hat,xi=kshape)
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
gev_k3_logf=function(params,x,kshape){
	l=params[1]
	s=pmax(params[2],.Machine$double.eps)
	logf=sum(extraDistr::dgev(x,mu=l,sigma=s,xi=kshape,log=TRUE))-log(s)
	return(logf)
}
#'  log-likelihood function
#' @inherit manloglik return
#' @inheritParams manf
gev_k3_loglik=function(vv,x,kshape){
	n=length(x)
	loglik=sum(extraDistr::dgev(x,mu=vv[1],sigma=max(vv[2],.Machine$double.eps),xi=kshape,log=TRUE))
	return(loglik)
}
#' MLE and RHP means
#' @inherit manmeans return
#' @inheritParams manf
gev_k3_means=function(means,ml_params,lddi,lddd,lambdad_rhp,nx,dim=2,kshape){

	if(means){
# intro
		eulerconstant=0.57721566490153286060651209008240243104215933593992

		if(kshape==0){
# gumbel case
			means=gumbel_means(means,ml_params,lddi,lddd,lambdad_rhp,nx,dim=2)
			ml_mean=means$ml_mean
			rh_mean=means$rh_mean
		} else{
# non-gumbel case

# -ml
			g1=gamma(1-kshape)
			ml_mean=ml_params[1]+ml_params[2]*(g1-1)/kshape
# -rhp
			meand1=array(0,c(2,1))
			meand1[1,1]=1
			meand1[2,1]=(g1-1)/kshape
			meand2=array(0,c(2,2,1)) #and they are all zero it seems
			dmean=dmgs(lddi,lddd,meand1,lambdad_rhp,meand2,dim=2)
			rh_mean=ml_mean+dmean/nx
		}
	}else{
		ml_mean="means not selected"
		rh_mean="means not selected"
	}

	list(ml_mean=ml_mean,rh_mean=rh_mean)
}
#' Densities from MLE and RHP
#' @inherit mandsub return
#' @inheritParams manf
dgev_k3sub=function(x,y,kshape){

		nx=length(x)

		v1start=mean(x)
		v2start=sd(x)
		opt=optim(c(v1start,v2start),gev_k3_loglik,x=x,kshape=kshape,
			control=list(fnscale=-1))
		v1hat=opt$par[1]
		v2hat=opt$par[2]
		ml_params=c(v1hat,v2hat)

		y=fixgevrange(y,v1hat,v2hat,kshape)

# mle
		ml_pdf=extraDistr::dgev(y,mu=v1hat,sigma=v2hat,xi=kshape)
		ml_cdf=extraDistr::pgev(y,mu=v1hat,sigma=v2hat,xi=kshape)

# return
		list(	ml_params=ml_params,
					ml_pdf=ml_pdf,
					ml_cdf=ml_cdf)
}



