#' Cases
#'
#' @param model			which distribution to test. Possibles values are
#'	"gev",
#'	"gpd_k1",
#'	"gev_pred1".
#' @param nx				length of training data
#' @param params		model parameters
#'
#' @return Two integers
reltest2_cases=function(model="gev",nx=50,params){

	nmethods=5
	case=1
	if(	model=="ntt_ppm"||
			model=="gamma_k1_ppm"||
			model=="gamma_ppm"||
			model=="gev_k12_ppm"){
		nmethods=4
		case=2
	}
	if( model=="gev_ppm"||
			model=="gpd_k1_ppm"||
			model=="gev_p1_ppm"){
		nmethods=6
		case=3
	}
	if(	model=="gev_mpd_ppm"){
		nmethods=5
		case=4
	}
  return(list(nmethods=nmethods,case=case))
}
##' Random training data from one model
#'
#' @param model			which distribution to test. Possibles values are
#'	"gev",
#'	"gpd_k1",
#'	"gev_pred1".
#' @param nx				the length of the training data.
#' @param tt				the predictor
#' @param params		values for the parameters for the specified distribution
#'
#' @return Vector
reltest2_simulate=function(model="gev",nx=50,tt,params){

# cp models have 5 outputs
# but some of the models lack rh_ml and cp
# so I need to notice that and set those results to zero
# which i do using flags, which are by default TRUE
	rh_ml_flag=TRUE
	cp_flag=TRUE

	if(model=="ntt_ppm"){
		xx=rnorm(nx,mean=params[1],sd=sqrt(params[1]))
		rh_ml_flag=FALSE
		cp_flag=FALSE
	}

	if(model=="gamma_k1_ppm"){
		xx=rgamma(nx,shape=params[2],scale=params[1])
		rh_ml_flag=FALSE
		cp_flag=FALSE
	}

	if(model=="gamma_ppm"){
		xx=rgamma(nx,shape=params[2],scale=params[1])
		rh_ml_flag=FALSE
		cp_flag=FALSE
	}

	if(model=="gev"||model=="gev_ppm"||model=="gev_mpd_ppm"){
		xihat=-10
		while(xihat<(-0.45)||(xihat>0.45)){ #0.46 also works...0.47 doesn't
			xx=rgev(nx,mu=params[1],sigma=params[2],xi=params[3])
			ics=matrix(0,3)
			ics=gev_setics(xx,ics)
			opt1=optim(ics,gev_loglik,x=xx,control=list(fnscale=-1))
			xihat=opt1$par[3]
		}
	}

	if(model=="gpd_k1"||model=="gpd_k1_ppm"){
		xihat=-10
		while(xihat<(-0.45)||(xihat>0.45)){ #negative values don't give sensible results
			xx=rgpd(nx,mu=params[1],sigma=params[2],xi=params[3])
			ics=matrix(0,2)
			ics=gpd_k1_setics(xx,ics)
			opt1=optim(ics,gpd_k1_loglik,x=xx,kloc=params[1],control=list(fnscale=-1))
			xihat=opt1$par[2]
		}
	}

	if(model=="gev_k12_ppm"){
		xihat=-10
		while(xihat<(-0.45)||(xihat>0.45)){ #0.46 also works...0.47 doesn't
			xx=rgev(nx,mu=params[1],sigma=params[2],xi=params[3])
			xihat=opt1$minimum
# optimize gives lots of warnings about replacing Inf with max value
# I think it's because the gev has a limited range
# and so sometimes, during the optimisation, the log density gets -Inf value
# when a value of x lies outside the range for the parameters being tested
# all the xihat estimates it returns are very reasonable and well within range
			}
		}
#
# slightly differently for p1 models
#
		if(model=="gev_p1"||model=="gev_p1_ppm"){
			xihat=-10
			while(xihat<(-0.45)||(xihat>0.45)){ #0.46 also works...0.47 doesn't
				xx=rgev(nx,mu=params[1]+params[2]*tt,		sigma=params[3],xi=params[4])
				ics=matrix(0,4)
				ics=gev_p1_setics(xx,tt,ics)
				opt1=optim(ics,gev_p1_loglik,x=xx,t=tt,control=list(fnscale=-1))
				xihat=opt1$par[4]
			}
			ml_params=matrix(0,4)
			d4=0.01
			for (i in 1:4){ml_params[i]=opt1$par[i]}
			gev_p1_checkmle(ml_params,minxi=-0.45,maxxi=0.45)
		}

	return(list(xx=xx,rh_ml_flag=rh_ml_flag,cp_flag=cp_flag))
}
#' Make prediction from one model
#' @param model			which distribution to test. Possibles values are
#'	"\code{exp}",
#'	"\code{pareto_k1}",
#'	"\code{halfnorm}",
#'	"\code{norm}",
#'	"\code{lnorm}",
#'	"\code{gumbel}",
#'	"\code{frechet_k1}",
#'	"\code{weibull}",
#'	"\code{gev_k3}",
#'	"\code{logis}",
#'	"\code{lst_k3}",
#'	"\code{cauchy}",
#'	"\code{norm_p1}",
#'	"\code{lnorm_p1}",
#'	"\code{logis_p1}",
#'	"\code{lst_k3p1}",
#'	"\code{gumbel_p1}",
#'	"\code{norm_p12}",
#'	"\code{gev}",
#'	"\code{gpd}",
#'	"\code{gev_p1}".
#' @param xx				training data
#' @param tt				predictor vector
#' @param n0				index for predictor vector
#' @param pp				probabilities to predict
#' @param params		model parameters
#' @param case			the case number: different models have different lists of methods
#' @param nmethods	the number of methods: different models have different numbers of methods
#' @return Vector
reltest2_predict=function(model="gev",xx,tt,n0,pp,params,case,nmethods){

###			if(model=="ntt_ppm")				pred0=qntt_ppm(xx,pp)
#			if(model=="gamma_k1_ppm")		pred0=qgamma_k1_ppm(xx,pp,kscale=params[1])
###			if(model=="gamma_k1_ppm")		pred0=qgamma_k1_ppm(xx,pp)
###			if(model=="gamma_ppm")			pred0=qgamma_ppm(xx,pp)

			if(model=="gev")						pred0=qgev_cp(xx,pp,extramodels=TRUE)
###			if(model=="gev_ppm")				pred0=qgev_ppm(xx,pp)

#			if(model=="gpd_k1")					pred0=qgpd_k1_cp(xx,pp,extramodels=TRUE,kloc=params[1])
			if(model=="gpd_k1")					pred0=qgpd_k1_cp(xx,pp,extramodels=TRUE)
#			if(model=="gpd_k1_ppm")			pred0=qgpd_k1_ppm(xx,pp,kloc=params[1])
###			if(model=="gpd_k1_ppm")			pred0=qgpd_k1_ppm(xx,pp)

#			if(model=="gev_k12_ppm")		pred0=qgev_k12_ppm(xx,pp,kloc=params[1],kscale=params[2])
###			if(model=="gev_k12_ppm")		pred0=qgev_k12_ppm(xx,pp)
###			if(model=="gev_mpd_ppm")		pred0=qgev_mpd_ppm(xx,pp)

			if(model=="gev_p1")					pred0=qgev_p1_cp(x=xx,t=tt,n0=n0,p=pp,extramodels=TRUE)
###			if(model=="gev_p1_ppm")			pred0=qgev_p1_ppm(x=xx,t=tt,n0=n0,p=pp)


# note that for some of the ppm models, pred[2,] and pred[5,] end up zero
# because those models don't support those kinds of predictions
#			message("flags=",rh_ml_flag,cp_flag)
			nalpha=length(pp)
			pred=matrix(0,nmethods,nalpha)
			if(case==1){
				pred[1,]=0 #pred0$flat_quantiles
				pred[2,]=0 #pred0$jp_quantiles
				pred[3,]=pred0$ml_quantiles
				pred[4,]=0 #pred0$rh_ml_quantiles
				pred[5,]=pred0$cp_quantiles #cp and ml last, so when it's plotted, it's on top
			}
			if(case==2){
				pred[1,]=0 #pred0$flat_quantiles
				pred[2,]=0 #pred0$jp_quantiles
				pred[3,]=pred0$ml_quantiles
				pred[4,]=pred0$mpd_quantiles
			}
			if(case==3){
				pred[1,]=0 #pred0$flat_quantiles
				pred[2,]=0 #pred0$jp_quantiles
				pred[3,]=pred0$ml_quantiles
				pred[4,]=0 #pred0$rh_ml_quantiles
				pred[5,]=pred0$cp_quantiles #cp and ml last, so when it's plotted, it's on top
				pred[6,]=pred0$mpd_quantiles
			}
			if(case==4){
				pred[1,]=0 #pred0$flat_quantiles
				pred[2,]=0 #pred0$jp_quantiles
				pred[3,]=pred0$ml_quantiles
				pred[4,]=pred0$rh_flat_quantiles #cp and ml last, so when it's plotted, it's on top
				pred[5,]=pred0$mpd_quantiles
			}

			return(pred)
}
#' Cases
#'
#' @param model			which distribution to test. Possibles values are
#'	"gev",
#'	"gpd_k1",
#'	"gev_pred1".
#' @param pred1			quantile predictions
#' @param tt0				value of predictor vector
#' @param	params		model parameters
#'
#' @return Vector
reltest2_makeep=function(model,pred1,tt0,params){

###			if(model=="ntt_ppm")				ep11=1-pnorm(pred1,mean=params[1],sd=sqrt(params[1]))
###			if(model=="gamma_k1_ppm")		ep11=1-pgamma(pred1,shape=params[2],scale=params[1])
###			if(model=="gamma_ppm")			ep11=1-pgamma(pred1,shape=params[2],scale=params[1])

			if(model=="gev")						ep11=1-pgev(pred1,mu=params[1],sigma=params[2],xi=params[3])
###			if(model=="gev_ppm")				ep11=1-pgev(pred1,mu=params[1],sigma=params[2],xi=params[3])

			if(model=="gpd_k1")					ep11=1-pgpd(pred1,mu=params[1],sigma=params[2],xi=params[3])
###			if(model=="gpd_k1_ppm")			ep11=1-pgpd(pred1,mu=params[1],sigma=params[2],xi=params[3])

###			if(model=="gev_k12_ppm")		ep11=1-pgev(pred1,mu=params[1],sigma=params[2],xi=params[3])
###			if(model=="gev_mpd_ppm")		ep11=1-pgev(pred1,mu=params[1],sigma=params[2],xi=params[3])

			if(model=="gev_p1")					ep11=1-pgev(pred1,mu=params[1]+params[2]*tt0,sigma=params[3],xi=params[4])
###			if(model=="gev_p1_ppm")			ep11=1-pgev(pred1,mu=params[1]+params[2]*tt0,sigma=params[3],xi=params[4])

			return(ep11)
}


