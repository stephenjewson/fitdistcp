#' Illustration of Model Selection Among 10 One Tail Distributions from the \code{fitdistcp} Package
#'
#' @description
#' Applies model selection using AIC, WAIC1, WAIC2 and leave-one-out logscore
#' to the input data \eqn{x},
#' for 10 one tailed models in the \code{fitdistcp} package
#' (although for the GPD, the logscore is NA for mathematical reasons).
#'
#' The code is straightforward, and the point is to illustrate what is
#' possible using the model selection outputs from the \code{fitdistcp} routines.
#'
#' The input data may be automatically shifted so that the minimum value is positive.
#'
#' For the Pareto, the data may be further shifted so that the minimum value is slightly greater than 1.
#'
#' @param x 								data vector
#' @param index							which data point to use for plotting positions
#' @param nyears						number of years for frequency calculations
#' @param plottype					What to plot? Possible values are 'both', 'empirical', 'cp'
#' @param plottingposition	Weibull or Hazen
#' @param quiet							logical for whether to print screen messages
#'
#' @details
#' The 10 models are:
#' \code{exp},
#' \code{pareto_k2},
#' \code{halfnorm},
#' \code{lnorm},
#' \code{frechet_k1},
#' \code{weibull},
#' \code{gamma},
#' \code{invgamma},
#' \code{invgauss} and
#' \code{gpd_k1}.

#' @return
#' Plots QQ plots to the screen, for each of the models,
#' and returns a data frame containing
#' \itemize{
#'  \item MLE parameter values
#'	\item AIC scores (times -0.5), AIC weights
#'	\item WAIC1 scores, WAIC1 weights
#'	\item WAIC2 scores, WAIC2 weights
#'	\item logscores, logscore weights
#'	\item maximum likelihood and calibrating prior means
#'	\item maximum likelihood and calibrating prior standard deviations
#' }
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@gmail.com}
#'
#' @example man/examples/example_004a_ms_flat_1tail.R
#'
#' @export
#'
ms_flat_1tail=function(x,index=1,nyears=10,
	plottype="empirical",plottingposition="Weibull",
	quiet=FALSE){
#
xx=x
#
# 1 set up labels for the graph
#
models=c(	"exp","pareto_k2","halfnorm","lnorm","frechet_k1",
						"weibull","gamma","invgamma","invgauss","gpd_k1")
nmodels=length(models)
#
# 2 shift the data if necessary
# -since all models require x>0, shift the data if there are negative or zero values
# -for the Pareto, shift the data so that x>1
# -and print a message to the screen to warn that the data has been shifted
#
dd1=0
if(min(xx)<0.001){
	dd1=0-min(xx)+0.001
	if(!quiet){
		message("Message from inside modelselection_flat_1tail:")
		message(" I'm increasing the data to eliminate negative or zero values")
		message("  adjustment=",dd1)
		message("  before adjustment:min(xx)     =",min(xx))
		message("  after  adjustment:min(xx+dd1) =",min(xx+dd1))
	}
}
dd2=0
if(min(xx)<1){
	dd2=1-min(xx)+0.001
	if(!quiet){
		message("Message from inside modelselection_flat_1tail:")
		message(" For the Pareto only, I increase the input data so that it's all >1")
		message(" I fit the model, and then adjust back")
		message("  adjustment=",dd2)
		message("  before Pareto adjustment:min(xx)     =",min(xx))
		message("  after  Pareto adjustment:min(xx) =",min(xx+dd2))
	}
}
#
# 3 fit the quantiles using q***
#
nx=length(xx)
if(plottingposition=="Weibull"){
	pp1=c(1:nx)/(nx+1)
} else if (plottingposition=="Hazen"){
	pp1=(c(1:nx)-0.5)/nx
} else{
	message("Plotting position not correctly specified, so stopping.\n")
	stop()
}
qq01	=qexp_cp				(xx+dd1,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq02	=qpareto_k2_cp	(xx+dd2,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq03	=qhalfnorm_cp		(xx+dd1,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq04	=qlnorm_cp			(xx+dd1,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq05	=qfrechet_k1_cp	(xx+dd1,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq06	=qweibull_cp		(xx+dd1,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq07	=qgamma_cp			(xx+dd1,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq08	=qinvgamma_cp		(xx+dd1,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq09	=qinvgauss_cp		(xx+dd1,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq10	=qgpd_k1_cp			(xx+dd1,	pp1,waicscores=TRUE,means=TRUE,kloc=0)
#
# 4 simulate randoms to calculate sd using r***
#
nran=100000
rr01 =rexp_cp					(nran,xx+dd1)
rr02 =rpareto_k2_cp		(nran,xx+dd2)
rr03 =rhalfnorm_cp		(nran,xx+dd1)
rr04 =rlnorm_cp				(nran,xx+dd1)
rr05 =rfrechet_k1_cp	(nran,xx+dd1)
rr06 =rweibull_cp			(nran,xx+dd1)
rr07 =rgamma_cp				(nran,xx+dd1)
rr08 =rinvgamma_cp		(nran,xx+dd1)
rr09 =rinvgauss_cp		(nran,xx+dd1)
rr10 =rgpd_k1_cp			(nran,xx+dd1)
#
# 5 make qq plots
# -horiz axis is the data we are trying to model
# -vertical axis in black are quantiles from fitting using maxlik
# -vertical axis in red are quantiles from fitting with parameter uncertainty
# -vertical axis is based on so-called 'Hazen plotting position'
# https://glossary.ametsoc.org/wiki/Plotting_position
#
oldpar=par(no.readonly = TRUE)
on.exit(par(oldpar))
par(mfrow=c(4,3))
sxx=sort(xx)
#
# 6 define the plotting function
#
plot1=function(x,y1,y2,n,plottype){

	if(plottype=="both"){
		ymin=min(y1,y2)
		ymax=max(y1,y2)
		plot(x,y2,main=models[n],
			ylim=c(ymin,ymax),
			xlab="data",ylab="model",
			col="black")
		points(x,y1,col="red")
		points(x,x,"l")
	} else if (plottype=="empirical"){
		ymin=min(y2)
		ymax=max(y2)
		plot(x,y2,main=models[n],
			ylim=c(ymin,ymax),
			xlab="data",ylab="model",
			col="black")
		points(x,x,"l")
	} else if (plottype=="cp"){
		ymin=min(y1)
		ymax=max(y1)
		plot(x,y1,main=models[n],
			ylim=c(ymin,ymax),
			xlab="data",ylab="model",
			col="red")
		points(x,x,"l")
	} else {
		message("plottype not defined, so stopping.")
		stop()
	}
}
#
# 7 plot the qq plots using the fitted quantiles
#
plot1(sxx,qq01	$cp_quantiles-dd1,	qq01	$ml_quantiles-dd1		,1, plottype)
plot1(sxx,qq02	$cp_quantiles-dd2,	qq02	$ml_quantiles-dd2		,2, plottype)
plot1(sxx,qq03	$cp_quantiles-dd1,	qq03	$ml_quantiles-dd1		,3, plottype)
plot1(sxx,qq04	$cp_quantiles-dd1,	qq04	$ml_quantiles-dd1		,4, plottype)
plot1(sxx,qq05	$cp_quantiles-dd1,	qq05	$ml_quantiles-dd1		,5, plottype)
plot1(sxx,qq05	$cp_quantiles-dd1,	qq06	$ml_quantiles-dd1		,6, plottype)
plot1(sxx,qq06	$cp_quantiles-dd1,	qq07	$ml_quantiles-dd1		,7, plottype)
plot1(sxx,qq08	$cp_quantiles-dd1,	qq08	$ml_quantiles-dd1		,8, plottype)
plot1(sxx,qq09	$cp_quantiles-dd1,	qq09	$ml_quantiles-dd1		,9, plottype)
plot1(sxx,qq10	$cp_quantiles-dd1,	qq10	$ml_quantiles-dd1		,10,plottype)
#
# 8 put mle loglik value and parameters into a data frame
#
mle=matrix(NA,nmodels,3)
mle[1,1:2]	=c(qq01$ml_value,qq01	$ml_params)
mle[2,1:2]	=c(qq02$ml_value,qq02	$ml_params)
mle[3,1:2]	=c(qq03$ml_value,qq03	$ml_params)
mle[4,1:3]	=c(qq04$ml_value,qq04	$ml_params)
mle[5,1:3]	=c(qq05$ml_value,qq05	$ml_params)
mle[6,1:3]	=c(qq06$ml_value,qq06	$ml_params)
mle[7,1:3]	=c(qq07$ml_value,qq07	$ml_params)
mle[8,1:3]	=c(qq08$ml_value,qq08	$ml_params)
mle[9,1:3]	=c(qq09$ml_value,qq09	$ml_params)
mle[10,1:3]	=c(qq10$ml_value,qq10	$ml_params)
mleparams_df=data.frame(models,mle)
colnames(mleparams_df)=c("models","maxlik value","param 1","param 2")
#
# 10 extract maic, waic and logscores from the qresults
#
maic=matrix(0,nmodels)
maic[ 1]=qq01	$maic[3]
maic[ 2]=qq02	$maic[3]
maic[ 3]=qq03	$maic[3]
maic[ 4]=qq04	$maic[3]
maic[ 5]=qq05	$maic[3]
maic[ 6]=qq06	$maic[3]
maic[ 7]=qq07	$maic[3]
maic[ 8]=qq08	$maic[3]
maic[ 9]=qq09	$maic[3]
maic[10]=qq10	$maic[3]
#
waic1=matrix(0,nmodels)
waic1[ 1]=qq01	$waic1[3]
waic1[ 2]=qq02	$waic1[3]
waic1[ 3]=qq03	$waic1[3]
waic1[ 4]=qq04	$waic1[3]
waic1[ 5]=qq05	$waic1[3]
waic1[ 6]=qq06	$waic1[3]
waic1[ 7]=qq07	$waic1[3]
waic1[ 8]=qq08	$waic1[3]
waic1[ 9]=qq09	$waic1[3]
waic1[10]=qq10$waic1[3]
#
waic2=matrix(0,nmodels)
waic2[ 1]=qq01	$waic2[3]
waic2[ 2]=qq02	$waic2[3]
waic2[ 3]=qq03	$waic2[3]
waic2[ 4]=qq04	$waic2[3]
waic2[ 5]=qq05	$waic2[3]
waic2[ 6]=qq06	$waic2[3]
waic2[ 7]=qq07	$waic2[3]
waic2[ 8]=qq08	$waic2[3]
waic2[ 9]=qq09	$waic2[3]
waic2[10]=qq10$waic2[3]
#
lsc=matrix(0,nmodels)
lsc[ 1]=qq01$cp_oos_logscore
lsc[ 2]=qq02$cp_oos_logscore
lsc[ 3]=qq03$cp_oos_logscore
lsc[ 4]=qq04$cp_oos_logscore
lsc[ 5]=qq05$cp_oos_logscore
lsc[ 6]=qq06$cp_oos_logscore
lsc[ 7]=qq07$cp_oos_logscore
lsc[ 8]=qq08$cp_oos_logscore
lsc[ 9]=qq09$cp_oos_logscore
lsc[10]=NA
#
# 11 calculate 3 sets of model weights (maic, waic1, waic2)
#
lscd=	exp(0.5*(max(lsc,na.rm=TRUE)+lsc))
maicd=exp(0.5*(max(maic)+maic))
waic1d=exp(0.5*(max(waic1)+waic1))
waic2d=exp(0.5*(max(waic2)+waic2))
logscore_weights=	round(100*lscd	/sum(lscd,na.rm=TRUE)	,digits=1)
maic_weights=			round(100*maicd	/sum(maicd),digits=1)
waic1_weights=		round(100*waic1d	/sum(waic1d),digits=1)
waic2_weights=		round(100*waic2d	/sum(waic2d),digits=1)
#
# 12 print weights table to screen
#
maic=maic
modelselection_df=data.frame(models,maic,maic_weights,lsc,logscore_weights,waic1,waic1_weights,waic2,waic2_weights)
colnames(modelselection_df)=c(	"model",
								"maic","weights1",
								"logscore","weights2",
								"waic1","weights3",
								"waic2","weights4")
#
# 13 extract means from the q results
# -except for lnorm cp mean, which has to be calculated from randoms
#
ml_means=matrix(0,nmodels)
ml_means[ 1]=qq01	$ml_mean
ml_means[ 2]=qq02	$ml_mean
ml_means[ 3]=qq03	$ml_mean
ml_means[ 4]=qq04	$ml_mean
ml_means[ 5]=qq05	$ml_mean
ml_means[ 6]=qq06	$ml_mean
ml_means[ 7]=qq07	$ml_mean
ml_means[ 8]=qq08	$ml_mean
ml_means[ 9]=qq09	$ml_mean
ml_means[10]=qq10	$ml_mean
#
cp_means=matrix(0,nmodels)
cp_means[ 1]=qq01	$cp_mean
cp_means[ 2]=qq02	$cp_mean
cp_means[ 3]=qq03	$cp_mean
cp_means[ 4]=mean(rr04	$cp_deviates) #for lnorm, have to calculate from randoms
cp_means[ 5]=qq05	$cp_mean
cp_means[ 6]=qq06	$cp_mean
cp_means[ 7]=qq07	$cp_mean
cp_means[ 8]=qq08  $cp_mean
cp_means[ 9]=qq09	$cp_mean
cp_means[10]=qq10	$cp_mean

# calculate sds from the random deviates
# -and set to inf if we know they are inf

ml_sds=matrix(0,nmodels)
ml_sds[ 1]=sd(rr01	$ml_deviates)
ml_sds[ 2]=ifelse(mle[2,1]>1,sd(rr02$ml_deviates),Inf)
ml_sds[ 3]=sd(rr03	$ml_deviates)
ml_sds[ 4]=sd(rr04	$ml_deviates)
ml_sds[ 5]=sd(rr05	$ml_deviates)
ml_sds[ 6]=sd(rr06	$ml_deviates)
ml_sds[ 7]=sd(rr07	$ml_deviates)
ml_sds[ 8]=sd(rr08	$ml_deviates)
ml_sds[ 9]=sd(rr09	$ml_deviates)
ml_sds[10]=sd(rr10$ml_deviates)

cp_sds=matrix(0,nmodels)
cp_sds[ 1]=sd(rr01	$cp_deviates)
cp_sds[ 2]=Inf
cp_sds[ 3]=sd(rr03	$cp_deviates)
cp_sds[ 4]=sd(rr04	$cp_deviates)
cp_sds[ 5]=Inf
cp_sds[ 6]=sd(rr06	$cp_deviates)
cp_sds[ 7]=sd(rr07	$cp_deviates)
cp_sds[ 8]=Inf
cp_sds[ 9]=sd(rr09	$cp_deviates)
cp_sds[10]=Inf
#
# 14 make data frames to return the results for means and sds
#
means_df=data.frame(models,ml_means,cp_means,ml_sds,cp_sds)
colnames(means_df)=c(	"models","ml_means","cp_means","ml_sds","cp_sds")
#
# 15 calculate the returnperiod for one particular value in the series
#
return_period=matrix(0,nmodels,3)
#
# method 1
# empirical, assuming the input data was sorted from low to high
ndatapoints=length(xx)
# first, calculate the frequency of any event at all per year
# which is the 'n' in the binomial np formula
f0=ndatapoints/nyears
# then calculate the conditional exceedance probability empirically
if(plottingposition=="Weibull"){
	prob=(ndatapoints+1-index)/(ndatapoints+1)
} else if (plottingposition=="Hazen"){
	prob=(ndatapoints+0.5-index)/ndatapoints
}
# then calculate the frequency of the index'th event
ff=f0*prob
# and the rp is the inverse o fthat
return_period[,1]=1/ff
#
# method 2
# using maxlik probabilities
returnlevel=xx[index]
return_period[2,2] =1/(f0*(1-pexp_cp				(xx+dd1,y=returnlevel+dd1)$ml_cdf))
return_period[2,2] =1/(f0*(1-ppareto_k2_cp	(xx+dd2,y=returnlevel+dd2)$ml_cdf))
return_period[3,2] =1/(f0*(1-phalfnorm_cp		(xx+dd1,y=returnlevel+dd1)$ml_cdf))
return_period[4,2] =1/(f0*(1-plnorm_cp			(xx+dd1,y=returnlevel+dd1)$ml_cdf))
return_period[5,2] =1/(f0*(1-pfrechet_k1_cp	(xx+dd1,y=returnlevel+dd1)$ml_cdf))
return_period[6,2] =1/(f0*(1-pweibull_cp		(xx+dd1,y=returnlevel+dd1)$ml_cdf))
return_period[7,2] =1/(f0*(1-pgamma_cp			(xx+dd1,y=returnlevel+dd1)$ml_cdf))
return_period[8,2] =1/(f0*(1-pinvgamma_cp		(xx+dd1,y=returnlevel+dd1)$ml_cdf))
return_period[9,2] =1/(f0*(1-pinvgauss_cp		(xx+dd1,y=returnlevel+dd1)$ml_cdf))
return_period[10,2]=1/(f0*(1-pgpd_k1_cp			(xx+dd1,y=returnlevel+dd1)$ml_cdf))
#
# method 3
# using cp probabilities
return_period[1,3] =1/(f0*(1-pexp_cp				(xx+dd1,y=returnlevel+dd1)$cp_cdf))
return_period[2,3] =1/(f0*(1-ppareto_k2_cp	(xx+dd2,y=returnlevel+dd2)$cp_cdf))
return_period[3,3] =1/(f0*(1-phalfnorm_cp		(xx+dd1,y=returnlevel+dd1)$cp_cdf))
return_period[4,3] =1/(f0*(1-plnorm_cp			(xx+dd1,y=returnlevel+dd1)$cp_cdf))
return_period[5,3] =1/(f0*(1-pfrechet_k1_cp	(xx+dd1,y=returnlevel+dd1)$cp_cdf))
return_period[6,3] =1/(f0*(1-pweibull_cp		(xx+dd1,y=returnlevel+dd1)$cp_cdf))
return_period[7,3] =1/(f0*(1-pgamma_cp			(xx+dd1,y=returnlevel+dd1)$cp_cdf))
return_period[8,3] =1/(f0*(1-pinvgamma_cp		(xx+dd1,y=returnlevel+dd1)$cp_cdf))
return_period[9,3] =1/(f0*(1-pinvgauss_cp		(xx+dd1,y=returnlevel+dd1)$cp_cdf))
return_period[10,3]=NA
#
return_period_df=data.frame(return_period)
colnames(return_period_df)=c(plottingposition,"maxlik","cp")
#
# 16 return model selection results, means and sds, return periods
#
return(list(mleparams_df=mleparams_df,
	modelselection_df=modelselection_df,
	means_df=means_df,
	return_period_df=return_period_df))

}
