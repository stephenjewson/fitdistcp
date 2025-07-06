#' Illustration of Model Selection Among 18 Distributions from the \code{fitdistcp} Package
#'
#' @description
#' Applies model selection using AIC, WAIC1, WAIC2 and leave-one-out logscore
#' to the input data \eqn{x},
#' for 7 two tailed models in the \code{fitdistcp} packages
#'
#' The code is straightforward, and the point is to illustrate what is
#' possible using the model selection outputs from the \code{fitdistcp} routines.
#'
#' @param x 			data vector
#'
#' @details
#' The 7 models are:
#' \code{norm},
#' \code{gnorm_k3},
#' \code{gumbel},
#' \code{logis},
#' \code{lst_k3},
#' \code{cauchy},
#' \code{gev}

#' @return
#' Plots QQ plots to the screen, for each of the models,
#' and returns a data frame containing
#' \itemize{
#'	\item AIC scores (times -0.5), AIC weights
#'	\item  WAIC1 scores, WAIC1 weights
#'	\item WAIC2 scores, WAIC2 weights
#'	\item logscores, logscore weights
#'	\item maximum likelihood and calibrating prior means
#'	\item maximum likelihood and calibrating prior standard deviations
#' }
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@gmail.com}
#'
#' @example man/examples/example_004b_ms_flat_2tail.R
#'
#' @export
#'
ms_flat_2tail=function(x){
#
xx=x
#
models=c(	"norm","gnorm_k3","gumbel","logis","lst_k3","cauchy","gev")
#
# fit the models
#
nx=length(xx)
pp1=(c(1:nx)-0.5)/nx
qq1=qnorm_cp				(xx,			pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq2=qgnorm_k3_cp		(xx,			pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq3=qgumbel_cp			(xx,			pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq4=qlogis_cp				(xx,			pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq5=qlst_k3_cp			(xx,			pp1,waicscores=TRUE,logscores=TRUE,means=TRUE,kdf=5)
qq6=qcauchy_cp			(xx,			pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq7=qgev_cp					(xx,			pp1,waicscores=TRUE,means=TRUE)
#
# simulate for means
#
nran=100000
rr1=rnorm_cp				(nran,xx)
rr2=rgnorm_k3_cp		(nran,xx)
rr3=rgumbel_cp			(nran,xx)
rr4=rlogis_cp				(nran,xx)
rr5=rlst_k3_cp			(nran,xx,kdf=5)
rr6=rcauchy_cp			(nran,xx)
rr7=rgev_cp					(nran,xx)
#
# qq plots
# -horiz axis is the data we are trying to model
# -vertical axis in black are quantiles from fitting using maxlik
# -vertical axis in red are quantiles from fitting with parameter uncertainty
# -vertical axis is based on so-called 'Hazen plotting position'
# https://glossary.ametsoc.org/wiki/Plotting_position
#
oldpar=par(no.readonly = TRUE)
on.exit(par(oldpar))
par(mfrow=c(3,3))
sxx=sort(xx)
plot1=function(x,y1,y2,n){
	ymin=min(y1,y2)
	ymax=max(y1,y2)
	plot(x,y2,main=models[n],
		ylim=c(ymin,ymax),
		xlab="data",ylab="model",
		col="black")
	points(x,y1,col="red")
	points(x,x,"l")
}
plot1(sxx,qq1 $cp_quantiles,qq1$ml_quantiles		,1)
plot1(sxx,qq2 $cp_quantiles,qq2$ml_quantiles		,2)
plot1(sxx,qq3 $cp_quantiles,qq3$ml_quantiles		,3)
plot1(sxx,qq4 $cp_quantiles,qq4$ml_quantiles		,4)
plot1(sxx,qq5 $cp_quantiles,qq5$ml_quantiles		,5)
plot1(sxx,qq6 $cp_quantiles,qq6$ml_quantiles		,6)
plot1(sxx,qq7 $cp_quantiles,qq7$ml_quantiles		,7)
#
# parameter values
#
nmodels=7
mle=matrix(0,nmodels,3)
mle[1,1:2]	=qq1$ml_params
mle[2,1:2]	=qq2$ml_params
mle[3,1:2]	=qq3$ml_params
mle[4,1:2]	=qq4$ml_params
mle[5,1:2]	=qq5$ml_params
mle[6,1:2]	=qq6$ml_params
mle[7,1:3]	=qq7$ml_params
mleparams_df=data.frame(models,mle)
colnames(mleparams_df)=c("models","param 1","param 2")
#
# extract maic, waic and logscores
#
maic=matrix(0,nmodels)
maic[ 1]=qq1$maic[3]
maic[ 2]=qq2$maic[3]
maic[ 3]=qq3$maic[3]
maic[ 4]=qq4$maic[3]
maic[ 5]=qq5$maic[3]
maic[ 6]=qq6$maic[3]
maic[ 7]=qq7$maic[3]

waic1=matrix(0,nmodels)
waic1[ 1]=qq1$waic1[3]
waic1[ 2]=qq2$waic1[3]
waic1[ 3]=qq3$waic1[3]
waic1[ 4]=qq4$waic1[3]
waic1[ 5]=qq5$waic1[3]
waic1[ 6]=qq6$waic1[3]
waic1[ 7]=qq7$waic1[3]

waic2=matrix(0,nmodels)
waic2[ 1]=qq1$waic2[3]
waic2[ 2]=qq2$waic2[3]
waic2[ 3]=qq3$waic2[3]
waic2[ 4]=qq4$waic2[3]
waic2[ 5]=qq5$waic2[3]
waic2[ 6]=qq6$waic2[3]
waic2[ 7]=qq7$waic2[3]

lsc=matrix(0,nmodels)
lsc[ 1]=qq1$cp_oos_logscore
lsc[ 2]=qq2$cp_oos_logscore
lsc[ 3]=qq3$cp_oos_logscore
lsc[ 4]=qq4$cp_oos_logscore
lsc[ 5]=qq5$cp_oos_logscore
lsc[ 6]=qq6$cp_oos_logscore
lsc[ 7]=NaN
#
# calculate weights
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
# print weights table
#
maic=maic
modelselection_df=data.frame(models,maic,maic_weights,lsc,logscore_weights,waic1,waic1_weights,waic2,waic2_weights)
colnames(modelselection_df)=c(	"models",
								"maic","weights",
								"logscore","weights",
								"waic1","weights",
								"waic2","weights")
#
# extract means from the q results
#
ml_means=matrix(0,nmodels)
ml_means[ 1]=qq1$ml_mean
ml_means[ 2]=qq2$ml_mean
ml_means[ 3]=qq3$ml_mean
ml_means[ 4]=qq4$ml_mean
ml_means[ 5]=qq5$ml_mean
ml_means[ 6]=qq6$ml_mean
ml_means[ 7]=qq7$ml_mean
#
cp_means=matrix(0,nmodels)
cp_means[ 1]=qq1$cp_mean
cp_means[ 2]=qq2$cp_mean
cp_means[ 3]=qq3$cp_mean
cp_means[ 4]=qq4$cp_mean
cp_means[ 5]=qq5$cp_mean
cp_means[ 6]=qq6$cp_mean
cp_means[ 7]=qq7$cp_mean
#
# calculate sds from the random deviates
# -and set to inf if we know they are inf, to avoid calculating nonsense
#
ml_sds=matrix(0,nmodels)
ml_sds[ 1]=sd(rr1$ml_deviates)
ml_sds[ 2]=sd(rr2$ml_deviates)
ml_sds[ 3]=sd(rr3$ml_deviates)
ml_sds[ 4]=sd(rr4$ml_deviates)
ml_sds[ 5]=sd(rr5$ml_deviates)
ml_sds[ 6]=Inf
ml_sds[ 7]=sd(rr7$ml_deviates)
#
cp_sds=matrix(0,nmodels)
cp_sds[ 1]=sd(rr1$cp_deviates)
cp_sds[ 2]=sd(rr2$cp_deviates)
cp_sds[ 3]=sd(rr3$cp_deviates)
cp_sds[ 4]=sd(rr4$cp_deviates)
cp_sds[ 5]=sd(rr5$cp_deviates)
cp_sds[ 6]=Inf
cp_sds[ 7]=Inf
#
means_df=data.frame(models,ml_means,cp_means,ml_sds,cp_sds)
colnames(means_df)=c(	"models","ml_means","cp_means","ml_sds","cp_sds")

return(list(mleparams_df=mleparams_df,
	modelselection_df=modelselection_df,
	means_df=means_df))

}
