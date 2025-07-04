#' Model Selection Among 6 Distributions with predictors from the \code{fitdistcp} Package
#'
#' @description
#' Applies model selection using AIC, WAIC1, WAIC2 and leave-one-out logscore
#' to the input data \eqn{x,t},
#' for 6 two tail models with predictors in the \code{fitdistcp} packages
#' (although for the GEV, the logscore is NA for mathematical reasons).
#'
#' The code is straightforward, and the point is to illustrate what is
#' possible using the model selection outputs from the \code{fitdistcp} routines.
#'
#' GEVD is temperamental in that
#' it doesn't work if the shape parameter is extreme.
#'
#' @param x 	data vector
#' @param t 	predictor vector
#'
#' @details
#' The 11 models are:
#' \code{norm_p1},
#' \code{gumbel_p1},
#' \code{logis_p1},
#' \code{lst_k3_p1},
#' \code{cauchy_p1} and
#' \code{gev_p1}.

#' @return
#' Plots QQ plots to the screen, for each of the 6 models,
#' and returns a data frame containing
#' \itemize{
#'	\item AIC scores, AIC weights
#'	\item WAIC1 scores, WAIC1 weights
#'	\item WAIC2 scores, WAIC2 weights
#'	\item logscores and logscore weights
#' }
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_03_ms_predictors_2tail.R
#'
#' @export
#'
ms_predictors_2tail=function(x,t){
#
xx=x
tt=t
#
models=c(	"norm_p1","gumbel_p1","logis_p1","lst_p1k3","cauchy_p1","gev_p1")
#
# fit the models
#
nx=length(xx)
nmodels=6
pp1=(c(1:nx)-0.5)/nx
qq1=	qnorm_p1_cp					(xx,			tt,n0=nx,p=pp1,waicscores=TRUE,logscores=TRUE)
qq2=	qgumbel_p1_cp				(xx,			tt,n0=nx,p=pp1,waicscores=TRUE,logscores=TRUE)
qq3=	qlogis_p1_cp				(xx,			tt,n0=nx,p=pp1,waicscores=TRUE,logscores=TRUE)
qq4=	qlst_p1k3_cp				(xx,			tt,n0=nx,p=pp1,kdf=3,waicscores=TRUE,logscores=TRUE)
qq5=	qcauchy_p1_cp				(xx,			tt,n0=nx,p=pp1,waicscores=TRUE,logscores=TRUE)
qq6=	qgev_p1_cp					(xx,			tt,n0=nx,p=pp1,waicscores=TRUE)
#
# plot qq plots
# -horiz axis is the data we are trying to model
# --but now adjusted using the fitted trend, to the point 'nx'
# --i.e. "detrended" to the point 'nx'
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
		xlab="residuals",ylab="model",col="black")
	points(x,y1,col="red")
	points(x,x,"l")
}
plot1(sort(qq1$adjustedx),	qq1$cp_quantiles	,qq1$ml_quantiles			,1)
plot1(sort(qq2$adjustedx),	qq2$cp_quantiles	,qq2$ml_quantiles			,2)
plot1(sort(qq3$adjustedx),	qq3$cp_quantiles	,qq3$ml_quantiles			,3)
plot1(sort(qq4$adjustedx),	qq4$cp_quantiles	,qq4$ml_quantiles			,4)
plot1(sort(qq5$adjustedx),	qq5$cp_quantiles	,qq5$ml_quantiles			,5)
plot1(sort(qq6$adjustedx),	qq6$cp_quantiles	,qq6$ml_quantiles			,6)
#
# print parameter values
#
mle=matrix(0,nmodels,4)
mle[1,1:3]	=qq1$ml_params
mle[2,1:3]	=qq2$ml_params
mle[3,1:3]	=qq3$ml_params
mle[4,1:3]	=qq4$ml_params
mle[5,1:3]	=qq5$ml_params
mle[6,1:4]	=qq6$ml_params
mleparams_df=data.frame(models,mle)
colnames(mleparams_df)=c("models","param 1","param 2","param 3","param 4")
#
# extract maic, waic and logscores
#
maic=matrix(0,nmodels)
maic[ 1]	=qq1$maic[3]
maic[ 2]	=qq2$maic[3]
maic[ 3]	=qq3$maic[3]
maic[ 4]	=qq4$maic[3]
maic[ 5]	=qq5$maic[3]
maic[ 6]	=qq6$maic[3]
#
waic1=matrix(0,nmodels)
waic1[ 1]	=qq1$waic1[3]
waic1[ 2]	=qq2$waic1[3]
waic1[ 3]	=qq3$waic1[3]
waic1[ 4]	=qq4$waic1[3]
waic1[ 5]	=qq5$waic1[3]
waic1[ 6]	=qq6$waic1[3]
#
waic2=matrix(0,nmodels)
waic2[ 1]	=qq1$waic2[3]
waic2[ 2]	=qq2$waic2[3]
waic2[ 3]	=qq3$waic2[3]
waic2[ 4]	=qq4$waic2[3]
waic2[ 5]	=qq5$waic2[3]
waic2[ 6]	=qq6$waic2[3]
#
lsc=matrix(0,nmodels)
lsc[ 1]	=qq1$cp_oos_logscore
lsc[ 2]	=qq2$cp_oos_logscore
lsc[ 3]	=qq3$cp_oos_logscore
lsc[ 4]	=qq4$cp_oos_logscore
lsc[ 5]	=qq5$cp_oos_logscore
lsc[ 6]	=NaN
#
# calculate weights
#
lscd=exp(0.5*(max(lsc,na.rm=TRUE)+lsc))
maicd=exp(0.5*(max(maic)+maic))
waic1d=exp(0.5*(max(waic1)+waic1))
waic2d=exp(0.5*(max(waic2)+waic2))

logscore_weights=round(90*lscd/sum(lscd,na.rm=TRUE),digits=1)
maic_weights=round(90*maicd/sum(maicd),digits=1)
waic1_weights=round(90*waic1d/sum(waic1d),digits=1)
waic2_weights=round(90*waic2d/sum(waic2d),digits=1)
#
# print weights table
#
maic=maic
df=data.frame(models,maic,maic_weights,lsc,logscore_weights,waic1,waic1_weights,waic2,waic2_weights)
colnames(df)=c(	"models",
								"maic","weights",
								"logscore","weights",
								"waic1","weights",
								"waic2","weights")

return(df)

}
