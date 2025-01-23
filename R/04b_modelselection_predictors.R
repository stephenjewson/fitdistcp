#' Model Selection Among 6 Distributions with predictors from the \code{fitdistcp} Package
#'
#' @description
#' Applies model selection using AIC, WAIC1, WAIC2 and leave-one-out logscore
#' to the input data \eqn{x,t},
#' for 11 models with predictors in the \code{fitdistcp} packages
#' (although for the GEV and GPD, the logscore is NA).
#'
#' The code is straightforward, and the point is to illustrate what is
#' possible using the outputs from the \code{fitdistcp} routines.
#'
#' The input data is automatically shifted so that the minimum value is positive,
#' to avoid errors from the models that require positive random variable.
#'
#' For the Pareto, the data is further shifted so that the minimum value is slightly greater than 1.
#'
#' GEVD is temperamental in that
#' it doesn't work if the shape parameter is extreme.
#'
#' @param x 	data vector
#' @param t 	predictor vector
#'
#' @details
#' The 11 models are:
#' \code{exp_p1}
#' \code{pareto_p1k2}
#' \code{norm_p1},
#' \code{lnorm_p1},
#' \code{gumbel_p1},
#' \code{frechet_p2k1},
#' \code{weibull_p2},
#' \code{logis_p1},
#' \code{lst_k3_p1},
#' \code{cauchy_p1} and
#' \code{gev_p1}.

#' @returns
#' Plots QQ plots to the screen, for each of the 11 models,
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
#' @example man/examples/example_03_modelselection_predictors.R
#'
#' @export
#'
modelselection_predictors=function(x,t){
#
xx=x
tt=t
#
models=c(	"exp_p1","pareto_p1k2","norm_p1","lnorm_p1",
					"gumbel_p1","frechet_p2k1","weibull_p2",
					"logis_p1","lst_p1k3","cauchy_p1","gev_p1")

cat("\n")
#
# for the models that require x>0, shift the data if there are negative values
# for the pareto, shift the data so that x>1
#
dd1=0
cat("Note that I adjust the input data, to make sure all values are positive.\n")
dd1=0
dd2=0
if(min(xx)<0)	dd1=0-min(xx)+0.0001
if(min(xx)<1)	dd2=1-min(xx)+0.0001
cat("min before adjustment: min(xx)    =",min(xx),"\n")
cat("min after  adjustment: min(xx+dd1)=",min(xx+dd1),"\n")
#
# fit the models
#
nx=length(xx)
nmodels=11
pp1=(c(1:nx)-0.5)/nx
qq1=	qexp_p1_cp					(xx+dd1,	tt,n0=nx,p=pp1,waicscores=TRUE,logscores=TRUE)
qq2=	qpareto_p1k2_cp			(xx+dd2,	tt,n0=nx,p=pp1,waicscores=TRUE,logscores=TRUE)
qq3=	qnorm_p1_cp					(xx,			tt,n0=nx,p=pp1,waicscores=TRUE,logscores=TRUE)
qq4=	qlnorm_p1_cp				(xx+dd1,	tt,n0=nx,p=pp1,waicscores=TRUE,logscores=TRUE)
qq5=	qgumbel_p1_cp				(xx,			tt,n0=nx,p=pp1,waicscores=TRUE,logscores=TRUE)
qq6=	qfrechet_p2k1_cp		(xx+dd1,	tt,n0=nx,p=pp1,waicscores=TRUE,logscores=TRUE)
qq7=	qweibull_p2_cp			(xx+dd1,	tt,n0=nx,p=pp1,waicscores=TRUE,logscores=TRUE)
qq8=	qlogis_p1_cp				(xx,			tt,n0=nx,p=pp1,waicscores=TRUE,logscores=TRUE)
qq9=	qlst_p1k3_cp				(xx,			tt,n0=nx,p=pp1,kdf=3,waicscores=TRUE,logscores=TRUE)
qq10=	qcauchy_p1_cp				(xx,			tt,n0=nx,p=pp1,waicscores=TRUE,logscores=TRUE)
qq11=	qgev_p1_cp					(xx,			tt,n0=nx,p=pp1,waicscores=TRUE)
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
par(mfrow=c(3,4))
sxx=sort(xx)
plot1=function(x,y1,y2,n){
	ymin=min(y1,y2)
	ymax=max(y1,y2)
#	cat("model=",models[n],"\n")
#	cat("x=",x,"\n")
#	cat("y1=",y1,"\n")
#	cat("y2=",y2,"\n")
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
plot1(sort(qq7$adjustedx),	qq7$cp_quantiles	,qq7$ml_quantiles			,7)
plot1(sort(qq8$adjustedx),	qq8$cp_quantiles	,qq8$ml_quantiles			,8)
plot1(sort(qq9$adjustedx),	qq9$cp_quantiles	,qq9$ml_quantiles	,9)
plot1(sort(qq10$adjustedx),	qq10$cp_quantiles	,qq10$ml_quantiles	,10)
plot1(sort(qq11$adjustedx),	qq11$cp_quantiles	,qq11$ml_quantiles	,11)
#
# print parameter values
#
mle=matrix(0,nmodels,4)
mle[1,1:2]	=qq1$ml_params
mle[2,1:2]	=qq2$ml_params
mle[3,1:3]	=qq3$ml_params
mle[4,1:3]	=qq4$ml_params
mle[5,1:3]	=qq5$ml_params
mle[6,1:3]	=qq6$ml_params
mle[7,1:3]	=qq7$ml_params
mle[8,1:3]	=qq8$ml_params
mle[9,1:3]	=qq9$ml_params
mle[10,1:3]	=qq10$ml_params
mle[11,1:4]	=qq11$ml_params
mledf=data.frame(models,mle)
print(mledf)
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
maic[ 7]	=qq7$maic[3]
maic[ 8]	=qq8$maic[3]
maic[ 9]	=qq9$maic[3]
maic[ 10]	=qq10$maic[3]
maic[ 11]	=qq11$maic[3]
waic1=matrix(0,nmodels)
waic1[ 1]	=qq1$waic1[3]
waic1[ 2]	=qq2$waic1[3]
waic1[ 3]	=qq3$waic1[3]
waic1[ 4]	=qq4$waic1[3]
waic1[ 5]	=qq5$waic1[3]
waic1[ 6]	=qq6$waic1[3]
#waic1[ 7]	=qq7$waic1[3]
waic1[ 7]	=qq7$waic1[3]
waic1[ 8]	=qq8$waic1[3]
waic1[ 9]=qq9$waic1[3]
waic1[ 10]=qq10$waic1[3]
waic1[ 11]=qq11$waic1[3]
waic2=matrix(0,nmodels)
waic2[ 1]	=qq1$waic2[3]
waic2[ 2]	=qq2$waic2[3]
waic2[ 3]	=qq3$waic2[3]
waic2[ 4]	=qq4$waic2[3]
waic2[ 5]	=qq5$waic2[3]
waic2[ 6]	=qq6$waic2[3]
#waic2[ 7]	=qq7$waic2[3]
waic2[ 7]	=qq7$waic2[3]
waic2[ 8]	=qq8$waic2[3]
waic2[ 9]=qq9$waic2[3]
waic2[ 10]=qq10$waic2[3]
waic2[ 11]=qq11$waic2[3]
lsc=matrix(0,nmodels)
lsc[ 1]	=qq1$cp_oos_logscore
lsc[ 2]	=qq2$cp_oos_logscore
lsc[ 3]	=qq3$cp_oos_logscore
lsc[ 4]	=qq4$cp_oos_logscore
lsc[ 5]	=qq5$cp_oos_logscore
lsc[ 6]	=qq6$cp_oos_logscore
lsc[ 7]	=qq7$cp_oos_logscore
lsc[ 8]	=qq8$cp_oos_logscore
lsc[ 9]=qq9$cp_oos_logscore
lsc[ 10]=qq10$cp_oos_logscore
lsc[ 11]=NA
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
cat("\n")
maic=maic
df=data.frame(models,maic,maic_weights,lsc,logscore_weights,waic1,waic1_weights,waic2,waic2_weights)
colnames(df)=c(	"models",
								"maic","weights",
								"logscore","weights",
								"waic1","weights",
								"waic2","weights")
cat("\n")



return(df)

}
