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
#' @param x 			data vector
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
#' @example man/examples/example_03_ms_flat_1tail.R
#'
#' @export
#'
ms_flat_1tail=function(x){
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
if(min(xx)<0){
	dd1=0-min(xx)+0.0001
	message("Message from inside modelselection_flat_1tail:")
	message(" I'm increasing the data to eliminate negative or zero values")
	message("  adjustment=",dd1)
	message("  before adjustment:min(xx)     =",min(xx))
	message("  after  adjustment:min(xx+dd1) =",min(xx+dd1))
}
dd2=0
if(min(xx)<1){
	dd2=1-min(xx)+0.0001
	message("Message from inside modelselection_flat_1tail:")
	message(" For the Pareto only, I increase the input data so that it's all >1")
	message(" I fit the model, and then adjust back")
	message("  adjustment=",dd2)
	message("  before Pareto adjustment:min(xx)     =",min(xx))
	message("  after  Pareto adjustment:min(xx) =",min(xx+dd2))
}
#
# 3 fit the quantiles using q***
#
nx=length(xx)
pp1=(c(1:nx)-0.5)/nx
qq1	=qexp_cp				(xx+dd1,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq2	=qpareto_k2_cp	(xx+dd2,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq3	=qhalfnorm_cp		(xx+dd1,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq4	=qlnorm_cp			(xx+dd1,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq5	=qfrechet_k1_cp	(xx+dd1,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq6	=qweibull_cp		(xx+dd1,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq7	=qgamma_cp			(xx,			pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq8	=qinvgamma_cp		(xx,			pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq9	=qinvgauss_cp		(xx,			pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq10=qgpd_k1_cp			(xx,			pp1,waicscores=TRUE,means=TRUE,kloc=0)
#
# 4 simulate randoms to calculate sd using r***
#
nran=100000
rr1 =rexp_cp				(nran,xx+dd1)
rr2 =rpareto_k2_cp	(nran,xx+dd2)
rr3 =rhalfnorm_cp		(nran,xx+dd1)
rr4 =rlnorm_cp			(nran,xx+dd1)
rr5 =rfrechet_k1_cp	(nran,xx+dd1)
rr6 =rweibull_cp		(nran,xx+dd1)
rr7 =rgamma_cp			(nran,xx)
rr8 =rinvgamma_cp		(nran,xx)
rr9 =rinvgauss_cp		(nran,xx)
rr10=rgpd_k1_cp			(nran,xx)
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
#
# 7 plot the qq plots using the fitted quantiles
#
plot1(sxx,qq1	$cp_quantiles-dd1,	qq1	$ml_quantiles-dd1		,1)
plot1(sxx,qq2	$cp_quantiles-dd2,	qq2	$ml_quantiles-dd2		,2)
plot1(sxx,qq3	$cp_quantiles-dd1,	qq3	$ml_quantiles-dd1		,3)
plot1(sxx,qq4	$cp_quantiles,			qq4	$ml_quantiles				,4)
plot1(sxx,qq5	$cp_quantiles-dd1,	qq5	$ml_quantiles-dd1		,5)
plot1(sxx,qq5	$cp_quantiles-dd1,	qq6	$ml_quantiles-dd1		,6)
plot1(sxx,qq6	$cp_quantiles,			qq7	$ml_quantiles				,7)
plot1(sxx,qq8	$cp_quantiles,			qq8	$ml_quantiles				,8)
plot1(sxx,qq9	$cp_quantiles,			qq9	$ml_quantiles				,9)
plot1(sxx,qq10$cp_quantiles,			qq10$ml_quantiles				,10)
#
# 8 put mle parameters into a data frame
#
mle=matrix(NA,nmodels,2)
mle[1,1:1]	=qq1	$ml_params
mle[2,1:1]	=qq2	$ml_params
mle[3,1:1]	=qq3	$ml_params
mle[4,1:2]	=qq4	$ml_params
mle[5,1:2]	=qq5	$ml_params
mle[6,1:2]	=qq6	$ml_params
mle[7,1:2]	=qq7	$ml_params
mle[8,1:2]	=qq8	$ml_params
mle[9,1:2]	=qq9	$ml_params
mle[10,1:2]	=qq10	$ml_params
mleparams_df=data.frame(models,mle)
colnames(mleparams_df)=c("models","param 1","param 2")
#
# 10 extract maic, waic and logscores from the qresults
#
maic=matrix(0,nmodels)
maic[ 1]=qq1	$maic[3]
maic[ 2]=qq2	$maic[3]
maic[ 3]=qq3	$maic[3]
maic[ 4]=qq4	$maic[3]
maic[ 5]=qq5	$maic[3]
maic[ 6]=qq6	$maic[3]
maic[ 7]=qq7	$maic[3]
maic[ 8]=qq8	$maic[3]
maic[ 9]=qq9	$maic[3]
maic[10]=qq10	$maic[3]
#
waic1=matrix(0,nmodels)
waic1[ 1]=qq1	$waic1[3]
waic1[ 2]=qq2	$waic1[3]
waic1[ 3]=qq3	$waic1[3]
waic1[ 4]=qq4	$waic1[3]
waic1[ 5]=qq5	$waic1[3]
waic1[ 6]=qq6	$waic1[3]
waic1[ 7]=qq7	$waic1[3]
waic1[ 8]=qq8	$waic1[3]
waic1[ 9]=qq9	$waic1[3]
waic1[10]=qq10$waic1[3]
#
waic2=matrix(0,nmodels)
waic2[ 1]=qq1	$waic2[3]
waic2[ 2]=qq2	$waic2[3]
waic2[ 3]=qq3	$waic2[3]
waic2[ 4]=qq4	$waic2[3]
waic2[ 5]=qq5	$waic2[3]
waic2[ 6]=qq6	$waic2[3]
waic2[ 7]=qq7	$waic2[3]
waic2[ 8]=qq8	$waic2[3]
waic2[ 9]=qq9	$waic2[3]
waic2[10]=qq10$waic2[3]
#
lsc=matrix(0,nmodels)
lsc[ 1]=qq1$cp_oos_logscore
lsc[ 2]=qq2$cp_oos_logscore
lsc[ 3]=qq3$cp_oos_logscore
lsc[ 4]=qq4$cp_oos_logscore
lsc[ 5]=qq5$cp_oos_logscore
lsc[ 6]=qq6$cp_oos_logscore
lsc[ 7]=qq7$cp_oos_logscore
lsc[ 8]=qq8$cp_oos_logscore
lsc[ 9]=qq9$cp_oos_logscore
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
colnames(modelselection_df)=c(	"models",
								"maic","weights",
								"logscore","weights",
								"waic1","weights",
								"waic2","weights")
#
# 13 extract means from the q results
# -except for lnorm cp mean, which has to be calculated from randoms
#
ml_means=matrix(0,nmodels)
ml_means[ 1]=qq1	$ml_mean
ml_means[ 2]=qq2	$ml_mean
ml_means[ 3]=qq3	$ml_mean
ml_means[ 4]=qq4	$ml_mean
ml_means[ 5]=qq5	$ml_mean
ml_means[ 6]=qq6	$ml_mean
ml_means[ 7]=qq7	$ml_mean
ml_means[ 8]=qq8	$ml_mean
ml_means[ 9]=qq9	$ml_mean
ml_means[10]=qq10	$ml_mean
#
cp_means=matrix(0,nmodels)
cp_means[ 1]=qq1	$cp_mean
cp_means[ 2]=qq2	$cp_mean
cp_means[ 3]=qq3	$cp_mean
cp_means[ 4]=mean(rr4	$cp_deviates) #for lnorm, have to calculate from randoms
cp_means[ 5]=qq5	$cp_mean
cp_means[ 6]=qq6	$cp_mean
cp_means[ 7]=qq7	$cp_mean
cp_means[ 8]=qq8  $cp_mean
cp_means[ 9]=qq9	$cp_mean
cp_means[10]=qq10	$cp_mean

# calculate sds from the random deviates
# -and set to inf if we know they are inf

ml_sds=matrix(0,nmodels)
ml_sds[ 1]=sd(rr1	$ml_deviates)
ml_sds[ 2]=ifelse(mle[2,1]>1,sd(rr2$ml_deviates),Inf)
ml_sds[ 3]=sd(rr3	$ml_deviates)
ml_sds[ 4]=sd(rr4	$ml_deviates)
ml_sds[ 5]=sd(rr5	$ml_deviates)
ml_sds[ 6]=sd(rr6	$ml_deviates)
ml_sds[ 7]=sd(rr7	$ml_deviates)
ml_sds[ 8]=sd(rr8	$ml_deviates)
ml_sds[ 9]=sd(rr9	$ml_deviates)
ml_sds[10]=sd(rr10$ml_deviates)

cp_sds=matrix(0,nmodels)
cp_sds[ 1]=sd(rr1	$cp_deviates)
cp_sds[ 2]=Inf
cp_sds[ 3]=sd(rr3	$cp_deviates)
cp_sds[ 4]=sd(rr4	$cp_deviates)
cp_sds[ 5]=Inf
cp_sds[ 6]=sd(rr6	$cp_deviates)
cp_sds[ 7]=sd(rr7	$cp_deviates)
cp_sds[ 8]=Inf
cp_sds[ 9]=sd(rr9	$cp_deviates)
cp_sds[10]=Inf
#
# 14 make data frames to return the results for means and sds
#
means_df=data.frame(models,ml_means,cp_means,ml_sds,cp_sds)
colnames(means_df)=c(	"models","ml_means","cp_means","ml_sds","cp_sds")
#
# 15 return model selection results, means and sds
#
return(list(mleparams_df=mleparams_df,
	modelselection_df=modelselection_df,
	means_df=means_df))

}
