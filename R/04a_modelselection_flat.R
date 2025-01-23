#' Illustration of Model Selection Among 18 Distributions from the \code{fitdistcp} Package
#'
#' @description
#' Applies model selection using AIC, WAIC1, WAIC2 and leave-one-out logscore
#' to the input data \eqn{x},
#' for 18 models in the \code{fitdistcp} packages
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
#' GEVD and GPD are not included by default, because they are temperamental in that
#' they do not work if the shape parameter is extreme.
#'
#' @param x 			data vector
#' @param evd		a flag to include GEVD and GPD from the analysis
#'
#' @details
#' The 18 models are:
#' \code{exp},
#' \code{pareto_k2},
#' \code{halfnorm},
#' \code{norm},
#' \code{lnorm},
#' \code{gnorm_k3},
#' \code{gumbel},
#' \code{frechet_k1},
#' \code{weibull},
#' \code{gev_k3},
#' \code{logis},
#' \code{lst_k3},
#' \code{cauchy},
#' \code{gamma},
#' \code{invgamma},
#' \code{invgauss},
#' \code{gev} and
#' \code{gpd_k1}.

#' @returns
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
#' @example man/examples/example_03_modelselection_flat.R
#'
#' @export
#'
modelselection_flat=function(x,evd=FALSE){
#
if(evd)cat("evd is selected\n")

	xx=x
#
if(evd){
	models=c(	"exp","pareto_k2","halfnorm","norm",
						"lnorm","gnorm_k3","gumbel","frechet_k1",
						"weibull","gev_k3","logis","lst_k3",
						"cauchy","gamma","invgamma","invgauss","gev","gpd_k1")
} else {
	models=c(	"exp","pareto_k2","halfnorm","norm",
						"lnorm","gnorm_k3","gumbel","frechet_k1",
						"weibull","gev_k3","logis","lst_k3",
						"cauchy","gamma","invgamma","invgauss")
}
cat("\n")
#
# for the models that require x>0, shift the data if there are negative values
# for the pareto, shift the data so that x>1
#
dd1=0
dd2=0
if(min(xx)<0)	dd1=0-min(xx)+0.0001
if(min(xx)<1)	dd2=1-min(xx)+0.0001
cat("Note that to make this code work, I adjust the input data, to make sure all values are positive.\n")
cat(" before adjustment  : min(xx)    =",min(xx),"\n")
cat(" after  adjustment 1: min(xx+dd1)=",min(xx+dd1),"\n")
cat(" after  adjustment 2: min(xx+dd2)=",min(xx+dd2),"\n")
#
# fit the models
#
nx=length(xx)
pp1=(c(1:nx)-0.5)/nx
qq1=qexp_cp					(xx+dd1,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq2=qpareto_k2_cp		(xx+dd2,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq3=qhalfnorm_cp		(xx+dd1,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq4=qnorm_cp				(xx,			pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq5=qlnorm_cp				(xx+dd1,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq6=qgnorm_k3_cp		(xx,			pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq7=qgumbel_cp			(xx,			pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq8=qfrechet_k1_cp	(xx+dd1,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq9=qweibull_cp			(xx+dd1,	pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq10=qgev_k3_cp			(xx,			pp1,waicscores=TRUE,means=TRUE,kshape=0.0)
qq11=qlogis_cp			(xx,			pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq12=qlst_k3_cp			(xx,			pp1,waicscores=TRUE,logscores=TRUE,means=TRUE,kdf=5)
qq13=qcauchy_cp			(xx,			pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq14=qgamma_cp			(xx,			pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq15=qinvgamma_cp		(xx,			pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
qq16=qinvgauss_cp		(xx,			pp1,waicscores=TRUE,logscores=TRUE,means=TRUE)
if(evd)qq17=qgev_cp				(xx,			pp1,waicscores=TRUE,means=TRUE)
if(evd)qq18=qgpd_k1_cp		(xx,			pp1,waicscores=TRUE,means=TRUE,kloc=0)
#
# simulate for means
#
nran=100000
rr1=rexp_cp					(nran,xx+dd1)
rr2=rpareto_k2_cp		(nran,xx+dd2)
rr3=rhalfnorm_cp		(nran,xx+dd1)
rr4=rnorm_cp				(nran,xx)
rr5=rgnorm_k3_cp		(nran,xx)
rr6=rlnorm_cp				(nran,xx+dd1)
rr7=rgumbel_cp			(nran,xx)
rr8=rfrechet_k1_cp	(nran,xx+dd1)
rr9=rweibull_cp			(nran,xx+dd1)
rr10=rgev_k3_cp			(nran,xx,kshape=0.0)
rr11=rlogis_cp			(nran,xx)
rr12=rlst_k3_cp			(nran,xx,kdf=5)
rr13=rcauchy_cp			(nran,xx)
rr14=rgamma_cp			(nran,xx)
rr15=rinvgamma_cp		(nran,xx)
rr16=rinvgauss_cp		(nran,xx)
if(evd)rr17=rgev_cp				(nran,xx)
if(evd)rr18=rgpd_k1_cp		(nran,xx)
#
# qq plots
# -horiz axis is the data we are trying to model
# -vertical axis in black are quantiles from fitting using maxlik
# -vertical axis in red are quantiles from fitting with parameter uncertainty
# -vertical axis is based on so-called 'Hazen plotting position'
# https://glossary.ametsoc.org/wiki/Plotting_position
#
par(mfrow=c(4,5))
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
plot1(sxx,qq1 $cp_quantiles-dd1,qq1$ml_quantiles-dd1		,1)
plot1(sxx,qq2 $cp_quantiles-dd2,qq2$ml_quantiles-dd2		,2)
plot1(sxx,qq3 $cp_quantiles-dd1,qq3$ml_quantiles-dd1		,3)
plot1(sxx,qq4 $cp_quantiles,		qq4$ml_quantiles				,4)
plot1(sxx,qq5 $cp_quantiles,		qq5$ml_quantiles				,5)
plot1(sxx,qq6 $cp_quantiles-dd1,qq6$ml_quantiles-dd1		,6)
plot1(sxx,qq7 $cp_quantiles,		qq7$ml_quantiles				,7)
plot1(sxx,qq8 $cp_quantiles-dd1,qq8$ml_quantiles-dd1		,8)
plot1(sxx,qq9 $cp_quantiles-dd1,qq9$ml_quantiles-dd1		,9)
plot1(sxx,qq10 $cp_quantiles,		qq10$ml_quantiles				,10)
plot1(sxx,qq11$cp_quantiles,		qq11$ml_quantiles				,11)
plot1(sxx,qq12$cp_quantiles,		qq12$ml_quantiles				,12)
plot1(sxx,qq13$cp_quantiles,		qq13$ml_quantiles				,13)
plot1(sxx,qq14$cp_quantiles,		qq14$ml_quantiles				,14)
plot1(sxx,qq15$cp_quantiles,		qq15$ml_quantiles				,15)
plot1(sxx,qq16$cp_quantiles,		qq16$ml_quantiles				,16)
if(evd)plot1(sxx,qq17$cp_quantiles,		qq17$ml_quantiles				,17)
if(evd)plot1(sxx,qq18$cp_quantiles,		qq18$ml_quantiles				,18)
#
# print parameter values
#
if(evd){
	nmodels=18
} else {
	nmodels=16
}
mle=matrix(0,nmodels,3)
mle[1,1:1]	=qq1$ml_params
mle[2,1:1]	=qq2$ml_params
mle[3,1:1]	=qq3$ml_params
mle[4,1:2]	=qq4$ml_params
mle[5,1:2]	=qq5$ml_params
mle[6,1:2]	=qq6$ml_params
mle[7,1:2]	=qq7$ml_params
mle[8,1:2]	=qq8$ml_params
mle[9,1:2]	=qq9$ml_params
mle[10,1:2]	=qq10$ml_params
mle[11,1:2]	=qq11$ml_params
mle[12,1:2]	=qq12$ml_params
mle[13,1:2]	=qq13$ml_params
mle[14,1:2]	=qq13$ml_params
mle[15,1:2]	=qq13$ml_params
mle[16,1:2]	=qq13$ml_params
if(evd)mle[17,1:3]	=qq17$ml_params
if(evd)mle[18,1:2]	=qq18$ml_params
mledf=data.frame(models,mle)
print(mledf)
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
maic[ 8]=qq8$maic[3]
maic[ 9]=qq9$maic[3]
maic[10]=qq10$maic[3]
maic[11]=qq11$maic[3]
maic[12]=qq12$maic[3]
maic[13]=qq13$maic[3]
maic[14]=qq14$maic[3]
maic[15]=qq15$maic[3]
maic[16]=qq16$maic[3]
if(evd)maic[17]=qq17$maic[3]
if(evd)maic[18]=qq18$maic[3]
waic1=matrix(0,nmodels)
waic1[ 1]=qq1$waic1[3]
waic1[ 2]=qq2$waic1[3]
waic1[ 3]=qq3$waic1[3]
waic1[ 4]=qq4$waic1[3]
waic1[ 5]=qq5$waic1[3]
waic1[ 6]=qq6$waic1[3]
waic1[ 7]=qq7$waic1[3]
waic1[ 8]=qq8$waic1[3]
waic1[ 9]=qq9$waic1[3]
waic1[10]=qq10$waic1[3]
waic1[11]=qq11$waic1[3]
waic1[12]=qq12$waic1[3]
waic1[13]=qq13$waic1[3]
waic1[14]=qq14$waic1[3]
waic1[15]=qq15$waic1[3]
waic1[16]=qq16$waic1[3]
if(evd)waic1[17]=qq17$waic1[3]
if(evd)waic1[18]=qq18$waic1[3]
waic2=matrix(0,nmodels)
waic2[ 1]=qq1$waic2[3]
waic2[ 2]=qq2$waic2[3]
waic2[ 3]=qq3$waic2[3]
waic2[ 4]=qq4$waic2[3]
waic2[ 5]=qq5$waic2[3]
waic2[ 6]=qq6$waic2[3]
waic2[ 7]=qq7$waic2[3]
waic2[ 8]=qq8$waic2[3]
waic2[ 9]=qq9$waic2[3]
waic2[10]=qq10$waic2[3]
waic2[11]=qq11$waic2[3]
waic2[12]=qq12$waic2[3]
waic2[13]=qq13$waic2[3]
waic2[14]=qq14$waic2[3]
waic2[15]=qq15$waic2[3]
waic2[16]=qq16$waic2[3]
if(evd)waic2[17]=qq17$waic2[3]
if(evd)waic2[18]=qq18$waic2[3]
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
lsc[11]=qq11$cp_oos_logscore
lsc[12]=qq12$cp_oos_logscore
lsc[13]=qq13$cp_oos_logscore
lsc[14]=qq14$cp_oos_logscore
lsc[15]=qq15$cp_oos_logscore
lsc[16]=qq16$cp_oos_logscore
if(evd)lsc[17]=NA
if(evd)lsc[18]=NA
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
cat("\n")
maic=maic
modelselection_df=data.frame(models,maic,maic_weights,lsc,logscore_weights,waic1,waic1_weights,waic2,waic2_weights)
colnames(modelselection_df)=c(	"models",
								"maic","weights",
								"logscore","weights",
								"waic1","weights",
								"waic2","weights")
cat("\n")
#
ml_means=matrix(0,nmodels)
ml_means[ 1]=qq1$ml_mean
ml_means[ 2]=Inf
ml_means[ 3]=qq3$ml_mean
ml_means[ 4]=qq4$ml_mean
ml_means[ 5]=qq5$ml_mean
ml_means[ 6]=qq5$ml_mean
ml_means[ 7]=qq6$ml_mean
ml_means[ 8]=Inf
ml_means[ 9]=qq9$ml_mean
ml_means[10]=qq10$ml_mean
ml_means[11]=qq11$ml_mean
ml_means[12]=qq12$ml_mean
ml_means[13]=qq13$ml_mean
ml_means[14]=qq14$ml_mean
ml_means[15]=qq15$ml_mean
ml_means[16]=qq16$ml_mean
if(evd)ml_means[17]=qq17$ml_mean
if(evd)ml_means[18]=qq18$ml_mean

cp_means=matrix(0,nmodels)
cp_means[ 1]=qq1$cp_mean
cp_means[ 2]=qq2$cp_mean
cp_means[ 3]=qq3$cp_mean
cp_means[ 4]=qq4$cp_mean
cp_means[ 5]=qq5$cp_mean
cp_means[ 6]=mean(rr5$cp_deviates) #lognormal special case
cp_means[ 7]=qq7$cp_mean
cp_means[ 8]=qq8$cp_mean
cp_means[ 9]=qq9$cp_mean
cp_means[10]=qq10$cp_mean
cp_means[11]=qq11$cp_mean
cp_means[12]=qq12$cp_mean
cp_means[13]=qq13$cp_mean
cp_means[14]=qq14$cp_mean
cp_means[15]=qq15$cp_mean
cp_means[16]=qq16$cp_mean
if(evd)cp_means[17]=qq17$cp_mean
if(evd)cp_means[18]=qq18$cp_mean

ml_sds=matrix(0,nmodels)
ml_sds[ 1]=sd(rr1$ml_deviates)
ml_sds[ 2]=Inf
ml_sds[ 3]=sd(rr3$ml_deviates)
ml_sds[ 4]=sd(rr4$ml_deviates)
ml_sds[ 5]=sd(rr5$ml_deviates)
ml_sds[ 6]=sd(rr6$ml_deviates)
ml_sds[ 7]=sd(rr7$ml_deviates)
ml_sds[ 8]=Inf
ml_sds[ 9]=sd(rr9$ml_deviates)
ml_sds[10]=sd(rr10$ml_deviates)
ml_sds[11]=sd(rr11$ml_deviates)
ml_sds[12]=sd(rr12$ml_deviates)
ml_sds[13]=sd(rr13$ml_deviates)
ml_sds[14]=sd(rr14$ml_deviates)
ml_sds[15]=sd(rr15$ml_deviates)
ml_sds[16]=sd(rr16$ml_deviates)
if(evd)ml_sds[17]=sd(rr17$ml_deviates)
if(evd)ml_sds[18]=sd(rr18$ml_deviates)

cp_sds=matrix(0,nmodels)
cp_sds[ 1]=sd(rr1$cp_deviates)
cp_sds[ 2]=Inf
cp_sds[ 3]=sd(rr3$cp_deviates)
cp_sds[ 4]=sd(rr4$cp_deviates)
cp_sds[ 5]=sd(rr5$cp_deviates)
cp_sds[ 6]=sd(rr6$cp_deviates)
cp_sds[ 7]=sd(rr7$cp_deviates)
cp_sds[ 8]=Inf
cp_sds[ 9]=sd(rr9$cp_deviates)
cp_sds[10]=sd(rr10$cp_deviates)
cp_sds[11]=sd(rr11$cp_deviates)
cp_sds[12]=sd(rr12$cp_deviates)
cp_sds[13]=sd(rr13$cp_deviates)
cp_sds[14]=sd(rr14$cp_deviates)
cp_sds[15]=sd(rr15$cp_deviates)
cp_sds[16]=sd(rr16$cp_deviates)
if(evd)cp_sds[17]=sd(rr17$cp_deviates)
if(evd)cp_sds[18]=sd(rr18$cp_deviates)

means_df=data.frame(models,ml_means,cp_means,ml_sds,cp_sds)
colnames(means_df)=c(	"models","ml_means","cp_means","ml_sds","cp_sds")

return(list(modelselection_df=modelselection_df,means_df=means_df))

}
