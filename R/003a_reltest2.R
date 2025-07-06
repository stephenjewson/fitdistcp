#' Evaluation of Reliability for Certain Additional Models in the \code{fitdistcp} Package
#'
#' @inherit man references seealso
#'
#' @description
#' This routine is mainly for reproducing certain results in Jewson et al. (2025),
#' and not of general interest.
#'
#' It uses simulations to evaluate the reliability of
#' the predictive quantiles produced by the
#' \code{qgev_cp}, \code{ggpd_cp} and \code{qgev_p1_cp}
#' routines in the \code{fitdistcp} package.
#' For each model, results for 5 models are calculated.
#' This is to illustrate that the calibrating prior predictions dominate
#' the \code{ml}, \code{flat}, \code{crhp_ml} and \code{jp} predictions,
#' in terms of reliability.
#'
#' @param model			which distribution to test. Possibles values are
#'	"\code{gev}",
#'	"\code{gpd_k1}",
#'	"\code{gev_p1}".
#' @param ntrials		the number of trials to run. 5000 typically gives good results.
#' @param nrepeats  the number of entire repeats of the test to run, to check for convergence. 3 is a good choice.
#' @param nx				the length of the training data.
#' @param params		values for the parameters for the specified distribution
#' @param alpha			the alpha values at which to test
#' @param plotflag	logical to turn the plotting on and off
#' @param verbose		logical to turn loop counting on and off
#'
#' @return
#' A plot showing 9 different reliability checks, and a list containing
#' various outputs, including the
#' probabilities shown in the plot.
#'
#' @details
#'
#' The maximum likelihood quantiles (plotted in blue) do not give good reliability. They typically underestimate the tails (see panel (f)).
#'
#' The \code{cp} predictive quantiles generally give reasonably good reliability, especially for sample sizes of ~100.
#' The other predictions generally give poor reliability.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_002_reltest2.R
#'
#' @export
#'
reltest2=function(model="gev",ntrials=100,nrepeats=3,nx=50,
	params=c(0,1,0),alpha=seq(0.005,0.995,0.005),plotflag=TRUE,verbose=TRUE){
#
# intro
#
	ny=nx*nx #increasing y doesn't seem to slow things down, so why not do 100
	pp=1-alpha
	nalpha=length(alpha)
	tt=c(1:nx)/nx
	n0=1
	tt0=tt[n0]
#
# define nmethods and case by model
	cases=reltest2_cases(model="gev",nx,params)
	nmethods=cases$nmethods
	case=cases$case
	if(verbose)message("nmethods=",nmethods)
	if(verbose)message("case=",case)
#
# the big testing loop
#
	ep1=array(0,c(nmethods,nalpha))
	epsum=array(0,c(nmethods,nrepeats,nalpha))
	if(verbose)message("model=",model)
	for (ir in 1:nrepeats){
		if(verbose)message(" \nrepeat:",ir)
		for (it in 1:ntrials){
			if(verbose)message(it)

# make xx and flags
			simop=reltest2_simulate(model,nx,tt,params)
			xx=simop$xx
			rh_ml_flag=simop$rh_ml_flag
			cp_flag=simop$cp_flag
#
# make predictions
# -params is only passed in to provide the various known parameters
			pred=reltest2_predict(model,xx,tt,n0,pp,params,case,nmethods)
#
# make ep
#
			for (im in 1:nmethods){
				ep1[im,]=reltest2_makeep(model,pred[im,],tt0,params)
				epsum[im,ir,]=epsum[im,ir,]+ep1[im,] #summing over trials
			}

	} #end of trials loop
} #end of repeats loop
	ep=epsum/ntrials
#
# 8 new metric: REF: ratio of exceedance frequencies
#
	freqexceeded=array(0,c(nmethods,nrepeats,nalpha))
 	for (ip in 1:nmethods){
 		for (ir in 1:nrepeats){
 			freqexceeded[ip,ir,]=(ep[ip,ir,])
 		}
 	}
#
# plotting
#
	if(plotflag)reltest2_plot(model,ntrials,nrepeats,nx,params,
		nmethods,alpha,freqexceeded,case)

	return(list(nmethods=nmethods,alpha=alpha,freqexceeded=freqexceeded))

}


