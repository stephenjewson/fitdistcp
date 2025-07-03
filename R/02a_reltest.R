#' Evaluation of Reliability for Models in the \code{fitdistcp} Package
#'
#' @inherit man references seealso
#'
#' @description
#' Uses simulations to evaluate the reliability of
#' the predictive quantiles produced by the \code{q****_cp} routines in the \code{fitdistcp} package.
#'
#' @param model			which distribution to test. Possibles values are
#'	"\code{exp}",
#'	"\code{pareto_k1}",
#'	"\code{halfnorm}",
#'	"\code{unif}",
#'	"\code{norm}",
#'	"\code{norm_dmgs}",
#'	"\code{gnorm_k3}",
#'	"\code{lnorm}",
#'	"\code{lnorm_dmgs}",
#'	"\code{logis}",
#'	"\code{lst_k3}",
#'	"\code{cauchy}",
#'	"\code{gumbel}",
#'	"\code{frechet_k1}",
#'	"\code{weibull}",
#'	"\code{gev_k3}",
#'	"\code{exp_p1}",
#'	"\code{pareto_p1k3}",
#'	"\code{norm_p1}",
#'	"\code{lnorm_p1}",
#'	"\code{logis_p1}",
#'	"\code{lst_p1k4}",
#'	"\code{cauchy_p1}",
#'	"\code{gumbel_p1}",
#'	"\code{frechet_p2k1}",
#'	"\code{weibull_p2}",
#'	"\code{gev_p1k4}",
#'	"\code{norm_p12}",
#'	"\code{lst_p12k5}",
#'	"\code{gamma}",
#'	"\code{invgamma}",
#'	"\code{invgauss}",
#'	"\code{gev}",
#'	"\code{gpd_k1}",
#'	"\code{gev_p1}".
#'	"\code{gev_p12}".
#'	"\code{gev_p123}".
#' @param ntrials		the number of trials to run. 5000 typically gives good results.
#' @param nrepeats  the number of entire repeats of the test to run, to check for convergence. 3 is a good choice.
#' @param nx				the length of the training data to use.
#' @param params		values for the parameters for the specified distribution
#' @param alpha			the exceedance probability values at which to test
#' @param plotflag	logical to turn the plotting on and off
#' @param verbose		logical to turn loop counting on and off
#' @param dmgs			logical to turn DMGS calculations on and off (to optimize speed for maxlik only calculations)
#' @param	debug			logical for turning debug messages on and off
#' @param aderivs		logical for whether to use analytic derivatives (instead of numerical)
#' @param unbiasedv logical for whether to use the unbiased variance instead of maxlik (for the normal)
#' @param pwm				logical for whether to use PWM instead of maxlik (for the GEV)
#' @param minxi			minimum value for EVT shape parameter
#' @param maxxi			maximum value for EVT shape parameter
#'
#' @return
#' A plot showing 9 different reliability checks, and a list containing
#' various outputs, including the
#' probabilities shown in the plot.
#'
#' @details
#'
#' The maximum likelihood quantiles (plotted in blue) do not give good reliability.
#' They typically underestimate the tails (see panel (f)).
#'
#' For
#' "\code{exp}",
#' "\code{pareto_k1}",
#' "\code{unif}",
#' "\code{norm}",
#' "\code{lnorm}",
#' "\code{norm_p1}" and
#' "\code{lnorm_p1}",
#' the calibrating prior quantiles are calculated using the right Haar prior
#' and an exact solution for the Bayesian prediction integral.
#' They will converge towards exact reliability with a large enough number of trials,
#' for any sample size.
#'
#' For
#' "\code{halfnorm}",
#' "\code{norm_dmgs}",
#' "\code{lnorm_dmgs}",
#' "\code{gnorm_k3}",
#' "\code{logis}",
#' "\code{lst_k3}",
#' "\code{cauchy}",
#' "\code{gumbel}",
#' "\code{frechet_k1}",
#' "\code{weibull}",
#' "\code{gev_k3}",
#' "\code{exp_p1}",
#' "\code{pareto_p1k3}",
#' "\code{gumbel_p1}",
#' "\code{logis_p1}" and
#' "\code{lst_p1k4}"
#' "\code{cauchy_p1}",
#' "\code{gumbel_p1}",
#' "\code{frechet_p2k1}",
#' "\code{weibull_p2}",
#' "\code{gev_p1k4}",
#' "\code{norm_p12}",
#' "\code{lst_p12k5}"
#' the calibrating prior quantiles are calculated using the right Haar prior,
#' with the DMGS asymptotic solution for the Bayesian prediction integral.
#' They will converge towards good reliability with a large enough number of trials,
#' with the only deviation from exact reliability being due to the neglect of
#' higher order terms in the asymptotic expansion.
#' They will converge towards exact reliability with a large enough number of trials
#' and a large enough sample size.
#'
#' For
#' "\code{gamma}",
#' "\code{invgamma}",
#' "\code{invgauss}",
#' "\code{gev}",
#' "\code{gpd_k1}" and
#' "\code{gev_p1}",
#' "\code{gev_p12}",
#' "\code{gev_p123}",
#' the calibrating prior quantiles are calculated using the "\code{fitdistcp}"
#' recommended calibrating priors,
#' with the DMGS asymptotic solution for the Bayesian prediction integral.
#' The chosen priors give reasonably good reliability with a
#' large enough number of trials,
#' and for large sample sizes, but may give poor reliability for small
#' sample sizes (e.g., n<20).
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_01_reltest.R
#'
#' @export
#'
reltest=function(model="exp",ntrials=1000,nrepeats=3,nx=20,params=c(1),
	alpha=seq(0.005,0.995,0.005),
	plotflag=TRUE,verbose=TRUE,dmgs=TRUE,debug=FALSE,aderivs=TRUE,
	unbiasedv=FALSE,pwm=FALSE,minxi=-10,maxxi=10){
#
# intro
#
	ny=nx*nx #increasing y doesn't seem to slow things down, so why not do 100
#	alpha=c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.25,0.2,0.15,0.1,1/20,1/30,1/40,1/50,1/60,1/70,1/80,1/90,1/100,1/110,1/120,1/130,1/140,1/150,1/160,1/170,1/180,1/190,1/200)
	pp=1-alpha
	rp=1/alpha
	nalpha=length(alpha)
	tt=c(1:nx)/nx
	tt=tt-mean(tt)
	n0=nx
	tt1=tt
	tt2=tt
	tt3=tt
# maybe the following don't work for gev_p12 and gev_p123 actually
#	tt2=log(runif(nx,min=0.9,max=1.1)) #this is for gev sd, which is always exp(a+bt)
#	tt3=runif(nx,min=-0.1,max=0.1) #this is for gev shape, which is always a+bt
	tt0=tt[n0]
	tt10=tt1[n0]
	tt20=tt2[n0]
	tt30=tt3[n0]
	known_scale=1 #pareto
	known_loc=0 #frechet
#
# the big testingloop
#
	if(dmgs==FALSE)	nmethods=1
	if(dmgs==TRUE)	nmethods=2
	ep1=array(0,c(nmethods,nalpha))
	epsum=array(0,c(nmethods,nrepeats,nalpha))
	maxepsum=array(0,c(nrepeats))
	if(verbose)message("\nmodel=",model)
	for (ir in 1:nrepeats){
		if(verbose)message(" \nrepeat:",ir)
		for (it in 1:ntrials){
			if(verbose)message(it)
#
# make random training and testing data
#
			xx=reltest_simulate(model,nx,tt,tt1,tt2,tt3,params,minxi=minxi,maxxi=maxxi)
#
# make predictions
# -params is only passed in to provide the various known parameters
			pred=reltest_predict(model,xx,tt,tt1,tt2,tt3,n0,n10=n0,n20=n0,n30=n0,pp,params,dmgs=dmgs,
				debug=debug,aderivs=aderivs,unbiasedv=unbiasedv,pwm=pwm,minxi=minxi,maxxi=maxxi)
#
# make ep
#
			for (im in 1:nmethods){
				ep1[im,]=reltest_makeep(model,pred$pred[im,],tt0,tt10,tt20,tt30,params)
				epsum[im,ir,]=epsum[im,ir,]+ep1[im,] #summing over trials
			}
#
# for evt models, make ep for the max
			if((model=="gev")||(model=="gpd_k1")||(model=="gev_p1")||(model=="gev_p12")){
				maxep1=reltest_makemaxep(model,pred$ml_max,tt0,tt10,tt20,tt30,params)
				maxepsum[ir]=maxepsum[ir]+maxep1 #summing over trials
			}

		} #end of trials loop
	} #end of repeats loop
#
# normalize counts and calculate differences versus expected
#
	ep=epsum/ntrials
	maxep=maxepsum/ntrials
#
# 8 copy for no good reason
#
	freqexceeded=array(0,c(nmethods,nrepeats,nalpha))
	maxfreqexceeded=array(0,c(nrepeats))
	for (ir in 1:nrepeats){
	 	for (ip in 1:nmethods){
 			freqexceeded[ip,ir,]=(ep[ip,ir,])
 		}
# and for evt models, the max
		if((model=="gev")||(model=="gpd_k1")||(model=="gev_p1")||(model=="gev_p12")){
			maxfreqexceeded[ir]=(maxep[ir])
		}
 	}
#
# plotting
#
	if(plotflag)testppm_plot(model,ntrials,nrepeats,nx,params,nmethods,alpha,freqexceeded)

	return(list(nmethods=nmethods,alpha=alpha,freqexceeded=freqexceeded,
		maxfreqexceeded=maxfreqexceeded))
}


