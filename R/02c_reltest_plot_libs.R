#' Plotting routine for testppm
#'
#' @description
#' Plots 9 diagnostics related to predictive probability matching.
#'
#' @param model					which distribution to test. Possibles values are
#' @param ntrials				the number of trials to run. 5000 typically gives good results.
#' @param nrepeats			the number of entire repeats of the test to run, to check for convergence
#' @param nx						the length of the training data.
#' @param params				values for the parameters for the specified distribution
#' @param nmethods			the number of methods being tested
#' @param alpha					the values of alpha being tested
#' @param freqexceeded	the exceedance counts
testppm_plot=function(model,ntrials,nrepeats,nx,params,nmethods,alpha,freqexceeded){
#
nalpha=length(alpha)
rp=1/alpha
actualrp=array(0,c(nmethods,nrepeats,nalpha))
diffs=array(0,c(nmethods,nrepeats,nalpha))
pcdiffs=array(0,c(nmethods,nrepeats,nalpha))
for (ip in 1:nmethods){
	for (ir in 1:nrepeats){
		actualrp[ip,ir,]=1/freqexceeded[ip,ir,]
		diffs[ip,ir,]=alpha[]-freqexceeded[ip,ir,] #this sign convention means tail too thin -> negative values
		pcdiffs[ip,ir,]=100*(diffs[ip,ir,]/alpha[]) #old method
	}
}
#
# graphics parameters
#
legcex=0.8
textcex=0.9
par(mfrow=c(3,3))
names=c(	"ML","CP")
cols=c("blue","red")
diag=c(-100,500)
dummy=pcdiffs[1,1,]
basex=c(-1000:1000)
#
# a) exceedance frequency vs alpha, over the whole range
#
# set up axes but don't plot anything
	plot(alpha,dummy,xlim=c(1,0),type="n",ylim=c(100,0),
		main="(a) Exceedance Freq. vs Alpha",
		ylab="Freq of exceedances (%)",
		xlab="Alpha",
		col="black",lwd=2)
# plot the repeats
	for (ir in 1:nrepeats){
		lines(alpha,100*freqexceeded[1,ir,],col="blue",lwd=2,lty=1)
		lines(alpha,100*freqexceeded[2,ir,],col="red",lwd=2,lty=1)
	}
	lines(diag,100*diag,lty=1)
	legend(1,0,names,col=cols,lty=1,lwd=2,cex=legcex)
	rect(0.1,10,0,0)
# dots on x axis
	points(alpha,rep(100,nalpha),cex=0.1)
# top left
	text(1.07,50,"too few exceedances",pos=4,cex=textcex)
	text(1.07,60,"tail too fat",pos=4,cex=textcex)
# bottom right
	text(0.6,85,"too many exceedances",pos=4,cex=textcex)
	text(0.6,95,"tail too thin",pos=4,cex=textcex)
#
# b) same as (a), but zoomed into the upper tail
#
# set up axes but don't plot anything
	plot(alpha,dummy,xlim=c(0.1,0),type="n",ylim=c(10,0),
		main="(b) Same but zoomed",
		ylab="Freq of exceedances (%)",
		xlab="Alpha",
		col="black",lwd=2)
# plot the repeats
	for (ir in 1:nrepeats){
		lines(alpha,100*freqexceeded[1,ir,],col="blue",lwd=2,lty=1)
		lines(alpha,100*freqexceeded[2,ir,],col="red",lwd=2,lty=1)
	}
	lines(diag,100*diag,col="black",lty=1)
	legend(0.1,0,names,col=cols,lty=1,lwd=2,cex=legcex)
	points(alpha,rep(10,nalpha),cex=0.1)
	text(0.04,9.5,"too thin",pos=4,cex=textcex)
	text(0.1,5.0,"too fat",pos=4,cex=textcex)
#
# (c) difference between actual frequency and alpha (whole range)
#
	y0=100*diffs[1,1,]
	ymin=100*min(diffs)
	ymax=100*max(diffs)
	dy=ymax-ymin
	ddy=0.1*dy
	plot(alpha,y0,xlim=c(1,0),type="n",ylim=c(ymin-ddy,ymax+ddy),
		main="(c) Difference (Freq-Alpha)",
		ylab="Probability Diff x 100",
		xlab="Alpha",
		col="black",lwd=2)
	for (ir in 1:nrepeats){
		y1i=100*diffs[1,ir,]
		y2i=100*diffs[2,ir,]
		lines(alpha,y1i,col="blue",lwd=2,lty=1)
		lines(alpha,y2i,col="red",lwd=2,lty=1)
	}
	lines(basex,basex*0,lty=1)
	legend(1,ymin+1.10*dy,names,col=cols,lty=1,lwd=2,cex=legcex)
	points(alpha,rep(-1,nalpha),cex=0.1)
	text(0.5,ymin-0.05*dy,"too thin",pos=4,cex=textcex)
	text(0.5,ymin+1.05*dy,"too fat",pos=4,cex=textcex)
#
# (d) same as (c) but zoomed into upper tail
#
	y0=100*diffs[1,1,]
	ymin=100*min(diffs)
	ymax=100*max(diffs)
	dy=ymax-ymin
	ddy=0.1*dy
	plot(alpha,y0,xlim=c(0.1,0),type="n",ylim=c(ymin-ddy,ymax+ddy),
		main="(d) Same but zoomed",
		ylab="Probability Diff x 100",
		xlab="Alpha",
		col="black",lwd=2)
	for (ir in 1:nrepeats){
		y1i=100*diffs[1,ir,]
		y2i=100*diffs[2,ir,]
		lines(alpha,y1i,col="blue",lwd=2,lty=1)
		lines(alpha,y2i,col="red",lwd=2,lty=1)
	}
	lines(basex,basex*0,lty=1)
	legend(0.1,ymin+1.10*dy,names,col=cols,lty=1,lwd=2,cex=legcex)
	points(alpha,rep(-1,nalpha),cex=0.1)
	text(0.05,ymin-0.05*dy,"too thin",pos=4,cex=textcex)
	text(0.05,ymin+1.05*dy,"too fat",pos=4,cex=textcex)
#
# (e) Same as (c), but vs rp
#
	y0=100*diffs[1,1,]
	ymin=100*min(diffs)
	ymax=100*max(diffs)
	dy=ymax-ymin
	ddy=0.1*dy
	plot(rp,y0,xlim=c(1,200),type="n",ylim=c(ymin-ddy,ymax+ddy),
		main="(e) Difference vs RP",
		ylab="Probability Diff x 100",
		xlab="RP",
		col="black",lwd=2)
	for (ir in 1:nrepeats){
		y1i=100*diffs[1,ir,]
		y2i=100*diffs[2,ir,]
		lines(rp,y1i,col="blue",lwd=2,lty=1)
		lines(rp,y2i,col="red",lwd=2,lty=1)
	}
	lines(basex,basex*0,lty=1)
	legend(1,ymin+1.10*dy,names,col=cols,lty=1,lwd=2,cex=legcex)
	points(1/alpha,rep(-1,nalpha),cex=0.1)
	text(80,ymin-0.05*dy,"too thin",pos=4,cex=textcex)
	text(80,ymin+1.05*dy,"too fat",pos=4,cex=textcex)
#
# f) RP 1-250
#
	plot(rp,pcdiffs[1,1,],xlim=c(1,250),type="n",ylim=c(1,250),
		main="(f) RP vs RP (0-250)",
		ylab="Actual RP of predictions on average",
		xlab="intended RP",
		col="black",lwd=2)
	for (ir in 1:nrepeats){
#		lines(rp,rp,col="grey",lwd=2)
		lines(rp,actualrp[1,ir,],col="blue",lwd=2,lty=1)
		lines(rp,actualrp[2,ir,],col="red",lwd=2,lty=1)
	}
	lines(diag,diag,col="black",lty=1)
	legend(1,250,names,col=cols,lty=1,lwd=2,cex=legcex)
	points((1/alpha),rep(0,nalpha),cex=0.1)
	text(150,20,"too thin",pos=4,cex=textcex)
	text(0,130,"too fat",pos=4,cex=textcex)
#
# g) RP 1-50 thick
#
	plot(rp,pcdiffs[1,1,],xlim=c(1,50),type="n",ylim=c(1,50),
		main="(g) RP vs RP (zoom: 0-50)",
		ylab="Actual RP of predictions on average",
		xlab="intended RP",
		col="black",lwd=2)
	for (ir in 1:nrepeats){
		lines(rp,actualrp[1,ir,],col="blue",lwd=2,lty=1)
		lines(rp,actualrp[2,ir,],col="red",lwd=2,lty=1)
	}
	lines(diag,diag,col="grey",lty=1)
	legend(1,50,names,col=cols,lty=1,lwd=2,cex=legcex)
	points(1/alpha,rep(0,nalpha),cex=0.1)
	text(25,2,"too thin",pos=4,cex=textcex)
	text(0,25,"too fat",pos=4,cex=textcex)
#
# h) pc probability bias (vs alpha zoom)
#
	y0=pcdiffs[1,1,]
	ymin=min(pcdiffs)
	ymax=max(pcdiffs)
	dy=ymax-ymin
	ddy=0.1*dy
	plot(alpha,y0,xlim=c(0.1,0),type="n",ylim=c(ymin-ddy,ymax+ddy),
		main="(h) % Difference zoomed",
		ylab="Probability Diff in %%",
		col="black",lwd=2)
	for (ir in 1:nrepeats){
# plot the 3 different models from within each routine (typically ml, ml+, ob)
		y1i=pcdiffs[1,ir,]
		y2i=pcdiffs[2,ir,]
		lines(alpha,y1i,col="blue",lwd=2,lty=1)
		lines(alpha,y2i,col="red",lwd=2,lty=1)
	}
	lines(basex,basex*0,lty=1)
	legend(0.1,ymin+0.3*dy,names,col=cols,lty=1,lwd=2,cex=legcex)
	points(alpha,rep(-100,nalpha),cex=0.1)
#
# i) pc probability bias (vs RP)
#
	y0=pcdiffs[1,1,]
	ymin=min(pcdiffs)
	ymax=max(pcdiffs)
	dy=ymax-ymin
	ddy=0.1*dy
	plot(rp,y0,xlim=c(1,200),type="n",ylim=c(ymin-ddy,ymax+ddy),
		main="(i) % Difference vs RP",
		ylab="Probability Diff in %%",
		col="black",lwd=2)
	for (ir in 1:nrepeats){
# plot the 3 different models from within each routine (typically ml, ml+, ob)
		y1i=pcdiffs[1,ir,]
		y2i=pcdiffs[2,ir,]
		lines(rp,y1i,col="blue",lwd=2,lty=1)
		lines(rp,y2i,col="red",lwd=2,lty=1)

#		lines(rp,pcdiffs[nmethods,ir,],col="grey",lwd=2)
	}
	lines(basex,basex*0,lty=1)
	legend(0,ymin+0.3*dy,names,col=cols,lty=1,lwd=2,cex=legcex)
	points((1/alpha),rep(-100,nalpha),cex=0.1)
#
# title
#
}
