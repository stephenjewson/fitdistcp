#' Plotting routine for ppm_testppm
#'
#' @description
#' Plots 9 diagnostics related to predictive probability matching.
#'
#' @param model			which distribution to test. Possibles values are
#'	"gev",
#'	"gpd",
#'	"gev_p1".
#' @param ntrials				the number of trials o run. 5000 typically gives good results.
#' @param nrepeats			the number of entire repeats of the test to run, to check for convergence
#' @param nx						the length of the training data.
#' @param params				values for the parameters for the specified distribution
#' @param nmethods			the number of methods being tested
#' @param alpha					the values of alpha being tested
#' @param freqexceeded	the exceedance counts
#' @param case				  there are 3 cases (must be set to case=1 except for my testing)
reltest2_plot=function(model,ntrials,nrepeats,nx,params,
	nmethods,alpha,freqexceeded,case){
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
legcex=0.6;
textcex=0.9;
par(mfrow=c(3,3))
if(case==1)names=c("Flat","JP","ML","RHP-ML","RHP-Flat")
if(case==2)names=c("Flat","JP","ML","MPD")
if(case==3)names=c("Flat","JP","ML","RHP-ML","RHP-Flat","MPD")
if(case==4)names=c("Flat","JP","ML","RHP-Flat","MPD")
if(case==1)cols=c("purple","orange","grey","blue","red")
if(case==2)cols=c("purple","orange","grey","black")
if(case==3)cols=c("purple","orange","grey","blue","red","black")
if(case==4)cols=c("purple","orange","grey","red","black")
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
		for (ip in 1:nmethods){
			lines(alpha,100*freqexceeded[ip,ir,],col=cols[ip],lwd=2,lty=1)
		}
	}
	lines(diag,100*diag,lty=1)
	legend(1,0,names,col=cols,lty=1,lwd=2,cex=legcex)
	rect(0.1,10,0,0)
# dots on x axis
	points(alpha,rep(100,nalpha),cex=0.1)
# top left
	text(1.07,60,"too few exceedances",pos=4,cex=textcex)
	text(1.07,70,"tail too fat",pos=4,cex=textcex)
# bottom right
	text(0.8,85,"too many exceedances",pos=4,cex=textcex)
	text(0.8,95,"tail too thin",pos=4,cex=textcex)
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
		for (ip in 1:nmethods){
			lines(alpha,100*freqexceeded[ip,ir,],col=cols[ip],lwd=2,lty=1)
		}
	}
	lines(diag,100*diag,col="black",lty=1)
	legend(0.1,0,names,col=cols,lty=1,lwd=2,cex=legcex)
	points(alpha,rep(10,nalpha),cex=0.1)
	text(0.04,9.5,"too thin",pos=4,cex=textcex)
	text(0.1,6.0,"too fat",pos=4,cex=textcex)
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
		for (ip in 1:nmethods){
			lines(alpha,100*diffs[ip,ir,],col=cols[ip],lwd=2,lty=1)
		}
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
	plot(alpha,-100*diffs[1,1,],xlim=c(0.1,0),type="n",ylim=c(ymin-ddy,ymax+ddy),
		main="(d) Same but zoomed",
		ylab="Probability Diff x 100",
		xlab="Alpha",
		col="black",lwd=2)
	for (ir in 1:nrepeats){
		for (ip in 1:nmethods){
			lines(alpha,100*diffs[ip,ir,],col=cols[ip],lwd=2,lty=1)
		}
	}
	lines(basex,basex*0,lty=1)
	legend(1,ymin+1.10*dy,names,col=cols,lty=1,lwd=2,cex=legcex)
	points(alpha,rep(-1,nalpha),cex=0.1)
	text(0.5,ymin-0.05*dy,"too thin",pos=4,cex=textcex)
	text(0.5,ymin+1.05*dy,"too fat",pos=4,cex=textcex)
#
# (e) Same as (c), but vs rp
#
	y0=100*diffs[1,1,]
	ymin=100*min(diffs)
	ymax=100*max(diffs)
	dy=ymax-ymin
	ddy=0.1*dy
	plot(rp,-100*diffs[1,1,],xlim=c(1,200),type="n",ylim=c(ymin-ddy,ymax+ddy),
		main="(e) Difference vs RP",
		ylab="Probability Diff x 100",
		xlab="RP",
		col="black",lwd=2)
	for (ir in 1:nrepeats){
		for (ip in 1:nmethods){
			lines(rp,100*diffs[ip,ir,],col=cols[ip],lwd=2,lty=1)
		}
	}
	lines(basex,basex*0,lty=1)
	legend(1,ymin+1.10*dy,names,col=cols,lty=1,lwd=2,cex=legcex)
	points(1/alpha,rep(-1,nalpha),cex=0.1)
	text(100,ymin-0.05*dy,"too thin",pos=4,cex=textcex)
	text(100,ymin+1.05*dy,"too fat",pos=4,cex=textcex)
#
# f) RP 1-250
#
	plot(rp,pcdiffs[1,1,],xlim=c(1,250),type="n",ylim=c(1,250),
		main="(f) RP vs RP (0-250)",
		ylab="Actual RP of predictions on average",
		xlab="intended RP",
		col="black",lwd=2)
	for (ir in 1:nrepeats){
		for (ip in 1:nmethods){
			lines(rp,actualrp[ip,ir,],col=cols[ip],lwd=2,lty=1)
		}
	}
	lines(diag,diag,col="black",lty=1)
	legend(1,250,names,col=cols,lty=1,lwd=2,cex=legcex)
	points((1/alpha),rep(0,nalpha),cex=0.1)
	text(150,20,"too thin",pos=4,cex=textcex)
	text(0,100,"too fat",pos=4,cex=textcex)
#
# g) RP 1-50 thick
#
	plot(rp,pcdiffs[1,1,],xlim=c(1,50),type="n",ylim=c(1,50),
		main="(g) RP vs RP (zoom: 0-50)",
		ylab="Actual RP of predictions on average",
		xlab="intended RP",
		col="black",lwd=2)
	for (ir in 1:nrepeats){
		for (ip in 1:nmethods){
			lines(rp,actualrp[ip,ir,],col=cols[ip],lwd=2,lty=1)
		}
	}
	lines(diag,diag,col="grey",lty=1)
	legend(1,50,names,col=cols,lty=1,lwd=2,cex=legcex)
	points(1/alpha,rep(0,nalpha),cex=0.1)
	text(25,2,"too thin",pos=4,cex=textcex)
	text(0,20,"too fat",pos=4,cex=textcex)
#
# h) pc probability bias (vs alpha zoom)
#
	y0=pcdiffs[1,1,]
	ymin=min(pcdiffs)
	ymax=max(pcdiffs)
	dy=ymax-ymin
	ddy=0.1*dy
	plot(alpha,pcdiffs[1,1,],xlim=c(0.1,0),type="n",ylim=c(ymin-ddy,ymax+ddy),
		main="(h) % Difference zoomed",
		ylab="Probability Diff in %%",
		col="black",lwd=2)
	for (ir in 1:nrepeats){
		for (ip in 1:nmethods){
			lines(alpha,pcdiffs[ip,ir,],col=cols[ip],lwd=2,lty=1)
		}
	}
	lines(basex,basex*0,lty=1)
	legend(0.1,ymin+1.10*dy,names,col=cols,lty=1,lwd=2,cex=legcex)
	points(alpha,rep(-100,nalpha),cex=0.1)
#
# i) pc probability bias (vs RP)
#
	y0=pcdiffs[1,1,]
	ymin=min(pcdiffs)
	ymax=max(pcdiffs)
	dy=ymax-ymin
	ddy=0.1*dy
	plot(rp,pcdiffs[1,1,],xlim=c(1,200),type="n",ylim=c(ymin-ddy,ymax+ddy),
		main="(i) % Difference vs RP",
		ylab="Probability Diff in %%",
		col="black",lwd=2)
	for (ir in 1:nrepeats){
		for (ip in 1:nmethods){
			lines(rp,pcdiffs[ip,ir,],col=cols[ip],lwd=2,lty=1)
		}
	}
	lines(basex,basex*0,lty=1)
	legend(1,ymin+1.10*dy,names,col=cols,lty=1,lwd=2,cex=legcex)
	points((1/alpha),rep(-100,nalpha),cex=0.1)
#
# title
#
}
