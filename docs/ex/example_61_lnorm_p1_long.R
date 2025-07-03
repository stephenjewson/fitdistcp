par(mfrow=c(2,3))
#
# example 1
x=fitdistcp::d61lnorm_p1_example_data_v1_x
tt=fitdistcp::d61lnorm_p1_example_data_v1_t
p=c(1:9)/10
n0=10
q=qlnorm_p1_cp(x,tt,n0=n0,p=p,rust=TRUE,nrust=1000)
xmin=min(q$ml_quantiles,q$cp_quantiles);
xmax=max(q$ml_quantiles,q$cp_quantiles);
plot(q$ml_quantiles,p,xlab="quantile estimates",xlim=c(xmin,xmax),
	sub="(from qlnorm_p1_cp)",
	main="Log-Normal w/ p1: quantile estimates");
points(q$cp_quantiles,p,col="red",lwd=2)
points(q$ru_quantiles,p,col="blue")
#
# example 2: focus on upper tail
invp=c(2:20)
p=1-1/invp
n0=10
q=qlnorm_p1_cp(x,tt,n0=n0,p=p,rust=TRUE,nrust=1000)
ymin=min(q$ml_quantiles,q$cp_quantiles);
ymax=max(q$ml_quantiles,q$cp_quantiles);
plot(invp,q$ml_quantiles,ylab="quantile estimates",ylim=c(ymin,ymax),
	sub="(from qlnorm_p1_cp)",
	main="Quantile estimates: tail");
points(invp,q$cp_quantiles,col="red",lwd=2)
points(invp,q$ru_quantiles,col="blue")
legend(x=8,y=ymin+0.25*(ymax-ymin),c("maxlik","cp analytic","cp rust"),
	col=c("black","red","blue"),pch=1)
#
# example 3: random number generation from the predictive distribution
n0=3
r=rlnorm_p1_cp(1000,x,tt,n0=n0,rust=TRUE)
mle=pmin(r$ml_deviates,25) #cut off the long tail
rhp=pmin(r$cp_deviates,25)
plot(density(mle),xlim=c(0,20),
	sub="(from rlnorm_p1_cp)",
	main="Density (from random deviates)")
lines(density(rhp),col="red",lwd=2)
lines(density(r$ru_deviates),col="blue")
#
# example 4: plot pdf and cdf
n0=3
xax=c(1:100)/40
pdfs=dlnorm_p1_cp(x=x,t=tt,t0=NA,n0=n0,y=xax,rust=TRUE)
cdfs=plnorm_p1_cp(x=x,t=tt,t0=NA,n0=n0,y=xax,rust=TRUE)
#
plot(xax,pdfs$ml_pdf,"l",
	sub="(from dlnorm_p1_cp)",
	main="PDF")
lines(xax,pdfs$cp_pdf,col="red",lwd=2)
lines(xax,pdfs$ru_pdf,col="blue")
#
plot(xax,pmax(exp(-10),pdfs$ml_pdf),"l",log="y",
	sub="(from dlnorm_p1_cp)",
	main="Log-PDF")
lines(xax,pmax(exp(-10),pdfs$cp_pdf),col="red",lwd=2)
lines(xax,pmax(exp(-10),pdfs$ru_pdf),col="blue")
#
plot(xax,cdfs$ml_cdf,"l",
	sub="(from plnorm_p1_cp)",
	main="CDF")
lines(xax,cdfs$cp_cdf,col="red",lwd=2)
lines(xax,cdfs$ru_cdf,col="blue")
