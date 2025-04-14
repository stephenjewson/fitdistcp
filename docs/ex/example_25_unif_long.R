par(mfrow=c(2,3))
#
# example 1
x=fitdistcp::d25unif_example_data_v1
cat("length(x)=",length(x),"\n")
p=c(1:9)/10
q=qunif_cp(x,p)
xmin=min(q$ml_quantiles,q$cp_quantiles);
xmax=max(q$ml_quantiles,q$cp_quantiles);
plot(q$ml_quantiles,p,xlab="quantile estimates",xlim=c(xmin,xmax),
	sub="(from qunif_cp)",
	main="unif: quantile estimates");
points(q$cp_quantiles,p,col="red",lwd=2)
#
# example 2: focus on upper tail
invp=c(2:20)
p=1-1/invp
q=qunif_cp(x,p)
ymin=min(q$ml_quantiles,q$cp_quantiles);
ymax=max(q$ml_quantiles,q$cp_quantiles);
plot(invp,q$ml_quantiles,ylab="quantile estimates",ylim=c(ymin,ymax),
	sub="(from qunif_cp)",
	main="Quantile estimates: tail");
points(invp,q$cp_quantiles,col="red",lwd=2)
legend(x=10,y=ymin+0.25*(ymax-ymin),c("maxlik","cp analytic"),col=c("black","red"),pch=1)
#
# example 3: random number generation from the predictive distribution
r=runif_cp(100000,x)
plot(density(r$ml_deviates),
	sub="(from runif_cp)",
	main="Density (from random deviates)")
lines(density(r$cp_deviates),col="red",lwd=2)
#
# example 4: plot pdf and cdf
xax=c(-100:200)/100
pdfs=dunif_cp(x,xax)
cdfs=punif_cp(x,xax)
#
plot(xax,pdfs$ml_pdf,"l",
	sub="(from dunif_cp)",
	main="PDF")
lines(xax,pdfs$cp_pdf,col="red",lwd=2)
#
plot(xax,pmax(exp(-10),pdfs$ml_pdf),"l",log="y",
	sub="(from dunif_cp)",
	main="Log-PDF")
lines(xax,pmax(exp(-10),pdfs$cp_pdf),col="red",lwd=2)
#
plot(xax,cdfs$ml_cdf,"l",
	sub="(from punif_cp)",
	main="CDF")
lines(xax,cdfs$cp_cdf,col="red",lwd=2)

