par(mfrow=c(2,3))
#
# example 1
x=fitdistcp::d32gnorm_k3_example_data_v1
p=c(1:9)/10
q=qgnorm_k3_cp(x,p,kbeta=4,rust=TRUE,nrust=1000)
xmin=min(q$ml_quantiles,q$cp_quantiles);
xmax=max(q$ml_quantiles,q$cp_quantiles);
plot(q$ml_quantiles,p,xlab="quantile estimates",xlim=c(xmin,xmax),
	sub="(from qgnorm_k3_cp)",
	main="gnorm: quantile estimates");
points(q$cp_quantiles,p,col="red",lwd=2)
points(q$ru_quantiles,p,col="blue")
#
# example 2: focus on upper tail
invp=c(2:20)
p=1-1/invp
ics=c(0,0,0)
q=qgnorm_k3_cp(x,p,kbeta=4,rust=TRUE,nrust=1000)
ymin=min(q$ml_quantiles,q$cp_quantiles);
ymax=max(q$ml_quantiles,q$cp_quantiles);
plot(invp,q$ml_quantiles,ylab="quantile estimates",ylim=c(ymin,ymax),
	sub="(from qgnorm_k3_cp)",
	main="Quantile estimates: tail");
points(invp,q$cp_quantiles,col="red",lwd=2)
points(invp,q$ru_quantiles,col="blue")
legend(x=10,y=ymin+0.25*(ymax-ymin),c("maxlik","cp dmgs","cp rust"),
	col=c("black","red","blue"),pch=1)
#
# example 3: random number generation from the predictive distribution
r=rgnorm_k3_cp(100000,x,kbeta=4,rust=TRUE)
plot(density(r$ml_deviates),
	sub="(from rgnorm_k3_cp)",
	main="Density (from random deviates)")
lines(density(r$cp_deviates),col="red",lwd=2)
lines(density(r$ru_deviates),col="blue")
#
# example 4: plot pdf and cdf
xax=c(-100:100)/20
pdfs=dgnorm_k3_cp(x,xax,kbeta=4,rust=TRUE)
cdfs=pgnorm_k3_cp(x,xax,kbeta=4,rust=TRUE)
#
plot(xax,pdfs$ml_pdf,"l",
	sub="(from dgnorm_k3_cp)",
	main="PDF")
lines(xax,pdfs$cp_pdf,col="red",lwd=2)
lines(xax,pdfs$ru_pdf,col="blue")
#
plot(xax,pmax(exp(-10),pdfs$ml_pdf),"l",log="y",
	sub="(from dgnorm_k3_cp)",
	main="Log-PDF")
lines(xax,pmax(exp(-10),pdfs$cp_pdf),col="red",lwd=2)
lines(xax,pmax(exp(-10),pdfs$ru_pdf),col="blue")
#
plot(xax,cdfs$ml_cdf,"l",
	sub="(from pgnorm_k3_cp)",
	main="CDF")
lines(xax,cdfs$cp_cdf,col="red",lwd=2)
lines(xax,cdfs$ru_cdf,col="blue")

