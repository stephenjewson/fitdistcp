#
# example 1
x=fitdistcp::d11pareto_k2_example_data_v1
p=c(1:9)/10
q=qpareto_k2_cp(x,p,rust=TRUE,nrust=1000)
xmin=min(q$ml_quantiles,q$cp_quantiles)
xmax=max(q$ml_quantiles,q$cp_quantiles)
plot(q$ml_quantiles,p,xlab="quantile estimates",xlim=c(xmin,xmax),
	sub="(from qpareto_k2_cp)",
	main="Pareto: quantile estimates")
points(q$cp_quantiles,p,col="red",lwd=2)
points(q$ru_quantiles,p,col="blue")
#
# example 2: focus on upper tail
invp=c(2:20)
p=1-1/invp
q=qpareto_k2_cp(x,p,rust=TRUE,nrust=1000)
ymin=min(q$ml_quantiles,q$cp_quantiles)
ymax=max(q$ml_quantiles,q$cp_quantiles)
plot(invp,q$ml_quantiles,ylab="quantile estimates",ylim=c(ymin,ymax),
	sub="(from qpareto_k2_cp)",
	main="Quantile estimates: tail")
points(invp,q$cp_quantiles,col="red",lwd=2)
points(invp,q$ru_quantiles,col="blue")
legend(x=8,y=ymin+0.25*(ymax-ymin),c("maxlik","cp analytic","cp rust"),
	col=c("black","red","blue"),pch=1)
#
# example 3: random number generation from the predictive distribution
r=rpareto_k2_cp(100000,x,rust=TRUE)
mle=pmin(r$ml_deviates,25) # cut the long tail
rhp=pmin(r$cp_deviates,25)
plot(density(mle),xlim=c(1,5),
	sub="(from rpareto_k2_cp)",
	main="Density (from random deviates)")
lines(density(rhp),col="red",lwd=2)
lines(density(r$ru_deviates),col="blue")

# example 4: plot pdf and cdf
xax=c(10:100)/10
pdfs=dpareto_k2_cp(x,xax,rust=TRUE)
cdfs=ppareto_k2_cp(x,xax,rust=TRUE)
#
plot(xax,pdfs$ml_pdf,"l",
	sub="(from dpareto_k2_cp)",
	main="PDF")
lines(xax,pdfs$cp_pdf,col="red",lwd=2)
lines(xax,pdfs$ru_pdf,col="blue")
#
plot(xax,pmax(exp(-10),pdfs$ml_pdf),"l",log="y",
	sub="(from dpareto_k2_cp)",
	main="Log-PDF")
lines(xax,pmax(exp(-10),pdfs$cp_pdf),col="red",lwd=2)
lines(xax,pmax(exp(-10),pdfs$ru_pdf),col="blue")
#
plot(xax,cdfs$ml_cdf,"l",
	sub="(from ppareto_k2_cp)",
	main="CDF")
lines(xax,cdfs$cp_cdf,col="red",lwd=2)
lines(xax,cdfs$ru_cdf,col="blue")

