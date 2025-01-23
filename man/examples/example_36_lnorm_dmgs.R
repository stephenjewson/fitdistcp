#
# example 1
x=fitdistcp::d35lnorm_example_data_v1
p=c(1:9)/10
q=qlnorm_dmgs_cp(x,p)
xmin=min(q$ml_quantiles,q$cp_quantiles);
xmax=max(q$ml_quantiles,q$cp_quantiles);
plot(q$ml_quantiles,p,xlab="quantile estimates",xlim=c(xmin,xmax),
	sub="(from qlnorm_dmgs_cp)",
	main="Log-normal_DMGS: quantile estimates");
points(q$cp_quantiles,p,col="red",lwd=2)
#
# example 2: focus on upper tail
invp=c(2:20)
p=1-1/invp
q=qlnorm_dmgs_cp(x,p)
ymin=min(q$ml_quantiles,q$cp_quantiles);
ymax=max(q$ml_quantiles,q$cp_quantiles);
plot(invp,q$ml_quantiles,ylab="quantile estimates",ylim=c(ymin,ymax),
	sub="(from qlnorm_dmgs_cp)",
	main="Quantile estimates: tail");
points(invp,q$cp_quantiles,col="red",lwd=2)
legend(x=8,y=ymin+0.25*(ymax-ymin),c("maxlik","cp dmgs"),col=c("black","red"),pch=1)
#
# example 3: random number generation from the predictive distribution
r=rlnorm_dmgs_cp(1000,x)
mle=pmin(r$ml_deviates,25) # cut the long tail
rhp=pmin(r$cp_deviates,25)
plot(density(mle),xlim=c(0,20),
	sub="(from rlnorm_dmgs_cp)",
	main="Density (from random deviates)")
lines(density(rhp),col="red",lwd=2)
#
# example 4: plot pdf and cdf
xax=c(1:100)/20
pdfs=dlnorm_dmgs_cp(x,xax)
cdfs=plnorm_dmgs_cp(x,xax)
#
plot(xax,pdfs$ml_pdf,"l",
	sub="(from dlnorm_dmgs_cp)",
	main="PDF")
lines(xax,pdfs$cp_pdf,col="red",lwd=2)
#
plot(xax,pmax(exp(-10),pdfs$ml_pdf),"l",log="y",
	sub="(from dlnorm_dmgs_cp)",
	main="Log-PDF")
lines(xax,pmax(exp(-10),pdfs$cp_pdf),col="red",lwd=2)
#
plot(xax,cdfs$ml_cdf,"l",
	sub="(from plnorm_dmgs_cp)",
	main="CDF")
lines(xax,cdfs$cp_cdf,col="red",lwd=2)

