#
# example 1
shape=-0.4
x=fitdistcp::d110gev_example_data_v1
p=c(1:9)/10
q=qgev_cp(x,p,rust=TRUE,nrust=1000)
xmin=min(q$ml_quantiles,q$cp_quantiles);
xmax=max(q$ml_quantiles,q$cp_quantiles);
plot(q$ml_quantiles,p,xlab="quantile estimates",xlim=c(xmin,xmax),
	sub="(from qgev_cp)",
	main="GEVD: quantile estimates");
points(q$cp_quantiles,p,col="red",lwd=2)
points(q$ru_quantiles,p,col="blue",lwd=2)
cat(" ml_params=",q$ml_params,"\n")
#
# example 2: focus on upper tail
invp=c(2:20)
p=1-1/invp
q=qgev_cp(x,p,rust=TRUE)
ymin=min(q$ml_quantiles,q$cp_quantiles);
ymax=max(q$ml_quantiles,q$cp_quantiles);
plot(invp,q$ml_quantiles,ylab="quantile estimates",ylim=c(ymin,ymax),
	sub="(from qgev_cp)",
	main="Quantile estimates: tail");
points(invp,q$cp_quantiles,col="red",lwd=2)
points(invp,q$ru_quantiles,col="blue",lwd=2)
legend(x=5,y=ymin+0.2*(ymax-ymin),
	c("maxlik","cp dmgs","cp rust"),
	col=c("black","red","blue"),pch=1)
#
# example 3: random number generation from the predictive distribution
r=rgev_cp(1000,x,rust=TRUE)
plot(density(r$ml_deviates),
	sub="(from rgev_cp)",
	main="Density (from random deviates)")
lines(density(r$cp_deviates),col="red",lwd=2)
lines(density(r$ru_deviates),col="blue",lwd=2)
#
# example 4: plot pdf and cdf
xax=c(-100:100)/20
np=10000
pp=(c(1:np)-0.5)/np
pdfs_qq=qgev_cp(x,p=pp,pdf=TRUE,rust=TRUE)
pdfs_dd=dgev_cp(x,xax,rust=TRUE)
cdfs=pgev_cp(x,xax,rust=TRUE)
#
plot(xax,pdfs_dd$ml_pdf,"l",
	sub="(from dgev_cp, qgev_cp)",
	main="PDF")
lines(xax,pdfs_dd$ru_pdf,col="blue",lwd=2)
lines(pdfs_qq$cp_quantiles,pdfs_qq$cp_pdf,col="orange")
legend(x=xax[1],y=0.99*max(pdfs_dd$ml_pdf),
	c("maxlik","cp rust","cp from dmgs q"),
	col=c("black","blue","orange"),pch=1,cex=0.7)
#
plot(xax,pmax(exp(-10),pdfs_dd$ml_pdf),"l",log="y",
	sub="(from dgev_cp, qgev_cp)",
	main="Log-PDF")
lines(xax,pmax(exp(-10),pdfs_dd$ru_pdf),col="blue",lwd=2)
lines(pdfs_qq$cp_quantiles,pmax(exp(-10),pdfs_qq$cp_pdf),col="orange")
legend(x=-5,y=ymin+0.25*(ymax-ymin),c("maxlik","cp rust","cp from dmgs q"),
	col=c("black","blue","orange"),pch=1,cex=0.7)
#
plot(xax,cdfs$ml_cdf,"l",
	sub="(from pgev_cp)",
	main="CDF")
lines(xax,cdfs$ru_cdf,col="blue",lwd=2)
lines(pdfs_qq$cp_quantiles,pp,col="orange")
#
