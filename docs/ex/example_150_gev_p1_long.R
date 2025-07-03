par(mfrow=c(2,3))
#
# example 1
x=fitdistcp::d150gev_p1_example_data_v1_x
tt=fitdistcp::d150gev_p1_example_data_v1_t
p=c(1:9)/10
n0=10
q=qgev_p1_cp(x=x,t=tt,n0=n0,t0=NA,p=p,rust=TRUE,nrust=1000)
xmin=min(q$ml_quantiles,q$cp_quantiles);
xmax=max(q$ml_quantiles,q$cp_quantiles);
plot(q$ml_quantiles,p,xlab="quantile estimates",xlim=c(xmin,xmax),
	sub="(from qgev_p1_cp)",
	main="GEVD w/ p1: quantile estimates");
points(q$cp_quantiles,p,col="red",lwd=2)
points(q$ru_quantiles,p,col="blue",lwd=2)
cat(" ml_params=",q$ml_params,"\n")
#
# example 2: focus on upper tail
invp=c(2:20)
p=1-1/invp
n0=10
q=qgev_p1_cp(x=x,t=tt,n0=n0,t0=NA,p=p,rust=TRUE,nrust=1000)
ymin=min(q$ml_quantiles,q$cp_quantiles);
ymax=max(q$ml_quantiles,q$cp_quantiles);
plot(invp,q$ml_quantiles,ylab="quantile estimates",ylim=c(ymin,ymax),
	sub="(from qgev_p1_cp)",
	main="Quantile estimates: tail");
points(invp,q$cp_quantiles,col="red",lwd=2)
points(invp,q$ru_quantiles,col="blue",lwd=2)
legend(x=10,y=ymin+0.25*(ymax-ymin),c("maxlik","cp dmgs","cp rust"),
	col=c("black","red","blue"),pch=1)
#
# example 3: random number generation from the predictive distribution
n0=3
r=rgev_p1_cp(n=10000,x=x,t=tt,n0=n0,t0=NA,rust=TRUE)
plot(density(r$ml_deviates),
	sub="(from rgev_p1_cp)",
	main="Density (from random deviates)")
lines(density(r$cp_deviates),col="red",lwd=2)
lines(density(r$ru_deviates),col="blue",lwd=2)
#
# example 4: plot pdf and cdf
#
nx=50
n0=3
xax=c(-75:100)/20
np=10000
pp=(c(1:np)-0.5)/np
pdfs_qq=qgev_p1_cp(x=x,t=tt,n0=n0,p=pp,pdf=TRUE,rust=TRUE)
pdfs=dgev_p1_cp(x=x,t=tt,n0=n0,y=xax,rust=TRUE)
cdfs=pgev_p1_cp(x=x,t=tt,n0=n0,y=xax,rust=TRUE)
#
# panel 4
#
plot(xax,pdfs$ml_pdf,"l",
	sub="(from dgev_p1_cp, qgev_p1_cp)",
	main="PDF")
lines(xax,pdfs$ru_pdf,col="blue",lwd=2)
lines(pdfs_qq$cp_quantiles,pdfs_qq$cp_pdf,col="orange")
legend(x=xax[1],y=0.99*max(pdfs$ml_pdf),
	c("maxlik","cp rust","cp from dmgs q"),
	col=c("black","blue","orange"),pch=1,cex=0.7)
#
#
# panel 5
#
plot(xax,pmax(exp(-10),pdfs$ml_pdf),"l",log="y",
	sub="(from dgev_p1_cp, qgev_p1_cp)",
	main="Log-PDF")
lines(xax,pmax(exp(-10),pdfs$ru_pdf),col="blue",lwd=2)
lines(pdfs_qq$cp_quantiles,pmax(exp(-10),pdfs_qq$cp_pdf),col="orange")
legend(x=xax[1],y=0.99*max(pdfs$ml_pdf),
	c("maxlik","cp rust","cp from dmgs q"),
	col=c("black","blue","orange"),pch=1,cex=0.7)
#
#
# panel 6 (note that I don't currently produce a cdf from qq, just a pdf)
#
plot(xax,cdfs$ml_cdf,"l",
	sub="(from pgev_p1_cp)",
	main="CDF")
lines(xax,cdfs$ru_cdf,col="blue",lwd=2)
lines(pdfs_qq$cp_quantiles,pp,col="orange")
legend(x=xax[1],y=1,
	c("maxlik","cp rust","cp from dmgs q"),
	col=c("black","blue","orange"),pch=1,cex=0.7)

