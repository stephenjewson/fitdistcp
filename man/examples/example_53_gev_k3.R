#
# example 1
kshape=-0.4
x=fitdistcp::d53gev_k3_example_data_v1
p=c(1:9)/10
q=qgev_k3_cp(x,p,kshape=kshape,rust=TRUE,nrust=1000)
xmin=min(q$ml_quantiles,q$cp_quantiles);
xmax=max(q$ml_quantiles,q$cp_quantiles);
plot(q$ml_quantiles,p,xlab="quantile estimates",xlim=c(xmin,xmax),
	sub="(from qgev_k3_cp)",
	main="GEV: quantile estimates");
points(q$cp_quantiles,p,col="red",lwd=2)
points(q$ru_quantiles,p,col="blue")
muhat=q$ml_params[1]
sghat=q$ml_params[2]
xi=kshape
qmax=ifelse(xi<0,muhat-sghat/xi,Inf)
cat(" ml_params=",q$ml_params,",")
cat(" qmax=",qmax,"\n")
#
# example 2: focus on upper tail
invp=c(2:20)
p=1-1/invp
ics=c(0,0,0)
q=qgev_k3_cp(x,p,kshape=kshape,rust=TRUE,nrust=1000)
ymin=min(q$ml_quantiles,q$cp_quantiles);
ymax=max(q$ml_quantiles,q$cp_quantiles);
plot(invp,q$ml_quantiles,ylab="quantile estimates",ylim=c(ymin,ymax),
	sub="(from qgev_k3_cp)",
	main="Quantile estimates: tail","p");
points(invp,q$cp_quantiles,col="red",lwd=2)
points(invp,q$ru_quantiles,col="blue")
# add marker for max value
lines(c(0,10000),rep(qmax,2),lty=2)
legend(x=10,y=ymin+0.25*(ymax-ymin),c("maxlik","cp dmgs","cp rust"),
	col=c("black","red","blue"),pch=1)
#
# example 3: random number generation from the predictive distribution
r=rgev_k3_cp(1000,x,kshape=kshape,rust=TRUE)
plot(density(r$cp_deviates),col="red",lwd=2,
	sub="(from rgev_k3_cp)",
	main="Density (from random deviates)")
lines(density(r$ml_deviates),col="black")
lines(density(r$ru_deviates),col="blue")
# add marker for max value
points(qmax,0,pch=4)

#
# example 4: plot pdf and cdf
xax=c(-100:100)/20
np=10000
pp=(c(1:np)-0.5)/np
pdfs_qq=qgev_k3_cp(x,p=pp,pdf=TRUE,kshape=kshape)
pdfs=dgev_k3_cp(x,xax,kshape=kshape,d1=0.01,fd2=0.01,rust=TRUE)
cdfs=pgev_k3_cp(x,xax,kshape=kshape,rust=TRUE)
#
plot(xax,pdfs$ml_pdf,"l",
	sub="(from dgev_k3_cp, qgev_k3_cp)",
	main="PDF")
lines(xax,pdfs$ru_pdf,col="blue")
points(qmax,0,pch=4)
lines(pdfs_qq$cp_quantiles,pdfs_qq$cp_pdf,col="orange")
legend(x=xax[1],y=0.99*max(pdfs$ml_pdf),
	c("maxlik","cp rust","cp from dmgs q"),
	col=c("black","blue","orange"),pch=1,cex=0.7)
#
plot(xax,pmax(exp(-10),pdfs$ml_pdf),"l",log="y",
	sub="(from dgev_k3_cp, qgev_k3_cp)",
	main="Log-PDF")
lines(xax,pmax(exp(-10),pdfs$ru_pdf),col="blue")
lines(pdfs_qq$cp_quantiles,pmax(exp(-10),pdfs_qq$cp_pdf),col="orange")
legend(x=xax[1],y=0.99*max(pdfs$ml_pdf),
	c("maxlik","cp rust","cp from dmgs q"),
	col=c("black","blue","orange"),pch=1,cex=0.7)
#
plot(xax,cdfs$ml_cdf,"l",
	sub="(from pgev_k3_cp)",
	main="CDF")
lines(xax,cdfs$ru_cdf,col="blue")
lines(pdfs_qq$cp_quantiles,pp,col="orange")
points(qmax,0,pch=4)
#

