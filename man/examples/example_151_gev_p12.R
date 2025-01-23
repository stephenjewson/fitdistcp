#
# example 1
x=fitdistcp::d151gev_p12_example_data_v1_x
tt=fitdistcp::d151gev_p12_example_data_v1_t
t1=tt[,1]
t2=tt[,2]
p=c(1:9)/10
n01=10
n02=10
q=qgev_p12_cp(x=x,t1=t1,t2=t2,n01=n01,n02=n02,t01=NA,t02=NA,p=p,rust=TRUE,nrust=1000)
xmin=min(q$ml_quantiles,q$cp_quantiles);
xmax=max(q$ml_quantiles,q$cp_quantiles);
plot(q$ml_quantiles,p,xlab="quantile estimates",xlim=c(xmin,xmax),
	sub="(from qgev_p12_cp)",
	main="GEVD w/ p12: quantile estimates");
points(q$cp_quantiles,p,col="red",lwd=2)
points(q$ru_quantiles,p,col="blue",lwd=2)
cat(" ml_params=",q$ml_params,"\n")
#
# example 2: focus on upper tail
invp=c(2:20)
p=1-1/invp
n01=10
n02=10
q=qgev_p12_cp(x=x,t1=t1,t2=t2,n01=n01,n02=n02,t01=NA,t02=NA,p=p,rust=TRUE,nrust=1000)
ymin=min(q$ml_quantiles,q$cp_quantiles);
ymax=max(q$ml_quantiles,q$cp_quantiles);
plot(invp,q$ml_quantiles,ylab="quantile estimates",ylim=c(ymin,ymax),
	sub="(from qgev_p12_cp)",
	main="Quantile estimates: tail");
points(invp,q$cp_quantiles,col="red",lwd=2)
points(invp,q$ru_quantiles,col="blue",lwd=2)
legend(x=10,y=ymin+0.25*(ymax-ymin),c("maxlik","cp dmgs","cp rust"),
	col=c("black","red","blue"),pch=1)
#
# example 3: random number generation from the predictive distribution
n01=3
n02=3
r=rgev_p12_cp(n=10000,x=x,t1=t1,t2=t2,n01=n01,n02=n02,t01=NA,t02=NA,rust=TRUE)
plot(density(r$ml_deviates),
	sub="(from rgev_p12_cp)",
	main="Density (from random deviates)")
lines(density(r$cp_deviates),col="red",lwd=2)
lines(density(r$ru_deviates),col="blue",lwd=2)
#
# example 4: plot pdf and cdf
#
nx=50
n01=3
n02=3
xax=c(-75:100)/20
np=10000
pp=(c(1:np)-0.5)/np
pdfs_qq=qgev_p12_cp(x=x,t1=t1,t2=t2,n01=n01,n02=n02,p=pp,pdf=TRUE,rust=FALSE)
pdfs=dgev_p12_cp(x=x,t1=t1,t2=t2,n01=n01,n02=n02,y=xax,rust=FALSE)
cdfs=pgev_p12_cp(x=x,t1=t1,t2=t2,n01=n01,n02=n02,y=xax,rust=FALSE)
#
plot(xax,pdfs$ml_pdf,"l",
	sub="(from dgev_p12_cp, qgev_p12_cp)",
	main="PDF")
lines(pdfs_qq$cp_quantiles,pdfs_qq$cp_pdf,col="orange")
legend(x=xax[1],y=0.99*max(pdfs$ml_pdf),
	c("maxlik","cp from dmgs q"),
	col=c("black","orange"),pch=1,cex=0.7)
#
# panel 5
#
plot(xax,pmax(exp(-10),pdfs$ml_pdf),"l",log="y",
	sub="(from dgev_p12_cp, qgev_p12_cp)",
	main="Log-PDF")
lines(pdfs_qq$cp_quantiles,pmax(exp(-10),pdfs_qq$cp_pdf),col="orange")
legend(x=xax[1],y=0.99*max(pdfs$ml_pdf),
	c("maxlik","cp from dmgs q"),
	col=c("black","orange"),pch=1,cex=0.7)
#
#
# panel 6
#
plot(xax,cdfs$ml_cdf,"l",
	sub="(from pgev_p12_cp)",
	main="CDF")
lines(pdfs_qq$cp_quantiles,pp,col="orange")
legend(x=xax[1],y=1,
	c("maxlik","cp from dmgs q"),
	col=c("black","orange"),pch=1,cex=0.7)

