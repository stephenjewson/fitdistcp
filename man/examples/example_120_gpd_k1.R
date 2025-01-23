#
# example 1
x=fitdistcp::d120gpd_k1_example_data_v1
p=c(1:9)/10
q=qgpd_k1_cp(x,p,rust=TRUE,nrust=1000)
xmin=min(q$ml_quantiles,q$cp_quantiles);
xmax=max(q$ml_quantiles,q$cp_quantiles);
plot(q$ml_quantiles,p,xlab="quantile estimates",xlim=c(xmin,xmax),
	sub="(from qgpd_k1_cp)",
	main="GPD: quantile estimates");
points(q$cp_quantiles,p,col="red",lwd=2)
points(q$ru_quantiles,p,col="blue",lwd=2)
cat(" ml_params=",q$ml_params,"\n")
#
# example 2: focus on upper tail
invp=c(2:20)
p=1-1/invp
q=qgpd_k1_cp(x,p,rust=TRUE,nrust=1000)
ymin=min(q$ml_quantiles,q$cp_quantiles);
ymax=max(q$ml_quantiles,q$cp_quantiles);
plot(invp,q$ml_quantiles,ylab="quantile estimates",ylim=c(ymin,ymax),
	sub="(from qgpd_k1_cp)",
	main="Quantile estimates: tail");
points(invp,q$cp_quantiles,col="red",lwd=2)
points(invp,q$ru_quantiles,col="blue",lwd=2)
legend(x=5,y=ymin+0.2*(ymax-ymin),
	c("maxlik","cp dmgs","cp rust"),
	col=c("black","red","blue"),pch=1)
#
# example 3: random number generation from the predictive distribution
nn=1000
#setting nn too large (i.e. 100,000)
#tends to break the density plotting routine
r=rgpd_k1_cp(nn,x,rust=TRUE)
dens_mle=density(r$ml_deviates)
dens_cp=density(r$cp_deviates)
dens_ru=density(r$ru_deviates)
ymax=max(dens_mle$y,dens_cp$y)
plot(dens_mle,xlim=c(0,10),ylim=c(0,ymax),
	sub="(from rgpd_k1_cp)",
	main="Density (from random deviates)")
lines(dens_cp,col="red",lwd=2)
lines(dens_ru,col="blue",lwd=2)
#
# example 4: plot pdf and cdf
xax=c(1:200)/20
np=10000
pp=(c(1:np)-0.5)/np
pdfs_qq=qgev_cp(x,p=pp,pdf=TRUE,rust=TRUE)
pdfs_dd=dgpd_k1_cp(x,xax,rust=TRUE)
cdfs=pgpd_k1_cp(x,xax,rust=TRUE)
#
pdf_mle=pdfs_dd$ml_pdf
pdf_cp=pdfs_dd$cp_pdf
pdf_ru=pdfs_dd$ru_pdf
ymax=max(pdf_mle,pdf_cp)
plot(xax,pdf_mle,"l",ylim=c(0,ymax),
	sub="(from dgpd_k1_cp, qgpd_k1_cp)",
	main="PDF")
lines(xax,pdf_ru,col="blue",lwd=2)
lines(pdfs_qq$cp_quantiles,pdfs_qq$cp_pdf,col="orange")
legend(x=3,y=ymin+0.25*(ymax-ymin),c("maxlik","cp rust","cp from dmgs q"),
	col=c("black","blue","orange"),pch=1,cex=0.7)
#
plot(xax,pmax(exp(-10),pdfs_dd$ml_pdf),"l",log="y",
	sub="(from dgpd_k1_cp, qgpd_k1_cp)",
	main="Log-PDF")
lines(xax,pmax(exp(-10),pdfs_dd$ru_pdf),col="blue",lwd=2)
lines(pdfs_qq$cp_quantiles,pmax(exp(-10),pdfs_qq$cp_pdf),col="orange")
legend(x=4,y=0.99*max(pdf_mle),
	c("maxlik","cp rust","cp from dmgs q"),
	col=c("black","blue","orange"),pch=1,cex=0.7)
#
plot(xax,cdfs$ml_cdf,"l",
	sub="(from pgpd_k1_cp)",
	main="CDF")
lines(xax,cdfs$ru_cdf,col="blue",lwd=2)
lines(pdfs_qq$cp_quantiles,pp,col="orange")
