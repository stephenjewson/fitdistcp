#
# example 1
x=fitdistcp::d060norm_p1_example_data_v1_x
tt=fitdistcp::d060norm_p1_example_data_v1_t
p=c(1:9)/10
n0=10
q=qnorm_p12_cp(x,t1=tt,t2=tt,n10=n0,n20=n0,p=p,rust=TRUE,nrust=1000)
xmin=min(q$ml_quantiles,q$cp_quantiles);
xmax=max(q$ml_quantiles,q$cp_quantiles);
plot(q$ml_quantiles,p,xlab="quantile estimates",xlim=c(xmin,xmax),
	sub="(from qnorm_p12_cp)",
	main="Normal w/ p1: quantile estimates");
points(q$cp_quantiles,p,col="red",lwd=2)
points(q$ru_quantiles,p,col="blue")
