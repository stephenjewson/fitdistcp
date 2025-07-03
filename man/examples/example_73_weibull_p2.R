#
# example 1
x=fitdistcp::d73weibull_p2_example_data_v1_x
tt=fitdistcp::d73weibull_p2_example_data_v1_t
p=c(1:9)/10
n0=10
q=qweibull_p2_cp(x,tt,n0=n0,p=p,rust=TRUE,nrust=1000)
xmin=min(q$ml_quantiles,q$cp_quantiles);
xmax=max(q$ml_quantiles,q$cp_quantiles);
plot(q$ml_quantiles,p,xlab="quantile estimates",xlim=c(xmin,xmax),
	sub="(from qweibull_p2_cp)",
	main="Weibull w/ p2: quantile estimates");
points(q$cp_quantiles,p,col="red",lwd=2)
points(q$ru_quantiles,p,col="blue")
