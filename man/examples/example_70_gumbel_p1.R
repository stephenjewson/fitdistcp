#
# example 1
x=fitdistcp::d70gumbel_p1_example_data_v1_x
tt=fitdistcp::d70gumbel_p1_example_data_v1_t
p=c(1:9)/10
n0=10
q=qgumbel_p1_cp(x,tt,n0=n0,p=p,rust=TRUE,nrust=1000)
xmin=min(q$ml_quantiles,q$cp_quantiles);
xmax=max(q$ml_quantiles,q$cp_quantiles);
plot(q$ml_quantiles,p,xlab="quantile estimates",xlim=c(xmin,xmax),
	sub="(from qgumbel_p1_cp)",
	main="Gumbel w/ p1: quantile estimates");
points(q$cp_quantiles,p,col="red",lwd=2)
points(q$ru_quantiles,p,col="blue")
