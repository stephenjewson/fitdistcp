#
# example 1
x=fitdistcp::d56pareto_p1k2_example_data_v1_x
tt=fitdistcp::d56pareto_p1k2_example_data_v1_t
p=c(1:9)/10
n0=10
q=qpareto_p1k2_cp(x,tt,n0=n0,p=p,rust=TRUE,nrust=1000)
xmin=min(q$ml_quantiles,q$cp_quantiles);
xmax=max(q$ml_quantiles,q$cp_quantiles);
plot(q$ml_quantiles,p,xlab="quantile estimates",xlim=c(xmin,xmax),
	sub="(from qpareto_p1k2_cp)",
	main="Pareto w/ p2: quantile estimates");
points(q$cp_quantiles,p,col="red",lwd=2)
points(q$ru_quantiles,p,col="blue")
