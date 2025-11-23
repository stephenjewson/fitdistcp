# example 1
x=fitdistcp::d150gev_p1_example_data_v1_x
t1=fitdistcp::d150gev_p1_example_data_v1_t
t2=sample(t1)
t=cbind(t1,t2)
p=c(1:9)/10
n0=c(10,10)
q=qgev_p1n_cp(x=x,t=t,n0=n0,t0=NA,p=p,rust=TRUE,nrust=1000)
xmin=min(q$ml_quantiles,q$cp_quantiles);
xmax=max(q$ml_quantiles,q$cp_quantiles);
plot(q$ml_quantiles,p,xlab="quantile estimates",xlim=c(xmin,xmax),
	sub="(from qgev_p1n_cp)",
	main="GEVD w/ p1: quantile estimates");
points(q$cp_quantiles,p,col="red",lwd=2)
points(q$ru_quantiles,p,col="blue",lwd=2)
cat(" ml_params=",q$ml_params,"\n")
