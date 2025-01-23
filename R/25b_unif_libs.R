#' Predictive PDFs
#' @inheritParams manf
dunif_formula=function(x,y){
	nx=length(x)
	ny=length(y)
	ml=numeric(ny)
	rh=numeric(ny)
	a0=max(x)
	b0=min(x)
	range=a0-b0
	for (iy in 1:ny){
		aa=max(y[iy],a0)
		bb=min(y[iy],b0)
		amb=aa-bb
		a0b=a0-b0
		ml[iy]=ifelse((y[iy]>a0)|(y[iy]<b0),0,1/range)
		t1=a0b^(nx-1)
		t2=amb^nx
		t3=(nx-1)/(nx+1)
		rh[iy]=(t1/t2)*t3
	}
	return(list(ml_pdf=ml,rh_pdf=rh))
}
#' Predictive CDFs
#' @inheritParams manf
punif_formula=function(x,y){
	nx=length(x)
	ny=length(y)
	ml=numeric(ny)
	rh=numeric(ny)
	a0=max(x)
	b0=min(x)
	range=a0-b0
	fact1=(a0-b0)^(nx-1)
	fact2=(nx-1)/(nx+1)
	fact=fact1*fact2
	for (iy in 1:ny){
		aa=max(y[iy],a0)
		bb=min(y[iy],b0)
		amb=aa-bb
		a0b=a0-b0
		if(y[iy]<=b0){

# left section
			ml[iy]=0
			t1=1/(nx-1)
			t2=(a0-y[iy])^(nx-1)
			rh[iy]=fact*t1/t2
		} else if ((y[iy]>b0)&(y[iy]<=a0)){

# middle section

			ml[iy]=(y[iy]-b0)/(a0-b0)
# -the first bit, which is now fixed
			t1=1/(nx-1)
			t2=(a0-b0)^(nx-1)
# -the second bit, which varies
			t3=y[iy]-b0
			t4=(a0-b0)^nx
			rh[iy]=fact*(t1/t2+t3/t4)
		} else{

# right section
			ml[iy]=1
# -the first bit
			t2=(a0-b0)^(nx-1)
# -the second bit, which is now fixed
# -the third bit, which now varies
			t3=(y[iy]-b0)^(nx-1)
			rh[iy]=fact*((nx+1)/t2-1/t3)/(nx-1)
		} #end of if
	} #end of y loop
	return(list(ml_cdf=ml,rh_cdf=rh))
}
#' Predictive Quantiles
#' @inheritParams manf
qunif_formula=function(x,p){
	nx=length(x)
	np=length(p)
	ml=numeric(np)
	rh=numeric(np)
	a0=max(x)
	b0=min(x)
	range=a0-b0
	fact1=(a0-b0)^(nx-1)
	fact2=(nx-1)/(nx+1)
	fact=fact1*fact2
	ifact=1/fact
	t1=1/(nx-1)
	t2=(a0-b0)^(nx-1)
	t4=(a0-b0)^(nx-1)
	pp1=fact*t1/t2
	pp2=fact*(t1/t2+1/t4)
#	cat("pp1,pp2=",pp1,pp2,"\n")
#	cat("a0,b0=",a0,b0,"\n")
	for (ip in 1:np){
		if(p[ip]<=pp1){

# left section
			ml[ip]=b0+p[ip]*(a0-b0)
			rh[ip]=a0-1/(ifact*p[ip]*(nx-1))^(1/(nx-1))

		} else if ((p[ip]>pp1)&(p[ip]<=pp2)){
# middle section
			ml[ip]=b0+p[ip]*(a0-b0)
			t1=(a0-b0)^(nx-1)
			t2=(nx-1)*t1
			rh[ip]=b0+(a0-b0)^nx*(ifact*p[ip]-1/t2)

		} else{
# right section
			ml[ip]=b0+p[ip]*(a0-b0)
			t1=(a0-b0)^(nx-1)
			tt=(nx+1)/t1-(nx-1)*(ifact*p[ip])
			rh[ip]=b0+1/(tt^(1/(nx-1)))

		} #end of if
	} #end of y loop
	mean=0.5*(a0+b0)
	return(list(ml_quantiles=ml,rh_quantiles=rh,mean=mean,ml_params=c(a0,b0)))
}
