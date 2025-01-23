######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
exp_fd=function (x, v1) 
{
    .e1 <- v1 * x
    (1 - .e1) * exp(-.e1)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
exp_fdd=function (x, v1) 
{
    .e1 <- v1 * x
    -(x * (2 - .e1) * exp(-.e1))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
exp_pd=function (x, v1) 
x * exp(-(v1 * x))
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
exp_pdd=function (x, v1) 
-(x^2 * exp(-(v1 * x)))
######################################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
exp_logfdd=function (x, v1) 
-(1/v1^2)
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
exp_logfddd=function (x, v1) 
2/v1^3
############################################################
#' The first derivative of the density
#' @inheritParams manf
exp_f1fa=function(x,v1){
	nx=length(x)
	f1=matrix(0,1,nx)
	vf=Vectorize(exp_fd)
	f1[1,]=vf(x,v1)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @inheritParams manf
exp_f2fa=function(x,v1){
	nx=length(x)
	f2=array(0,c(1,1,nx))
	vf=Vectorize(exp_fdd)
	f2[1,1,]=vf(x,v1)
	return(f2)
}
############################################################
#' The first derivative of the cdf
#' @inheritParams manf
exp_p1fa=function(x,v1){
	nx=length(x)
	p1=matrix(0,1,nx)
	vf=Vectorize(exp_pd)
	p1[1,]=vf(x,v1)
	return(p1)
}
############################################################
#' The second derivative of the cdf
#' @inheritParams manf
exp_p2fa=function(x,v1){
	nx=length(x)
	p2=array(0,c(1,1,nx))
	vf=Vectorize(exp_pdd)
	p2[1,1,]=vf(x,v1)
	return(p2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @inheritParams manf
exp_ldda=function(x,v1){
	nx=length(x)
	ldd=matrix(0,1,1)
	vf=Vectorize(exp_logfdd)
	ldd[1,1]=sum(vf(x,v1))/nx
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @inheritParams manf
exp_lddda=function(x,v1){
	nx=length(x)
	lddd=array(0,c(1,1,1))
	vf=Vectorize(exp_logfddd)
	lddd[1,1,1]=sum(vf(x,v1))/nx
	return(lddd)
}
