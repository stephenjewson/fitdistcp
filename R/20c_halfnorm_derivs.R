######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
halfnorm_fd=function (x, v1) 
{
    .e3 <- v1^2 * x^2/pi
    (2 - 4 * .e3) * exp(-.e3)/pi
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
halfnorm_fdd=function (x, v1) 
{
    .e1 <- x^2
    .e4 <- v1^2 * .e1/pi
    -(v1 * .e1 * (2 * (2 - 4 * .e4) + 8) * exp(-.e4)/pi^2)
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
halfnorm_logfdd=function (x, v1) 
-(1/v1^2 + 2 * (x^2/pi))
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
halfnorm_logfddd=function (x, v1) 
2/v1^3
############################################################
#' The first derivative of the density
#' @returns Vector
#' @inheritParams manf
halfnorm_f1fa=function(x,v1){
	nx=length(x)
	f1=matrix(0,1,nx)
	vf=Vectorize(halfnorm_fd)
	f1[1,]=vf(x,v1)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @returns Matrix
#' @inheritParams manf
halfnorm_f2fa=function(x,v1){
	nx=length(x)
	f2=array(0,c(1,1,nx))
	vf=Vectorize(halfnorm_fdd)
	f2[1,1,]=vf(x,v1)
	return(f2)
}
############################################################
############################################################
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
halfnorm_ldda=function(x,v1){
	nx=length(x)
	ldd=matrix(0,1,1)
	vf=Vectorize(halfnorm_logfdd)
	ldd[1,1]=sum(vf(x,v1))/nx
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
halfnorm_lddda=function(x,v1){
	nx=length(x)
	lddd=array(0,c(1,1,1))
	vf=Vectorize(halfnorm_logfddd)
	lddd[1,1,1]=sum(vf(x,v1))/nx
	return(lddd)
}
