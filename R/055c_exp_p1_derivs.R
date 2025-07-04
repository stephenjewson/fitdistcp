######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
exp_p1_fd=function (x, t, v1, v2) 
{
    .e2 <- exp(-(t * v2 + v1))
    .e3 <- x * .e2
    .e5 <- exp(-.e3)
    .e6 <- .e3 - 1
    c(v1 = .e2 * .e5 * .e6, v2 = t * .e2 * .e5 * .e6)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
exp_p1_fdd=function (x, t, v1, v2) 
{
    .e2 <- exp(-(t * v2 + v1))
    .e3 <- x * .e2
    .e6 <- (.e3 - 1)^2 - .e3
    .e7 <- exp(-.e3)
    .e10 <- t * .e6 * .e2 * .e7
    c(v1 = c(v1 = .e6 * .e2 * .e7, v2 = .e10), v2 = c(v1 = .e10, 
        v2 = t^2 * .e6 * .e2 * .e7))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
exp_p1_pd=function (x, t, v1, v2) 
{
    .e2 <- exp(-(t * v2 + v1))
    .e3 <- x * .e2
    .e5 <- exp(-.e3)
    c(v1 = -(.e3 * .e5), v2 = -(t * x * .e2 * .e5))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
exp_p1_pdd=function (x, t, v1, v2) 
{
    .e2 <- exp(-(t * v2 + v1))
    .e3 <- x * .e2
    .e5 <- exp(-.e3)
    .e6 <- .e3 - 1
    .e7 <- -(t * x * .e2 * .e5 * .e6)
    c(v1 = c(v1 = -(.e3 * .e5 * .e6), v2 = .e7), v2 = c(v1 = .e7, 
        v2 = -(t^2 * x * .e2 * .e5 * .e6)))
}
######################################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
exp_p1_logfdd=function (x, t, v1, v2) 
{
    .e2 <- exp(-(t * v2 + v1))
    .e3 <- -(t * x * .e2)
    c(v1 = c(v1 = -(x * .e2), v2 = .e3), v2 = c(v1 = .e3, v2 = -(t^2 * 
        x * .e2)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
exp_p1_logfddd=function (x, t, v1, v2) 
{
    .e2 <- exp(-(t * v2 + v1))
    .e4 <- t * x * .e2
    .e7 <- t^2 * x * .e2
    .e8 <- c(v1 = .e4, v2 = .e7)
    c(v1 = c(v1 = c(v1 = x * .e2, v2 = .e4), v2 = .e8), v2 = c(v1 = .e8, 
        v2 = c(v1 = .e7, v2 = t^3 * x * .e2)))
}
############################################################
#' The first derivative of the density
#' @returns Vector
#' @inheritParams manf
exp_p1_f1fa=function(x,t,v1,v2){
	vf=Vectorize(exp_p1_fd,"x")
	f1=vf(x,t,v1,v2)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @returns Matrix
#' @inheritParams manf
exp_p1_f2fa=function(x,t,v1,v2){
	nx=length(x)
	vf=Vectorize(exp_p1_fdd,"x")
	temp1=vf(x,t,v1,v2)
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}
############################################################
#' The first derivative of the cdf
#' @returns Vector
#' @inheritParams manf
exp_p1_p1fa=function(x,t,v1,v2){
	vf=Vectorize(exp_p1_pd,"x")
	p1=vf(x,t,v1,v2)
	return(p1)
}
############################################################
#' The second derivative of the cdf
#' @returns Matrix
#' @inheritParams manf
exp_p1_p2fa=function(x,t,v1,v2){
	nx=length(x)
	p2=array(0,c(2,2,nx))
	vf=Vectorize(exp_p1_pdd,"x")
	temp1=vf(x,t,v1,v2)
	p2=deriv_copyfdd(temp1,nx,dim=2)
	return(p2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @returns Vector
#' @inheritParams manf
exp_p1_mu1fa=function(alpha,t,v1,v2){
	x=qexp((1-alpha),rate=exp(-v1-v2*t))
	vf=Vectorize(exp_p1_pd,"x")
	mu1=-vf(x,t,v1,v2)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @returns Matrix
#' @inheritParams manf
exp_p1_mu2fa=function(alpha,t,v1,v2){
	x=qexp((1-alpha),rate=exp(-v1-v2*t))
	nalpha=length(alpha)
	vf=Vectorize(exp_p1_pdd,"x")
	temp1=vf(x,t,v1,v2)
	mu2=-deriv_copyfdd(temp1,nalpha,dim=2)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
exp_p1_ldda=function(x,t,v1,v2){
	nx=length(x)
	ldd=matrix(0,2,2)
	vf=Vectorize(exp_p1_logfdd,"x")
	temp1=vf(x,t,v1,v2)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
exp_p1_lddda=function(x,t,v1,v2){
	nx=length(x)
	lddd=array(0,c(2,2,2))
	vf=Vectorize(exp_p1_logfddd,"x")
	temp1=vf(x,t,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}
