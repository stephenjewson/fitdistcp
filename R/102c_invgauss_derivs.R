######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
invgauss_fd=function (x, v1, v2) 
{
    .e1 <- v1^2
    .e2 <- x - v1
    .e3 <- .e1 * x
    .e4 <- .e2^2
    .e6 <- 2 * .e3
    .e7 <- sqrt(v2/(2 * (pi * x^3)))
    .e9 <- exp(-(v2 * .e4/.e6))
    c(v1 = v2 * (1/.e3 + 4 * (v1 * x * .e2/.e6^2)) * .e9 * .e7 * 
        .e2, v2 = (1/(4 * (pi * x^2 * .e7)) - .e7 * .e4/(2 * 
        .e1)) * .e9/x)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
invgauss_fdd=function (x, v1, v2) 
{
    .e1 <- v1^2
    .e2 <- x - v1
    .e3 <- .e1 * x
    .e4 <- 2 * .e3
    .e5 <- pi * x^3
    .e7 <- sqrt(v2/(2 * .e5))
    .e8 <- .e2^2
    .e9 <- .e4^2
    .e10 <- x^2
    .e12 <- v2 * .e8/.e4
    .e14 <- 1/.e3
    .e15 <- 2 * .e1
    .e16 <- exp(-.e12)
    .e17 <- .e14 + 4 * (v1 * x * .e2/.e9)
    .e18 <- 4 * (pi * .e10 * .e7)
    .e20 <- 1/.e18 - .e7 * .e8/.e15
    c(v1 = c(v1 = v2 * ((v2 * .e17^2 * .e2 + x * (4 * ((x - v1 * 
        (16 * (v1^3 * .e10 * .e2/.e9) + 2))/.e9) - v1 * (2/.e3^2 + 
        4/.e9))) * .e2 - .e14) * .e16 * .e7, v2 = ((1/.e1 + 4 * 
        (v1 * .e2/.e15^2)) * .e7 + v2 * .e20 * .e17) * .e16 * 
        .e2/x), v2 = c(v1 = ((1 - .e12) * .e7 + v2/(4 * (.e5 * 
        .e7))) * .e17 * .e16 * .e2, v2 = -(((.e8/(8 * (pi * .e1 * 
        .e10)) + 1/.e18^2)/.e7 + .e20 * .e8/.e15) * .e16/.e10)))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
invgauss_logfdd=function (x, v1, v2) 
{
    .e2 <- v1^2 * x
    .e3 <- x - v1
    .e4 <- (2 * .e2)^2
    .e5 <- 1/.e2
    .e6 <- (.e5 + 4 * (v1 * x * .e3/.e4)) * .e3
    c(v1 = c(v1 = v2 * (x * (4 * ((x - v1 * (16 * (v1^3 * x^2 * 
        .e3/.e4) + 2))/.e4) - v1 * (2/.e2^2 + 4/.e4)) * .e3 - 
        .e5), v2 = .e6), v2 = c(v1 = .e6, v2 = -(0.5/v2^2)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
invgauss_logfddd=function (x, v1, v2) 
{
    .e2 <- v1^2 * x
    .e3 <- 2 * .e2
    .e4 <- .e3^2
    .e5 <- x - v1
    .e6 <- x^2
    .e8 <- v1^3 * .e6
    .e9 <- 16 * (.e8 * .e5/.e4)
    .e10 <- .e2^2
    .e13 <- x - v1 * (.e9 + 2)
    .e15 <- 2/.e10
    .e16 <- 4 * (.e13/.e4)
    .e17 <- 4/.e4
    .e23 <- x * (.e16 - v1 * (.e15 + .e17)) * .e5 - 1/.e2
    .e24 <- c(v1 = .e23, v2 = 0)
    c(v1 = c(v1 = c(v1 = v2 * x * ((v1^4 * .e6 * (64/.e3^4 + 
        8/.e2^4) - ((4 + 4 * (2 + .e8 * (16 * (3 * .e5 - v1 * 
        (1 + .e9)) + 16 * .e13 + 16 * .e5)/.e4))/.e4 + .e15)) * 
        .e5 + v1 * (.e17 + 4/.e10) - .e16), v2 = .e23), v2 = .e24), 
        v2 = c(v1 = .e24, v2 = c(v1 = 0, v2 = 1/v2^3)))
}
############################################################
#' The first derivative of the density
#' @inheritParams manf
invgauss_f1fa=function(x,v1,v2){
	vf=Vectorize(invgauss_fd)
	f1=vf(x,v1,v2)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @inheritParams manf
invgauss_f2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(invgauss_fdd)
	temp1=vf(x,v1,v2)
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}
############################################################
############################################################
############################################################
#' The second derivative of the normalized log-likelihood
#' @inheritParams manf
invgauss_ldda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(invgauss_logfdd)
	temp1=vf(x,v1,v2)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @inheritParams manf
invgauss_lddda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(invgauss_logfddd)
	temp1=vf(x,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}
