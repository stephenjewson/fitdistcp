######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
cauchy_fd=function (x, v1, v2) 
{
    .e1 <- x - v1
    .e2 <- .e1^2
    .e4 <- .e2/v2^2 + 1
    c(v1 = 2 * (.e1/(pi * v2^3 * .e4^2)), v2 = (2 * (.e2/(pi * 
        v2^4 * .e4)) - pi/(pi * v2)^2)/.e4)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
cauchy_fdd=function (x, v1, v2) 
{
    .e1 <- x - v1
    .e2 <- .e1^2
    .e3 <- v2^2
    .e5 <- .e2/.e3 + 1
    .e6 <- pi * v2
    .e8 <- pi * v2^4 * .e5
    .e9 <- v2^3
    .e12 <- pi * .e9 * .e5^2
    .e13 <- .e3 * .e5
    .e15 <- .e12^2
    .e16 <- .e8^2
    .e19 <- 2 * (.e2/.e8) - pi/.e6^2
    c(v1 = c(v1 = 2 * (4 * (.e6 * .e5 * .e2/.e15) - 1/.e12), 
        v2 = (2 * (.e19/.e13) + 2 * (2 * (pi * .e3 * .e2/.e16) - 
            2/.e8)) * .e1/.e5), v2 = c(v1 = -(2 * (pi * .e5 * 
        (3 * .e13 - 4 * .e2) * .e1/.e15)), v2 = (2 * (.e19 * 
        .e2/(.e9 * .e5)) + .e6 * (2 * (pi^2/.e6^4) - 2 * ((4 * 
        .e13 - 2 * .e2) * .e2/.e16)))/.e5))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
cauchy_logfdd=function (x, v1, v2) 
{
    .e1 <- x - v1
    .e3 <- (.e1/v2)^2 + 1
    .e5 <- v2^2 * .e3
    .e6 <- .e1^2
    .e7 <- .e5^2
    .e10 <- 2 * (.e6/.e7)
    .e13 <- 2 * (v2 * .e3) - 2 * (.e6/v2)
    c(v1 = c(v1 = 2 * (.e10 - 1/.e5), v2 = 2 * ((.e10 - 2/.e5) * 
        .e1/v2)), v2 = c(v1 = -(2 * (.e13 * .e1/.e7)), v2 = -(((2 * 
        (.e6/.e5) - 1)/v2 + 2 * (.e13 * .e6/.e7))/v2)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
cauchy_logfddd=function (x, v1, v2) 
{
    .e1 <- x - v1
    .e3 <- (.e1/v2)^2 + 1
    .e4 <- v2^2
    .e5 <- .e4 * .e3
    .e6 <- .e1^2
    .e7 <- .e5^2
    .e11 <- 2 * (v2 * .e3) - 2 * (.e6/v2)
    .e12 <- 4 * (.e5 * .e6/.e7)
    .e14 <- 2 * (.e6/.e4)
    .e15 <- 2/.e5
    .e16 <- .e12 - 2
    .e18 <- (2 * (.e6/.e7) - .e15)/v2
    .e22 <- 2 * (.e3 - .e14) + .e14 - 2 * (.e5 * .e11^2/.e7)
    .e23 <- 2 * .e16
    c(v1 = c(v1 = c(v1 = 2 * ((.e23 - 2) * .e1/.e7), v2 = 2 * 
        (((.e23 - 6) * .e6/.e7 + .e15)/v2)), v2 = c(v1 = -(2 * 
        (.e11 * (.e12 - 1)/.e7)), v2 = -((2 * .e18 + 2 * (.e11 * 
        .e16/.e7)) * .e1/v2))), v2 = c(v1 = c(v1 = 2 * ((1 - 
        .e12) * .e11/.e7), v2 = 2 * (((2 - .e12) * .e11/.e7 - 
        .e18) * .e1/v2)), v2 = c(v1 = -(2 * (.e22 * .e1/.e7)), 
        v2 = -((2 * (.e22 * .e6/.e7) - 2 * (((2 * (.e6/.e5) - 
            1)/v2 + 2 * (.e11 * .e6/.e7))/v2))/v2))))
}
############################################################
#' The first derivative of the density
#' @returns Vector
#' @inheritParams manf
cauchy_f1fa=function(x,v1,v2){
	vf=Vectorize(cauchy_fd)
	f1=vf(x,v1,v2)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @returns Matrix
#' @inheritParams manf
cauchy_f2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(cauchy_fdd)
	temp1=vf(x,v1,v2)
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
cauchy_ldda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(cauchy_logfdd)
	temp1=vf(x,v1,v2)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
cauchy_lddda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(cauchy_logfddd)
	temp1=vf(x,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}
