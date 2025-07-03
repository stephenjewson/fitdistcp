######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
logis_fd=function (x, v1, v2) 
{
    .e1 <- x - v1
    .e3 <- exp(-(.e1/v2))
    .e4 <- 1 + .e3
    .e6 <- 1 - 2 * (.e3/.e4)
    .e8 <- v2^2 * .e4^2
    c(v1 = .e6 * .e3/.e8, v2 = (.e6 * .e1/v2 - 1) * .e3/.e8)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
logis_fdd=function (x, v1, v2) 
{
    .e1 <- x - v1
    .e3 <- exp(-(.e1/v2))
    .e4 <- 1 + .e3
    .e5 <- .e3/.e4
    .e6 <- .e4^2
    .e7 <- 1 - 2 * .e5
    .e9 <- .e7 * .e1/v2
    .e10 <- (v2^2 * .e6)^2
    .e11 <- 1 - .e5
    .e12 <- .e9 - ((2 * (.e11 * .e1/v2) - 2) * .e3/.e4 + 2)
    .e13 <- .e9 - 1
    .e16 <- 1 - (2 + 2 * .e11) * .e3/.e4
    .e18 <- 2 * (.e3 * .e1) + 2 * (v2 * .e4)
    .e20 <- v2^3 * .e6
    .e22 <- v2^4 * .e6
    c(v1 = c(v1 = (.e16/.e20 - 2 * (v2 * .e7 * .e4 * .e3/.e10)) * 
        .e3, v2 = (.e12/.e20 - 2 * (v2 * .e13 * .e4 * .e3/.e10)) * 
        .e3), v2 = c(v1 = (.e16 * .e1/.e22 - .e7 * .e4 * .e18/.e10) * 
        .e3, v2 = (.e12 * .e1/.e22 - .e13 * .e4 * .e18/.e10) * 
        .e3))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
logis_pd=function (x, v1, v2) 
{
    .e1 <- x - v1
    .e3 <- exp(-(.e1/v2))
    .e4 <- (1 + .e3)^2
    c(v1 = -(.e3/(v2 * .e4)), v2 = -(.e3 * .e1/(v2^2 * .e4)))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
logis_pdd=function (x, v1, v2) 
{
    .e1 <- x - v1
    .e2 <- .e1/v2
    .e4 <- exp(-.e2)
    .e5 <- 1 + .e4
    .e6 <- .e5^2
    .e8 <- v2^2 * .e6
    .e9 <- (v2 * .e6)^2
    .e10 <- .e8^2
    .e11 <- v2 * .e5
    c(v1 = c(v1 = -((1/.e8 - 2 * (.e5 * .e4/.e9)) * .e4), v2 = -(((.e2 - 
        1)/.e8 - 2 * (.e11 * .e4 * .e1/.e10)) * .e4)), v2 = c(v1 = -((.e1/(v2^3 * 
        .e6) - ((1 + 2 * .e2) * .e4 + 1) * .e5/.e9) * .e4), v2 = -((.e1/(v2^4 * 
        .e6) - .e5 * (2 * (.e4 * .e1) + 2 * .e11)/.e10) * .e4 * 
        .e1)))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
logis_logfdd=function (x, v1, v2) 
{
    .e1 <- x - v1
    .e3 <- exp(-(.e1/v2))
    .e4 <- 1 + .e3
    .e5 <- 1 - .e3/.e4
    .e6 <- v2^2
    .e9 <- 2 * (.e5 * .e1/v2)
    .e10 <- -(((.e9 - 2) * .e3/.e4 + 1)/.e6)
    c(v1 = c(v1 = -(2 * (.e5 * .e3/(.e6 * .e4))), v2 = .e10), 
        v2 = c(v1 = .e10, v2 = -((((.e9 - 4) * .e3/.e4 + 2) * 
            .e1/v2 - 1)/.e6)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
logis_logfddd=function (x, v1, v2) 
{
    .e1 <- x - v1
    .e2 <- .e1/v2
    .e4 <- exp(-.e2)
    .e5 <- 1 + .e4
    .e6 <- 1 - .e4/.e5
    .e7 <- .e6 * .e1
    .e9 <- 2 * (.e7/v2)
    .e10 <- v2 * .e5
    .e11 <- .e4 * .e1
    .e12 <- v2^3
    .e14 <- (.e9 - 2) * .e4/.e5
    .e16 <- 2 * .e2
    .e17 <- 2 * (1 + .e11/.e10)
    .e18 <- (.e6 * (.e16 - .e17) - (.e14 + 2)) * .e4
    .e20 <- (.e9 - 4) * .e4/.e5
    .e21 <- .e12 * .e5
    .e22 <- -((.e18 * .e1/.e10 - 2 * (.e14 + 1))/.e12)
    .e23 <- -(.e18/.e21)
    .e26 <- ((.e6 * (.e16 - (2 + .e17)) - (.e20 + 4)) * .e1/v2 + 
        4) * .e4/.e5 - 2
    .e27 <- (v2^2 * .e5)^2
    c(v1 = c(v1 = c(v1 = -(2 * ((.e6/.e21 - v2 * .e4/.e27) * 
        .e6 * .e4)), v2 = .e23), v2 = c(v1 = .e23, v2 = -(.e26/.e12))), 
        v2 = c(v1 = c(v1 = -(2 * ((.e7/(v2^4 * .e5) - (2 * .e10 + 
            .e11)/.e27) * .e6 * .e4)), v2 = .e22), v2 = c(v1 = .e22, 
            v2 = -((.e26 * .e1/v2 - 2 * ((.e20 + 2) * .e1/v2 - 
                1))/.e12))))
}
############################################################
#' The first derivative of the density
#' @returns Vector
#' @inheritParams manf
logis_f1fa=function(x,v1,v2){
	vf=Vectorize(logis_fd)
	f1=vf(x,v1,v2)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @returns Matrix
#' @inheritParams manf
logis_f2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(logis_fdd)
	temp1=vf(x,v1,v2)
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}
############################################################
#' The first derivative of the cdf
#' @returns Vector
#' @inheritParams manf
logis_p1fa=function(x,v1,v2){
	vf=Vectorize(logis_pd)
	p1=vf(x,v1,v2)
	return(p1)
}
############################################################
#' The second derivative of the cdf
#' @returns Vector
#' @inheritParams manf
logis_p2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(logis_pdd)
	temp1=vf(x,v1,v2)
	p2=deriv_copyfdd(temp1,nx,dim=2)
	return(p2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @returns Vector
#' @inheritParams manf
logis_mu1fa=function(alpha,v1,v2){
	x=qlogis((1-alpha),location=v1,scale=v2)
	vf=Vectorize(logis_pd)
	mu1=-vf(x,v1,v2)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @returns Matrix
#' @inheritParams manf
logis_mu2fa=function(alpha,v1,v2){
	x=qlogis((1-alpha),location=v1,scale=v2)
	nx=length(x)
	vf=Vectorize(logis_pdd)
	temp1=vf(x,v1,v2)
	mu2=-deriv_copyfdd(temp1,nx,dim=2)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
logis_ldda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(logis_logfdd)
	temp1=vf(x,v1,v2)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
logis_lddda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(logis_logfddd)
	temp1=vf(x,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}
