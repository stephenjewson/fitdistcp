######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gev_k3_fd=function (x, v1, v2, v3) 
{
    .e1 <- 1/v3
    .e2 <- x - v1
    .e3 <- 1 + v3 * .e2/v2
    .e4 <- 1 + .e1
    .e9 <- exp(-.e3^-.e1)
    .e10 <- v2^2
    .e13 <- v3 * .e4/.e3^(.e1 + 2) - 1/.e3^(2 * .e4)
    c(v1 = .e9 * .e13/.e10, v2 = (.e13 * .e2/v2 - 1/.e3^.e4) * 
        .e9/.e10)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_k3_fdd=function (x, v1, v2, v3) 
{
    .e1 <- 1/v3
    .e2 <- x - v1
    .e3 <- 1 + v3 * .e2/v2
    .e4 <- 1 + .e1
    .e5 <- .e1 + 2
    .e6 <- 2 * .e4
    .e7 <- .e3^.e5
    .e9 <- 1/.e3^.e6
    .e10 <- v3 * .e4
    .e11 <- .e3^.e4
    .e13 <- .e10/.e7 - .e9
    .e18 <- exp(-.e3^-.e1)
    .e19 <- v2^3
    .e22 <- v3 * .e3^(.e4 - 2 * .e5) * .e5 - 2/.e3^(1 + .e6)
    .e25 <- .e13 * .e2/v2 - 1/.e11
    .e34 <- .e9 + v3 * (.e22 * .e2/v2 - (.e3^(.e1 - .e6) + 1/.e7)) * 
        .e4 - .e25/.e11
    .e36 <- .e10 * .e22 - .e13/.e11
    c(v1 = c(v1 = .e18 * .e36/.e19, v2 = .e34 * .e18/.e19), v2 = c(v1 = (.e36 * 
        .e2/v2 - 2 * .e13) * .e18/.e19, v2 = (.e34 * .e2/v2 - 
        2 * .e25) * .e18/.e19))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gev_k3_pd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- 1 + v3 * .e1/v2
    .e3 <- 1/v3
    .e5 <- .e2^(1 + .e3)
    .e6 <- exp(-.e2^-.e3)
    c(v1 = -(.e6/(v2 * .e5)), v2 = -(.e6 * .e1/(v2^2 * .e5)))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_k3_pdd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- 1/v3
    .e3 <- 1 + v3 * .e1/v2
    .e4 <- 1 + .e2
    .e5 <- .e3^.e4
    .e7 <- .e3^.e2
    .e8 <- exp(-.e3^-.e2)
    .e9 <- v2 * .e5
    .e10 <- v2^2
    .e11 <- .e3^(2 * .e4)
    .e12 <- .e10 * .e5
    .e14 <- v3 * .e4 * .e7
    .e15 <- .e9^2
    .e16 <- .e12^2
    .e17 <- .e14 * .e1
    c(v1 = c(v1 = -(.e8 * (.e14/.e15 - 1/(.e10 * .e11))), v2 = -(.e8 * 
        (v2 * v3 * .e4 * .e7 * .e1/.e16 - (.e1/.e9 + 1)/.e12))), 
        v2 = c(v1 = ((.e5 - .e17/v2)/.e15 + .e1/(v2^3 * .e11)) * 
            .e8, v2 = ((2 * .e9 - .e17)/.e16 + .e1/(v2^4 * .e11)) * 
            .e8 * .e1))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_k3_logfdd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- 1 + v3 * .e1/v2
    .e3 <- 1/v3
    .e4 <- 1 + .e3
    .e6 <- 1/.e2^.e3
    .e7 <- v2 * .e2
    .e8 <- v3 * .e4
    .e9 <- .e8 - .e6
    .e10 <- .e7^2
    .e11 <- .e2^(.e3 + 2)
    .e13 <- .e9/.e10 + .e1/(v2^3 * .e11)
    .e14 <- v3 * .e9
    c(v1 = c(v1 = .e14/.e10 - 1/(v2^2 * .e11), v2 = (.e14 * .e1/.e10 - 
        (.e1/(v2 * .e2^.e4) + .e8 - .e6)/.e7)/v2), v2 = c(v1 = -.e13, 
        v2 = -(((.e9 * .e1/.e7 - 1)/v2 + .e13 * .e1)/v2)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
gev_k3_logfddd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- 1 + v3 * .e1/v2
    .e3 <- 1/v3
    .e4 <- 1 + .e3
    .e5 <- v2 * .e2
    .e6 <- .e2^.e4
    .e7 <- .e5^2
    .e8 <- .e2^.e3
    .e9 <- .e3 + 2
    .e10 <- v3 * .e4
    .e11 <- 1/.e8
    .e12 <- .e10 - .e11
    .e13 <- .e2^.e9
    .e14 <- v2 * .e6
    .e15 <- v2^2
    .e17 <- v2^3 * .e13
    .e18 <- .e5 * .e12
    .e19 <- .e15 * .e6
    .e20 <- .e17^2
    .e21 <- .e1/.e14
    .e23 <- .e1/.e19 + 2 * (.e18/.e7)
    .e24 <- (2 * (v2 * v3 * .e2 * .e12/.e7) - 1/.e14)/.e7
    .e26 <- .e21 + .e10 - .e11
    .e27 <- v2 * .e13
    .e30 <- v3 * .e6 * .e9 * .e1
    .e33 <- .e23/.e7 + v2 * (3 * .e27 - .e30) * .e1/.e20
    .e34 <- .e24 + .e15 * v3 * .e6 * .e9 * .e1/.e20
    .e35 <- .e14^2
    .e36 <- (.e15 * .e13)^2
    .e37 <- .e12/.e7
    .e38 <- (v3 * .e12 * .e1/.e7 - .e26/.e5)/v2
    .e40 <- v3 * .e23/.e7
    .e42 <- .e10 * .e8 * .e1
    c(v1 = c(v1 = c(v1 = v3 * (.e24 - .e14 * .e9/.e36), v2 = (v3 * 
        (2/.e8 + v3 * (2 * (.e18 * .e1/.e7) - (2 + 2/v3)) - 2 * 
            .e21)/.e7 - (.e42/.e35 - 2/.e14)/.e5)/v2), v2 = c(v1 = -(.e34 - 
        1/.e17), v2 = -(((.e34 - 2/.e17) * .e1 + .e38 - .e37)/v2))), 
        v2 = c(v1 = c(v1 = (2 * .e27 - .e30)/.e36 - .e40, v2 = ((((.e6 - 
            .e42/v2)/.e35 + 1/.e19)/.e5 - .e40) * .e1 + .e26/.e7 - 
            .e38)/v2), v2 = c(v1 = .e33, v2 = (.e33 * .e1 + 2 * 
            (((.e12 * .e1/.e5 - 1)/v2 + (.e37 + .e1/.e17) * .e1)/v2))/v2)))
}
############################################################
#' The first derivative of the density
#' @returns Vector
#' @inheritParams manf
gev_k3_f1fa=function(x,v1,v2,kshape){
	kshape=movexiawayfromzero(kshape)
	vf=Vectorize(gev_k3_fd)
	f1=vf(x,v1,v2,kshape)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @returns Matrix
#' @inheritParams manf
gev_k3_f2fa=function(x,v1,v2,kshape){
	nx=length(x)

	kshape=movexiawayfromzero(kshape)

	vf=Vectorize(gev_k3_fdd)
	temp1=vf(x,v1,v2,kshape)
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @returns Vector
#' @inheritParams manf
gev_k3_mu1fa=function(alpha,v1,v2,kshape){
	x=extraDistr::qgev((1-alpha),mu=v1,sigma=v2,xi=kshape)
	kshape=movexiawayfromzero(kshape)
	vf=Vectorize(gev_k3_pd)
	mu1=-vf(x,v1,v2,kshape)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @returns Matrix
#' @inheritParams manf
gev_k3_mu2fa=function(alpha,v1,v2,kshape){
	x=extraDistr::qgev((1-alpha),mu=v1,sigma=v2,xi=kshape)
	nx=length(x)

	kshape=movexiawayfromzero(kshape)

	vf=Vectorize(gev_k3_pdd)
	temp1=vf(x,v1,v2,kshape)
	mu2=-deriv_copyfdd(temp1,nx,dim=2)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
gev_k3_ldda=function(x,v1,v2,kshape){
	nx=length(x)

	kshape=movexiawayfromzero(kshape)

	vf=Vectorize(gev_k3_logfdd)
	temp1=vf(x,v1,v2,kshape)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
gev_k3_lddda=function(x,v1,v2,kshape){
	nx=length(x)
	vf=Vectorize(gev_k3_logfddd)

	kshape=movexiawayfromzero(kshape)

	temp1=vf(x,v1,v2,kshape)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}
