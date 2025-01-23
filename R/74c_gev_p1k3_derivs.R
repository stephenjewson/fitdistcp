######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
gev_p1k3_fd=function (x, t, v1, v2, v3, v4) 
{
    .e1 <- 1/v4
    .e4 <- x - (t * v2 + v1)
    .e5 <- 1 + v4 * .e4/v3
    .e6 <- 1 + .e1
    .e11 <- exp(-.e5^-.e1)
    .e12 <- v3^2
    .e15 <- v4 * .e6/.e5^(.e1 + 2) - 1/.e5^(2 * .e6)
    c(v1 = .e11 * .e15/.e12, v2 = t * .e11 * .e15/.e12, v3 = (.e15 * 
        .e4/v3 - 1/.e5^.e6) * .e11/.e12)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
gev_p1k3_fdd=function (x, t, v1, v2, v3, v4) 
{
    .e1 <- 1/v4
    .e4 <- x - (t * v2 + v1)
    .e5 <- 1 + v4 * .e4/v3
    .e6 <- 1 + .e1
    .e7 <- .e1 + 2
    .e8 <- 2 * .e6
    .e9 <- v4 * .e6
    .e10 <- .e5^.e7
    .e12 <- 1/.e5^.e8
    .e13 <- .e5^.e6
    .e15 <- .e9/.e10 - .e12
    .e20 <- exp(-.e5^-.e1)
    .e21 <- v3^3
    .e24 <- v4 * .e5^(.e6 - 2 * .e7) * .e7 - 2/.e5^(1 + .e8)
    .e27 <- .e9 * .e24 - .e15/.e13
    .e30 <- .e15 * .e4/v3 - 1/.e13
    .e38 <- .e12 + v4 * (.e24 * .e4/v3 - (.e5^(.e1 - .e8) + 1/.e10)) * 
        .e6 - .e30/.e13
    .e41 <- .e27 * .e4/v3 - 2 * .e15
    .e44 <- t * .e20 * .e27/.e21
    c(v1 = c(v1 = .e20 * .e27/.e21, v2 = .e44, v3 = .e38 * .e20/.e21), 
        v2 = c(v1 = .e44, v2 = t^2 * .e20 * .e27/.e21, v3 = t * 
            .e38 * .e20/.e21), v3 = c(v1 = .e41 * .e20/.e21, 
            v2 = t * .e41 * .e20/.e21, v3 = (.e38 * .e4/v3 - 
                2 * .e30) * .e20/.e21))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
gev_p1k3_pd=function (x, t, v1, v2, v3, v4) 
{
    .e3 <- x - (t * v2 + v1)
    .e4 <- 1 + v4 * .e3/v3
    .e5 <- 1/v4
    .e7 <- .e4^(1 + .e5)
    .e8 <- exp(-.e4^-.e5)
    .e9 <- v3 * .e7
    c(v1 = -(.e8/.e9), v2 = -(t * .e8/.e9), v3 = -(.e8 * .e3/(v3^2 * 
        .e7)))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
gev_p1k3_pdd=function (x, t, v1, v2, v3, v4) 
{
    .e3 <- x - (t * v2 + v1)
    .e4 <- 1/v4
    .e5 <- 1 + v4 * .e3/v3
    .e6 <- 1 + .e4
    .e7 <- .e5^.e6
    .e9 <- .e5^.e4
    .e10 <- exp(-.e5^-.e4)
    .e11 <- v3 * .e7
    .e12 <- v3^2
    .e13 <- .e5^(2 * .e6)
    .e15 <- v4 * .e6 * .e9
    .e16 <- .e11^2
    .e17 <- .e12 * .e7
    .e20 <- .e15/.e16 - 1/(.e12 * .e13)
    .e21 <- .e17^2
    .e22 <- t * .e10
    .e23 <- .e15 * .e3
    .e24 <- -(.e22 * .e20)
    .e26 <- (.e7 - .e23/v3)/.e16 + .e3/(v3^3 * .e13)
    .e33 <- v3 * v4 * .e6 * .e9 * .e3/.e21 - (.e3/.e11 + 1)/.e17
    c(v1 = c(v1 = -(.e10 * .e20), v2 = .e24, v3 = -(.e10 * .e33)), 
        v2 = c(v1 = .e24, v2 = -(t^2 * .e10 * .e20), v3 = -(.e22 * 
            .e33)), v3 = c(v1 = .e26 * .e10, v2 = t * .e26 * 
            .e10, v3 = ((2 * .e11 - .e23)/.e21 + .e3/(v3^4 * 
            .e13)) * .e10 * .e3))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
gev_p1k3_logfdd=function (x, t, v1, v2, v3, v4) 
{
    .e3 <- x - (t * v2 + v1)
    .e4 <- 1 + v4 * .e3/v3
    .e5 <- 1/v4
    .e6 <- 1 + .e5
    .e8 <- 1/.e4^.e5
    .e9 <- v3 * .e4
    .e10 <- v4 * .e6
    .e11 <- .e10 - .e8
    .e12 <- .e9^2
    .e13 <- .e4^(.e5 + 2)
    .e14 <- v4 * .e11
    .e17 <- .e14/.e12 - 1/(v3^2 * .e13)
    .e19 <- .e11/.e12 + .e3/(v3^3 * .e13)
    .e21 <- t * .e17
    .e24 <- .e14 * .e3/.e12 - (.e3/(v3 * .e4^.e6) + .e10 - .e8)/.e9
    c(v1 = c(v1 = .e17, v2 = .e21, v3 = .e24/v3), v2 = c(v1 = .e21, 
        v2 = t^2 * .e17, v3 = t * .e24/v3), v3 = c(v1 = -.e19, 
        v2 = -(t * .e19), v3 = -(((.e11 * .e3/.e9 - 1)/v3 + .e19 * 
            .e3)/v3)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
gev_p1k3_logfddd=function (x, t, v1, v2, v3, v4) 
{
    .e3 <- x - (t * v2 + v1)
    .e4 <- 1 + v4 * .e3/v3
    .e5 <- 1/v4
    .e6 <- 1 + .e5
    .e7 <- v3 * .e4
    .e8 <- .e4^.e6
    .e9 <- .e7^2
    .e10 <- .e5 + 2
    .e11 <- .e4^.e5
    .e12 <- v4 * .e6
    .e13 <- v3 * .e8
    .e14 <- 1/.e11
    .e15 <- .e4^.e10
    .e16 <- .e12 - .e14
    .e17 <- v3^2
    .e19 <- v3^3 * .e15
    .e20 <- (2 * (v3 * v4 * .e4 * .e16/.e9) - 1/.e13)/.e9
    .e21 <- .e7 * .e16
    .e22 <- (.e17 * .e15)^2
    .e23 <- .e17 * .e8
    .e24 <- .e3/.e13
    .e25 <- .e19^2
    .e27 <- .e3/.e23 + 2 * (.e21/.e9)
    .e28 <- .e20 - .e13 * .e10/.e22
    .e29 <- v3 * .e15
    .e32 <- v4 * .e8 * .e10 * .e3
    .e33 <- .e20 + .e17 * v4 * .e8 * .e10 * .e3/.e25
    .e34 <- .e13^2
    .e36 <- .e24 + .e12 - .e14
    .e37 <- t^2
    .e39 <- v4 * .e27/.e9
    .e41 <- .e12 * .e11 * .e3
    .e44 <- (2 * .e29 - .e32)/.e22 - .e39
    .e45 <- .e33 - 1/.e19
    .e47 <- (v4 * .e16 * .e3/.e9 - .e36/.e7)/v3
    .e56 <- v4 * (2/.e11 + v4 * (2 * (.e21 * .e3/.e9) - (2 + 
        2/v4)) - 2 * .e24)/.e9 - (.e41/.e34 - 2/.e13)/.e7
    .e58 <- .e27/.e9 + v3 * (3 * .e29 - .e32) * .e3/.e25
    .e59 <- .e16/.e9
    .e61 <- t * v4 * .e28
    .e63 <- .e37 * v4 * .e28
    .e64 <- -(t * .e45)
    .e67 <- (((.e8 - .e41/v3)/.e34 + 1/.e23)/.e7 - .e39) * .e3 + 
        .e36/.e9 - .e47
    .e70 <- (.e33 - 2/.e19) * .e3 + .e47 - .e59
    .e71 <- c(v1 = .e61, v2 = .e63, v3 = t * .e56/v3)
    .e72 <- t * .e44
    c(v1 = c(v1 = c(v1 = v4 * .e28, v2 = .e61, v3 = .e56/v3), 
        v2 = .e71, v3 = c(v1 = -.e45, v2 = .e64, v3 = -(.e70/v3))), 
        v2 = c(v1 = .e71, v2 = c(v1 = .e63, v2 = t^3 * v4 * .e28, 
            v3 = .e37 * .e56/v3), v3 = c(v1 = .e64, v2 = -(.e37 * 
            .e45), v3 = -(t * .e70/v3))), v3 = c(v1 = c(v1 = .e44, 
            v2 = .e72, v3 = .e67/v3), v2 = c(v1 = .e72, v2 = .e37 * 
            .e44, v3 = t * .e67/v3), v3 = c(v1 = .e58, v2 = t * 
            .e58, v3 = (.e58 * .e3 + 2 * (((.e16 * .e3/.e7 - 
            1)/v3 + (.e59 + .e3/.e19) * .e3)/v3))/v3)))
}
############################################################
#' The first derivative of the density
#' @inheritParams manf
gev_p1k3_f1fa=function(x,t,v1,v2,v3,kshape){
	kshape=movexiawayfromzero(kshape)
	vf=Vectorize(gev_p1k3_fd)
	f1=vf(x,t,v1,v2,v3,kshape)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @inheritParams manf
gev_p1k3_f2fa=function(x,t,v1,v2,v3,kshape){
	nx=length(x)

	kshape=movexiawayfromzero(kshape)

	vf=Vectorize(gev_p1k3_fdd)
	temp1=vf(x,t,v1,v2,v3,kshape)
	f2=deriv_copyfdd(temp1,nx,dim=3)
	return(f2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @inheritParams manf
gev_p1k3_mu1fa=function(alpha,t,v1,v2,v3,kshape){
	x=extraDistr::qgev((1-alpha),mu=v1+v2*t,sigma=v3,xi=kshape)
	kshape=movexiawayfromzero(kshape)
	vf=Vectorize(gev_p1k3_pd)
	mu1=-vf(x,t,v1,v2,v3,kshape)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @inheritParams manf
gev_p1k3_mu2fa=function(alpha,t,v1,v2,v3,kshape){
	x=extraDistr::qgev((1-alpha),mu=v1+v2*t,sigma=v3,xi=kshape)
	nx=length(x)

	kshape=movexiawayfromzero(kshape)

	vf=Vectorize(gev_p1k3_pdd)
	temp1=vf(x,t,v1,v2,v3,kshape)
	mu2=-deriv_copyfdd(temp1,nx,dim=3)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @inheritParams manf
gev_p1k3_ldda=function(x,t,v1,v2,v3,kshape){
	nx=length(x)

	kshape=movexiawayfromzero(kshape)

	vf=Vectorize(gev_p1k3_logfdd)
	temp1=vf(x,t,v1,v2,v3,kshape)
	ldd=deriv_copyldd(temp1,nx,dim=3)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @inheritParams manf
gev_p1k3_lddda=function(x,t,v1,v2,v3,kshape){
	nx=length(x)
	vf=Vectorize(gev_p1k3_logfddd)

	kshape=movexiawayfromzero(kshape)

	temp1=vf(x,t,v1,v2,v3,kshape)
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}
