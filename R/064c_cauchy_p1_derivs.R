######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
cauchy_p1_fd=function (x, t, v1, v2, v3) 
{
    .e3 <- x - (t * v2 + v1)
    .e4 <- .e3^2
    .e6 <- .e4/v3^2 + 1
    .e9 <- pi * v3^3 * .e6^2
    c(v1 = 2 * (.e3/.e9), v2 = 2 * (t * .e3/.e9), v3 = (2 * (.e4/(pi * 
        v3^4 * .e6)) - pi/(pi * v3)^2)/.e6)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @returns Matrix
#' @inheritParams manf
cauchy_p1_fdd=function (x, t, v1, v2, v3) 
{
    .e3 <- x - (t * v2 + v1)
    .e4 <- .e3^2
    .e5 <- v3^2
    .e7 <- .e4/.e5 + 1
    .e8 <- v3^3
    .e11 <- pi * .e8 * .e7^2
    .e12 <- pi * v3
    .e14 <- pi * v3^4 * .e7
    .e15 <- .e11^2
    .e16 <- .e5 * .e7
    .e19 <- 4 * (.e12 * .e7 * .e4/.e15) - 1/.e11
    .e21 <- .e14^2
    .e24 <- 2 * (.e4/.e14) - pi/.e12^2
    .e27 <- 2 * (.e24/.e16) + 2 * (2 * (pi * .e5 * .e4/.e21) - 
        2/.e14)
    .e28 <- 2 * (t * .e19)
    .e30 <- 3 * .e16 - 4 * .e4
    c(v1 = c(v1 = 2 * .e19, v2 = .e28, v3 = .e27 * .e3/.e7), 
        v2 = c(v1 = .e28, v2 = 2 * (t^2 * .e19), v3 = t * .e27 * 
            .e3/.e7), v3 = c(v1 = -(2 * (pi * .e7 * .e30 * .e3/.e15)), 
            v2 = -(2 * (pi * t * .e7 * .e30 * .e3/.e15)), v3 = (2 * 
                (.e24 * .e4/(.e8 * .e7)) + .e12 * (2 * (pi^2/.e12^4) - 
                2 * ((4 * .e16 - 2 * .e4) * .e4/.e21)))/.e7))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
cauchy_p1_logfdd=function (x, t, v1, v2, v3) 
{
    .e3 <- x - (t * v2 + v1)
    .e5 <- (.e3/v3)^2 + 1
    .e7 <- v3^2 * .e5
    .e8 <- .e3^2
    .e9 <- .e7^2
    .e11 <- 2 * (.e8/.e9)
    .e13 <- .e11 - 1/.e7
    .e17 <- 2 * (v3 * .e5) - 2 * (.e8/v3)
    .e18 <- .e11 - 2/.e7
    .e19 <- 2 * (t * .e13)
    c(v1 = c(v1 = 2 * .e13, v2 = .e19, v3 = 2 * (.e18 * .e3/v3)), 
        v2 = c(v1 = .e19, v2 = 2 * (t^2 * .e13), v3 = 2 * (t * 
            .e18 * .e3/v3)), v3 = c(v1 = -(2 * (.e17 * .e3/.e9)), 
            v2 = -(2 * (t * .e17 * .e3/.e9)), v3 = -(((2 * (.e8/.e7) - 
                1)/v3 + 2 * (.e17 * .e8/.e9))/v3)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
cauchy_p1_logfddd=function (x, t, v1, v2, v3) 
{
    .e3 <- x - (t * v2 + v1)
    .e5 <- (.e3/v3)^2 + 1
    .e6 <- v3^2
    .e7 <- .e6 * .e5
    .e8 <- .e7^2
    .e9 <- .e3^2
    .e10 <- 4 * (.e7 * .e9/.e8)
    .e14 <- 2 * (v3 * .e5) - 2 * (.e9/v3)
    .e15 <- .e10 - 2
    .e16 <- 2 * .e15
    .e17 <- .e16 - 2
    .e18 <- 2/.e7
    .e20 <- 2 * (.e9/.e6)
    .e21 <- t^2
    .e22 <- (2 * (.e9/.e8) - .e18)/v3
    .e25 <- (.e16 - 6) * .e9/.e8 + .e18
    .e26 <- 1 - .e10
    .e27 <- .e10 - 1
    .e32 <- 2 * (.e5 - .e20) + .e20 - 2 * (.e7 * .e14^2/.e8)
    .e33 <- 2 * (t * .e17 * .e3/.e8)
    .e34 <- 2 * (.e21 * .e17 * .e3/.e8)
    .e35 <- -(2 * (t * .e14 * .e27/.e8))
    .e38 <- (2 - .e10) * .e14/.e8 - .e22
    .e42 <- 2 * .e22 + 2 * (.e14 * .e15/.e8)
    .e44 <- 2 * (t * .e26 * .e14/.e8)
    .e45 <- c(v1 = .e33, v2 = .e34, v3 = 2 * (t * .e25/v3))
    c(v1 = c(v1 = c(v1 = 2 * (.e17 * .e3/.e8), v2 = .e33, v3 = 2 * 
        (.e25/v3)), v2 = .e45, v3 = c(v1 = -(2 * (.e14 * .e27/.e8)), 
        v2 = .e35, v3 = -(.e42 * .e3/v3))), v2 = c(v1 = .e45, 
        v2 = c(v1 = .e34, v2 = 2 * (t^3 * .e17 * .e3/.e8), v3 = 2 * 
            (.e21 * .e25/v3)), v3 = c(v1 = .e35, v2 = -(2 * (.e21 * 
            .e14 * .e27/.e8)), v3 = -(t * .e42 * .e3/v3))), v3 = c(v1 = c(v1 = 2 * 
        (.e26 * .e14/.e8), v2 = .e44, v3 = 2 * (.e38 * .e3/v3)), 
        v2 = c(v1 = .e44, v2 = 2 * (.e21 * .e26 * .e14/.e8), 
            v3 = 2 * (t * .e38 * .e3/v3)), v3 = c(v1 = -(2 * 
            (.e32 * .e3/.e8)), v2 = -(2 * (t * .e32 * .e3/.e8)), 
            v3 = -((2 * (.e32 * .e9/.e8) - 2 * (((2 * (.e9/.e7) - 
                1)/v3 + 2 * (.e14 * .e9/.e8))/v3))/v3))))
}
############################################################
#' The first derivative of the density for DMGS
#' @returns Vector
#' @inheritParams manf
cauchy_p1_f1fa=function(x,t0,v1,v2,v3){
	vf=Vectorize(cauchy_p1_fd,"x")
	f1=vf(x,t0,v1,v2,v3)
	return(f1)
}
############################################################
#' The first derivative of the density for WAIC
#' @returns Vector
#' @inheritParams manf
cauchy_p1_f1fw=function(x,t,v1,v2,v3){
	vf=Vectorize(cauchy_p1_fd,c("x","t"))
	f1=vf(x,t,v1,v2,v3)
	return(f1)
}
############################################################
#' The second derivative of the density for DMGS
#' @returns Matrix
#' @inheritParams manf
cauchy_p1_f2fa=function(x,t0,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(cauchy_p1_fdd,"x")
	temp1=vf(x,t0,v1,v2,v3)
	f2=deriv_copyfdd(temp1,nx,dim=3)
	return(f2)
}
############################################################
#' The second derivative of the density for WAIC
#' @returns Matrix
#' @inheritParams manf
cauchy_p1_f2fw=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(cauchy_p1_fdd,c("x","t"))
	temp1=vf(x,t,v1,v2,v3)
	f2=deriv_copyfdd(temp1,nx,dim=3)
	return(f2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
cauchy_p1_ldda=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(cauchy_p1_logfdd,c("x","t"))
	temp1=vf(x,t,v1,v2,v3)
	ldd=deriv_copyldd(temp1,nx,dim=3)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
cauchy_p1_lddda=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(cauchy_p1_logfddd,c("x","t"))
	temp1=vf(x,t,v1,v2,v3)
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}
