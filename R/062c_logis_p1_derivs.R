######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
logis_p1_fd=function (x, t, v1, v2, v3) 
{
    .e3 <- x - (t * v2 + v1)
    .e5 <- exp(-(.e3/v3))
    .e6 <- 1 + .e5
    .e8 <- 1 - 2 * (.e5/.e6)
    .e10 <- v3^2 * .e6^2
    c(v1 = .e8 * .e5/.e10, v2 = t * .e8 * .e5/.e10, v3 = (.e8 * 
        .e3/v3 - 1) * .e5/.e10)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
logis_p1_fdd=function (x, t, v1, v2, v3) 
{
    .e3 <- x - (t * v2 + v1)
    .e5 <- exp(-(.e3/v3))
    .e6 <- 1 + .e5
    .e7 <- .e5/.e6
    .e8 <- .e6^2
    .e9 <- 1 - 2 * .e7
    .e10 <- (v3^2 * .e8)^2
    .e11 <- 1 - .e7
    .e13 <- .e9 * .e3/v3
    .e16 <- 1 - (2 + 2 * .e11) * .e5/.e6
    .e18 <- v3^3 * .e8
    .e20 <- .e16/.e18 - 2 * (v3 * .e9 * .e6 * .e5/.e10)
    .e21 <- .e13 - ((2 * (.e11 * .e3/v3) - 2) * .e5/.e6 + 2)
    .e22 <- .e13 - 1
    .e24 <- 2 * (.e5 * .e3) + 2 * (v3 * .e6)
    .e26 <- v3^4 * .e8
    .e28 <- .e21/.e18 - 2 * (v3 * .e22 * .e6 * .e5/.e10)
    .e31 <- .e16 * .e3/.e26 - .e9 * .e6 * .e24/.e10
    .e33 <- t * .e20 * .e5
    c(v1 = c(v1 = .e20 * .e5, v2 = .e33, v3 = .e28 * .e5), v2 = c(v1 = .e33, 
        v2 = t^2 * .e20 * .e5, v3 = t * .e28 * .e5), v3 = c(v1 = .e31 * 
        .e5, v2 = t * .e31 * .e5, v3 = (.e21 * .e3/.e26 - .e22 * 
        .e6 * .e24/.e10) * .e5))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
logis_p1_pd=function (x, t, v1, v2, v3) 
{
    .e3 <- x - (t * v2 + v1)
    .e5 <- exp(-(.e3/v3))
    .e6 <- (1 + .e5)^2
    .e7 <- v3 * .e6
    c(v1 = -(.e5/.e7), v2 = -(t * .e5/.e7), v3 = -(.e5 * .e3/(v3^2 * 
        .e6)))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
logis_p1_pdd=function (x, t, v1, v2, v3) 
{
    .e3 <- x - (t * v2 + v1)
    .e4 <- .e3/v3
    .e6 <- exp(-.e4)
    .e7 <- 1 + .e6
    .e8 <- .e7^2
    .e10 <- v3^2 * .e8
    .e11 <- (v3 * .e8)^2
    .e15 <- 1/.e10 - 2 * (.e7 * .e6/.e11)
    .e16 <- .e10^2
    .e17 <- v3 * .e7
    .e18 <- -(t * .e15 * .e6)
    .e22 <- (.e4 - 1)/.e10 - 2 * (.e17 * .e6 * .e3/.e16)
    .e24 <- .e3/(v3^3 * .e8) - ((1 + 2 * .e4) * .e6 + 1) * .e7/.e11
    c(v1 = c(v1 = -(.e15 * .e6), v2 = .e18, v3 = -(.e22 * .e6)), 
        v2 = c(v1 = .e18, v2 = -(t^2 * .e15 * .e6), v3 = -(t * 
            .e22 * .e6)), v3 = c(v1 = -(.e24 * .e6), v2 = -(t * 
            .e24 * .e6), v3 = -((.e3/(v3^4 * .e8) - .e7 * (2 * 
            (.e6 * .e3) + 2 * .e17)/.e16) * .e6 * .e3)))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
logis_p1_logfdd=function (x, t, v1, v2, v3) 
{
    .e3 <- x - (t * v2 + v1)
    .e5 <- exp(-(.e3/v3))
    .e6 <- 1 + .e5
    .e7 <- 1 - .e5/.e6
    .e8 <- v3^2
    .e11 <- 2 * (.e7 * .e3/v3)
    .e14 <- (.e11 - 2) * .e5/.e6 + 1
    .e15 <- .e8 * .e6
    .e16 <- -(.e14/.e8)
    .e17 <- -(2 * (t * .e7 * .e5/.e15))
    .e18 <- -(t * .e14/.e8)
    c(v1 = c(v1 = -(2 * (.e7 * .e5/.e15)), v2 = .e17, v3 = .e16), 
        v2 = c(v1 = .e17, v2 = -(2 * (t^2 * .e7 * .e5/.e15)), 
            v3 = .e18), v3 = c(v1 = .e16, v2 = .e18, v3 = -((((.e11 - 
            4) * .e5/.e6 + 2) * .e3/v3 - 1)/.e8)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
logis_p1_logfddd=function (x, t, v1, v2, v3) 
{
    .e3 <- x - (t * v2 + v1)
    .e4 <- .e3/v3
    .e6 <- exp(-.e4)
    .e7 <- 1 + .e6
    .e8 <- 1 - .e6/.e7
    .e9 <- .e8 * .e3
    .e10 <- v3 * .e7
    .e11 <- v3^3
    .e13 <- 2 * (.e9/v3)
    .e14 <- .e6 * .e3
    .e16 <- (.e13 - 2) * .e6/.e7
    .e17 <- .e11 * .e7
    .e19 <- 2 * .e4
    .e20 <- 2 * (1 + .e14/.e10)
    .e22 <- .e8 * (.e19 - .e20) - (.e16 + 2)
    .e23 <- (v3^2 * .e7)^2
    .e25 <- .e8/.e17 - v3 * .e6/.e23
    .e26 <- .e22 * .e6
    .e27 <- t^2
    .e28 <- -(t * .e22 * .e6/.e17)
    .e31 <- .e26 * .e3/.e10 - 2 * (.e16 + 1)
    .e33 <- .e9/(v3^4 * .e7) - (2 * .e10 + .e14)/.e23
    .e35 <- (.e13 - 4) * .e6/.e7
    .e36 <- -(2 * (t * .e25 * .e8 * .e6))
    .e37 <- -(2 * (.e27 * .e25 * .e8 * .e6))
    .e40 <- ((.e8 * (.e19 - (2 + .e20)) - (.e35 + 4)) * .e3/v3 + 
        4) * .e6/.e7 - 2
    .e41 <- -(.e31/.e11)
    .e42 <- -(.e26/.e17)
    .e43 <- -(2 * (t * .e33 * .e8 * .e6))
    .e44 <- -(t * .e31/.e11)
    .e45 <- -(.e27 * .e22 * .e6/.e17)
    .e46 <- c(v1 = .e36, v2 = .e37, v3 = .e28)
    c(v1 = c(v1 = c(v1 = -(2 * (.e25 * .e8 * .e6)), v2 = .e36, 
        v3 = .e42), v2 = .e46, v3 = c(v1 = .e42, v2 = .e28, v3 = -(.e40/.e11))), 
        v2 = c(v1 = .e46, v2 = c(v1 = .e37, v2 = -(2 * (t^3 * 
            .e25 * .e8 * .e6)), v3 = .e45), v3 = c(v1 = .e28, 
            v2 = .e45, v3 = -(t * .e40/.e11))), v3 = c(v1 = c(v1 = -(2 * 
            (.e33 * .e8 * .e6)), v2 = .e43, v3 = .e41), v2 = c(v1 = .e43, 
            v2 = -(2 * (.e27 * .e33 * .e8 * .e6)), v3 = .e44), 
            v3 = c(v1 = .e41, v2 = .e44, v3 = -((.e40 * .e3/v3 - 
                2 * ((.e35 + 2) * .e3/v3 - 1))/.e11))))
}
############################################################
#' The first derivative of the density for DMGS
#' @returns Vector
#' @inheritParams manf
logis_p1_f1fa=function(x,t0,v1,v2,v3){
	vf=Vectorize(logis_p1_fd,"x")
	f1=vf(x,t0,v1,v2,v3)
	return(f1)
}
############################################################
#' The first derivative of the density for WAIC
#' @returns Vector
#' @inheritParams manf
logis_p1_f1fw=function(x,t,v1,v2,v3){
	vf=Vectorize(logis_p1_fd,c("x","t"))
	f1=vf(x,t,v1,v2,v3)
	return(f1)
}
############################################################
#' The second derivative of the density for DMGS
#' @returns Matrix
#' @inheritParams manf
logis_p1_f2fa=function(x,t0,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(logis_p1_fdd,"x")
	temp1=vf(x,t0,v1,v2,v3)
	f2=deriv_copyfdd(temp1,nx,dim=3)
	return(f2)
}
############################################################
#' The second derivative of the density for WAIC
#' @returns Matrix
#' @inheritParams manf
logis_p1_f2fw=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(logis_p1_fdd,c("x","t"))
	temp1=vf(x,t,v1,v2,v3)
	f2=deriv_copyfdd(temp1,nx,dim=3)
	return(f2)
}
############################################################
#' The first derivative of the cdf
#' @returns Vector
#' @inheritParams manf
logis_p1_p1fa=function(x,t0,v1,v2,v3){
	vf=Vectorize(logis_p1_pd,"x")
	p1=vf(x,t0,v1,v2,v3)
	return(p1)
}
############################################################
#' The second derivative of the cdf
#' @returns Matrix
#' @inheritParams manf
logis_p1_p2fa=function(x,t0,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(logis_p1_pdd,"x")
	temp1=vf(x,t0,v1,v2,v3)
	p2=deriv_copyfdd(temp1,nx,dim=3)
	return(p2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @returns Vector
#' @inheritParams manf
logis_p1_mu1fa=function(alpha,t0,v1,v2,v3){
	x=qlogis((1-alpha),location=v1+v2*t0,scale=v3)
	vf=Vectorize(logis_p1_pd,"x")
	mu1=-vf(x,t0,v1,v2,v3)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @returns Matrix
#' @inheritParams manf
logis_p1_mu2fa=function(alpha,t0,v1,v2,v3){
	x=qlogis((1-alpha),location=v1+v2*t0,scale=v3)
	nx=length(x)
	vf=Vectorize(logis_p1_pdd,"x")
	temp1=vf(x,t0,v1,v2,v3)
	mu2=-deriv_copyfdd(temp1,nx,dim=3)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
logis_p1_ldda=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(logis_p1_logfdd,c("x","t"))
	temp1=vf(x,t,v1,v2,v3)
	ldd=deriv_copyldd(temp1,nx,dim=3)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
logis_p1_lddda=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(logis_p1_logfddd,c("x","t"))
	temp1=vf(x,t,v1,v2,v3)
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}
