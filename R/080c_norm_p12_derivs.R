######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
norm_p12_fd=function (x, t1, t2, v1, v2, v3, v4)
{
    .e1 <- t2 * v4
    .e4 <- exp(2 * .e1 + 2 * v3)
    .e7 <- x - (t1 * v2 + v1)
    .e8 <- .e7^2
    .e9 <- 2 * .e4
    .e12 <- exp(-(.e8/.e9))
    .e13 <- exp(.e1 + v3)
    .e14 <- sqrt(2 * pi)
    .e17 <- 4 * (.e4 * .e8/.e9^2) - 1
    .e19 <- .e4 * .e13 * .e14
    .e20 <- .e13 * .e14
    c(v1 = .e12 * .e7/.e19, v2 = t1 * .e12 * .e7/.e19, v3 = .e17 *
        .e12/.e20, v4 = t2 * .e17 * .e12/.e20)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
norm_p12_fdd=function (x, t1, t2, v1, v2, v3, v4)
{
    .e1 <- t2 * v4
    .e4 <- exp(2 * .e1 + 2 * v3)
    .e7 <- x - (t1 * v2 + v1)
    .e8 <- 2 * .e4
    .e9 <- .e7^2
    .e11 <- exp(.e1 + v3)
    .e12 <- sqrt(2 * pi)
    .e13 <- .e8^2
    .e15 <- exp(-(.e9/.e8))
    .e17 <- 4 * (.e4 * .e9/.e13) - 1
    .e19 <- .e4 * .e11 * .e12
    .e21 <- .e13 * .e11 * .e12
    .e22 <- .e11 * .e12
    .e26 <- (4 * (2 - 16 * (.e4^2/.e13)) + 4 * .e17) * .e4 *
        .e9/.e21 - .e17 * .e11 * .e12/.e22^2
    .e28 <- .e17/.e4 - 8 * (.e4/.e13)
    .e32 <- .e9/.e4 - 1
    .e35 <- 4 * (.e9/.e21) - 3 * (.e19/.e19^2)
    .e38 <- t1 * .e32 * .e15/.e19
    .e39 <- t1 * t2
    .e41 <- t2 * .e26 * .e15
    c(v1 = c(v1 = .e32 * .e15/.e19, v2 = .e38, v3 = .e28 * .e15 *
        .e7/.e22, v4 = t2 * .e28 * .e15 * .e7/.e22), v2 = c(v1 = .e38,
        v2 = t1^2 * .e32 * .e15/.e19, v3 = t1 * .e28 * .e15 *
            .e7/.e22, v4 = .e39 * .e28 * .e15 * .e7/.e22), v3 = c(v1 = .e35 *
        .e15 * .e7, v2 = t1 * .e35 * .e15 * .e7, v3 = .e26 *
        .e15, v4 = .e41), v4 = c(v1 = t2 * .e35 * .e15 * .e7,
        v2 = .e39 * .e35 * .e15 * .e7, v3 = .e41, v4 = t2^2 *
            .e26 * .e15))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
norm_p12_pd=function (x, t1, t2, v1, v2, v3, v4)
{
    .e1 <- exp(t2 * v4 + v3)
    .e4 <- x - (t1 * v2 + v1)
    .e6 <- dnorm(.e4/.e1, 0, 1)
    c(v1 = -(.e6/.e1), v2 = -(t1 * .e6/.e1), v3 = -(.e6 * .e4/.e1),
        v4 = -(t2 * .e6 * .e4/.e1))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
norm_p12_pdd=function (x, t1, t2, v1, v2, v3, v4)
{
    .e1 <- exp(t2 * v4 + v3)
    .e4 <- x - (t1 * v2 + v1)
    .e6 <- dnorm(.e4/.e1, 0, 1)
    .e9 <- .e4^2/.e1^2 - 1
    .e10 <- .e1^3
    .e12 <- t2 * .e9 * .e6
    .e13 <- .e9 * .e6
    .e14 <- -(.e13/.e1)
    .e15 <- -(t1 * .e9 * .e6/.e1)
    .e16 <- -(t1 * .e6 * .e4/.e10)
    .e17 <- -(t1 * t2 * .e9 * .e6/.e1)
    .e18 <- -(.e12 * .e4/.e1)
    .e19 <- -(.e12/.e1)
    c(v1 = c(v1 = -(.e6 * .e4/.e10), v2 = .e16, v3 = .e14, v4 = .e19),
        v2 = c(v1 = .e16, v2 = -(t1^2 * .e6 * .e4/.e10), v3 = .e15,
            v4 = .e17), v3 = c(v1 = .e14, v2 = .e15, v3 = -(.e13 *
            .e4/.e1), v4 = .e18), v4 = c(v1 = .e19, v2 = .e17,
            v3 = .e18, v4 = -(t2^2 * .e9 * .e6 * .e4/.e1)))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
norm_p12_logfdd=function (x, t1, t2, v1, v2, v3, v4)
{
    .e3 <- exp(2 * (t2 * v4) + 2 * v3)
    .e4 <- (2 * .e3)^2
    .e7 <- x - (t1 * v2 + v1)
    .e8 <- .e7^2
    .e10 <- 2 - 16 * (.e3^2/.e4)
    .e11 <- -(t1/.e3)
    .e12 <- 4 * (t2 * .e10 * .e3 * .e8/.e4)
    .e13 <- t1 * t2
    c(v1 = c(v1 = -(1/.e3), v2 = .e11, v3 = -(8 * (.e3 * .e7/.e4)),
        v4 = -(8 * (t2 * .e3 * .e7/.e4))), v2 = c(v1 = .e11,
        v2 = -(t1^2/.e3), v3 = -(8 * (t1 * .e3 * .e7/.e4)), v4 = -(8 *
            (.e13 * .e3 * .e7/.e4))), v3 = c(v1 = -(2 * (.e7/.e3)),
        v2 = -(2 * (t1 * .e7/.e3)), v3 = 4 * (.e10 * .e3 * .e8/.e4),
        v4 = .e12), v4 = c(v1 = -(2 * (t2 * .e7/.e3)), v2 = -(2 *
        (.e13 * .e7/.e3)), v3 = .e12, v4 = 4 * (t2^2 * .e10 *
        .e3 * .e8/.e4)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
norm_p12_logfddd=function (x, t1, t2, v1, v2, v3, v4)
{
    .e3 <- exp(2 * (t2 * v4) + 2 * v3)
    .e4 <- (2 * .e3)^2
    .e5 <- .e3^2
    .e6 <- 16 * (.e5/.e4)
    .e9 <- x - (t1 * v2 + v1)
    .e10 <- 2 - .e6
    .e11 <- t1 * t2
    .e12 <- t2^2
    .e14 <- .e10^2 - 16 * ((4 - .e6) * .e5/.e4)
    .e15 <- .e9^2
    .e16 <- t1^2
    .e17 <- -(8 * (.e11 * .e10 * .e3 * .e9/.e4))
    .e18 <- -(8 * (t2 * .e10 * .e3 * .e9/.e4))
    .e19 <- 2 * (.e11/.e3)
    .e20 <- 2 * (t1/.e3)
    .e21 <- 4 * (t2 * .e14 * .e3 * .e15/.e4)
    .e22 <- 4 * (.e12 * .e14 * .e3 * .e15/.e4)
    .e23 <- t1 * .e12
    .e24 <- .e16 * t2
    .e37 <- c(v1 = 0, v2 = 0, v3 = 8 * (t1 * .e3/.e4), v4 = 8 *
        (.e11 * .e3/.e4))
    .e38 <- c(v1 = .e19, v2 = 2 * (.e24/.e3), v3 = .e17, v4 = -(8 *
        (.e23 * .e10 * .e3 * .e9/.e4)))
    .e39 <- c(v1 = .e20, v2 = 2 * (.e16/.e3), v3 = -(8 * (t1 *
        .e10 * .e3 * .e9/.e4)), v4 = .e17)
    .e40 <- c(v1 = 2 * (t2/.e3), v2 = .e19, v3 = .e18, v4 = -(8 *
        (.e12 * .e10 * .e3 * .e9/.e4)))
    .e41 <- c(v1 = 2/.e3, v2 = .e20, v3 = -(8 * (.e10 * .e3 *
        .e9/.e4)), v4 = .e18)
    .e42 <- c(v1 = 4 * (t2 * .e9/.e3), v2 = 4 * (.e11 * .e9/.e3),
        v3 = .e21, v4 = .e22)
    c(v1 = c(v1 = c(v1 = 0, v2 = 0, v3 = 8 * (.e3/.e4), v4 = 8 *
        (t2 * .e3/.e4)), v2 = .e37, v3 = .e41, v4 = .e40), v2 = c(v1 = .e37,
        v2 = c(v1 = 0, v2 = 0, v3 = 8 * (.e16 * .e3/.e4), v4 = 8 *
            (.e24 * .e3/.e4)), v3 = .e39, v4 = .e38), v3 = c(v1 = .e41,
        v2 = .e39, v3 = c(v1 = 4 * (.e9/.e3), v2 = 4 * (t1 *
            .e9/.e3), v3 = 4 * (.e14 * .e3 * .e15/.e4), v4 = .e21),
        v4 = .e42), v4 = c(v1 = .e40, v2 = .e38, v3 = .e42, v4 = c(v1 = 4 *
        (.e12 * .e9/.e3), v2 = 4 * (.e23 * .e9/.e3), v3 = .e22,
        v4 = 4 * (t2^3 * .e14 * .e3 * .e15/.e4))))
}
############################################################
#' The first derivative of the density for DMGS
#' @returns Vector
#' @inheritParams manf
norm_p12_f1fa=function(x,t01,t02,v1,v2,v3,v4){
	vf=Vectorize(norm_p12_fd,"x")
	f1=vf(x,t01,t02,v1,v2,v3,v4)
	return(f1)
}
############################################################
#' The first derivative of the density for WAIC
#' @returns Vector
#' @inheritParams manf
norm_p12_f1fw=function(x,t1,t2,v1,v2,v3,v4){
	vf=Vectorize(norm_p12_fd,c("x","t1","t2"))
	f1=vf(x,t1,t2,v1,v2,v3,v4)
	return(f1)
}
############################################################
#' The second derivative of the density for DMGS
#' @returns Matrix
#' @inheritParams manf
norm_p12_f2fa=function(x,t01,t02,v1,v2,v3,v4){
	nx=length(x)
	vf=Vectorize(norm_p12_fdd,"x")
	temp1=vf(x,t01,t02,v1,v2,v3,v4)
	f2=deriv_copyfdd(temp1,nx,dim=4)
	return(f2)
}
############################################################
#' The second derivative of the density for WAIC
#' @returns Matrix
#' @inheritParams manf
norm_p12_f2fw=function(x,t1,t2,v1,v2,v3,v4){
	nx=length(x)
	vf=Vectorize(norm_p12_fdd,c("x","t1","t2"))
	temp1=vf(x,t1,t2,v1,v2,v3,v4)
	f2=deriv_copyfdd(temp1,nx,dim=4)
	return(f2)
}
############################################################
#' The first derivative of the cdf
#' @returns Vector
#' @inheritParams manf
norm_p12_p1fa=function(x,t01,t02,v1,v2,v3,v4){
	vf=Vectorize(norm_p12_pd,"x")
	p1=vf(x,t01,t02,v1,v2,v3,v4)
	return(p1)
}
############################################################
#' The second derivative of the cdf
#' @returns Matrix
#' @inheritParams manf
norm_p12_p2fa=function(x,t01,t02,v1,v2,v3,v4){
	nx=length(x)
	vf=Vectorize(norm_p12_pdd,"x")
	temp1=vf(x,t01,t02,v1,v2,v3,v4)
	p2=deriv_copyfdd(temp1,nx,dim=4)
	return(p2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @returns Vector
#' @inheritParams manf
norm_p12_mu1fa=function(alpha,t01,t02,v1,v2,v3,v4){
	x=qnorm((1-alpha),mean=v1+v2*t01,sd=exp(v3+v4*t02))
	vf=Vectorize(norm_p12_pd,"x")
	mu1=-vf(x,t01,t02,v1,v2,v3,v4)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @returns Matrix
#' @inheritParams manf
norm_p12_mu2fa=function(alpha,t01,t02,v1,v2,v3,v4){
	x=qnorm((1-alpha),mean=v1+v2*t01,sd=exp(v3+v4*t02))
	nx=length(x)
	vf=Vectorize(norm_p12_pdd,"x")
	temp1=vf(x,t01,t02,v1,v2,v3,v4)
	mu2=-deriv_copyfdd(temp1,nx,dim=4)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
norm_p12_ldda=function(x,t1,t2,v1,v2,v3,v4){
	nx=length(x)
	vf=Vectorize(norm_p12_logfdd,c("x","t1","t2"))
	temp1=vf(x,t1,t2,v1,v2,v3,v4)
	ldd=deriv_copyldd(temp1,nx,dim=4)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
norm_p12_lddda=function(x,t1,t2,v1,v2,v3,v4){
	nx=length(x)
	vf=Vectorize(norm_p12_logfddd,c("x","t1","t2"))
	temp1=vf(x,t1,t2,v1,v2,v3,v4)
	lddd=deriv_copylddd(temp1,nx,dim=4)
	return(lddd)
}
