######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gev_p12k3_fd=function (x, t1, t2, v1, v2, v3, v4, v5) 
{
    .e2 <- exp(-(t2 * v4 + v3))
    .e3 <- 1/v5
    .e6 <- x - (t1 * v2 + v1)
    .e7 <- 1 + v5 * .e2 * .e6
    .e8 <- 1 + .e3
    .e13 <- exp(-.e7^-.e3)
    .e16 <- v5 * .e8/.e7^(.e3 + 2) - 1/.e7^(2 * .e8)
    .e21 <- .e2 * .e16 * .e6 - 1/.e7^.e8
    .e22 <- .e2^2
    c(v1 = .e13 * .e22 * .e16, v2 = t1 * .e13 * .e22 * .e16, 
        v3 = .e13 * .e2 * .e21, v4 = t2 * .e13 * .e2 * .e21)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_p12k3_fdd=function (x, t1, t2, v1, v2, v3, v4, v5) 
{
    .e1 <- 1/v5
    .e3 <- exp(-(t2 * v4 + v3))
    .e6 <- x - (t1 * v2 + v1)
    .e7 <- 1 + v5 * .e3 * .e6
    .e8 <- 1 + .e1
    .e9 <- .e1 + 2
    .e10 <- 2 * .e8
    .e11 <- v5 * .e8
    .e12 <- .e7^.e8
    .e13 <- .e7^.e9
    .e15 <- 1/.e7^.e10
    .e20 <- exp(-.e7^-.e1)
    .e22 <- .e11/.e13 - .e15
    .e25 <- v5 * .e7^(.e8 - 2 * .e9) * .e9 - 2/.e7^(1 + .e10)
    .e29 <- .e15 + .e11 * (.e3 * .e25 * .e6 - (.e7^(.e1 - .e10) + 
        1/.e13))
    .e32 <- .e3 * .e22 * .e6 - 1/.e12
    .e34 <- .e3 * .e6/.e12
    .e35 <- .e3^2
    .e39 <- .e29 * .e3 * .e6 - (1 + .e34) * .e32
    .e43 <- .e29 - .e32/.e12
    .e44 <- .e3^3
    .e46 <- .e11 * .e25 - .e22/.e12
    .e50 <- .e11 * .e3 * .e25 * .e6 - (2 + .e34) * .e22
    .e51 <- t1 * .e20
    .e53 <- .e51 * .e44 * .e46
    .e54 <- t1 * t2
    .e57 <- t2 * .e39 * .e20 * .e3
    c(v1 = c(v1 = .e20 * .e44 * .e46, v2 = .e53, v3 = .e43 * 
        .e20 * .e35, v4 = t2 * .e43 * .e20 * .e35), v2 = c(v1 = .e53, 
        v2 = t1^2 * .e20 * .e44 * .e46, v3 = t1 * .e43 * .e20 * 
            .e35, v4 = .e54 * .e43 * .e20 * .e35), v3 = c(v1 = .e20 * 
        .e35 * .e50, v2 = .e51 * .e35 * .e50, v3 = .e39 * .e20 * 
        .e3, v4 = .e57), v4 = c(v1 = t2 * .e20 * .e35 * .e50, 
        v2 = .e54 * .e20 * .e35 * .e50, v3 = .e57, v4 = t2^2 * 
            .e39 * .e20 * .e3))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gev_p12k3_pd=function (x, t1, t2, v1, v2, v3, v4, v5) 
{
    .e2 <- exp(-(t2 * v4 + v3))
    .e5 <- x - (t1 * v2 + v1)
    .e6 <- 1 + v5 * .e2 * .e5
    .e7 <- 1/v5
    .e9 <- .e6^(1 + .e7)
    .e10 <- exp(-.e6^-.e7)
    .e11 <- .e10 * .e2
    c(v1 = -(.e11/.e9), v2 = -(t1 * .e10 * .e2/.e9), v3 = -(.e11 * 
        .e5/.e9), v4 = -(t2 * .e10 * .e2 * .e5/.e9))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_p12k3_pdd=function (x, t1, t2, v1, v2, v3, v4, v5) 
{
    .e2 <- exp(-(t2 * v4 + v3))
    .e3 <- 1/v5
    .e6 <- x - (t1 * v2 + v1)
    .e7 <- 1 + .e3
    .e8 <- 1 + v5 * .e2 * .e6
    .e10 <- 2 * .e7
    .e13 <- exp(-.e8^-.e3)
    .e15 <- v5 * .e7 * .e8^(.e3 - .e10)
    .e19 <- .e15 * .e2 * .e6 - (1 + .e2 * .e6/.e8^.e7)/.e8^.e7
    .e22 <- .e2^2
    .e23 <- t1 * .e13
    .e26 <- t2 * .e13 * .e2 * .e19
    .e27 <- .e15 - 1/.e8^.e10
    .e29 <- .e13 * .e2 * .e19
    .e30 <- -.e29
    .e31 <- -(.e23 * .e2 * .e19)
    .e32 <- -(.e23 * .e22 * .e27)
    .e33 <- -(t1 * t2 * .e13 * .e2 * .e19)
    .e34 <- -(.e26 * .e6)
    .e35 <- -.e26
    c(v1 = c(v1 = -(.e13 * .e22 * .e27), v2 = .e32, v3 = .e30, 
        v4 = .e35), v2 = c(v1 = .e32, v2 = -(t1^2 * .e13 * .e22 * 
        .e27), v3 = .e31, v4 = .e33), v3 = c(v1 = .e30, v2 = .e31, 
        v3 = -(.e29 * .e6), v4 = .e34), v4 = c(v1 = .e35, v2 = .e33, 
        v3 = .e34, v4 = -(t2^2 * .e13 * .e2 * .e19 * .e6)))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_p12k3_logfdd=function (x, t1, t2, v1, v2, v3, v4, v5) 
{
    .e2 <- exp(-(t2 * v4 + v3))
    .e5 <- x - (t1 * v2 + v1)
    .e6 <- 1 + v5 * .e2 * .e5
    .e7 <- 1/v5
    .e9 <- 1/.e6^.e7
    .e11 <- v5 * (1 + .e7)
    .e14 <- v5 * (.e11 - .e9) - .e9
    .e17 <- (1/.e6^(.e7 - 1) + .e2 * .e14 * .e5)/.e6 - .e11
    .e18 <- .e6^2
    .e19 <- .e2^2
    .e21 <- t2 * .e17 * .e2
    .e22 <- .e17 * .e2
    .e23 <- .e22/.e6
    .e26 <- t1 * .e17 * .e2/.e6
    .e29 <- t1 * .e19 * .e14/.e18
    .e33 <- t1 * t2 * .e17 * .e2/.e6
    .e35 <- .e21 * .e5/.e6
    .e36 <- .e21/.e6
    c(v1 = c(v1 = .e19 * .e14/.e18, v2 = .e29, v3 = .e23, v4 = .e36), 
        v2 = c(v1 = .e29, v2 = t1^2 * .e19 * .e14/.e18, v3 = .e26, 
            v4 = .e33), v3 = c(v1 = .e23, v2 = .e26, v3 = .e22 * 
            .e5/.e6, v4 = .e35), v4 = c(v1 = .e36, v2 = .e33, 
            v3 = .e35, v4 = t2^2 * .e17 * .e2 * .e5/.e6))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
gev_p12k3_logfddd=function (x, t1, t2, v1, v2, v3, v4, v5) 
{
    .e2 <- exp(-(t2 * v4 + v3))
    .e5 <- x - (t1 * v2 + v1)
    .e6 <- 1/v5
    .e7 <- 1 + v5 * .e2 * .e5
    .e8 <- .e7^.e6
    .e9 <- 1 + .e6
    .e10 <- 1/.e8
    .e11 <- .e6 - 1
    .e12 <- v5 * .e9
    .e13 <- .e7^.e11
    .e14 <- 1/.e13
    .e17 <- v5 * (.e12 - .e10) - .e10
    .e20 <- .e2 * .e17 * .e5
    .e22 <- .e7^(.e6 - (2 + 2 * .e11)) * .e11
    .e24 <- (1/.e7^.e9 + v5/.e7^.e9) * .e2 * .e5
    .e25 <- 2 * .e12
    .e26 <- 2/.e13
    .e33 <- ((2/.e8 + v5 * (.e22 + (.e14 + .e26 + .e20)/.e7 - 
        .e25) - .e24) * .e2 * .e5 - .e14)/.e7 + v5 * (((.e14 + 
        .e20)/.e7 - .e12) * .e2 * .e5/.e7 + 1 + .e6)
    .e34 <- .e7^2
    .e35 <- 2 * .e17
    .e36 <- .e2^2
    .e40 <- .e10 + v5 * (.e22 + (.e14 + 2 * .e20 + .e26)/.e7 - 
        .e25) - .e24
    .e43 <- v5 * (.e35 - .e10) - .e10
    .e44 <- t1 * t2
    .e45 <- t1^2
    .e47 <- t2 * .e33 * .e2
    .e48 <- t2^2
    .e49 <- .e7^3
    .e53 <- .e2 * .e43 * .e5/.e7 - .e35
    .e54 <- .e2^3
    .e57 <- .e44 * .e33 * .e2/.e7
    .e58 <- .e47/.e7
    .e60 <- .e48 * .e33 * .e2
    .e61 <- .e33 * .e2
    .e64 <- t1 * .e40 * .e36/.e34
    .e67 <- .e44 * .e40 * .e36/.e34
    .e68 <- .e61/.e7
    .e71 <- t1 * .e33 * .e2/.e7
    .e74 <- t1 * .e54 * .e43/.e49
    .e78 <- t1 * .e48 * .e33 * .e2/.e7
    .e81 <- .e45 * .e54 * .e43/.e49
    .e82 <- .e45 * t2
    .e84 <- .e47 * .e5/.e7
    .e86 <- .e60 * .e5/.e7
    .e87 <- .e60/.e7
    .e89 <- .e40 * .e36/.e34
    .e90 <- c(v1 = .e74, v2 = .e81, v3 = .e64, v4 = .e67)
    .e91 <- c(v1 = .e58, v2 = .e57, v3 = .e84, v4 = .e86)
    .e94 <- t1 * .e36 * .e53/.e34
    .e97 <- .e44 * .e36 * .e53/.e34
    .e100 <- .e45 * .e40 * .e36/.e34
    .e103 <- .e82 * .e40 * .e36/.e34
    .e106 <- t2 * .e40 * .e36/.e34
    c(v1 = c(v1 = c(v1 = .e54 * .e43/.e49, v2 = .e74, v3 = .e89, 
        v4 = .e106), v2 = .e90, v3 = c(v1 = .e89, v2 = .e64, 
        v3 = .e68, v4 = .e58), v4 = c(v1 = .e106, v2 = .e67, 
        v3 = .e58, v4 = .e87)), v2 = c(v1 = .e90, v2 = c(v1 = .e81, 
        v2 = t1^3 * .e54 * .e43/.e49, v3 = .e100, v4 = .e103), 
        v3 = c(v1 = .e64, v2 = .e100, v3 = .e71, v4 = .e57), 
        v4 = c(v1 = .e67, v2 = .e103, v3 = .e57, v4 = .e78)), 
        v3 = c(v1 = c(v1 = .e36 * .e53/.e34, v2 = .e94, v3 = .e68, 
            v4 = .e58), v2 = c(v1 = .e94, v2 = .e45 * .e36 * 
            .e53/.e34, v3 = .e71, v4 = .e57), v3 = c(v1 = .e68, 
            v2 = .e71, v3 = .e61 * .e5/.e7, v4 = .e84), v4 = .e91), 
        v4 = c(v1 = c(v1 = t2 * .e36 * .e53/.e34, v2 = .e97, 
            v3 = .e58, v4 = .e87), v2 = c(v1 = .e97, v2 = .e82 * 
            .e36 * .e53/.e34, v3 = .e57, v4 = .e78), v3 = .e91, 
            v4 = c(v1 = .e87, v2 = .e78, v3 = .e86, v4 = t2^3 * 
                .e33 * .e2 * .e5/.e7)))
}
############################################################
#' The first derivative of the density
#' @returns Vector
#' @inheritParams manf
gev_p12k3_f1fa=function(x,t,v1,v2,v3,v4,kshape){
	kshape=movexiawayfromzero(kshape)
	vf=Vectorize(gev_p12k3_fd,"x")
	f1=vf(x,t[,1],t[,2],v1,v2,v3,v4,kshape)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @returns Matrix
#' @inheritParams manf
gev_p12k3_f2fa=function(x,t,v1,v2,v3,v4,kshape){
	nx=length(x)

	kshape=movexiawayfromzero(kshape)

	vf=Vectorize(gev_p12k3_fdd,"x")
	temp1=vf(x,t[,1],t[,2],v1,v2,v3,v4,kshape)
	f2=deriv_copyfdd(temp1,nx,dim=4)
	return(f2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @returns Vector
#' @inheritParams manf
gev_p12k3_mu1fa=function(alpha,t,v1,v2,v3,v4,kshape){
	x=extraDistr::qgev((1-alpha),mu=v1+v2*t[,1],sigma=exp(v3+v4*t[,2]),xi=kshape)
	kshape=movexiawayfromzero(kshape)
	vf=Vectorize(gev_p12k3_pd,"x")
	mu1=-vf(x,t[,1],t[,2],v1,v2,v3,v4,kshape)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @returns Matrix
#' @inheritParams manf
gev_p12k3_mu2fa=function(alpha,t,v1,v2,v3,v4,kshape){
	x=extraDistr::qgev((1-alpha),mu=v1+v2*t[,1],sigma=exp(v3+v4*t[,2]),xi=kshape)
	nx=length(x)

	kshape=movexiawayfromzero(kshape)

	vf=Vectorize(gev_p12k3_pdd,"x")
	temp1=vf(x,t[,1],t[,2],v1,v2,v3,v4,kshape)
	mu2=-deriv_copyfdd(temp1,nx,dim=4)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
gev_p12k3_ldda=function(x,t,v1,v2,v3,v4,kshape){
	nx=length(x)

	kshape=movexiawayfromzero(kshape)

	vf=Vectorize(gev_p12k3_logfdd,"x")
	temp1=vf(x,t[,1],t[,2],v1,v2,v3,v4,kshape)
	ldd=deriv_copyldd(temp1,nx,dim=4)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
gev_p12k3_lddda=function(x,t,v1,v2,v3,v4,kshape){
	nx=length(x)
	vf=Vectorize(gev_p12k3_logfddd,"x")

	kshape=movexiawayfromzero(kshape)

	temp1=vf(x,t[,1],t[,2],v1,v2,v3,v4,kshape)
	lddd=deriv_copylddd(temp1,nx,dim=4)
	return(lddd)
}
