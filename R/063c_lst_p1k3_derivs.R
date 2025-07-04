######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
lst_p1k3_fd=function (x, t, v1, v2, v3, v4) 
{
    .e1 <- 1 + v4
    .e2 <- .e1/2
    .e5 <- x - (t * v2 + v1)
    .e6 <- v3^2
    .e7 <- .e5^2
    .e8 <- .e6 * v4
    .e10 <- .e7/.e8 + 1
    .e11 <- .e10^(.e2 + 1)
    .e12 <- gamma(.e2)
    .e13 <- gamma(v4/2)
    .e15 <- sqrt(pi * v4)
    .e20 <- v3^3 * v4 * .e11 * .e13 * .e15
    c(v1 = .e1 * .e12 * .e5/.e20, v2 = t * .e1 * .e12 * .e5/.e20, 
        v3 = .e12 * (v4 * .e1 * .e7/(.e11 * .e8^2) - 1/(.e6 * 
            .e10^.e2))/(.e13 * .e15))
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
lst_p1k3_fdd=function (x, t, v1, v2, v3, v4) 
{
    .e1 <- 1 + v4
    .e2 <- v3^2
    .e3 <- .e1/2
    .e6 <- .e2 * v4
    .e7 <- x - (t * v2 + v1)
    .e8 <- .e7^2
    .e10 <- .e8/.e6 + 1
    .e11 <- .e3 + 1
    .e12 <- gamma(v4/2)
    .e14 <- sqrt(pi * v4)
    .e15 <- .e10^.e11
    .e16 <- .e10^.e3
    .e17 <- .e6^2
    .e22 <- v3^3 * v4 * .e15 * .e12 * .e14
    .e23 <- gamma(.e3)
    .e24 <- .e22^2
    .e25 <- .e15 * .e17
    .e28 <- 2 * (v3 * .e11 * .e16 * .e12 * .e14 * .e8/.e24) - 
        1/.e22
    .e29 <- .e25^2
    .e30 <- .e11 * .e16
    .e31 <- .e10^(.e3 - 1)
    .e32 <- (.e2 * .e16)^2
    .e33 <- .e12 * .e14
    .e34 <- t * .e1
    .e41 <- 2 * (.e30 * .e17 * .e8/(.e2 * .e29)) - (.e31/(v4 * 
        .e32) + 2 * (v4/.e25))
    .e44 <- 3 * .e15 - 2 * (.e6 * .e11 * .e16 * .e8/.e17)
    .e46 <- .e34 * .e28 * .e23
    c(v1 = c(v1 = .e1 * .e28 * .e23, v2 = .e46, v3 = .e1 * .e41 * 
        .e23 * .e7/.e33), v2 = c(v1 = .e46, v2 = t^2 * .e1 * 
        .e28 * .e23, v3 = .e34 * .e41 * .e23 * .e7/.e33), v3 = c(v1 = -(.e6 * 
        .e1 * .e44 * .e23 * .e12 * .e14 * .e7/.e24), v2 = -(t * 
        .e2 * v4 * .e1 * .e44 * .e23 * .e12 * .e14 * .e7/.e24), 
        v3 = v3 * ((2 * .e16 - .e6 * .e31 * .e1 * .e8/.e17)/.e32 - 
            v4^2 * .e1 * (4 * (.e6 * .e15) - 2 * (.e30 * .e8)) * 
                .e8/.e29) * .e23/.e33))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
lst_p1k3_logfdd=function (x, t, v1, v2, v3, v4) 
{
    .e1 <- v3^2
    .e2 <- .e1 * v4
    .e5 <- x - (t * v2 + v1)
    .e6 <- .e5^2
    .e8 <- .e6/.e2 + 1
    .e9 <- .e2 * .e8
    .e10 <- .e2^2
    .e11 <- 1 + v4
    .e12 <- .e8 * .e10
    .e13 <- .e9^2
    .e17 <- 2 * (.e6/.e13) - 1/.e9
    .e18 <- .e12^2
    .e19 <- t * .e11
    .e20 <- v3 * v4
    .e24 <- 2 * (.e10 * .e6/(v3 * .e18)) - 2 * (.e20/.e12)
    .e26 <- 2 * .e8 - 2 * (.e2 * .e6/.e10)
    .e27 <- .e19 * .e17
    c(v1 = c(v1 = .e11 * .e17, v2 = .e27, v3 = .e11 * .e24 * 
        .e5), v2 = c(v1 = .e27, v2 = t^2 * .e11 * .e17, v3 = .e19 * 
        .e24 * .e5), v3 = c(v1 = -(.e20 * .e11 * .e26 * .e5/.e13), 
        v2 = -(t * v3 * v4 * .e11 * .e26 * .e5/.e13), v3 = 1/.e1 + 
            v4 * .e11 * (1/.e12 - .e2 * (4 * .e9 - 2 * .e6)/.e18) * 
                .e6))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
lst_p1k3_logfddd=function (x, t, v1, v2, v3, v4) 
{
    .e1 <- v3^2
    .e2 <- .e1 * v4
    .e5 <- x - (t * v2 + v1)
    .e6 <- .e5^2
    .e8 <- .e6/.e2 + 1
    .e9 <- .e2^2
    .e10 <- .e2 * .e8
    .e11 <- .e8 * .e9
    .e12 <- .e10^2
    .e13 <- .e11^2
    .e14 <- 1 + v4
    .e15 <- 4 * (.e10 * .e6/.e12)
    .e17 <- 2 * .e8 - 2 * (.e2 * .e6/.e9)
    .e18 <- 2 * .e6
    .e20 <- 4 * .e10 - .e18
    .e22 <- 2 * (.e15 - 2) - 2
    .e23 <- v4^2
    .e24 <- (v3 * .e13)^2
    .e25 <- t^2
    .e26 <- v3 * v4
    .e28 <- .e8 * (4 * (.e2 * .e17 * .e6/.e12) - 2) + (6 * (.e2/.e9) - 
        4/.e2) * .e6
    .e35 <- (2 * (4 * (.e8 * .e2^4 * .e6/(v4 * .e24)) - 2/.e13) - 
        6/.e13) * .e9 * .e6/v3 + 2 * (.e26/.e11)
    .e36 <- 1 - .e15
    .e38 <- 1/.e11 - .e2 * .e20/.e13
    .e39 <- t * .e14
    .e41 <- t * v3 * v4
    .e43 <- v3^4 * .e23
    .e46 <- .e39 * .e22 * .e5/.e12
    .e50 <- .e25 * .e14 * .e22 * .e5/.e12
    .e52 <- .e10 * .e20 * .e9
    .e53 <- -(.e41 * .e28 * .e14/.e12)
    .e58 <- ((2/.e2 - 4 * (.e8 * .e20 * .e9/.e13)) * .e9 + 4 * 
        .e2) * .e6/.e13 - 2 * .e38
    .e60 <- .e8 * (2 - 2 * (.e43 * .e17^2/.e12)) - .e2 * (2 * 
        (2 - 4 * (.e43/.e9)) + 6) * .e6/.e9
    .e63 <- 2 * ((4 * (.e1 * .e23/.e13) - (.e13 + 2 * .e52) * 
        .e9/.e24) * .e6) - 2 * (v4 * .e38)
    .e64 <- c(v1 = .e46, v2 = .e50, v3 = t * .e35 * .e14)
    .e68 <- .e41 * .e36 * .e14 * .e17/.e12
    .e69 <- t * v4
    .e71 <- .e25 * v3 * v4
    c(v1 = c(v1 = c(v1 = .e14 * .e22 * .e5/.e12, v2 = .e46, v3 = .e35 * 
        .e14), v2 = .e64, v3 = c(v1 = -(.e26 * .e28 * .e14/.e12), 
        v2 = .e53, v3 = v4 * .e58 * .e14 * .e5)), v2 = c(v1 = .e64, 
        v2 = c(v1 = .e50, v2 = t^3 * .e14 * .e22 * .e5/.e12, 
            v3 = .e25 * .e35 * .e14), v3 = c(v1 = .e53, v2 = -(.e71 * 
            .e28 * .e14/.e12), v3 = .e69 * .e58 * .e14 * .e5)), 
        v3 = c(v1 = c(v1 = .e26 * .e36 * .e14 * .e17/.e12, v2 = .e68, 
            v3 = .e14 * .e63 * .e5), v2 = c(v1 = .e68, v2 = .e71 * 
            .e36 * .e14 * .e17/.e12, v3 = .e39 * .e63 * .e5), 
            v3 = c(v1 = -(v4 * .e60 * .e14 * .e5/.e12), v2 = -(.e69 * 
                .e60 * .e14 * .e5/.e12), v3 = -(2/v3^3 + v3 * 
                .e23 * ((2 - 2 * (.e52/.e13)) * .e20 + .e2 * 
                (4 * .e8 + 4 * .e17) - .e18) * .e14 * .e6/.e13))))
}
############################################################
#' The first derivative of the density
#' @returns Vector
#' @inheritParams manf
lst_p1k3_f1fa=function(x,t,v1,v2,v3,kdf){
	vf=Vectorize(lst_p1k3_fd,"x")
	f1=vf(x,t,v1,v2,v3,kdf)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @returns Matrix
#' @inheritParams manf
lst_p1k3_f2fa=function(x,t,v1,v2,v3,kdf){
	nx=length(x)
	vf=Vectorize(lst_p1k3_fdd,"x")
	temp1=vf(x,t,v1,v2,v3,kdf)
	f2=deriv_copyfdd(temp1,nx,dim=3)
	return(f2)
}
############################################################
############################################################
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
lst_p1k3_ldda=function(x,t,v1,v2,v3,kdf){
	nx=length(x)
	vf=Vectorize(lst_p1k3_logfdd,"x")
	temp1=vf(x,t,v1,v2,v3,kdf)
	ldd=deriv_copyldd(temp1,nx,dim=3)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
lst_p1k3_lddda=function(x,t,v1,v2,v3,kdf){
	nx=length(x)
	vf=Vectorize(lst_p1k3_logfddd,"x")
	temp1=vf(x,t,v1,v2,v3,kdf)
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}
