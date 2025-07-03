######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
lst_k3_fd=function (x, v1, v2, v3) 
{
    .e1 <- 1 + v3
    .e2 <- .e1/2
    .e3 <- v2^2
    .e4 <- x - v1
    .e5 <- .e4^2
    .e6 <- .e3 * v3
    .e8 <- .e5/.e6 + 1
    .e9 <- .e8^(.e2 + 1)
    .e10 <- gamma(.e2)
    .e11 <- gamma(v3/2)
    .e13 <- sqrt(pi * v3)
    c(v1 = .e1 * .e10 * .e4/(v2^3 * v3 * .e9 * .e11 * .e13), 
        v2 = .e10 * (v3 * .e1 * .e5/(.e9 * .e6^2) - 1/(.e3 * 
            .e8^.e2))/(.e11 * .e13))
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
lst_k3_fdd=function (x, v1, v2, v3) 
{
    .e1 <- 1 + v3
    .e2 <- v2^2
    .e3 <- .e2 * v3
    .e4 <- .e1/2
    .e5 <- x - v1
    .e6 <- .e5^2
    .e8 <- .e6/.e3 + 1
    .e9 <- .e4 + 1
    .e10 <- .e8^.e9
    .e11 <- .e8^.e4
    .e12 <- gamma(v3/2)
    .e14 <- sqrt(pi * v3)
    .e15 <- .e3^2
    .e16 <- gamma(.e4)
    .e17 <- .e10 * .e15
    .e22 <- v2^3 * v3 * .e10 * .e12 * .e14
    .e23 <- .e17^2
    .e24 <- .e9 * .e11
    .e25 <- .e8^(.e4 - 1)
    .e26 <- (.e2 * .e11)^2
    .e27 <- .e22^2
    .e28 <- .e12 * .e14
    c(v1 = c(v1 = .e1 * (2 * (v2 * .e9 * .e11 * .e12 * .e14 * 
        .e6/.e27) - 1/.e22) * .e16, v2 = .e1 * (2 * (.e24 * .e15 * 
        .e6/(.e2 * .e23)) - (.e25/(v3 * .e26) + 2 * (v3/.e17))) * 
        .e16 * .e5/.e28), v2 = c(v1 = -(.e3 * .e1 * (3 * .e10 - 
        2 * (.e3 * .e9 * .e11 * .e6/.e15)) * .e16 * .e12 * .e14 * 
        .e5/.e27), v2 = v2 * ((2 * .e11 - .e3 * .e25 * .e1 * 
        .e6/.e15)/.e26 - v3^2 * .e1 * (4 * (.e3 * .e10) - 2 * 
        (.e24 * .e6)) * .e6/.e23) * .e16/.e28))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
lst_k3_logfdd=function (x, v1, v2, v3) 
{
    .e1 <- v2^2
    .e2 <- .e1 * v3
    .e3 <- x - v1
    .e4 <- .e3^2
    .e6 <- .e4/.e2 + 1
    .e7 <- .e2^2
    .e8 <- .e6 * .e7
    .e9 <- 1 + v3
    .e10 <- .e2 * .e6
    .e11 <- .e8^2
    .e12 <- .e10^2
    .e13 <- v2 * v3
    c(v1 = c(v1 = .e9 * (2 * (.e4/.e12) - 1/.e10), v2 = .e9 * 
        (2 * (.e7 * .e4/(v2 * .e11)) - 2 * (.e13/.e8)) * .e3), 
        v2 = c(v1 = -(.e13 * .e9 * (2 * .e6 - 2 * (.e2 * .e4/.e7)) * 
            .e3/.e12), v2 = 1/.e1 + v3 * .e9 * (1/.e8 - .e2 * 
            (4 * .e10 - 2 * .e4)/.e11) * .e4))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
lst_k3_logfddd=function (x, v1, v2, v3) 
{
    .e1 <- v2^2
    .e2 <- .e1 * v3
    .e3 <- x - v1
    .e4 <- .e3^2
    .e6 <- .e4/.e2 + 1
    .e7 <- .e2^2
    .e8 <- .e2 * .e6
    .e9 <- .e6 * .e7
    .e10 <- .e9^2
    .e11 <- .e8^2
    .e12 <- 1 + v3
    .e13 <- 2 * .e4
    .e15 <- 4 * .e8 - .e13
    .e17 <- 2 * .e6 - 2 * (.e2 * .e4/.e7)
    .e18 <- v3^2
    .e19 <- v2 * v3
    .e20 <- (v2 * .e10)^2
    .e22 <- 1/.e9 - .e2 * .e15/.e10
    .e23 <- 4 * (.e8 * .e4/.e11)
    .e25 <- .e8 * .e15 * .e7
    .e27 <- v2^4 * .e18
    c(v1 = c(v1 = c(v1 = .e12 * (2 * (.e23 - 2) - 2) * .e3/.e11, 
        v2 = ((2 * (4 * (.e6 * .e2^4 * .e4/(v3 * .e20)) - 2/.e10) - 
            6/.e10) * .e7 * .e4/v2 + 2 * (.e19/.e9)) * .e12), 
        v2 = c(v1 = -(.e19 * (.e6 * (4 * (.e2 * .e17 * .e4/.e11) - 
            2) + (6 * (.e2/.e7) - 4/.e2) * .e4) * .e12/.e11), 
            v2 = v3 * (((2/.e2 - 4 * (.e6 * .e15 * .e7/.e10)) * 
                .e7 + 4 * .e2) * .e4/.e10 - 2 * .e22) * .e12 * 
                .e3)), v2 = c(v1 = c(v1 = .e19 * (1 - .e23) * 
        .e12 * .e17/.e11, v2 = .e12 * (2 * ((4 * (.e1 * .e18/.e10) - 
        (.e10 + 2 * .e25) * .e7/.e20) * .e4) - 2 * (v3 * .e22)) * 
        .e3), v2 = c(v1 = -(v3 * (.e6 * (2 - 2 * (.e27 * .e17^2/.e11)) - 
        .e2 * (2 * (2 - 4 * (.e27/.e7)) + 6) * .e4/.e7) * .e12 * 
        .e3/.e11), v2 = -(2/v2^3 + v2 * .e18 * ((2 - 2 * (.e25/.e10)) * 
        .e15 + .e2 * (4 * .e6 + 4 * .e17) - .e13) * .e12 * .e4/.e10))))
}
############################################################
#' The first derivative of the density
#' @returns Vector
#' @inheritParams manf
lst_k3_f1fa=function(x,v1,v2,kdf){
	vf=Vectorize(lst_k3_fd)
	f1=vf(x,v1,v2,kdf)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @returns Matrix
#' @inheritParams manf
lst_k3_f2fa=function(x,v1,v2,kdf){
	nx=length(x)
	vf=Vectorize(lst_k3_fdd)
	temp1=vf(x,v1,v2,kdf)
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}
############################################################
############################################################
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
lst_k3_ldda=function(x,v1,v2,kdf){
	nx=length(x)
	vf=Vectorize(lst_k3_logfdd)
	temp1=vf(x,v1,v2,kdf)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
lst_k3_lddda=function(x,v1,v2,kdf){
	nx=length(x)
	vf=Vectorize(lst_k3_logfddd)
	temp1=vf(x,v1,v2,kdf)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}
