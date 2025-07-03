######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
frechet_k1_fd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- .e1/v2
    .e3 <- -1 - v3
    .e4 <- .e2^.e3
    .e6 <- 2 * v3
    .e7 <- exp(-.e2^-v3)
    c(v2 = -(v3 * ((.e3/.e2^(2 + v3) + v3/.e2^(2 + .e6)) * .e1/v2 + 
        .e4) * .e7/v2^2), v3 = (.e4 + v3 * (1/.e2^(1 + .e6) - 
        .e4) * (log(.e1) - log(v2))) * .e7/v2)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
frechet_k1_fdd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- .e1/v2
    .e3 <- -1 - v3
    .e4 <- 2 * v3
    .e5 <- .e2^.e3
    .e6 <- 2 + v3
    .e7 <- 1 + .e4
    .e8 <- 2 + .e4
    .e9 <- .e2^.e6
    .e12 <- log(.e1) - log(v2)
    .e13 <- .e3/.e9
    .e14 <- .e2^.e8
    .e16 <- 1/.e2^.e7
    .e21 <- (.e13 + v3/.e14) * .e1/v2 + .e5
    .e22 <- exp(-.e2^-v3)
    .e23 <- 1 + v3
    .e24 <- .e16 - .e5
    .e25 <- .e5 + v3 * .e24 * .e12
    .e26 <- .e2^.e23
    .e27 <- .e2^v3
    .e28 <- 1/.e14
    .e29 <- 2 * .e7
    .e30 <- 2 * .e8
    .e31 <- 2 * .e6
    .e32 <- v2^2
    c(v2 = c(v2 = -(v3 * (((.e3 * .e2^(.e23 - .e31) * .e6 + v3 * 
        .e2^(.e7 - .e30) * .e8) * .e1/v2 - (2 * .e13 + v3 * (.e21/.e26 + 
        .e28))) * .e1/v2 - 2 * .e21) * .e22/v2^3), v3 = .e22 * 
        (v3 * (((.e13 + .e2^(.e4 - .e29) * .e7) * .e12 - .e25/.e26) * 
            .e1/v2 + .e5 - (.e24 * .e12 + .e16)) - (.e3 * .e1/(v2 * 
            .e9) + .e5))/.e32), v3 = c(v2 = -((.e21 + v3 * ((.e21/.e27 - 
        .e5) * .e12 + (.e28 - ((.e3 * .e2^(.e6 - .e31) + 2 * 
        (v3 * .e2^(.e8 - .e30))) * .e12 + 1/.e9)) * .e1/v2)) * 
        .e22/.e32), v3 = (.e25/.e27 + .e16 + v3 * (.e5 - 2 * 
        .e2^(.e7 - .e29)) * .e12 - 2 * .e5) * .e22 * .e12/v2))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
frechet_k1_pd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- .e1/v2
    .e4 <- exp(-.e2^-v3)
    c(v2 = -(v3 * .e4 * .e1/(v2^2 * .e2^(1 + v3))), v3 = .e4 * 
        (log(.e1) - log(v2))/.e2^v3)
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
frechet_k1_pdd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- .e1/v2
    .e3 <- 1 + v3
    .e4 <- .e2^.e3
    .e7 <- log(.e1) - log(v2)
    .e9 <- .e2^v3
    .e10 <- exp(-.e2^-v3)
    .e11 <- v2^2
    .e12 <- v2 * .e4
    .e13 <- .e11 * .e4
    .e14 <- v3 * .e7
    .e15 <- .e13^2
    .e17 <- .e14 * .e1/.e12
    c(v2 = c(v2 = v3 * ((2 * .e12 - .e9 * .e3 * .e1)/.e15 + v3 * 
        .e1/(v2^4 * .e2^(2 * .e3))) * .e10 * .e1, v3 = .e10 * 
        (.e17 - (1 + .e17)/.e9)/v2), v3 = c(v2 = -(((1 + .e14/.e9)/.e13 - 
        .e11 * v3 * .e4 * .e7/.e15) * .e10 * .e1), v3 = (1/.e2^(2 * 
        v3) - 1/.e9) * .e10 * .e7^2))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
frechet_k1_logfdd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- .e1/v2
    .e3 <- 1 + v3
    .e4 <- .e2^.e3
    .e5 <- v2 * .e4
    .e6 <- .e2^v3
    .e9 <- log(.e1) - log(v2)
    .e10 <- .e5^2
    c(v2 = c(v2 = v3 * ((.e4 - .e6 * .e3 * .e1/v2) * .e1/.e10 - 
        (1 - .e1/.e5)/v2)/v2, v3 = (1 + v3 * .e9 * .e1/.e5 - 
        1/.e6)/v2), v3 = c(v2 = ((v2 * v3 * .e4 * .e9/.e10 - 
        1/.e5) * .e1 + 1)/v2, v3 = -(.e9^2/.e6 + 1/v3^2)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
frechet_k1_logfddd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- .e1/v2
    .e3 <- 1 + v3
    .e4 <- .e2^.e3
    .e5 <- v2 * .e4
    .e6 <- .e2^v3
    .e9 <- log(.e1) - log(v2)
    .e10 <- .e5^2
    .e11 <- .e6 * .e3
    .e13 <- .e11 * .e1/v2
    .e14 <- .e4 - .e13
    .e15 <- v2^2
    .e16 <- 2 * .e3
    .e17 <- .e2^.e16
    .e18 <- (1 - .e1/.e5)/v2
    .e19 <- 1/.e6
    .e20 <- 1/.e5
    .e24 <- v2 * v3 * .e4 * .e9/.e10
    .e27 <- v3 * .e9 * .e1/.e5
    c(v2 = c(v2 = c(v2 = v3 * (((.e6 + v3 * .e2^(v3 - 1) * .e1/v2 - 
        .e6) * .e3 * .e1/.e15 - 2 * (v2 * .e14^2 * .e4/.e10)) * 
        .e1/.e10 - 2 * ((.e14 * .e1/.e10 - .e18)/v2))/v2, v3 = -(((1 + 
        .e27 - .e19)/v2 + v3 * (.e14 * .e9/.e10 + 2/(.e15 * .e4)) * 
        .e1)/v2)), v3 = c(v2 = ((.e4 + v3 * (.e14 * (1 - 2 * 
        (.e15 * .e17/.e10)) * .e9 - .e4) - .e13) * .e1/.e10 - 
        ((.e24 - .e20) * .e1 + 1)/v2)/v2, v3 = -(.e9 * (.e27 - 
        2/.e6)/v2))), v3 = c(v2 = c(v2 = ((.e4 + v3 * ((.e4 - 
        (.e4 + 2 * (.e15 * .e14 * .e17/.e10))) * .e9 - (.e6 + 
        .e11 * .e9) * .e1/v2) - .e13) * .e1/.e10 - .e18)/v2, 
        v3 = ((.e20 - .e24) * .e1 + .e19) * .e9/v2), v3 = c(v2 = (2 * 
        .e4 + v3 * (.e4 - 2 * (.e15 * .e2^(1 + .e16 + v3)/.e10)) * 
        .e9) * .e9 * .e1/.e10, v3 = .e9^3/.e6 + 2/v3^3)))
}
############################################################
#' The first derivative of the density
#' @returns Vector
#' @inheritParams manf
frechet_k1_f1fa=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	vf=Vectorize(frechet_k1_fd)
	f1=vf(x,kloc,v1,v2)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @returns Matrix
#' @inheritParams manf
frechet_k1_f2fa=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	nx=length(x)
	vf=Vectorize(frechet_k1_fdd)
	temp1=vf(x,kloc,v1,v2)
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}
############################################################
#' The first derivative of the cdf
#' @returns Vector
#' @inheritParams manf
frechet_k1_p1fa=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	vf=Vectorize(frechet_k1_pd)
	p1=vf(x,kloc,v1,v2)
	return(p1)
}
############################################################
#' The second derivative of the cdf
#' @returns Matrix
#' @inheritParams manf
frechet_k1_p2fa=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	nx=length(x)
	vf=Vectorize(frechet_k1_pdd)
	temp1=vf(x,kloc,v1,v2)
	p2=deriv_copyfdd(temp1,nx,dim=2)
	return(p2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @returns Vector
#' @inheritParams manf
frechet_k1_mu1fa=function(alpha,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	x=qfrechet((1-alpha),mu=kloc,sigma=v1,lambda=v2)
	vf=Vectorize(frechet_k1_pd)
	mu1=-vf(x,kloc,v1,v2)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @returns Matrix
#' @inheritParams manf
frechet_k1_mu2fa=function(alpha,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	x=qfrechet((1-alpha),mu=kloc,sigma=v1,lambda=v2)
	nx=length(x)
	vf=Vectorize(frechet_k1_pdd)
	temp1=vf(x,kloc,v1,v2)
	mu2=-deriv_copyfdd(temp1,nx,dim=2)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
frechet_k1_ldda=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	nx=length(x)
	vf=Vectorize(frechet_k1_logfdd)
	temp1=vf(x,kloc,v1,v2)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
frechet_k1_lddda=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	nx=length(x)
	vf=Vectorize(frechet_k1_logfddd)
	temp1=vf(x,kloc,v1,v2) #these are in mu, sigma, lambda order
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}
