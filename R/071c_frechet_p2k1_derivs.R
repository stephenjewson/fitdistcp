######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
frechet_p2k1_fd=function (x, t, v1, v2, v3, v4) 
{
    .e2 <- t * v3 + v2
    .e3 <- exp(.e2)
    .e4 <- x - v1
    .e5 <- .e4/.e3
    .e6 <- -1 - v4
    .e7 <- .e5^.e6
    .e9 <- 2 * v4
    .e10 <- exp(-.e5^-v4)
    .e15 <- (.e6/.e5^(2 + v4) + v4/.e5^(2 + .e9)) * .e4/.e3 + 
        .e7
    c(v2 = -(v4 * .e15 * .e10/.e3), v3 = -(t * v4 * .e15 * .e10/.e3), 
        v4 = (.e7 + v4 * (1/.e5^(1 + .e9) - .e7) * (log(.e4) - 
            .e2)) * .e10/.e3)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
frechet_p2k1_fdd=function (x, t, v1, v2, v3, v4) 
{
    .e2 <- t * v3 + v2
    .e3 <- x - v1
    .e4 <- exp(.e2)
    .e5 <- .e3/.e4
    .e6 <- -1 - v4
    .e7 <- 2 * v4
    .e8 <- 2 + v4
    .e9 <- 2 + .e7
    .e10 <- .e5^.e6
    .e11 <- .e5^.e8
    .e12 <- .e5^.e9
    .e13 <- 1 + .e7
    .e15 <- log(.e3) - .e2
    .e16 <- .e6/.e11
    .e17 <- 1 + v4
    .e18 <- 1/.e12
    .e20 <- exp(-.e5^-v4)
    .e24 <- (.e16 + v4/.e12) * .e3/.e4 + .e10
    .e26 <- 1/.e5^.e13
    .e27 <- .e5^.e17
    .e28 <- 1/.e11
    .e29 <- 2 * .e9
    .e30 <- 2 * .e8
    .e31 <- .e26 - .e10
    .e42 <- ((.e6 * .e5^(.e17 - .e30) * .e8 + v4 * .e5^(.e13 - 
        .e29) * .e9) * .e3/.e4 - (.e6 * (.e28 + 2/.e11) + v4 * 
        (.e24/.e27 + .e18 + .e18))) * .e3/.e4 - .e10
    .e43 <- .e10 + v4 * .e31 * .e15
    .e44 <- .e5^v4
    .e45 <- 2 * .e13
    .e46 <- -(t * v4 * .e42 * .e20/.e4)
    .e57 <- .e24 + v4 * ((.e24/.e44 - .e10) * .e15 + (.e18 - 
        ((.e6 * .e5^(.e8 - .e30) + 2 * (v4 * .e5^(.e9 - .e29))) * 
            .e15 + .e28)) * .e3/.e4)
    .e63 <- v4 * (((.e16 + .e5^(.e7 - .e45) * .e13) * .e15 - 
        .e43/.e27) * .e3/.e4 + .e10 - (.e31 * .e15 + .e26)) - 
        (.e6 * .e3/(.e11 * .e4) + .e10)
    c(v2 = c(v2 = -(v4 * .e42 * .e20/.e4), v3 = .e46, v4 = .e20 * 
        .e63/.e4), v3 = c(v2 = .e46, v3 = -(t^2 * v4 * .e42 * 
        .e20/.e4), v4 = t * .e20 * .e63/.e4), v4 = c(v2 = -(.e57 * 
        .e20/.e4), v3 = -(t * .e57 * .e20/.e4), v4 = (.e43/.e44 + 
        .e26 + v4 * (.e10 - 2 * .e5^(.e13 - .e45)) * .e15 - 2 * 
        .e10) * .e20 * .e15/.e4))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
frechet_p2k1_pd=function (x, t, v1, v2, v3, v4) 
{
    .e2 <- t * v3 + v2
    .e3 <- x - v1
    .e4 <- exp(.e2)
    .e5 <- .e3/.e4
    .e7 <- exp(-.e5^-v4)
    .e9 <- .e5^(1 + v4) * .e4
    c(v2 = -(v4 * .e7 * .e3/.e9), v3 = -(t * v4 * .e7 * .e3/.e9), 
        v4 = .e7 * (log(.e3) - .e2)/.e5^v4)
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
frechet_p2k1_pdd=function (x, t, v1, v2, v3, v4) 
{
    .e2 <- t * v3 + v2
    .e3 <- x - v1
    .e4 <- exp(.e2)
    .e5 <- .e3/.e4
    .e6 <- 1 + v4
    .e7 <- .e5^.e6
    .e8 <- .e7 * .e4
    .e10 <- .e5^v4
    .e11 <- exp(-.e5^-v4)
    .e13 <- log(.e3) - .e2
    .e14 <- .e8^2
    .e15 <- v4 * .e13
    .e17 <- (.e8 - .e10 * .e6 * .e3)/.e14 + v4 * .e3/(.e5^(2 * 
        .e6) * .e4^2)
    .e19 <- .e15 * .e3/.e8
    .e22 <- (1 + .e15/.e10)/.e8 - v4 * .e7 * .e4 * .e13/.e14
    .e26 <- t * v4 * .e17 * .e11 * .e3
    .e27 <- .e19 - (1 + .e19)/.e10
    c(v2 = c(v2 = v4 * .e17 * .e11 * .e3, v3 = .e26, v4 = .e11 * 
        .e27), v3 = c(v2 = .e26, v3 = t^2 * v4 * .e17 * .e11 * 
        .e3, v4 = t * .e11 * .e27), v4 = c(v2 = -(.e22 * .e11 * 
        .e3), v3 = -(t * .e22 * .e11 * .e3), v4 = (1/.e5^(2 * 
        v4) - 1/.e10) * .e11 * .e13^2))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
frechet_p2k1_logfdd=function (x, t, v1, v2, v3, v4) 
{
    .e2 <- t * v3 + v2
    .e3 <- x - v1
    .e5 <- exp(-.e2)
    .e6 <- 1 + v4
    .e7 <- .e5 * .e3
    .e8 <- .e7^.e6
    .e9 <- 1/.e8
    .e10 <- 2 * .e6
    .e12 <- log(.e3) - .e2
    .e16 <- .e6 * .e5 * .e7^(v4 - .e10) * .e3 - .e9
    .e17 <- .e7^v4
    .e18 <- -(t * v4 * .e16 * .e5 * .e3)
    .e20 <- 1 + .e5 * (v4 * .e7^(.e6 - .e10) * .e12 - .e9) * 
        .e3
    .e22 <- 1 + v4 * .e5 * .e12 * .e3/.e8 - 1/.e17
    c(v2 = c(v2 = -(v4 * .e16 * .e5 * .e3), v3 = .e18, v4 = .e22), 
        v3 = c(v2 = .e18, v3 = -(t^2 * v4 * .e16 * .e5 * .e3), 
            v4 = t * .e22), v4 = c(v2 = .e20, v3 = t * .e20, 
            v4 = -(.e12^2/.e17 + 1/v4^2)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
frechet_p2k1_logfddd=function (x, t, v1, v2, v3, v4) 
{
    .e2 <- t * v3 + v2
    .e3 <- x - v1
    .e5 <- exp(-.e2)
    .e6 <- 1 + v4
    .e7 <- .e5 * .e3
    .e8 <- 2 * .e6
    .e9 <- v4 - .e8
    .e10 <- .e7^.e9
    .e11 <- log(.e3)
    .e12 <- .e11 - .e2
    .e13 <- .e7^.e6
    .e14 <- .e6 - .e8
    .e15 <- 1/.e13
    .e16 <- .e7^.e14
    .e18 <- .e6 * .e5 * .e10
    .e22 <- (.e10 + 2 * .e10 + .e5 * .e7^(v4 - (1 + .e8)) * .e9 * 
        .e3) * .e6 * .e5 * .e3 - .e15
    .e23 <- .e18 * .e3
    .e24 <- t^2
    .e25 <- .e7^v4
    .e26 <- t * v4
    .e36 <- .e18 * .e12 * .e3 - ((1 + .e11 - .e2)/.e13 + .e15)
    .e38 <- .e23 + v4 * ((.e14 * .e5 * .e10 * .e3 + .e16) * .e12 + 
        .e16) - .e15
    .e40 <- .e23 + v4 * ((.e10 - .e6 * .e10 * .e12) * .e5 * .e3 + 
        .e16 * .e12) - .e15
    .e41 <- .e24 * v4
    .e43 <- v4 * .e16 * .e12
    .e46 <- .e26 * .e22 * .e5 * .e3
    .e49 <- .e41 * .e22 * .e5 * .e3
    .e50 <- -(t * .e38 * .e5 * .e3)
    .e51 <- -(t * .e40 * .e5 * .e3)
    .e54 <- (.e15 - .e43) * .e5 * .e3 + 1/.e25
    .e56 <- 2 * .e16 - .e43
    .e58 <- c(v2 = .e46, v3 = .e49, v4 = .e26 * .e36 * .e5 * 
        .e3)
    .e63 <- v4 * .e5 * .e12 * .e3/.e13 - 2/.e25
    c(v2 = c(v2 = c(v2 = v4 * .e22 * .e5 * .e3, v3 = .e46, v4 = v4 * 
        .e36 * .e5 * .e3), v3 = .e58, v4 = c(v2 = -(.e38 * .e5 * 
        .e3), v3 = .e50, v4 = -(.e12 * .e63))), v3 = c(v2 = .e58, 
        v3 = c(v2 = .e49, v3 = t^3 * v4 * .e22 * .e5 * .e3, v4 = .e41 * 
            .e36 * .e5 * .e3), v4 = c(v2 = .e50, v3 = -(.e24 * 
            .e38 * .e5 * .e3), v4 = -(t * .e12 * .e63))), v4 = c(v2 = c(v2 = -(.e40 * 
        .e5 * .e3), v3 = .e51, v4 = .e54 * .e12), v3 = c(v2 = .e51, 
        v3 = -(.e24 * .e40 * .e5 * .e3), v4 = t * .e54 * .e12), 
        v4 = c(v2 = .e56 * .e5 * .e12 * .e3, v3 = t * .e56 * 
            .e5 * .e12 * .e3, v4 = .e12^3/.e25 + 2/v4^3)))
}
############################################################
#' The first derivative of the density for DMGS
#' @returns Vector
#' @inheritParams manf
frechet_p2k1_f1fa=function(x,t0,v1,v2,v3,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	vf=Vectorize(frechet_p2k1_fd,"x")
	f1=vf(x,t0,kloc,v1,v2,v3)
	return(f1)
}
############################################################
#' The first derivative of the density for WAIC
#' @returns Vector
#' @inheritParams manf
frechet_p2k1_f1fw=function(x,t,v1,v2,v3,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	vf=Vectorize(frechet_p2k1_fd,c("x","t"))
	f1=vf(x,t,kloc,v1,v2,v3)
	return(f1)
}
############################################################
#' The second derivative of the density for DMGS
#' @returns Matrix
#' @inheritParams manf
frechet_p2k1_f2fa=function(x,t0,v1,v2,v3,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	nx=length(x)
	vf=Vectorize(frechet_p2k1_fdd,"x")
	temp1=vf(x,t0,kloc,v1,v2,v3)
	f2=deriv_copyfdd(temp1,nx,dim=3)
	return(f2)
}
############################################################
#' The second derivative of the density for WAIC
#' @returns Matrix
#' @inheritParams manf
frechet_p2k1_f2fw=function(x,t,v1,v2,v3,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	nx=length(x)
	vf=Vectorize(frechet_p2k1_fdd,c("x","t"))
	temp1=vf(x,t,kloc,v1,v2,v3)
	f2=deriv_copyfdd(temp1,nx,dim=3)
	return(f2)
}
############################################################
#' The first derivative of the cdf
#' @returns Vector
#' @inheritParams manf
frechet_p2k1_p1fa=function(x,t0,v1,v2,v3,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	vf=Vectorize(frechet_p2k1_pd,"x")
	p1=vf(x,t0,kloc,v1,v2,v3)
	return(p1)
}
############################################################
#' The second derivative of the cdf
#' @returns Matrix
#' @inheritParams manf
frechet_p2k1_p2fa=function(x,t0,v1,v2,v3,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	nx=length(x)
	vf=Vectorize(frechet_p2k1_pdd,"x")
	temp1=vf(x,t0,kloc,v1,v2,v3)
	p2=deriv_copyfdd(temp1,nx,dim=3)
	return(p2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @returns Vector
#' @inheritParams manf
frechet_p2k1_mu1fa=function(alpha,t0,v1,v2,v3,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	x=qfrechet((1-alpha),mu=kloc,sigma=exp(v1+v2*t0),lambda=v3)
#	x=qfrechet((1-alpha),mu=kloc,sigma=exp(v2+v3*t0),lambda=vkloc)
	vf=Vectorize(frechet_p2k1_pd,"x")
	mu1=-vf(x,t0,kloc,v1,v2,v3)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @returns Matrix
#' @inheritParams manf
frechet_p2k1_mu2fa=function(alpha,t0,v1,v2,v3,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	x=qfrechet((1-alpha),mu=kloc,sigma=exp(v1+v2*t0),lambda=v3)
	nx=length(x)
	vf=Vectorize(frechet_p2k1_pdd,"x")
	temp1=vf(x,t0,kloc,v1,v2,v3)
	mu2=-deriv_copyfdd(temp1,nx,dim=3)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
frechet_p2k1_ldda=function(x,t,v1,v2,v3,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	nx=length(x)
	vf=Vectorize(frechet_p2k1_logfdd,c("x","t"))
	temp1=vf(x,t,kloc,v1,v2,v3)
	ldd=deriv_copyldd(temp1,nx,dim=3)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
frechet_p2k1_lddda=function(x,t,v1,v2,v3,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	nx=length(x)
	vf=Vectorize(frechet_p2k1_logfddd,c("x","t"))
	temp1=vf(x,t,kloc,v1,v2,v3) #these are in mu, sigma, lambda order
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}
