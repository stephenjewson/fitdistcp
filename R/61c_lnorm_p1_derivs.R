######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
lnorm_p1_fd=function (x, t, v1, v2, v3) 
{
    .e2 <- log(x) - (t * v2 + v1)
    .e3 <- v3^2
    .e4 <- .e2^2
    .e5 <- 2 * .e3
    .e8 <- exp(-(.e4/.e5))
    .e9 <- sqrt(2 * pi)
    .e12 <- v3^3 * x * .e9
    c(v1 = .e8 * .e2/.e12, v2 = t * .e8 * .e2/.e12, v3 = (4 * 
        (.e4/.e5^2) - 1/.e3) * .e8/(x * .e9))
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
lnorm_p1_fdd=function (x, t, v1, v2, v3) 
{
    .e1 <- v3^2
    .e3 <- log(x) - (t * v2 + v1)
    .e4 <- .e3^2
    .e5 <- 2 * .e1
    .e7 <- sqrt(2 * pi)
    .e9 <- .e5^2
    .e10 <- exp(-(.e4/.e5))
    .e11 <- v3^3
    .e13 <- .e11 * x * .e7
    .e15 <- .e4/.e1 - 1
    .e16 <- .e1 * x
    .e20 <- 4 * (.e4/.e9) - 1/.e1
    .e21 <- x * .e7
    .e23 <- .e20/.e1 - 8/.e9
    .e28 <- 4 * (.e4/(.e16 * .e9 * .e7)) - 3 * (.e16 * .e7/.e13^2)
    .e31 <- t * .e15 * .e10/.e13
    c(v1 = c(v1 = .e15 * .e10/.e13, v2 = .e31, v3 = .e23 * .e10 * 
        .e3/.e21), v2 = c(v1 = .e31, v2 = t^2 * .e15 * .e10/.e13, 
        v3 = t * .e23 * .e10 * .e3/.e21), v3 = c(v1 = .e28 * 
        .e10 * .e3, v2 = t * .e28 * .e10 * .e3, v3 = (2/.e11 + 
        v3 * (4 * .e20 - 64 * (.e1/.e9)) * .e4/.e9) * .e10/.e21))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
lnorm_p1_pd=function (x, t, v1, v2, v3) 
{
    .e2 <- log(x) - (t * v2 + v1)
    .e4 <- dnorm(.e2/v3, 0, 1)
    c(v1 = -(.e4/v3), v2 = -(t * .e4/v3), v3 = -(.e4 * .e2/v3^2))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
lnorm_p1_pdd=function (x, t, v1, v2, v3) 
{
    .e2 <- log(x) - (t * v2 + v1)
    .e4 <- dnorm(.e2/v3, 0, 1)
    .e5 <- v3^2
    .e7 <- .e2^2/.e5
    .e8 <- v3^3
    .e9 <- .e7 - 1
    .e10 <- -(.e9 * .e4/.e5)
    .e11 <- -(t * .e9 * .e4/.e5)
    .e12 <- -(t * .e4 * .e2/.e8)
    c(v1 = c(v1 = -(.e4 * .e2/.e8), v2 = .e12, v3 = .e10), v2 = c(v1 = .e12, 
        v2 = -(t^2 * .e4 * .e2/.e8), v3 = .e11), v3 = c(v1 = .e10, 
        v2 = .e11, v3 = -((.e7 - 2) * .e4 * .e2/.e8)))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
lnorm_p1_logfdd=function (x, t, v1, v2, v3) 
{
    .e1 <- v3^2
    .e3 <- log(x) - (t * v2 + v1)
    .e4 <- (2 * .e1)^2
    .e5 <- -(t/.e1)
    .e6 <- 1/.e1
    .e7 <- v3^3
    c(v1 = c(v1 = -.e6, v2 = .e5, v3 = -(8 * (v3 * .e3/.e4))), 
        v2 = c(v1 = .e5, v2 = -(t^2/.e1), v3 = -(8 * (t * v3 * 
            .e3/.e4))), v3 = c(v1 = -(2 * (.e3/.e7)), v2 = -(2 * 
            (t * .e3/.e7)), v3 = .e6 + 4 * ((1 - 16 * (v3^4/.e4)) * 
            .e3^2/.e4)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
lnorm_p1_logfddd=function (x, t, v1, v2, v3) 
{
    .e1 <- 2 * v3^2
    .e2 <- .e1^2
    .e3 <- v3^3
    .e4 <- v3^4
    .e6 <- log(x) - (t * v2 + v1)
    .e7 <- 16 * (.e4/.e2)
    .e8 <- 1 - .e7
    .e9 <- 2 * (t/.e3)
    .e10 <- 2/.e3
    .e11 <- t^2
    .e16 <- c(v1 = 0, v2 = 0, v3 = 8 * (t * v3/.e2))
    .e17 <- c(v1 = .e9, v2 = 2 * (.e11/.e3), v3 = -(8 * (t * 
        .e8 * .e6/.e2)))
    .e18 <- c(v1 = .e10, v2 = .e9, v3 = -(8 * (.e8 * .e6/.e2)))
    c(v1 = c(v1 = c(v1 = 0, v2 = 0, v3 = 8 * (v3/.e2)), v2 = .e16, 
        v3 = .e18), v2 = c(v1 = .e16, v2 = c(v1 = 0, v2 = 0, 
        v3 = 8 * (.e11 * v3/.e2)), v3 = .e17), v3 = c(v1 = .e18, 
        v2 = .e17, v3 = c(v1 = 6 * (.e6/.e4), v2 = 6 * (t * .e6/.e4), 
            v3 = -(.e10 + 4 * (.e3 * (16 * .e8 + 16 * (4 - .e7)) * 
                .e6^2/.e1^4)))))
}
############################################################
#' The first derivative of the density
#' @returns Vector
#' @inheritParams manf
lnorm_p1_f1fa=function(x,t,v1,v2,v3){
	vf=Vectorize(lnorm_p1_fd)
	f1=vf(x,t,v1,v2,v3)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @returns Matrix
#' @inheritParams manf
lnorm_p1_f2fa=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(lnorm_p1_fdd)
	temp1=vf(x,t,v1,v2,v3)
	f2=deriv_copyfdd(temp1,nx,dim=3)
	return(f2)
}
############################################################
#' The first derivative of the cdf
#' @returns Vector
#' @inheritParams manf
lnorm_p1_p1fa=function(x,t,v1,v2,v3){
	vf=Vectorize(lnorm_p1_pd)
	p1=vf(x,t,v1,v2,v3)
	return(p1)
}
############################################################
#' The second derivative of the cdf
#' @returns Matrix
#' @inheritParams manf
lnorm_p1_p2fa=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(lnorm_p1_pdd)
	temp1=vf(x,t,v1,v2,v3)
	p2=deriv_copyfdd(temp1,nx,dim=3)
	return(p2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @returns Vector
#' @inheritParams manf
lnorm_p1_mu1fa=function(alpha,t,v1,v2,v3){
	x=qlnorm((1-alpha),meanlog=v1+v2*t,sdlog=v3)
	vf=Vectorize(lnorm_p1_pd)
	mu1=-vf(x,t,v1,v2,v3)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @returns Matrix
#' @inheritParams manf
lnorm_p1_mu2fa=function(alpha,t,v1,v2,v3){
	x=qlnorm((1-alpha),meanlog=v1+v2*t,sdlog=v3)
	nx=length(x)
	vf=Vectorize(lnorm_p1_pdd)
	temp1=vf(x,t,v1,v2,v3)
	mu2=-deriv_copyfdd(temp1,nx,dim=3)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
lnorm_p1_ldda=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(lnorm_p1_logfdd)
	temp1=vf(x,t,v1,v2,v3)
	ldd=deriv_copyldd(temp1,nx,dim=3)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
lnorm_p1_lddda=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(lnorm_p1_logfddd)
	temp1=vf(x,t,v1,v2,v3)
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}
