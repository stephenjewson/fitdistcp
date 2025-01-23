######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
norm_p1_fd=function (x, t, v1, v2, v3) 
{
    .e3 <- x - (t * v2 + v1)
    .e4 <- v3^2
    .e5 <- .e3^2
    .e6 <- 2 * .e4
    .e9 <- exp(-(.e5/.e6))
    .e10 <- sqrt(2 * pi)
    .e12 <- v3^3 * .e10
    c(v1 = .e9 * .e3/.e12, v2 = t * .e9 * .e3/.e12, v3 = (4 * 
        (.e5/.e6^2) - 1/.e4) * .e9/.e10)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
norm_p1_fdd=function (x, t, v1, v2, v3) 
{
    .e1 <- v3^2
    .e4 <- x - (t * v2 + v1)
    .e5 <- .e4^2
    .e6 <- 2 * .e1
    .e8 <- sqrt(2 * pi)
    .e10 <- .e6^2
    .e11 <- exp(-(.e5/.e6))
    .e12 <- v3^3
    .e13 <- .e12 * .e8
    .e15 <- .e5/.e1 - 1
    .e19 <- 4 * (.e5/.e10) - 1/.e1
    .e21 <- .e19/.e1 - 8/.e10
    .e26 <- 4 * (.e5/(.e1 * .e10 * .e8)) - 3 * (.e1 * .e8/.e13^2)
    .e29 <- t * .e15 * .e11/.e13
    c(v1 = c(v1 = .e15 * .e11/.e13, v2 = .e29, v3 = .e21 * .e11 * 
        .e4/.e8), v2 = c(v1 = .e29, v2 = t^2 * .e15 * .e11/.e13, 
        v3 = t * .e21 * .e11 * .e4/.e8), v3 = c(v1 = .e26 * .e11 * 
        .e4, v2 = t * .e26 * .e11 * .e4, v3 = (2/.e12 + v3 * 
        (4 * .e19 - 64 * (.e1/.e10)) * .e5/.e10) * .e11/.e8))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
norm_p1_pd=function (x, t, v1, v2, v3) 
{
    .e3 <- x - (t * v2 + v1)
    .e5 <- dnorm(.e3/v3, 0, 1)
    c(v1 = -(.e5/v3), v2 = -(t * .e5/v3), v3 = -(.e5 * .e3/v3^2))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
norm_p1_pdd=function (x, t, v1, v2, v3) 
{
    .e3 <- x - (t * v2 + v1)
    .e5 <- dnorm(.e3/v3, 0, 1)
    .e6 <- v3^2
    .e8 <- .e3^2/.e6
    .e9 <- v3^3
    .e10 <- .e8 - 1
    .e11 <- -(.e10 * .e5/.e6)
    .e12 <- -(t * .e10 * .e5/.e6)
    .e13 <- -(t * .e5 * .e3/.e9)
    c(v1 = c(v1 = -(.e5 * .e3/.e9), v2 = .e13, v3 = .e11), v2 = c(v1 = .e13, 
        v2 = -(t^2 * .e5 * .e3/.e9), v3 = .e12), v3 = c(v1 = .e11, 
        v2 = .e12, v3 = -((.e8 - 2) * .e5 * .e3/.e9)))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
norm_p1_logfdd=function (x, t, v1, v2, v3) 
{
    .e1 <- v3^2
    .e4 <- x - (t * v2 + v1)
    .e5 <- (2 * .e1)^2
    .e6 <- -(t/.e1)
    .e7 <- 1/.e1
    .e8 <- v3^3
    c(v1 = c(v1 = -.e7, v2 = .e6, v3 = -(8 * (v3 * .e4/.e5))), 
        v2 = c(v1 = .e6, v2 = -(t^2/.e1), v3 = -(8 * (t * v3 * 
            .e4/.e5))), v3 = c(v1 = -(2 * (.e4/.e8)), v2 = -(2 * 
            (t * .e4/.e8)), v3 = .e7 + 4 * ((1 - 16 * (v3^4/.e5)) * 
            .e4^2/.e5)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
norm_p1_logfddd=function (x, t, v1, v2, v3) 
{
    .e1 <- 2 * v3^2
    .e2 <- .e1^2
    .e3 <- v3^3
    .e4 <- v3^4
    .e7 <- x - (t * v2 + v1)
    .e8 <- 16 * (.e4/.e2)
    .e9 <- 1 - .e8
    .e10 <- 2 * (t/.e3)
    .e11 <- 2/.e3
    .e12 <- t^2
    .e17 <- c(v1 = 0, v2 = 0, v3 = 8 * (t * v3/.e2))
    .e18 <- c(v1 = .e10, v2 = 2 * (.e12/.e3), v3 = -(8 * (t * 
        .e9 * .e7/.e2)))
    .e19 <- c(v1 = .e11, v2 = .e10, v3 = -(8 * (.e9 * .e7/.e2)))
    c(v1 = c(v1 = c(v1 = 0, v2 = 0, v3 = 8 * (v3/.e2)), v2 = .e17, 
        v3 = .e19), v2 = c(v1 = .e17, v2 = c(v1 = 0, v2 = 0, 
        v3 = 8 * (.e12 * v3/.e2)), v3 = .e18), v3 = c(v1 = .e19, 
        v2 = .e18, v3 = c(v1 = 6 * (.e7/.e4), v2 = 6 * (t * .e7/.e4), 
            v3 = -(.e11 + 4 * (.e3 * (16 * .e9 + 16 * (4 - .e8)) * 
                .e7^2/.e1^4)))))
}
############################################################
#' The first derivative of the density
#' @inheritParams manf
norm_p1_f1fa=function(x,t,v1,v2,v3){
	vf=Vectorize(norm_p1_fd)
	f1=vf(x,t,v1,v2,v3)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @inheritParams manf
norm_p1_f2fa=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(norm_p1_fdd)
	temp1=vf(x,t,v1,v2,v3)
	f2=deriv_copyfdd(temp1,nx,dim=3)
	return(f2)
}
############################################################
#' The first derivative of the cdf
#' @inheritParams manf
norm_p1_p1fa=function(x,t,v1,v2,v3){
	vf=Vectorize(norm_p1_pd)
	p1=vf(x,t,v1,v2,v3)
	return(p1)
}
############################################################
#' The second derivative of the cdf
#' @inheritParams manf
norm_p1_p2fa=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(norm_p1_pdd)
	temp1=vf(x,t,v1,v2,v3)
	p2=deriv_copyfdd(temp1,nx,dim=3)
	return(p2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @inheritParams manf
norm_p1_mu1fa=function(alpha,t,v1,v2,v3){
	x=qnorm((1-alpha),mean=v1+v2*t,sd=v3)
	vf=Vectorize(norm_p1_pd)
	mu1=-vf(x,t,v1,v2,v3)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @inheritParams manf
norm_p1_mu2fa=function(alpha,t,v1,v2,v3){
	x=qnorm((1-alpha),mean=v1+v2*t,sd=v3)
	nx=length(x)
	vf=Vectorize(norm_p1_pdd)
	temp1=vf(x,t,v1,v2,v3)
	mu2=-deriv_copyfdd(temp1,nx,dim=3)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @inheritParams manf
norm_p1_ldda=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(norm_p1_logfdd)
	temp1=vf(x,t,v1,v2,v3)
	ldd=deriv_copyldd(temp1,nx,dim=3)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @inheritParams manf
norm_p1_lddda=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(norm_p1_logfddd)
	temp1=vf(x,t,v1,v2,v3)
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}
