######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
weibull_p2_fd=function (x, t, v1, v2, v3) 
{
    .e2 <- t * v3 + v2
    .e3 <- exp(.e2)
    .e4 <- x/.e3
    .e5 <- v1 - 1
    .e6 <- .e4^.e5
    .e8 <- exp(-.e4^v1)
    .e15 <- x * (v1 * .e4^(2 * .e5) - .e5 * .e4^(v1 - 2))/.e3 - 
        .e6
    c(v1 = (.e6 + v1 * (.e6 - .e4^(2 * v1 - 1)) * (log(x) - .e2)) * 
        .e8/.e3, v2 = v1 * .e8 * .e15/.e3, v3 = t * v1 * .e8 * 
        .e15/.e3)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
weibull_p2_fdd=function (x, t, v1, v2, v3) 
{
    .e2 <- t * v3 + v2
    .e3 <- exp(.e2)
    .e4 <- x/.e3
    .e5 <- v1 - 1
    .e6 <- .e4^.e5
    .e7 <- v1 - 2
    .e8 <- .e4^.e7
    .e9 <- 2 * .e5
    .e11 <- log(x) - .e2
    .e12 <- .e4^.e9
    .e13 <- 2 * v1
    .e14 <- .e4^v1
    .e15 <- .e13 - 1
    .e16 <- .e5 * .e8
    .e18 <- .e4^.e15
    .e19 <- exp(-.e14)
    .e20 <- v1 * .e12
    .e24 <- x * (.e20 - .e16)/.e3 - .e6
    .e25 <- .e6 - .e18
    .e28 <- .e6 + x * ((.e8 + 2 * .e8 + x * (.e7 * .e4^(v1 - 
        3) - 2 * (v1 * .e4^(.e9 - 1)))/.e3) * .e5 + v1 * (.e24 * 
        .e6 - 2 * .e12))/.e3
    .e29 <- .e6 + v1 * .e25 * .e11
    .e35 <- (1 - v1 * .e11 * .e14) * .e24 + v1 * (x * ((2 * .e20 - 
        .e16) * .e11 + .e12 - .e8)/.e3 - .e11 * .e6)
    .e42 <- t * v1 * .e28 * .e19/.e3
    .e44 <- v1 * (.e18 + x * ((.e15 * .e4^(.e13 - 2) - .e16) * 
        .e11 + .e29 * .e6)/.e3 - (.e25 * .e11 + .e6)) - (.e6 + 
        x * .e5 * .e8/.e3)
    c(v1 = c(v1 = (2 * .e6 + v1 * (.e6 - 2 * .e18) * .e11 - (.e29 * 
        .e14 + .e18)) * .e19 * .e11/.e3, v2 = .e35 * .e19/.e3, 
        v3 = t * .e35 * .e19/.e3), v2 = c(v1 = .e19 * .e44/.e3, 
        v2 = v1 * .e28 * .e19/.e3, v3 = .e42), v3 = c(v1 = t * 
        .e19 * .e44/.e3, v2 = .e42, v3 = t^2 * v1 * .e28 * .e19/.e3))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
weibull_p2_pd=function (x, t, v1, v2, v3) 
{
    .e2 <- t * v3 + v2
    .e4 <- exp(-.e2)
    .e5 <- x * .e4
    .e6 <- .e5^v1
    .e8 <- exp(-.e6)
    .e9 <- .e5^(v1 - 1)
    c(v1 = .e8 * (log(x) - .e2) * .e6, v2 = -(v1 * x * .e4 * 
        .e8 * .e9), v3 = -(t * v1 * x * .e4 * .e8 * .e9))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
weibull_p2_pdd=function (x, t, v1, v2, v3) 
{
    .e2 <- t * v3 + v2
    .e4 <- exp(-.e2)
    .e5 <- x * .e4
    .e6 <- v1 - 1
    .e7 <- .e5^.e6
    .e8 <- .e5^v1
    .e10 <- exp(-.e8)
    .e12 <- log(x) - .e2
    .e13 <- v1 * x
    .e14 <- .e13 * .e4
    .e16 <- (.e14 * .e7 - 1) * .e7 - .e5 * .e6 * .e5^(v1 - 2)
    .e17 <- v1 * .e12
    .e19 <- .e14 * .e12 * .e7
    .e20 <- -(t * v1 * x * .e16 * .e4 * .e10)
    .e22 <- (1 - .e17 * .e8) * .e7 + .e17 * .e7
    .e24 <- (.e19 - 1) * .e8 - .e19
    c(v1 = c(v1 = (.e8 - .e5^(2 * v1)) * .e10 * .e12^2, v2 = -(x * 
        .e22 * .e4 * .e10), v3 = -(t * x * .e22 * .e4 * .e10)), 
        v2 = c(v1 = .e24 * .e10, v2 = -(.e13 * .e16 * .e4 * .e10), 
            v3 = .e20), v3 = c(v1 = t * .e24 * .e10, v2 = .e20, 
            v3 = -(t^2 * v1 * x * .e16 * .e4 * .e10)))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
weibull_p2_logfdd=function (x, t, v1, v2, v3) 
{
    .e2 <- t * v3 + v2
    .e3 <- exp(.e2)
    .e4 <- x/.e3
    .e5 <- v1 - 1
    .e6 <- .e4^.e5
    .e8 <- log(x) - .e2
    .e9 <- .e6 + x * .e5 * .e4^(v1 - 2)/.e3
    .e10 <- .e4^v1
    .e11 <- v1 * x
    .e12 <- -(t * v1 * x * .e9/.e3)
    .e15 <- 1 - (.e10 + .e11 * .e8 * .e6/.e3)
    .e18 <- x * (.e6 + v1 * .e8 * .e6)/.e3 - 1
    c(v1 = c(v1 = -(.e8^2 * .e10 + 1/v1^2), v2 = .e18, v3 = t * 
        .e18), v2 = c(v1 = -.e15, v2 = -(.e11 * .e9/.e3), v3 = .e12), 
        v3 = c(v1 = -(t * .e15), v2 = .e12, v3 = -(t^2 * v1 * 
            x * .e9/.e3)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
weibull_p2_logfddd=function (x, t, v1, v2, v3) 
{
    .e2 <- t * v3 + v2
    .e3 <- exp(.e2)
    .e4 <- x/.e3
    .e5 <- v1 - 1
    .e6 <- .e4^.e5
    .e7 <- v1 - 2
    .e8 <- .e4^.e7
    .e10 <- log(x) - .e2
    .e13 <- x * .e5 * .e8/.e3
    .e15 <- (.e6 + .e13) * .e10 + .e6
    .e16 <- .e6 + x * (.e8 + 2 * .e8 + x * .e7 * .e4^(v1 - 3)/.e3) * 
        .e5/.e3
    .e17 <- t^2
    .e18 <- .e4^v1
    .e20 <- t * v1 * x
    .e21 <- t * x
    .e22 <- .e15 + .e6
    .e29 <- .e6 + v1 * .e15 + .e13
    .e31 <- .e6 + v1 * (.e10 * .e6 + x * (.e10 * .e5 * .e8 + 
        .e8)/.e3) + .e13
    .e33 <- .e17 * v1 * x
    .e35 <- v1 * .e10 * .e6
    .e36 <- v1 * x
    .e38 <- .e20 * .e16/.e3
    .e40 <- .e33 * .e16/.e3
    .e42 <- -(.e21 * .e29/.e3)
    .e43 <- -(.e21 * .e31/.e3)
    .e45 <- .e18 + x * (.e6 + .e35)/.e3
    .e47 <- 2 * .e6 + .e35
    .e49 <- 2 * .e18 + .e36 * .e10 * .e6/.e3
    .e50 <- c(v1 = -(.e20 * .e22/.e3), v2 = .e38, v3 = .e40)
    .e51 <- .e17 * x
    c(v1 = c(v1 = c(v1 = -(.e10^3 * .e18 - 2/v1^3), v2 = x * 
        .e47 * .e10/.e3, v3 = .e21 * .e47 * .e10/.e3), v2 = c(v1 = .e45 * 
        .e10, v2 = -(x * .e31/.e3), v3 = .e43), v3 = c(v1 = t * 
        .e45 * .e10, v2 = .e43, v3 = -(.e51 * .e31/.e3))), v2 = c(v1 = c(v1 = .e49 * 
        .e10, v2 = -(x * .e29/.e3), v3 = .e42), v2 = c(v1 = -(.e36 * 
        .e22/.e3), v2 = .e36 * .e16/.e3, v3 = .e38), v3 = .e50), 
        v3 = c(v1 = c(v1 = t * .e49 * .e10, v2 = .e42, v3 = -(.e51 * 
            .e29/.e3)), v2 = .e50, v3 = c(v1 = -(.e33 * .e22/.e3), 
            v2 = .e40, v3 = t^3 * v1 * x * .e16/.e3)))
}
############################################################
#' The first derivative of the density
#' @inheritParams manf
weibull_p2_f1fa=function(x,t,v1,v2,v3){
	vf=Vectorize(weibull_p2_fd)
	f1=vf(x,t,v1,v2,v3)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @inheritParams manf
weibull_p2_f2fa=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(weibull_p2_fdd)
	temp1=vf(x,t,v1,v2,v3)
	f2=deriv_copyfdd(temp1,nx,dim=3)
	return(f2)
}
############################################################
#' The first derivative of the cdf
#' @inheritParams manf
weibull_p2_p1fa=function(x,t,v1,v2,v3){
	vf=Vectorize(weibull_p2_pd)
	p1=vf(x,t,v1,v2,v3)
	return(p1)
}
############################################################
#' The second derivative of the cdf
#' @inheritParams manf
weibull_p2_p2fa=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(weibull_p2_pdd)
	temp1=vf(x,t,v1,v2,v3)
	p2=deriv_copyfdd(temp1,nx,dim=3)
	return(p2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @inheritParams manf
weibull_p2_mu1fa=function(alpha,t,v1,v2,v3){
	x=qweibull((1-alpha),shape=v1,scale=exp(v2+v3*t))
	vf=Vectorize(weibull_p2_pd)
	mu1=-vf(x,t,v1,v2,v3)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @inheritParams manf
weibull_p2_mu2fa=function(alpha,t,v1,v2,v3){
	x=qweibull((1-alpha),shape=v1,scale=exp(v2+v3*t))
	nx=length(x)
	vf=Vectorize(weibull_p2_pdd)
	temp1=vf(x,t,v1,v2,v3)
	mu2=-deriv_copyfdd(temp1,nx,dim=3)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @inheritParams manf
weibull_p2_ldda=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(weibull_p2_logfdd)
	temp1=vf(x,t,v1,v2,v3)
	ldd=deriv_copyldd(temp1,nx,dim=3)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @inheritParams manf
weibull_p2_lddda=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(weibull_p2_logfddd)
	temp1=vf(x,t,v1,v2,v3)
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}
