######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
weibull_fd=function (x, v1, v2) 
{
    .e1 <- x/v2
    .e2 <- v1 - 1
    .e3 <- .e1^.e2
    .e5 <- exp(-.e1^v1)
    c(v1 = (.e3 + v1 * (.e3 - .e1^(2 * v1 - 1)) * (log(x) - log(v2))) * 
        .e5/v2, v2 = v1 * .e5 * (x * (v1 * .e1^(2 * .e2) - .e2 * 
        .e1^(v1 - 2))/v2 - .e3)/v2^2)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
weibull_fdd=function (x, v1, v2) 
{
    .e1 <- x/v2
    .e2 <- v1 - 1
    .e3 <- .e1^.e2
    .e6 <- log(x) - log(v2)
    .e7 <- v1 - 2
    .e8 <- .e1^.e7
    .e9 <- 2 * v1
    .e10 <- 2 * .e2
    .e11 <- .e9 - 1
    .e12 <- .e1^.e10
    .e13 <- .e1^.e11
    .e14 <- .e1^v1
    .e15 <- .e2 * .e8
    .e17 <- exp(-.e14)
    .e18 <- v1 * .e12
    .e19 <- .e3 - .e13
    .e23 <- x * (.e18 - .e15)/v2 - .e3
    .e24 <- .e3 + v1 * .e19 * .e6
    .e25 <- v2^2
    c(v1 = c(v1 = (2 * .e3 + v1 * (.e3 - 2 * .e13) * .e6 - (.e24 * 
        .e14 + .e13)) * .e17 * .e6/v2, v2 = ((1 - v1 * .e6 * 
        .e14) * .e23 + v1 * (x * ((2 * .e18 - .e15) * .e6 + .e12 - 
        .e8)/v2 - .e6 * .e3)) * .e17/.e25), v2 = c(v1 = .e17 * 
        (v1 * (.e13 + x * ((.e11 * .e1^(.e9 - 2) - .e15) * .e6 + 
            .e24 * .e3)/v2 - (.e19 * .e6 + .e3)) - (.e3 + x * 
            .e2 * .e8/v2))/.e25, v2 = v1 * .e17 * (x * ((2 * 
        .e8 + x * (.e7 * .e1^(v1 - 3) - 2 * (v1 * .e1^(.e10 - 
        1)))/v2) * .e2 + v1 * (.e23 * .e3 - .e12))/v2 - 2 * .e23)/v2^3))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
weibull_pd=function (x, v1, v2) 
{
    .e1 <- x/v2
    .e2 <- .e1^v1
    .e4 <- exp(-.e2)
    c(v1 = .e4 * (log(x) - log(v2)) * .e2, v2 = -(v1 * x * .e4 * 
        .e1^(v1 - 1)/v2^2))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
weibull_pdd=function (x, v1, v2) 
{
    .e1 <- x/v2
    .e2 <- .e1^v1
    .e3 <- v1 - 1
    .e4 <- .e1^.e3
    .e7 <- log(x) - log(v2)
    .e9 <- exp(-.e2)
    .e10 <- v1 * x
    .e11 <- v1 * .e7
    .e14 <- .e10 * .e7 * .e4/v2
    c(v1 = c(v1 = (.e2 - .e1^(2 * v1)) * .e9 * .e7^2, v2 = -(x * 
        ((1 - .e11 * .e2) * .e4 + .e11 * .e4) * .e9/v2^2)), v2 = c(v1 = ((.e14 - 
        1) * .e2 - .e14) * .e9/v2, v2 = -(.e10 * .e9 * (x * (v1 * 
        .e1^(2 * .e3) - .e3 * .e1^(v1 - 2))/v2 - 2 * .e4)/v2^3)))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
weibull_logfdd=function (x, v1, v2) 
{
    .e1 <- x/v2
    .e2 <- v1 - 1
    .e3 <- .e1^.e2
    .e6 <- log(x) - log(v2)
    .e7 <- .e1^v1
    c(v1 = c(v1 = -(.e6^2 * .e7 + 1/v1^2), v2 = (x * (.e3 + v1 * 
        .e6 * .e3)/v2 - 1)/v2), v2 = c(v1 = -((1 - (.e7 + v1 * 
        x * .e6 * .e3/v2))/v2), v2 = -(v1 * (x * (2 * .e3 + x * 
        .e2 * .e1^(v1 - 2)/v2)/v2 - 1)/v2^2)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
weibull_logfddd=function (x, v1, v2) 
{
    .e1 <- x/v2
    .e2 <- v1 - 1
    .e3 <- .e1^.e2
    .e6 <- log(x) - log(v2)
    .e7 <- v1 - 2
    .e8 <- .e1^.e7
    .e9 <- 2 * .e3
    .e12 <- x * .e2 * .e8/v2
    .e13 <- .e1^v1
    .e14 <- v2^2
    .e15 <- .e9 + .e12
    .e17 <- .e15 * .e6 + .e3
    .e19 <- v1 * .e6 * .e3
    .e20 <- v1 * x
    c(v1 = c(v1 = c(v1 = -(.e6^3 * .e13 - 2/v1^3), v2 = x * (.e9 + 
        .e19) * .e6/.e14), v2 = c(v1 = (.e13 + x * (.e3 + .e19)/v2) * 
        .e6/v2, v2 = -((x * (.e9 + v1 * (2 * (.e6 * .e3) + x * 
        (.e6 * .e2 * .e8 + .e8)/v2) + .e12)/v2 - 1)/.e14))), 
        v2 = c(v1 = c(v1 = (2 * .e13 + .e20 * .e6 * .e3/v2) * 
            .e6/v2, v2 = -((x * (.e9 + v1 * .e17 + .e12)/v2 - 
            1)/.e14)), v2 = c(v1 = -((.e13 + .e20 * (.e17 + .e3)/v2 - 
            1)/.e14), v2 = v1 * (2 * (x * .e15/v2 - 1) + x * 
            (.e9 + x * (4 * .e8 + x * .e7 * .e1^(v1 - 3)/v2) * 
                .e2/v2)/v2)/v2^3)))
}
############################################################
#' The first derivative of the density
#' @returns Vector
#' @inheritParams manf
weibull_f1fa=function(x,v1,v2){
	vf=Vectorize(weibull_fd,"x")
	f1=vf(x,v1,v2)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @returns Matrix
#' @inheritParams manf
weibull_f2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(weibull_fdd,"x")
	temp1=vf(x,v1,v2)
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}
############################################################
#' The first derivative of the cdf
#' @returns Vector
#' @inheritParams manf
weibull_p1fa=function(x,v1,v2){
	vf=Vectorize(weibull_pd,"x")
	p1=vf(x,v1,v2)
	return(p1)
}
############################################################
#' The second derivative of the cdf
#' @returns Matrix
#' @inheritParams manf
weibull_p2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(weibull_pdd,"x")
	temp1=vf(x,v1,v2)
	p2=deriv_copyfdd(temp1,nx,dim=2)
	return(p2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @returns Vector
#' @inheritParams manf
weibull_mu1fa=function(alpha,v1,v2){
	x=qweibull((1-alpha),shape=v1,scale=v2)
	vf=Vectorize(weibull_pd,"x")
	mu1=-vf(x,v1,v2)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @returns Matrix
#' @inheritParams manf
weibull_mu2fa=function(alpha,v1,v2){
	x=qweibull((1-alpha),shape=v1,scale=v2)
	nx=length(x)
	vf=Vectorize(weibull_pdd,"x")
	temp1=vf(x,v1,v2)
	mu2=-deriv_copyfdd(temp1,nx,dim=2)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
weibull_ldda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(weibull_logfdd,"x")
	temp1=vf(x,v1,v2)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
weibull_lddda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(weibull_logfddd,"x")
	temp1=vf(x,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}
