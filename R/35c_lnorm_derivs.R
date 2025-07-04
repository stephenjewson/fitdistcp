######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
lnorm_fd=function (x, v1, v2) 
{
    .e2 <- log(x) - v1
    .e3 <- v2^2
    .e4 <- .e2^2
    .e5 <- 2 * .e3
    .e8 <- exp(-(.e4/.e5))
    .e9 <- sqrt(2 * pi)
    c(v1 = .e8 * .e2/(v2^3 * x * .e9), v2 = (4 * (.e4/.e5^2) - 
        1/.e3) * .e8/(x * .e9))
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
lnorm_fdd=function (x, v1, v2) 
{
    .e1 <- v2^2
    .e3 <- log(x) - v1
    .e4 <- 2 * .e1
    .e5 <- .e3^2
    .e6 <- .e4^2
    .e8 <- sqrt(2 * pi)
    .e10 <- exp(-(.e5/.e4))
    .e11 <- v2^3
    .e15 <- 4 * (.e5/.e6) - 1/.e1
    .e16 <- .e1 * x
    .e18 <- .e11 * x * .e8
    .e19 <- x * .e8
    c(v1 = c(v1 = (.e5/.e1 - 1) * .e10/.e18, v2 = (.e15/.e1 - 
        8/.e6) * .e10 * .e3/.e19), v2 = c(v1 = (4 * (.e5/(.e16 * 
        .e6 * .e8)) - 3 * (.e16 * .e8/.e18^2)) * .e10 * .e3, 
        v2 = (2/.e11 + v2 * (4 * .e15 - 64 * (.e1/.e6)) * .e5/.e6) * 
            .e10/.e19))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
lnorm_pd=function (x, v1, v2) 
{
    .e2 <- log(x) - v1
    .e4 <- dnorm(.e2/v2, 0, 1)
    c(v1 = -(.e4/v2), v2 = -(.e4 * .e2/v2^2))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
lnorm_pdd=function (x, v1, v2) 
{
    .e2 <- log(x) - v1
    .e3 <- v2^2
    .e5 <- dnorm(.e2/v2, 0, 1)
    .e7 <- .e2^2/.e3
    .e8 <- -((.e7 - 1) * .e5/.e3)
    .e9 <- v2^3
    c(v1 = c(v1 = -(.e5 * .e2/.e9), v2 = .e8), v2 = c(v1 = .e8, 
        v2 = -((.e7 - 2) * .e5 * .e2/.e9)))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
lnorm_logfdd=function (x, v1, v2) 
{
    .e1 <- v2^2
    .e2 <- (2 * .e1)^2
    .e4 <- log(x) - v1
    .e5 <- 1/.e1
    c(v1 = c(v1 = -.e5, v2 = -(8 * (v2 * .e4/.e2))), v2 = c(v1 = -(2 * 
        (.e4/v2^3)), v2 = .e5 + 4 * ((1 - 16 * (v2^4/.e2)) * 
        .e4^2/.e2)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
lnorm_logfddd=function (x, v1, v2) 
{
    .e1 <- 2 * v2^2
    .e2 <- .e1^2
    .e3 <- v2^4
    .e4 <- 16 * (.e3/.e2)
    .e6 <- log(x) - v1
    .e7 <- v2^3
    .e8 <- 1 - .e4
    .e9 <- 2/.e7
    .e11 <- c(v1 = .e9, v2 = -(8 * (.e8 * .e6/.e2)))
    c(v1 = c(v1 = c(v1 = 0, v2 = 8 * (v2/.e2)), v2 = .e11), v2 = c(v1 = .e11, 
        v2 = c(v1 = 6 * (.e6/.e3), v2 = -(.e9 + 4 * (.e7 * (16 * 
            .e8 + 16 * (4 - .e4)) * .e6^2/.e1^4)))))
}
############################################################
#' The first derivative of the density
#' @returns Vector
#' @inheritParams manf
lnorm_f1fa=function(x,v1,v2){
	vf=Vectorize(lnorm_fd,"x")
	f1=vf(x,v1,v2)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @returns Matrix
#' @inheritParams manf
lnorm_f2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(lnorm_fdd,"x")
	temp1=vf(x,v1,v2)
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}
############################################################
#' The first derivative of the cdf
#' @returns Vector
#' @inheritParams manf
lnorm_p1fa=function(x,v1,v2){
	vf=Vectorize(lnorm_pd,"x")
	p1=vf(x,v1,v2)
	return(p1)
}
############################################################
#' The second derivative of the cdf
#' @returns Matrix
#' @inheritParams manf
lnorm_p2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(lnorm_pdd,"x")
	temp1=vf(x,v1,v2)
	p2=deriv_copyfdd(temp1,nx,dim=2)
	return(p2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @returns Vector
#' @inheritParams manf
lnorm_mu1fa=function(alpha,v1,v2){
	x=qlnorm((1-alpha),meanlog=v1,sdlog=v2)
	vf=Vectorize(lnorm_pd,"x")
	mu1=-vf(x,v1,v2)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @returns Matrix
#' @inheritParams manf
lnorm_mu2fa=function(alpha,v1,v2){
	x=qlnorm((1-alpha),meanlog=v1,sdlog=v2)
	nx=length(x)
	vf=Vectorize(lnorm_pdd,"x")
	temp1=vf(x,v1,v2)
	mu2=-deriv_copyfdd(temp1,nx,dim=2)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
lnorm_ldda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(lnorm_logfdd,"x")
	temp1=vf(x,v1,v2)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
lnorm_lddda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(lnorm_logfddd,"x")
	temp1=vf(x,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}
