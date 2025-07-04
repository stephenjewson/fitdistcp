######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
norm_fd=function (x, v1, v2) 
{
    .e1 <- v2^2
    .e2 <- x - v1
    .e3 <- .e2^2
    .e4 <- 2 * .e1
    .e7 <- exp(-(.e3/.e4))
    .e8 <- sqrt(2 * pi)
    c(v1 = .e7 * .e2/(v2^3 * .e8), v2 = (4 * (.e3/.e4^2) - 1/.e1) * 
        .e7/.e8)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
norm_fdd=function (x, v1, v2) 
{
    .e1 <- v2^2
    .e2 <- x - v1
    .e3 <- 2 * .e1
    .e4 <- .e2^2
    .e5 <- .e3^2
    .e7 <- sqrt(2 * pi)
    .e9 <- exp(-(.e4/.e3))
    .e10 <- v2^3
    .e14 <- 4 * (.e4/.e5) - 1/.e1
    .e15 <- .e10 * .e7
    c(v1 = c(v1 = (.e4/.e1 - 1) * .e9/.e15, v2 = (.e14/.e1 - 
        8/.e5) * .e9 * .e2/.e7), v2 = c(v1 = (4 * (.e4/(.e1 * 
        .e5 * .e7)) - 3 * (.e1 * .e7/.e15^2)) * .e9 * .e2, v2 = (2/.e10 + 
        v2 * (4 * .e14 - 64 * (.e1/.e5)) * .e4/.e5) * .e9/.e7))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
norm_pd=function (x, v1, v2) 
{
    .e1 <- x - v1
    .e3 <- dnorm(.e1/v2, 0, 1)
    c(v1 = -(.e3/v2), v2 = -(.e3 * .e1/v2^2))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
norm_pdd=function (x, v1, v2) 
{
    .e1 <- x - v1
    .e2 <- v2^2
    .e4 <- dnorm(.e1/v2, 0, 1)
    .e6 <- .e1^2/.e2
    .e7 <- -((.e6 - 1) * .e4/.e2)
    .e8 <- v2^3
    cbind(v1 = cbind(v1 = -(.e4 * .e1/.e8), v2 = .e7), v2 = cbind(v1 = .e7, 
        v2 = -((.e6 - 2) * .e4 * .e1/.e8)))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
norm_logfdd=function (x, v1, v2) 
{
    .e1 <- v2^2
    .e2 <- (2 * .e1)^2
    .e3 <- x - v1
    .e4 <- 1/.e1
    c(v1 = c(v1 = -.e4, v2 = -(8 * (v2 * .e3/.e2))), v2 = c(v1 = -(2 * 
        (.e3/v2^3)), v2 = .e4 + 4 * ((1 - 16 * (v2^4/.e2)) * 
        .e3^2/.e2)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
norm_logfddd=function (x, v1, v2) 
{
    .e1 <- 2 * v2^2
    .e2 <- .e1^2
    .e3 <- v2^4
    .e4 <- 16 * (.e3/.e2)
    .e5 <- v2^3
    .e6 <- x - v1
    .e7 <- 1 - .e4
    .e8 <- 2/.e5
    .e10 <- c(v1 = .e8, v2 = -(8 * (.e7 * .e6/.e2)))
    c(v1 = c(v1 = c(v1 = 0, v2 = 8 * (v2/.e2)), v2 = .e10), v2 = c(v1 = .e10, 
        v2 = c(v1 = 6 * (.e6/.e3), v2 = -(.e8 + 4 * (.e5 * (16 * 
            .e7 + 16 * (4 - .e4)) * .e6^2/.e1^4)))))
}
############################################################
#' The first derivative of the density
#' @returns Vector
#' @inheritParams manf
norm_f1fa=function(x,v1,v2){
	vf=Vectorize(norm_fd,"x")
	f1=vf(x,v1,v2)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @returns Matrix
#' @inheritParams manf
norm_f2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(norm_fdd)
	temp1=vf(x,v1,v2)
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}
############################################################
#' The first derivative of the cdf
#' @returns Vector
#' @inheritParams manf
norm_p1fa=function(x,v1,v2){
	vf=Vectorize(norm_pd)
	p1=vf(x,v1,v2)
	return(p1)
}
############################################################
#' The second derivative of the cdf
#' @returns Matrix
#' @inheritParams manf
norm_p2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(norm_pdd)
	temp1=vf(x,v1,v2)
	p2=deriv_copyfdd(temp1,nx,dim=2)
	return(p2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @returns Vector
#' @inheritParams manf
norm_mu1fa=function(alpha,v1,v2){
	x=qnorm((1-alpha),mean=v1,sd=v2)
	vf=Vectorize(norm_pd)
	mu1=-vf(x,v1,v2)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @returns Matrix
#' @inheritParams manf
norm_mu2fa=function(alpha,v1,v2){
	x=qnorm((1-alpha),mean=v1,sd=v2)
	nx=length(x)
	vf=Vectorize(norm_pdd)
	temp1=vf(x,v1,v2)
	mu2=-deriv_copyfdd(temp1,nx,dim=2)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
norm_ldda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(norm_logfdd)
	temp1=vf(x,v1,v2)
#	temp1=norm_logfdd(x,v1,v2)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
norm_lddda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(norm_logfddd)
	temp1=vf(x,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}
