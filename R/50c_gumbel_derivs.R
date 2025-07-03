######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gumbel_fd=function (x, v1, v2) 
{
    .e1 <- x - v1
    .e2 <- .e1/v2
    .e4 <- exp(-.e2)
    .e6 <- exp(-(.e2 + .e4))
    .e7 <- .e4 - 1
    .e8 <- v2^2
    c(v1 = -(.e6 * .e7/.e8), v2 = -((.e7 * .e1/v2 + 1) * .e6/.e8))
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gumbel_fdd=function (x, v1, v2) 
{
    .e1 <- x - v1
    .e2 <- .e1/v2
    .e4 <- exp(-.e2)
    .e5 <- .e4 - 1
    .e7 <- exp(-(.e2 + .e4))
    .e8 <- v2^3
    .e11 <- .e5 * .e1/v2 + 1
    .e15 <- (.e2 - 1) * .e4 + 1 - .e11 * .e5
    .e17 <- .e4 - .e5^2
    c(v1 = c(v1 = -(.e7 * .e17/.e8), v2 = -(.e15 * .e7/.e8)), 
        v2 = c(v1 = -((.e17 * .e1/v2 - 2 * .e5) * .e7/.e8), v2 = -((.e15 * 
            .e1/v2 - 2 * .e11) * .e7/.e8)))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gumbel_pd=function (x, v1, v2) 
{
    .e1 <- x - v1
    .e3 <- exp(-(.e1/v2))
    .e5 <- .e3 * exp(-.e3)
    c(v1 = -(.e5/v2), v2 = -(.e5 * .e1/v2^2))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gumbel_pdd=function (x, v1, v2) 
{
    .e1 <- x - v1
    .e3 <- exp(-(.e1/v2))
    .e5 <- 1 - .e3
    .e6 <- exp(-.e3)
    .e8 <- .e5 * .e1/v2
    .e9 <- v2^2
    .e10 <- -((.e8 - 1) * .e3 * .e6/.e9)
    c(v1 = c(v1 = -(.e5 * .e3 * .e6/.e9), v2 = .e10), v2 = c(v1 = .e10, 
        v2 = -((.e8 - 2) * .e3 * .e6 * .e1/v2^3)))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gumbel_logfdd=function (x, v1, v2) 
{
    .e1 <- x - v1
    .e2 <- .e1/v2
    .e4 <- exp(-.e2)
    .e5 <- v2^2
    .e6 <- -(((.e2 - 1) * .e4 + 1)/.e5)
    c(v1 = c(v1 = -(.e4/.e5), v2 = .e6), v2 = c(v1 = .e6, v2 = -((((.e2 - 
        2) * .e4 + 2) * .e1/v2 - 1)/.e5)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
gumbel_logfddd=function (x, v1, v2) 
{
    .e1 <- x - v1
    .e2 <- .e1/v2
    .e4 <- exp(-.e2)
    .e5 <- v2^3
    .e6 <- (.e2 - 2) * .e4
    .e7 <- -(.e6/.e5)
    .e8 <- -((.e6 * .e1/v2 - 2 * ((.e2 - 1) * .e4 + 1))/.e5)
    .e10 <- ((.e2 - 4) * .e1/v2 + 2) * .e4 - 2
    c(v1 = c(v1 = c(v1 = -(.e4/.e5), v2 = .e7), v2 = c(v1 = .e7, 
        v2 = -(.e10/.e5))), v2 = c(v1 = c(v1 = .e7, v2 = .e8), 
        v2 = c(v1 = .e8, v2 = -((.e10 * .e1/v2 - 2 * ((.e6 + 
            2) * .e1/v2 - 1))/.e5))))
}
############################################################
#' The first derivative of the density
#' @returns Vector
#' @inheritParams manf
gumbel_f1fa=function(x,v1,v2){
	vf=Vectorize(gumbel_fd)
	f1=vf(x,v1,v2)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @returns Matrix
#' @inheritParams manf
gumbel_f2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(gumbel_fdd)
	temp1=vf(x,v1,v2)
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}
############################################################
#' The first derivative of the cdf
#' @returns Vector
#' @inheritParams manf
gumbel_p1fa=function(x,v1,v2){
	vf=Vectorize(gumbel_pd)
	p1=vf(x,v1,v2)
	return(p1)
}
############################################################
#' The second derivative of the cdf
#' @returns Matrix
#' @inheritParams manf
gumbel_p2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(gumbel_pdd)
	temp1=vf(x,v1,v2)
	p2=deriv_copyfdd(temp1,nx,dim=2)
	return(p2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @returns Vector
#' @inheritParams manf
gumbel_mu1fa=function(alpha,v1,v2){
	x=qgumbel((1-alpha),mu=v1,sigma=v2)
	vf=Vectorize(gumbel_pd)
	mu1=-vf(x,v1,v2)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @returns Matrix
#' @inheritParams manf
gumbel_mu2fa=function(alpha,v1,v2){
	x=qgumbel((1-alpha),mu=v1,sigma=v2)
	nx=length(x)
	vf=Vectorize(gumbel_pdd)
	temp1=vf(x,v1,v2)
	mu2=-deriv_copyfdd(temp1,nx,dim=2)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
gumbel_ldda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(gumbel_logfdd)
	temp1=vf(x,v1,v2)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
gumbel_lddda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(gumbel_logfddd)
	temp1=vf(x,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}
