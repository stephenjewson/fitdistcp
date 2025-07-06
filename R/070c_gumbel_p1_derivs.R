######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gumbel_p1_fd=function (x, t, v1, v2, v3) 
{
    .e3 <- x - (t * v2 + v1)
    .e4 <- .e3/v3
    .e6 <- exp(-.e4)
    .e8 <- exp(-(.e4 + .e6))
    .e9 <- .e6 - 1
    .e10 <- v3^2
    c(v1 = -(.e8 * .e9/.e10), v2 = -(t * .e8 * .e9/.e10), v3 = -((.e9 * 
        .e3/v3 + 1) * .e8/.e10))
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gumbel_p1_fdd=function (x, t, v1, v2, v3) 
{
    .e3 <- x - (t * v2 + v1)
    .e4 <- .e3/v3
    .e6 <- exp(-.e4)
    .e7 <- .e6 - 1
    .e9 <- exp(-(.e4 + .e6))
    .e10 <- v3^3
    .e12 <- .e6 - .e7^2
    .e15 <- .e7 * .e3/v3 + 1
    .e19 <- (.e4 - 1) * .e6 + 1 - .e15 * .e7
    .e20 <- -(t * .e9 * .e12/.e10)
    .e23 <- .e12 * .e3/v3 - 2 * .e7
    c(v1 = c(v1 = -(.e9 * .e12/.e10), v2 = .e20, v3 = -(.e19 * 
        .e9/.e10)), v2 = c(v1 = .e20, v2 = -(t^2 * .e9 * .e12/.e10), 
        v3 = -(t * .e19 * .e9/.e10)), v3 = c(v1 = -(.e23 * .e9/.e10), 
        v2 = -(t * .e23 * .e9/.e10), v3 = -((.e19 * .e3/v3 - 
            2 * .e15) * .e9/.e10)))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gumbel_p1_pd=function (x, t, v1, v2, v3) 
{
    .e3 <- x - (t * v2 + v1)
    .e5 <- exp(-(.e3/v3))
    .e7 <- exp(-.e5)
    .e8 <- .e5 * .e7
    c(v1 = -(.e8/v3), v2 = -(t * .e5 * .e7/v3), v3 = -(.e8 * 
        .e3/v3^2))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gumbel_p1_pdd=function (x, t, v1, v2, v3) 
{
    .e3 <- x - (t * v2 + v1)
    .e5 <- exp(-(.e3/v3))
    .e7 <- 1 - .e5
    .e8 <- exp(-.e5)
    .e9 <- v3^2
    .e11 <- .e7 * .e3/v3
    .e12 <- .e11 - 1
    .e13 <- -(.e12 * .e5 * .e8/.e9)
    .e14 <- -(t * .e12 * .e5 * .e8/.e9)
    .e15 <- -(t * .e7 * .e5 * .e8/.e9)
    c(v1 = c(v1 = -(.e7 * .e5 * .e8/.e9), v2 = .e15, v3 = .e13), 
        v2 = c(v1 = .e15, v2 = -(t^2 * .e7 * .e5 * .e8/.e9), 
            v3 = .e14), v3 = c(v1 = .e13, v2 = .e14, v3 = -((.e11 - 
            2) * .e5 * .e8 * .e3/v3^3)))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gumbel_p1_logfdd=function (x, t, v1, v2, v3) 
{
    .e3 <- x - (t * v2 + v1)
    .e4 <- .e3/v3
    .e6 <- exp(-.e4)
    .e7 <- v3^2
    .e9 <- (.e4 - 1) * .e6 + 1
    .e10 <- -(.e9/.e7)
    .e11 <- -(t * .e9/.e7)
    .e12 <- -(t * .e6/.e7)
    c(v1 = c(v1 = -(.e6/.e7), v2 = .e12, v3 = .e10), v2 = c(v1 = .e12, 
        v2 = -(t^2 * .e6/.e7), v3 = .e11), v3 = c(v1 = .e10, 
        v2 = .e11, v3 = -((((.e4 - 2) * .e6 + 2) * .e3/v3 - 1)/.e7)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
gumbel_p1_logfddd=function (x, t, v1, v2, v3) 
{
    .e3 <- x - (t * v2 + v1)
    .e4 <- .e3/v3
    .e6 <- exp(-.e4)
    .e7 <- v3^3
    .e8 <- .e4 - 2
    .e9 <- .e8 * .e6
    .e10 <- -(t * .e8 * .e6/.e7)
    .e11 <- t^2
    .e16 <- .e9 * .e3/v3 - 2 * ((.e4 - 1) * .e6 + 1)
    .e17 <- -(.e9/.e7)
    .e18 <- -(t * .e6/.e7)
    .e19 <- -(.e11 * .e8 * .e6/.e7)
    .e20 <- -(.e11 * .e6/.e7)
    .e22 <- ((.e4 - 4) * .e3/v3 + 2) * .e6 - 2
    .e23 <- -(.e16/.e7)
    .e24 <- -(t * .e16/.e7)
    .e25 <- c(v1 = .e18, v2 = .e20, v3 = .e10)
    c(v1 = c(v1 = c(v1 = -(.e6/.e7), v2 = .e18, v3 = .e17), v2 = .e25, 
        v3 = c(v1 = .e17, v2 = .e10, v3 = -(.e22/.e7))), v2 = c(v1 = .e25, 
        v2 = c(v1 = .e20, v2 = -(t^3 * .e6/.e7), v3 = .e19), 
        v3 = c(v1 = .e10, v2 = .e19, v3 = -(t * .e22/.e7))), 
        v3 = c(v1 = c(v1 = .e17, v2 = .e10, v3 = .e23), v2 = c(v1 = .e10, 
            v2 = .e19, v3 = .e24), v3 = c(v1 = .e23, v2 = .e24, 
            v3 = -((.e22 * .e3/v3 - 2 * ((.e9 + 2) * .e3/v3 - 
                1))/.e7))))
}
############################################################
#' The first derivative of the density for DMGS
#' @returns Vector
#' @inheritParams manf
gumbel_p1_f1fa=function(x,t0,v1,v2,v3){
	vf=Vectorize(gumbel_p1_fd,"x")
	f1=vf(x,t0,v1,v2,v3)
	return(f1)
}
############################################################
#' The first derivative of the density for WAIC
#' @returns Vector
#' @inheritParams manf
gumbel_p1_f1fw=function(x,t,v1,v2,v3){
	vf=Vectorize(gumbel_p1_fd,c("x","t"))
	f1=vf(x,t,v1,v2,v3)
	return(f1)
}
############################################################
#' The second derivative of the density for DMGS
#' @returns Matrix
#' @inheritParams manf
gumbel_p1_f2fa=function(x,t0,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(gumbel_p1_fdd,"x")
	temp1=vf(x,t0,v1,v2,v3)
	f2=deriv_copyfdd(temp1,nx,dim=3)
	return(f2)
}
############################################################
#' The second derivative of the density for WAIC
#' @returns Matrix
#' @inheritParams manf
gumbel_p1_f2fw=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(gumbel_p1_fdd,c("x","t"))
	temp1=vf(x,t,v1,v2,v3)
	f2=deriv_copyfdd(temp1,nx,dim=3)
	return(f2)
}
############################################################
#' The first derivative of the cdf
#' @returns Vector
#' @inheritParams manf
gumbel_p1_p1fa=function(x,t0,v1,v2,v3){
	vf=Vectorize(gumbel_p1_pd,"x")
	p1=vf(x,t0,v1,v2,v3)
	return(p1)
}
############################################################
#' The second derivative of the cdf
#' @returns Matrix
#' @inheritParams manf
gumbel_p1_p2fa=function(x,t0,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(gumbel_p1_pdd,"x")
	temp1=vf(x,t0,v1,v2,v3)
	p2=deriv_copyfdd(temp1,nx,dim=3)
	return(p2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @returns Vector
#' @inheritParams manf
gumbel_p1_mu1fa=function(alpha,t0,v1,v2,v3){
	x=qgumbel((1-alpha),mu=v1+v2*t0,sigma=v3)
	vf=Vectorize(gumbel_p1_pd,"x")
	mu1=-vf(x,t0,v1,v2,v3)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @returns Matrix
#' @inheritParams manf
gumbel_p1_mu2fa=function(alpha,t0,v1,v2,v3){
	x=qgumbel((1-alpha),mu=v1+v2*t0,sigma=v3)
	nx=length(x)
	vf=Vectorize(gumbel_p1_pdd,"x")
	temp1=vf(x,t0,v1,v2,v3)
	mu2=-deriv_copyfdd(temp1,nx,dim=3)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
gumbel_p1_ldda=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(gumbel_p1_logfdd,c("x","t"))
	temp1=vf(x,t,v1,v2,v3)
	ldd=deriv_copyldd(temp1,nx,dim=3)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
gumbel_p1_lddda=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(gumbel_p1_logfddd,c("x","t"))
	temp1=vf(x,t,v1,v2,v3)
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}
