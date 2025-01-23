######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
invgamma_fd=function (x, v1, v2) 
{
    .e1 <- v2/x
    .e2 <- .e1^v1
    .e4 <- exp(-.e1)
    .e5 <- gamma(v1)
    c(v1 = ((log(v2) - log(x)) * .e2 - digamma(v1) * .e2) * .e4/(x * 
        .e5), v2 = .e4 * (v1 * .e1^(v1 - 1) - .e2)/(x^2 * .e5))
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
invgamma_fdd=function (x, v1, v2) 
{
    .e1 <- v2/x
    .e2 <- .e1^v1
    .e3 <- gamma(v1)
    .e4 <- v1 - 1
    .e5 <- .e1^.e4
    .e6 <- digamma(v1)
    .e8 <- log(v2) - log(x)
    .e10 <- exp(-.e1)
    .e11 <- .e8 * .e2
    .e12 <- .e6 * .e2
    .e13 <- x * .e3
    .e14 <- x^2
    .e15 <- .e11 - .e12
    .e17 <- v1 * .e5 - .e2
    .e18 <- .e14 * .e3
    c(v1 = c(v1 = ((.e15 * .e8 - trigamma(v1) * .e2)/.e13 - x * 
        .e15 * .e6 * .e3/.e13^2) * .e10, v2 = ((.e8 * .e17 + 
        .e5)/.e18 - .e14 * .e6 * .e3 * .e17/.e18^2) * .e10), 
        v2 = c(v1 = ((.e12 + v1 * (.e8 * .e5 - .e6 * .e5) - .e11)/x + 
            .e2/v2) * .e10/.e13, v2 = (.e2 + v1 * (.e4 * .e1^(v1 - 
            2) - 2 * .e5)) * .e10/(x^3 * .e3)))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
invgamma_logfdd=function (x, v1, v2) 
{
    .e1 <- 1/v2
    c(v1 = c(v1 = -trigamma(v1), v2 = .e1), v2 = c(v1 = .e1, 
        v2 = -(v1/v2^2)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
invgamma_logfddd=function (x, v1, v2) 
{
    .e1 <- -(1/v2^2)
    .e2 <- c(v1 = 0, v2 = .e1)
    c(v1 = c(v1 = c(v1 = -psigamma(v1, 2L), v2 = 0), v2 = .e2), 
        v2 = c(v1 = .e2, v2 = c(v1 = .e1, v2 = 2 * (v1/v2^3))))
}
############################################################
#' The first derivative of the density
#' @inheritParams manf
invgamma_f1fa=function(x,v1,v2){
	vf=Vectorize(invgamma_fd)
	f1=vf(x,v1,v2)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @inheritParams manf
invgamma_f2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(invgamma_fdd)
	temp1=vf(x,v1,v2)
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}
############################################################
############################################################
############################################################
#' The second derivative of the normalized log-likelihood
#' @inheritParams manf
invgamma_ldda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(invgamma_logfdd)
	temp1=vf(x,v1,v2)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @inheritParams manf
invgamma_lddda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(invgamma_logfddd)
	temp1=vf(x,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}
