######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
gamma_fd=function (x, v1, v2) 
{
    .e2 <- v2^v1
    .e3 <- x^(v1 - 1)
    .e5 <- exp(-(x/v2))
    .e6 <- gamma(v1)
    c(v1 = .e5 * (.e3 * log(x)/.e2 - .e3 * (digamma(v1)/.e2 + 
        log(v2)/.e2))/.e6, v2 = .e5 * (x^v1/v2^(2 + v1) - v1 * 
        .e3/v2^(1 + v1))/.e6)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
gamma_fdd=function (x, v1, v2) 
{
    .e1 <- v1 - 1
    .e2 <- x^.e1
    .e3 <- v2^v1
    .e4 <- 1 + v1
    .e5 <- log(v2)
    .e6 <- 2 + v1
    .e7 <- digamma(v1)
    .e8 <- log(x)
    .e9 <- v2^.e4
    .e12 <- .e7/.e3 + .e5/.e3
    .e13 <- exp(-(x/v2))
    .e14 <- gamma(v1)
    .e15 <- x^v1
    .e17 <- v1 * .e2
    .e18 <- v2^.e6
    .e19 <- .e2 * .e12
    .e20 <- x^(v1 - 2)
    .e21 <- 2 * .e4
    .e22 <- 2 * .e6
    .e24 <- v2^2
    .e26 <- .e2 * .e8/.e3
    .e28 <- .e15/.e18 - .e17/.e9
    c(v1 = c(v1 = .e13 * (.e8 * (.e26 - (.e19 + .e2 * .e5/.e3)) - 
        (.e7 * (.e26 - .e19) + .e2 * (trigamma(v1)/.e3 - .e12 * 
            .e5)))/.e14, v2 = .e13 * (.e5 * (v1 * v2^(.e4 - .e21) * 
        .e2 - v2^(.e6 - .e22) * .e15) + .e15 * .e8/.e18 - ((.e17 * 
        .e8 + .e2)/.e9 + .e7 * .e28))/.e14), v2 = c(v1 = x * 
        ((.e8 * (.e2/.e3 - v1 * .e20/v2^.e1) - .e19)/.e24 - .e20 * 
            (1/.e9 - v1 * (.e7/.e9 + .e5/.e9))) * .e13/.e14, 
        v2 = x * (.e28/.e24 + v1 * v2^(v1 - .e21) * .e20 * .e4 - 
            v2^(.e4 - .e22) * .e2 * .e6) * .e13/.e14))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
gamma_logfdd=function (x, v1, v2) 
{
    .e1 <- -(1/v2)
    c(v1 = c(v1 = -trigamma(v1), v2 = .e1), v2 = c(v1 = .e1, 
        v2 = -((2 * (x/v2) - v1)/v2^2)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
gamma_logfddd=function (x, v1, v2) 
{
    .e1 <- 1/v2^2
    .e2 <- 2 * (x/v2)
    .e3 <- c(v1 = 0, v2 = .e1)
    c(v1 = c(v1 = c(v1 = -psigamma(v1, 2L), v2 = 0), v2 = .e3), 
        v2 = c(v1 = .e3, v2 = c(v1 = .e1, v2 = (2 * (.e2 - v1) + 
            .e2)/v2^3)))
}
############################################################
#' The first derivative of the density
#' @inheritParams manf
gamma_f1fa=function(x,v1,v2){
	vf=Vectorize(gamma_fd)
	f1=vf(x,v1,v2)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @inheritParams manf
gamma_f2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(gamma_fdd)
	temp1=vf(x,v1,v2)
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @inheritParams manf
gamma_ldda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(gamma_logfdd)
	temp1=vf(x,v1,v2)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @inheritParams manf
gamma_lddda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(gamma_logfddd)
	temp1=vf(x,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}
