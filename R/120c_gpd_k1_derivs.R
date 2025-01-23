######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
gpd_k1_fd=function (x, v1, v2, v3) 
{
    .e1 <- 1 + v3
    .e2 <- x - v1
    .e3 <- .e1/v3
    .e5 <- v3 * .e2/v2
    .e6 <- 1 + .e5
    .e8 <- .e6^.e3
    .e10 <- .e1 * .e2/(v2 * .e6^(.e3 + 1))
    c(v2 = (.e10 - 1/.e8)/v2^2, v3 = -(((1 - .e3) * log1p(.e5)/.e8 + 
        .e10)/(v2 * v3)))
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
gpd_k1_fdd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- 1 + v3
    .e3 <- .e2/v3
    .e5 <- v3 * .e1/v2
    .e6 <- 1 + .e5
    .e7 <- .e3 + 1
    .e8 <- .e6^.e7
    .e9 <- .e6^.e3
    .e10 <- log1p(.e5)
    .e11 <- v2 * .e8
    .e12 <- 1 - .e3
    .e13 <- .e11^2
    .e14 <- v2 * v3
    .e15 <- v2^2
    .e17 <- .e2 * .e1/.e11
    .e24 <- .e12 * .e9 * .e10 + .e2 * .e6^(.e3 - 1) * .e1/v2
    .e27 <- .e12 * .e10/.e9 + .e17
    .e28 <- .e8 - v3 * .e7 * .e9 * .e1/v2
    .e30 <- (1/.e11 - v2 * (.e7 * .e9 * .e1/v2 + .e12 * .e8 * 
        .e10/v3) * .e2/.e13) * .e1
    .e31 <- .e14^2
    .e32 <- v3 * .e6^(2 * .e3)
    c(v2 = c(v2 = -(((.e28/.e13 + 1/(.e15 * .e8)) * .e2 * .e1 + 
        2 * ((.e17 - 1/.e9)/v2))/.e15), v3 = -(((.e2 * .e10/.e8 - 
        v3/.e8) * .e12/.e15 - .e28 * .e2/.e13) * .e1/.e14 - v3 * 
        .e27/.e31)), v3 = c(v2 = (.e24/.e32 + .e30)/.e15, v3 = -((((.e1/(v2 * 
        .e6) - .e10/v3)/.e9 - .e24 * .e10/.e32) * .e12 + .e30)/.e14 - 
        v2 * .e27/.e31)))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
gpd_k1_pd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e3 <- v3 * .e1/v2
    .e4 <- 1 + .e3
    .e5 <- 1/v3
    .e6 <- .e4^(1 + .e5)
    c(v2 = -(.e1/(v2^2 * .e6)), v3 = -((log1p(.e3)/(v3 * .e4^.e5) - 
        .e1/(v2 * .e6))/v3))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
gpd_k1_pdd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e3 <- v3 * .e1/v2
    .e4 <- 1/v3
    .e5 <- 1 + .e3
    .e6 <- 1 + .e4
    .e7 <- .e5^.e6
    .e8 <- .e5^.e4
    .e9 <- log1p(.e3)
    .e10 <- v2 * .e7
    .e11 <- v2^2
    .e12 <- v3 * .e8
    .e16 <- .e6 * .e8 * .e1/v2 - .e7 * .e9/v3^2
    .e17 <- .e5^(.e4 - 1)
    .e18 <- .e10^2
    .e19 <- (.e11 * .e7)^2
    .e20 <- .e12^2
    .e23 <- v3 * .e6 * .e8 * .e1
    c(v2 = c(v2 = (2 * .e10 - .e23) * .e1/.e19, v3 = -(((.e7 - 
        .e23/v2)/.e18 + (v3 * .e17 * .e9/.e20 - 1/.e7)/.e11) * 
        .e1/v3)), v3 = c(v2 = .e11 * .e16 * .e1/.e19, v3 = -(((1/(v2 * 
        v3 * .e7) + v2 * .e16/.e18) * .e1 - ((.e17 * .e1/v2 + 
        .e8 - .e8 * .e9/v3) * .e9/.e20 + (.e9/.e12 - .e1/.e10)/v3))/v3)))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
gpd_k1_logfdd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e3 <- v3 * .e1/v2
    .e4 <- 1 + .e3
    .e5 <- 1 + v3
    .e6 <- v2 * .e4
    .e7 <- .e5 * .e1
    .e8 <- .e6^2
    .e9 <- .e7/.e8
    .e11 <- 1 - .e5/v3
    .e12 <- .e7/.e6
    .e13 <- (1/.e6 - .e9) * .e1
    .e14 <- log1p(.e3)
    c(v2 = c(v2 = -(((.e12 - 1)/v2 + .e9)/v2), v3 = (.e5/.e8 + 
        v3 * .e11/(v2^2 * .e4)) * .e1/v3), v3 = c(v2 = .e13/v2, 
        v3 = -(((.e1/.e6 - .e14/v3) * .e11 + .e13 - (.e11 * .e14 + 
            .e12)/v3)/v3)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
gpd_k1_logfddd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- v3 * .e1
    .e3 <- .e2/v2
    .e4 <- 1 + .e3
    .e5 <- v2 * .e4
    .e6 <- 1 + v3
    .e7 <- .e5^2
    .e9 <- 1 - .e6/v3
    .e11 <- v2 * .e6 * .e4
    .e12 <- .e11 * .e1
    .e13 <- .e6 * .e1
    .e14 <- 2 * (.e12/.e7)
    .e16 <- v2^2 * .e4
    .e17 <- .e13/.e7
    .e18 <- .e1^2
    .e20 <- 1/.e5 - .e17
    .e21 <- log1p(.e3)
    .e22 <- v3 * .e9
    .e23 <- (.e6/.e7 + .e22/.e16)/v3
    .e24 <- (1 - .e14)/.e7
    .e25 <- .e13/.e5
    .e26 <- .e20/v2
    .e27 <- (2 - .e14) * .e18
    .e28 <- (.e14 - 1)/.e7
    .e29 <- .e5^4
    .e30 <- .e16^2
    .e32 <- .e1/.e5 - .e21/v3
    c(v2 = c(v2 = c(v2 = (2 * (((.e25 - 1)/v2 + .e17)/v2) + 2 * 
        (.e12/.e29))/v2, v3 = -((2 * (.e11/.e29) + .e22 * (2 * 
        .e5 - .e2)/.e30) * .e1/v3)), v3 = c(v2 = (.e28 - .e26) * 
        .e1/v2, v3 = -((.e23 + .e9 * (1/.e16 - 1/.e7) + .e28) * 
        .e1/v3))), v3 = c(v2 = c(v2 = -((.e24 + .e26) * .e1/v2), 
        v3 = (.e24 - (.e23 + v2 * v3 * .e9 * .e1/.e30)) * .e1/v3), 
        v3 = c(v2 = -(.e27/(v2 * .e7)), v3 = ((.e18/.e7 + 2 * 
            (.e32/v3)) * .e9 + .e27/.e7 + 2 * ((.e32 * .e9 + 
            .e20 * .e1 - (.e9 * .e21 + .e25)/v3)/v3))/v3)))
}
############################################################
#' The first derivative of the density
#' @inheritParams manf
gpd_k1_f1fa=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my pu code
	vf=Vectorize(gpd_k1_fd)

	v2=movexiawayfromzero(v2)

	f1=vf(x,kloc,v1,v2) #these are in mu, sigma, xi order
	return(f1)
}
############################################################
#' The second derivative of the density
#' @inheritParams manf
gpd_k1_f2fa=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my pu code
	nx=length(x)
	vf=Vectorize(gpd_k1_fdd)

	v2=movexiawayfromzero(v2)

	temp1=vf(x,kloc,v1,v2) #these are in mu, sigma, xi order
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @inheritParams manf
gpd_k1_mu1fa=function(alpha,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my pu code
	x=extraDistr::qgpd((1-alpha),mu=kloc,sigma=v1,xi=v2)
	vf=Vectorize(gpd_k1_pd)

	v2=movexiawayfromzero(v2)

	mu1=-vf(x,kloc,v1,v2) #these are in mu, sigma, xi order
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @inheritParams manf
gpd_k1_mu2fa=function(alpha,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my pu code
	x=extraDistr::qgpd((1-alpha),mu=kloc,sigma=v1,xi=v2)
	nx=length(x)
	vf=Vectorize(gpd_k1_pdd)

	v2=movexiawayfromzero(v2)

	temp1=vf(x,kloc,v1,v2) #these are in mu, sigma, xi order
	mu2=-deriv_copyfdd(temp1,nx,dim=2)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @inheritParams manf
gpd_k1_ldda=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my pu code
	nx=length(x)
	vf=Vectorize(gpd_k1_logfdd)

	v2=movexiawayfromzero(v2)

	temp1=vf(x,kloc,v1,v2) #these are in mu, sigma, xi order
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @inheritParams manf
gpd_k1_lddda=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my pu code
# I have to switch because my pu code orders sigma and lambda differently
	nx=length(x)
	vf=Vectorize(gpd_k1_logfddd)

	v2=movexiawayfromzero(v2)

	temp1=vf(x,kloc,v1,v2) #these are in mu, sigma, xi order
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}
