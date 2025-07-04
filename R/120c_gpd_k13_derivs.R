######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gpd_k13_fd=function (x, v1, v2, v3) 
{
    .e1 <- 1 + v3
    .e2 <- x - v1
    .e3 <- .e1/v3
    .e4 <- 1 + v3 * .e2/v2
    (.e1 * .e2/(v2 * .e4^(.e3 + 1)) - 1/.e4^.e3)/v2^2
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gpd_k13_fdd=function (x, v1, v2, v3) 
{
    .e1 <- 1 + v3
    .e2 <- x - v1
    .e3 <- .e1/v3
    .e4 <- 1 + v3 * .e2/v2
    .e5 <- .e3 + 1
    .e6 <- .e4^.e5
    .e7 <- .e4^.e3
    .e8 <- v2 * .e6
    .e9 <- v2^2
    -((((.e6 - v3 * .e5 * .e7 * .e2/v2)/.e8^2 + 1/(.e9 * .e6)) * 
        .e1 * .e2 + 2 * ((.e1 * .e2/.e8 - 1/.e7)/v2))/.e9)
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gpd_k13_pd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    -(.e1/(v2^2 * (1 + v3 * .e1/v2)^(1 + 1/v3)))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gpd_k13_pdd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- 1/v3
    .e3 <- 1 + .e2
    .e4 <- 1 + v3 * .e1/v2
    .e5 <- .e4^.e3
    (2 * (v2 * .e5) - v3 * .e3 * .e4^.e2 * .e1) * .e1/(v2^2 * 
        .e5)^2
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gpd_k13_logfdd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- (1 + v3) * .e1
    .e4 <- v2 * (1 + v3 * .e1/v2)
    -(((.e2/.e4 - 1)/v2 + .e2/.e4^2)/v2)
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
gpd_k13_logfddd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- 1 + v3 * .e1/v2
    .e3 <- 1 + v3
    .e4 <- v2 * .e2
    .e5 <- .e3 * .e1
    (2 * (((.e5/.e4 - 1)/v2 + .e5/.e4^2)/v2) + 2 * (v2 * .e3 * 
        .e2 * .e1/.e4^4))/v2
}
############################################################
#' The first derivative of the density
#' @returns Vector
#' @inheritParams manf
gpd_k13_f1fa=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code

	nx=length(x)
	f1=matrix(0,1,nx)
	vf=Vectorize(gpd_k13_fd,"x")
	v2=movexiawayfromzero(v2)
	f1[1,]=vf(x,kloc,v1,v2)

	return(f1)
}
############################################################
#' The second derivative of the density
#' @returns Matrix
#' @inheritParams manf
gpd_k13_f2fa=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code

	nx=length(x)
	f2=array(0,c(1,1,nx))
	vf=Vectorize(gpd_k13_fdd,"x")
	v2=movexiawayfromzero(v2)
	f2[1,1,]=vf(x,kloc,v1,v2)

	return(f2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @returns Vector
#' @inheritParams manf
gpd_k13_mu1fa=function(alpha,v1,v2,kloc){
	x=extraDistr::qgpd((1-alpha),mu=kloc,sigma=v1,xi=v2)
	nx=length(x)
	mu1=array(0,c(1,nx))
	vf=Vectorize(gpd_k13_pd,"x")
	mu1[1,]=-vf(x,v1,v2,kloc)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @returns Matrix
#' @inheritParams manf
gpd_k13_mu2fa=function(alpha,v1,v2,kloc){
	x=extraDistr::qgpd((1-alpha),mu=kloc,sigma=v1,xi=v2)
	nx=length(x)
	mu2=array(0,c(1,1,nx))
	vf=Vectorize(gpd_k13_pdd,"x")
	mu2[1,1,]=-vf(x,v1,v2,kloc)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
gpd_k13_ldda=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code

	nx=length(x)
	ldd=matrix(0,1,1)
	vf=Vectorize(gpd_k13_logfdd,"x")
	v2=movexiawayfromzero(v2)
	ldd[1,1]=sum(vf(x,kloc,v1,v2))/nx

	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
gpd_k13_lddda=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch because my cp code orders sigma and lambda differently

	nx=length(x)
	lddd=array(0,c(1,1,1))
	vf=Vectorize(gpd_k13_logfddd,"x")
	v2=movexiawayfromzero(v2)
	lddd[1,1,1]=sum(vf(x,kloc,v1,v2))/nx

	return(lddd)
}
