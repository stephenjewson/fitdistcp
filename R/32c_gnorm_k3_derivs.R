######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gnorm_k3_fd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- abs(.e1)
    .e3 <- .e2/v2
    .e5 <- .e3^(v3 - 1)
    .e7 <- exp(-.e3^v3)
    .e8 <- gamma(1/v3)
    c(v1 = v3^2 * .e5 * .e7 * sign(.e1)/(2 * (v2^2 * .e8)), v2 = v3 * 
        .e7 * (v3 * .e2 * .e5/(2 * v2^3) - 2/(2 * v2)^2)/.e8)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gnorm_k3_fdd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- abs(.e1)
    .e3 <- .e2/v2
    .e4 <- v3 - 1
    .e5 <- .e3^.e4
    .e7 <- gamma(1/v3)
    .e9 <- .e3^(v3 - 2)
    .e10 <- exp(-.e3^v3)
    .e11 <- v2^2
    .e12 <- v2^3
    .e13 <- .e9 * .e4
    .e14 <- 2 * v2
    .e15 <- 2 * .e12
    .e16 <- sign(.e1)
    .e17 <- v3^2
    .e20 <- .e5 * (v3 * .e2 * .e5/.e15 - 2/.e14^2)
    .e22 <- v3 * .e3^(2 * .e4) - .e13
    c(v1 = c(v1 = .e17 * .e10 * .e16^2 * .e22/(2 * (.e12 * .e7)), 
        v2 = .e17 * (.e20 - (.e5 + .e2 * .e9 * .e4/v2)/(2 * .e11)) * 
            .e10 * .e16/(v2 * .e7)), v2 = c(v1 = .e17 * (.e2 * 
        .e22/(2 * (v2^4 * .e7)) - 4 * (v2 * .e5 * .e7/(2 * (.e11 * 
        .e7))^2)) * .e10 * .e16, v2 = v3 * (16 * (v2/.e14^4) + 
        v3 * ((.e20/.e11 - 6 * (.e11 * .e5/.e15^2)) * .e2 - .e13 * 
            .e1^2/(2 * v2^5))) * .e10/.e7))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gnorm_k3_logfdd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- abs(.e1)
    .e3 <- .e2/v2
    .e4 <- v3 - 1
    .e5 <- .e3^(v3 - 2)
    .e6 <- v2^2
    .e7 <- .e3^.e4
    .e8 <- sign(.e1)
    .e9 <- -(v3 * (.e7 + .e2 * .e5 * .e4/v2) * .e8/.e6)
    c(v1 = c(v1 = -(v3 * .e5 * .e8^2 * .e4/.e6), v2 = .e9), v2 = c(v1 = .e9, 
        v2 = -((v3 * (.e5 * .e4 * .e1^2/v2 + 2 * (.e2 * .e7))/v2 - 
            1)/.e6)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
gnorm_k3_logfddd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- abs(.e1)
    .e3 <- .e2/v2
    .e4 <- v3 - 2
    .e5 <- v3 - 1
    .e6 <- .e3^.e4
    .e7 <- .e3^(v3 - 3)
    .e8 <- sign(.e1)
    .e9 <- v2^3
    .e10 <- .e3^.e5
    .e11 <- .e2 * .e6
    .e12 <- .e1^2
    .e13 <- 2 * .e6
    .e16 <- .e2 * .e7 * .e4/v2
    .e17 <- .e10 + .e11 * .e5/v2
    .e23 <- v3 * (.e13 + .e16) * .e8^2 * .e5/.e9
    .e27 <- 2 * (.e2 * .e10)
    .e30 <- v3 * ((.e7 * .e4 * .e12/v2 + 2 * .e11) * .e5/v2 + 
        2 * .e17) * .e8/.e9
    c(v1 = c(v1 = c(v1 = v3 * .e7 * .e8^3 * .e5 * .e4/.e9, v2 = .e23), 
        v2 = c(v1 = .e23, v2 = v3 * ((.e7 * .e8 * .e4 * .e1/v2 + 
            .e13) * .e5 * .e1/v2 + 2 * (.e17 * .e8))/.e9)), v2 = c(v1 = c(v1 = .e23, 
        v2 = .e30), v2 = c(v1 = .e30, v2 = (2 * (v3 * (.e6 * 
        .e5 * .e12/v2 + .e27)/v2 - 1) + v3 * ((4 * .e6 + .e16) * 
        .e5 * .e12/v2 + .e27)/v2)/.e9)))
}
############################################################
#' The first derivative of the density
#' @returns Vector
#' @inheritParams manf
gnorm_k3_f1fa=function(x,v1,v2,kbeta){
	vf=Vectorize(gnorm_k3_fd,"x")
	f1=vf(x,v1,v2,kbeta)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @returns Matrix
#' @returns Matrix
#' @inheritParams manf
gnorm_k3_f2fa=function(x,v1,v2,kbeta){
	nx=length(x)
	vf=Vectorize(gnorm_k3_fdd,"x")
	temp1=vf(x,v1,v2,kbeta)
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}
############################################################
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
gnorm_k3_ldda=function(x,v1,v2,kbeta){
	nx=length(x)
	vf=Vectorize(gnorm_k3_logfdd,"x")
	temp1=vf(x,v1,v2,kbeta)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
gnorm_k3_lddda=function(x,v1,v2,kbeta){
	nx=length(x)
	vf=Vectorize(gnorm_k3_logfddd,"x")
	temp1=vf(x,v1,v2,kbeta)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}
