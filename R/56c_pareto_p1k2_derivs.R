######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
pareto_p1k2_fd=function (x, t, v1, v2, v3) 
{
    .e2 <- exp(-(t * v2 + v1))
    .e3 <- 1 + .e2
    .e4 <- v3^.e2
    .e11 <- .e4 * x^(.e3 - 2 * .e3) * .e2 * log(x) - (.e4 + .e4 * 
        .e2 * log(v3))/x^.e3
    c(v1 = .e2 * .e11, v2 = t * .e2 * .e11)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
pareto_p1k2_fdd=function (x, t, v1, v2, v3) 
{
    .e2 <- exp(-(t * v2 + v1))
    .e3 <- 1 + .e2
    .e4 <- v3^.e2
    .e5 <- log(v3)
    .e7 <- x^(.e3 - 2 * .e3)
    .e9 <- .e4 * .e2 * .e5
    .e12 <- .e4 + .e9
    .e13 <- x^.e3
    .e15 <- ((2 * .e4 + .e9) * .e5/.e13 + (.e2 * (.e4 * .e7 * 
        log(x) - .e4 * .e7 * .e5) - (2 * (.e4 * .e7) + .e7 * 
        .e12)) * log(x)) * .e2 + .e12/.e13
    .e17 <- t * .e15 * .e2
    c(v1 = c(v1 = .e15 * .e2, v2 = .e17), v2 = c(v1 = .e17, v2 = t^2 * 
        .e15 * .e2))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
pareto_p1k2_pd=function (x, t, v1, v2, v3) 
{
    .e2 <- exp(-(t * v2 + v1))
    .e3 <- (v3/x)^.e2
    .e5 <- log(v3) - log(x)
    c(v1 = .e2 * .e5 * .e3, v2 = t * .e2 * .e5 * .e3)
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
pareto_p1k2_pdd=function (x, t, v1, v2, v3) 
{
    .e2 <- exp(-(t * v2 + v1))
    .e5 <- log(v3) - log(x)
    .e6 <- (v3/x)^.e2 + .e2 * .e5 * (v3/x)^.e2
    .e7 <- -(t * .e6 * .e2 * .e5)
    c(v1 = c(v1 = -(.e6 * .e2 * .e5), v2 = .e7), v2 = c(v1 = .e7, 
        v2 = -(t^2 * .e6 * .e2 * .e5)))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
pareto_p1k2_logfdd=function (x, t, v1, v2, v3) 
{
    .e2 <- exp(-(t * v2 + v1))
    .e5 <- log(x) - log(v3)
    .e6 <- -(t * .e2 * .e5)
    c(v1 = c(v1 = -(.e2 * .e5), v2 = .e6), v2 = c(v1 = .e6, v2 = -(t^2 * 
        .e2 * .e5)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
pareto_p1k2_logfddd=function (x, t, v1, v2, v3) 
{
    .e2 <- exp(-(t * v2 + v1))
    .e5 <- log(x) - log(v3)
    .e7 <- t * .e2 * .e5
    .e10 <- t^2 * .e2 * .e5
    .e11 <- c(v1 = .e7, v2 = .e10)
    c(v1 = c(v1 = c(v1 = .e2 * .e5, v2 = .e7), v2 = .e11), v2 = c(v1 = .e11, 
        v2 = c(v1 = .e10, v2 = t^3 * .e2 * .e5)))
}
############################################################
#' The first derivative of the density
#' @inheritParams manf
pareto_p1k2_f1fa=function(x,t,v1,v2,kscale){
	vf=Vectorize(pareto_p1k2_fd)
	f1=vf(x,t,v1,v2,kscale)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @inheritParams manf
pareto_p1k2_f2fa=function(x,t,v1,v2,kscale){
	nx=length(x)
	vf=Vectorize(pareto_p1k2_fdd)
	temp1=vf(x,t,v1,v2,kscale)
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}
############################################################
#' The first derivative of the cdf
#' @inheritParams manf
pareto_p1k2_p1fa=function(x,t,v1,v2,kscale){
	vf=Vectorize(pareto_p1k2_pd)
	p1=vf(x,t,v1,v2,kscale)
	return(p1)
}
############################################################
#' The second derivative of the cdf
#' @inheritParams manf
pareto_p1k2_p2fa=function(x,t,v1,v2,kscale){
	nx=length(x)
	vf=Vectorize(pareto_p1k2_pdd)
	temp1=vf(x,t,v1,v2,kscale)
	p2=deriv_copyfdd(temp1,nx,dim=2)
	return(p2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @inheritParams manf
pareto_p1k2_mu1fa=function(alpha,t,v1,v2,kscale){
	x=extraDistr::qpareto((1-alpha),a=exp(-v1-v2*t),b=kscale)
	vf=Vectorize(pareto_p1k2_pd)
	mu1=-vf(x,t,v1,v2,kscale)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @inheritParams manf
pareto_p1k2_mu2fa=function(alpha,t,v1,v2,kscale){
	x=extraDistr::qpareto((1-alpha),a=exp(-v1-v2*t),b=kscale)
	nalpha=length(alpha)
	vf=Vectorize(pareto_p1k2_pdd)
	temp1=vf(x,t,v1,v2,kscale)
	mu2=-deriv_copyfdd(temp1,nalpha,dim=2)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @inheritParams manf
pareto_p1k2_ldda=function(x,t,v1,v2,kscale){
	nx=length(x)
	vf=Vectorize(pareto_p1k2_logfdd)
	temp1=vf(x,t,v1,v2,kscale)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @inheritParams manf
pareto_p1k2_lddda=function(x,t,v1,v2,kscale){
	nx=length(x)
	vf=Vectorize(pareto_p1k2_logfddd)
	temp1=vf(x,t,v1,v2,kscale)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}
