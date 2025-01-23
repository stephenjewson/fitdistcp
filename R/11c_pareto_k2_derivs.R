######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
pareto_k2_fd=function (x, v1, v2) 
{
    .e1 <- 1 + v1
    .e2 <- v2^v1
    .e3 <- v1 * .e2
    (.e3 * log(v2) + .e2)/x^.e1 - .e3 * x^(.e1 - 2 * .e1) * log(x)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
pareto_k2_fdd=function (x, v1, v2) 
{
    .e1 <- v2^v1
    .e2 <- 1 + v1
    .e3 <- log(v2)
    .e4 <- v1 * .e1
    .e6 <- log(x)
    .e8 <- .e4 * .e3 + .e1
    .e9 <- x^(.e2 - 2 * .e2)
    .e3 * (.e8 + .e1)/x^.e2 - (2 * (.e9 * .e8) - .e4 * .e9 * 
        .e6) * .e6
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
pareto_k2_pd=function (x, v1, v2) 
-((log(v2) - log(x)) * (v2/x)^v1)
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
pareto_k2_pdd=function (x, v1, v2) 
-((log(v2) - log(x))^2 * (v2/x)^v1)
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
pareto_k2_logfdd=function (x, v1, v2) 
-(1/v1^2)
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
pareto_k2_logfddd=function (x, v1, v2) 
2/v1^3
############################################################
#' The first derivative of the density
#' @inheritParams manf
pareto_k2_f1fa=function(x,v1,kscale){
	nx=length(x)
	f1=matrix(0,1,nx)
	vf=Vectorize(pareto_k2_fd)
	f1[1,]=vf(x,v1,v2=kscale)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @inheritParams manf
pareto_k2_f2fa=function(x,v1,kscale){
	nx=length(x)
	f2=array(0,c(1,1,nx))
	vf=Vectorize(pareto_k2_fdd)
	f2[1,1,]=vf(x,v1,v2=kscale)
	return(f2)
}
############################################################
#' The first derivative of the cdf
#' @inheritParams manf
pareto_k2_p1fa=function(x,v1,kscale){
	nx=length(x)
	p1=matrix(0,1,nx)
	vf=Vectorize(pareto_k2_pd)
	p1[1,]=vf(x,v1,v2=kscale)
	return(p1)
}
############################################################
#' The second derivative of the cdf
#' @inheritParams manf
pareto_k2_p2fa=function(x,v1,kscale){
	nx=length(x)
	p2=array(0,c(1,1,nx))
	vf=Vectorize(pareto_k2_pdd)
	p2[1,1,]=vf(x,v1,v2=kscale)
	return(p2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @inheritParams manf
pareto_k2_mu1fa=function(alpha,v1,kscale){
	x=extraDistr::qpareto((1-alpha),a=v1,b=kscale)
	vf=Vectorize(pareto_k2_pd)
	mu1=-vf(x,v1,kscale)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @inheritParams manf
pareto_k2_mu2fa=function(alpha,v1,kscale){
	x=qpareto((1-alpha),a=v1,b=kscale)
	nx=length(x)
	vf=Vectorize(pareto_k2_pdd)
	mu2=-vf(x,v1,kscale)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @inheritParams manf
pareto_k2_ldda=function(x,v1,kscale){
	nx=length(x)
	ldd=matrix(0,1,1)
	vf=Vectorize(pareto_k2_logfdd)
	ldd[1,1]=sum(vf(x,v1,v2=kscale))/nx
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @inheritParams manf
pareto_k2_lddda=function(x,v1,kscale){
	nx=length(x)
	lddd=array(0,c(1,1,1))
	temp=0
	vf=Vectorize(pareto_k2_logfddd)
	lddd[1,1,1]=sum(vf(x,v1,v2=kscale))/nx
	return(lddd)
}
