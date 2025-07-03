#' Jeffreys' Prior with two parameters
#' @param ggd 	gradient of the expected information matrix
#' @param detg	determinant of the expected information matrix
#' @param ggi 	inverse of the expected information matrix
#'
#' @return Vector of 2 values
jpf2p=function(ggd,detg,ggi){
	dgd1=ggd[1,,]
	dgd2=ggd[2,,]
	lambdad_jp=matrix(0,2)
	if(!is.na(detg)){
		if(detg>=0){
			lambdad_jp[1]=0.5*sqrt(max(0,detg))*sum(diag(ggi%*%dgd1))
			lambdad_jp[2]=0.5*sqrt(max(0,detg))*sum(diag(ggi%*%dgd2))
		}
	}
	return(lambdad_jp)
}
#' Jeffreys' Prior with three parameters
#' @param ggd 	gradient of the expected information matrix
#' @param detg	determinant of the expected information matrix
#' @param ggi 	inverse of the expected information matrix
#'
#' @return Vector of 3 values
jpf3p=function(ggd,detg,ggi){
	dgd1=ggd[1,,]
	dgd2=ggd[2,,]
	dgd3=ggd[3,,]
	lambdad_jp=matrix(0,3)
	if(!is.na(detg)){
		if(detg>=0){
			lambdad_jp[1]=0.5*sqrt(max(0,detg))*sum(diag(ggi%*%dgd1))
			lambdad_jp[2]=0.5*sqrt(max(0,detg))*sum(diag(ggi%*%dgd2))
			lambdad_jp[3]=0.5*sqrt(max(0,detg))*sum(diag(ggi%*%dgd3))
		}
	}
	return(lambdad_jp)
}
#' Jeffreys' Prior with four parameters
#' @param ggd 	gradient of the expected information matrix
#' @param detg	determinant of the expected information matrix
#' @param ggi 	inverse of the expected information matrix
#' @return Vector of 4 values
jpf4p=function(ggd,detg,ggi){
	dgd1=ggd[1,,]
	dgd2=ggd[2,,]
	dgd3=ggd[3,,]
	dgd4=ggd[4,,]
	lambdad_jp=matrix(0,4)
	if(!is.na(detg)){
		if(detg>=0){
			lambdad_jp[1]=0.5*sqrt(max(0,detg))*sum(diag(ggi%*%dgd1))
			lambdad_jp[2]=0.5*sqrt(max(0,detg))*sum(diag(ggi%*%dgd2))
			lambdad_jp[3]=0.5*sqrt(max(0,detg))*sum(diag(ggi%*%dgd3))
			lambdad_jp[4]=0.5*sqrt(max(0,detg))*sum(diag(ggi%*%dgd4))
		}
	}
	return(lambdad_jp)
}
