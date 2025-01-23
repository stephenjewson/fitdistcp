#' Evaluate DMGS equation 3.3
#' @param lddi inverse of second derivative of observed log-likelihood
#' @param lddd third derivative of observed log-likelihood
#' @param mu1 DMGS mu1 vector
#' @param pidopi1 first part of the prior term
#' @param pidopi2 second part of the prior term
#' @param mu2 DMGS mu2 matrix
#' @param dim number of parameters
bayesian_dq_4terms_v1=function(lddi,lddd,mu1,pidopi1,pidopi2,mu2,dim){
	nalpha=length(mu1[1,])
	beps1=rep(0,nalpha)
	beps2=rep(0,nalpha)
	beps3=rep(0,nalpha)
	beps4=rep(0,nalpha)
# beps1
	for (t in 1:dim){
		for (j in 1:dim){
			for (r in 1:dim){
				for (s in 1:dim){
					beps1=beps1+0.5*lddi[s,t]*lddi[j,r]*lddd[j,r,s]*mu1[t,]
				}
			}
		}
	}
# beps2
	for (t in 1:dim){
 	 for (s in 1:dim){
			beps2=beps2-lddi[s,t]*pidopi1[s]*mu1[t,] #this has to be negative because ldd is the opposite sign to c
		}
	}
# beps3
	for (t in 1:dim){
 	 for (s in 1:dim){
			beps3=beps3-lddi[s,t]*pidopi2[s]*mu1[t,] #this has to be negative because ldd is the opposite sign to c
		}
	}
# beps4
	for (j in 1:dim){
		for (r in 1:dim){
			beps4=beps4-0.5*lddi[j,r]*mu2[j,r,] #this has to be negative because ldd is the opposite sign to c
		}
	}
	beps=beps1+beps2+beps3+beps4
	return(beps)
}
#' Evaluate DMGS equation 3.3
#' @param lddi inverse of second derivative of observed log-likelihood
#' @param lddd third derivative of observed log-likelihood
#' @param mu1 DMGS mu1 vector
#' @param pidopi derivative of log prior
#' @param mu2 DMGS mu2 matrix
#' @param dim number of parameters
dmgs=function(lddi,lddd,mu1,pidopi,mu2,dim){
	pidopi1=2*pidopi
	pidopi2=-pidopi
	nalpha=length(mu1[1,])
	beps1=rep(0,nalpha)
	beps2=rep(0,nalpha)
	beps3=rep(0,nalpha)
	beps4=rep(0,nalpha)
# beps1
	for (t in 1:dim){
		for (j in 1:dim){
			for (r in 1:dim){
				for (s in 1:dim){
					beps1=beps1+0.5*lddi[s,t]*lddi[j,r]*lddd[j,r,s]*mu1[t,]
				}
			}
		}
	}
# beps2
	for (t in 1:dim){
 	 for (s in 1:dim){
			beps2=beps2-lddi[s,t]*pidopi1[s]*mu1[t,] #this has to be negative because ldd is the opposite sign to c
		}
	}
# beps3
	for (t in 1:dim){
 	 for (s in 1:dim){
			beps3=beps3-lddi[s,t]*pidopi2[s]*mu1[t,] #this has to be negative because ldd is the opposite sign to c
		}
	}
# beps4
	for (j in 1:dim){
		for (r in 1:dim){
			beps4=beps4-0.5*lddi[j,r]*mu2[j,r,] #this has to be negative because ldd is the opposite sign to c
		}
	}
	beps=beps1+beps2+beps3+beps4
	return(beps)
}

