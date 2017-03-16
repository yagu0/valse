#' Discard models which have the same relevant variables - for EMGLLF
#'
#' @param B1 array of relevant coefficients (of size p*m*length(gridlambda))
#' @param B2 array of irrelevant coefficients (of size p*m*length(gridlambda))
#' @param glambda grid of regularization parameters (vector)
#' @param rho covariance matrix (of size m*m*K*size(gridLambda))
#' @param pi weight parameters (of size K*size(gridLambda))
#'
#' @return a list with update B1, B2, glambda, rho and pi, and ind the vector of indices
#'	of selected models.
#' @export
discardSimilarModels_EMGLLF = function(B1,B2,glambda,rho,pi)
{
	ind = c()
	for (j in 1:length(glambda))
	{
		for (ll in 1:(l-1))
		{
			if(B1[,,l] == B1[,,ll])
				ind = c(ind, l)
		}
	}
	ind = unique(ind)
	B1 = B1[,,-ind]
	glambda = glambda[-ind]
	B2 = B2[,,-ind]
	rho = rho[,,,-ind] 
	pi = pi[,-ind]

	return (list("B1"=B1,"B2"=B2,"glambda"=glambda,"rho"=rho,"pi"=pi,"ind"=ind))
}

#' Discard models which have the same relevant variables
#'   - for Lasso-rank procedure (focus on columns)
#'
#' @param B1 array of relevant coefficients (of size p*m*length(gridlambda))
#' @param rho covariance matrix
#' @param pi weight parameters
#'
#' @return a list with B1, in, rho, pi
#' @export
discardSimilarModels_EMGrank = function(B1,rho,pi)
{
	ind = c()
	dim_B1 = dim(B1)
	B2 = array(0,dim=c(dim_B1[1],dim_B1[2],dim_B1[3]))
	sizeLambda=dim_B1[3]
	glambda = rep(0,sizeLambda)

	suppressmodel = discardSimilarModels_EMGLLF(B1,B2,glambda,rho,pi)
	return (list("B1" = suppressmodel$B1, "ind" = suppressmodel$ind,
		"rho" = suppressmodel$rho, "pi" = suppressmodel$pi))
}
