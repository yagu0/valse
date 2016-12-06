#' Discard models which have the same relevant variables
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
discardSimilarModels = function(B1,B2,glambda,rho,pi)
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

	return (list(B1=B1,B2=B2,glambda=glambda,rho=rho,pi=pi,ind=ind))
}
