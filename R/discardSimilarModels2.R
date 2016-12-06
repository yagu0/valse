#' Similar to discardSimilarModels, for Lasso-rank procedure (focus on columns)
#'
#' @param B1 array of relevant coefficients (of size p*m*length(gridlambda))
#' @param rho covariance matrix
#' @param pi weight parameters
#'
#' @return a list with B1, in, rho, pi
#' @export
#'
discardSimilarModels2 = function(B1,rho,pi)
{	ind = c()
	dim_B1 = dim(B1)
	B2 = array(0,dim=c(dim_B1[1],dim_B1[2],dim_B1[3]))
	sizeLambda=dim_B1[3]
	glambda = rep(0,sizeLambda)

	suppressmodel = discardSimilarModels(B1,B2,glambda,rho,pi)
	return (list(B1 = suppressmodel$B1, ind = suppressmodel$ind,
		rho = suppressmodel$rho, pi = suppressmodel$pi))
}
