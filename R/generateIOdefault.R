#' Generate a sample of (X,Y) of size n with default values
#' @param n sample size
#' @param p number of covariates
#' @param m size of the response
#' @param k number of clusters
#' @return list with X and Y
#' @export
#-----------------------------------------------------------------------
generateIOdefault = function(n, p, m, k)
{
	covX = diag(p)
	covY = array(0, dim=c(m,m,k))
	for(r in 1:k)
	{
		covY[,,r] = diag(m)
	}

	pi = rep(1./k,k)

	beta = array(0, dim=c(p,m,k))
	for(j in 1:p)
	{
		nonZeroCount = ceiling(m * runif(1))
		beta[j,1:nonZeroCount,] = matrix(runif(nonZeroCount*k), ncol=k)
	}

	sample_IO = generateIO(covX, covY, pi, beta, n)
	return (list(X=sample_IO$X,Y=sample_IO$Y))
}
