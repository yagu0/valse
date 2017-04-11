#' Generate a sample of (X,Y) of size n with default values
#'
#' @param n sample size
#' @param p number of covariates
#' @param m size of the response
#' @param k number of clusters
#'
#' @return list with X and Y
#'
generateXYdefault = function(n, p, m, k)
{
	meanX = rep(0, p)
	covX = diag(p)
	covY = array(dim=c(m,m,k))
	for(r in 1:k)
		covY[,,r] = diag(m)
	π = rep(1./k,k)
	#initialize beta to a random number of non-zero random value
	β = array(0, dim=c(p,m,k))
	for (j in 1:p)
	{
		nonZeroCount = sample(1:m, 1)
		β[j,1:nonZeroCount,] = matrix(runif(nonZeroCount*k), ncol=k)
	}

	sample_IO = generateXY(n, π, meanX, β, covX, covY)
	return (list(X=sample_IO$X,Y=sample_IO$Y))
}

#' Initialize the parameters in a basic way (zero for the conditional mean, uniform for
#' weights, identity for covariance matrices, and uniformly distributed for the
#' clustering)
#'
#' @param n sample size
#' @param p number of covariates
#' @param m size of the response
#' @param k number of clusters
#'
#' @return list with phiInit, rhoInit,piInit,gamInit
#'
basicInitParameters = function(n,p,m,k)
{
	phiInit = array(0, dim=c(p,m,k))

	piInit = (1./k)*rep(1,k)

	rhoInit = array(dim=c(m,m,k))
	for (i in 1:k)
		rhoInit[,,i] = diag(m)

	gamInit = 0.1 * matrix(1, nrow=n, ncol=k)
	R = sample(1:k, n, replace=TRUE)
	for (i in 1:n)
		gamInit[i,R[i]] = 0.9
	gamInit = gamInit/sum(gamInit[1,])

	return (list("phiInit"=phiInit, "rhoInit"=rhoInit, "piInit"=piInit, "gamInit"=gamInit))
}
