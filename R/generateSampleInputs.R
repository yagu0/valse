#' Generate a sample of (X,Y) of size n
#' @param meanX matrix of group means for covariates (of size p*K)
#' @param covX covariance for covariates (of size p*p*K)
#' @param covY covariance for the response vector (of size m*m*K)
#' @param pi	 proportion for each cluster
#' @param beta regression matrix
#' @param n		sample size
#'
#' @return list with X and Y
#' @export
generateXY = function(meanX, covX, covY, pi, beta, n)
{
	p = dim(covX)[1]
	m = dim(covY)[1]
	k = dim(covX)[3]

	X = matrix(nrow=n,ncol=p)
	Y = matrix(nrow=n,ncol=m)

	require(MASS) #simulate from a multivariate normal distribution
	for (i in 1:n)
	{
		class = sample(1:k, 1, prob=pi)
		X[i,] = mvrnorm(1, meanX[,class], covX[,,class])
		Y[i,] = mvrnorm(1, X[i,] %*% beta[,,class], covY[,,class])
	}

	return (list(X=X,Y=Y))
}

#' Generate a sample of (X,Y) of size n with default values
#' @param n sample size
#' @param p number of covariates
#' @param m size of the response
#' @param k number of clusters
#' @return list with X and Y
#' @export
generateXYdefault = function(n, p, m, k)
{
	rangeX = 100
	meanX = rangeX * matrix(1 - 2*runif(p*k), ncol=k)
	covX = array(dim=c(p,p,k))
	covY = array(dim=c(m,m,k))
	for(r in 1:k)
	{
		covX[,,r] = diag(p)
		covY[,,r] = diag(m)
	}
	pi = rep(1./k,k)
	#initialize beta to a random number of non-zero random value
	beta = array(0, dim=c(p,m,k))
	for (j in 1:p)
	{
		nonZeroCount = sample(1:m, 1)
		beta[j,1:nonZeroCount,] = matrix(runif(nonZeroCount*k), ncol=k)
	}

	sample_IO = generateXY(meanX, covX, covY, pi, beta, n)
	return (list(X=sample_IO$X,Y=sample_IO$Y))
}

#' Initialize the parameters in a basic way (zero for the conditional mean, uniform for weights,
#' identity for covariance matrices, and uniformly distributed for the clustering)
#' @param n sample size
#' @param p number of covariates
#' @param m size of the response
#' @param k number of clusters
#' @return list with phiInit, rhoInit,piInit,gamInit
#' @export
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

	return (list("phiInit" = phiInit, "rhoInit" = rhoInit, "piInit" = piInit, "gamInit" = gamInit))
}
