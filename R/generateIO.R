generateIO = function(covX, covY, pi, beta, n)
{
	size_covX = dim(covX)
	p = size_covX[1]
	k = size_covX[3]

	size_covY = dim(covY)
	m = size_covY[1]

	Y = matrix(0,n,m)
	BX = array(0, dim=c(n,m,k))

	require(MASS) #simulate from a multivariate normal distribution
	for (i in 1:n)
	{
		for (r in 1:k)
		{
			BXir = rep(0,m)
			for (mm in 1:m)
				Bxir[[mm]] = X[i,] %*% beta[,mm,r]
			Y[i,] = Y[i,] + pi[r] * mvrnorm(1,BXir, covY[,,r])
		}
	}

	return (list(X=X,Y=Y))
}
