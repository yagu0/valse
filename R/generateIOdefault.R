generateIOdefault = function(n, p, m, k)
{
	covX = array(0, dim=c(p,p,k))
	covY = array(0, dim=c(m,m,k))
	for(r in 1:k)
	{
		covX[,,r] = diag(p)
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
