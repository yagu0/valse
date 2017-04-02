#' generateXY
#'
#' Generate a sample of (X,Y) of size n
#'
#' @param n sample size
#' @param π proportion for each cluster
#' @param meanX matrix of group means for covariates (of size p)
#' @param covX covariance for covariates (of size p*p)
#' @param β regression matrix, of size p*m*k
#' @param covY covariance for the response vector (of size m*m*K)
#'
#' @return list with X and Y
#'
#' @export
generateXY = function(n, π, meanX, β, covX, covY)
{
	p <- dim(covX)[1]
	m <- dim(covY)[1]
	k <- dim(covY)[3]
	
	X <- matrix(nrow=0, ncol=p)
	Y <- matrix(nrow=0, ncol=m)

	#random generation of the size of each population in X~Y (unordered)
	sizePop <- rmultinom(1, n, pi)
	class <- c() #map i in 1:n --> index of class in 1:k

	for (i in 1:k)
	{
		class <- c(class, rep(i, sizePop[i]))
		newBlockX <- MASS::mvrnorm(sizePop[i], meanX, covX)
		X <- rbind( X, newBlockX )
		Y <- rbind( Y, apply( newBlockX, 1, function(row)
			mvrnorm(1, row %*% beta[,,i], covY[,,i]) ) )
	}

	shuffle = sample(n)
	list("X"=X[shuffle,], "Y"=Y[shuffle,], "class"=class[shuffle])
}