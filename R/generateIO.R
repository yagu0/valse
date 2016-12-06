#' Generate a sample of (X,Y) of size n
#' @param covX covariance for covariates (of size p*p*K)
#' @param covY covariance for the response vector (of size m*m*K)
#' @param pi   proportion for each cluster
#' @param beta regression matrix
#' @param n    sample size
#' 
#' @return list with X and Y
#' @export
#-----------------------------------------------------------------------
generateIO = function(covX, covY, pi, beta, n)
{
  p = dim(covX)[1]
  
  m = dim(covY)[1]
  k = dim(covY)[3]
  
  Y = matrix(0,n,m)
  require(mvtnorm)
  X = rmvnorm(n, mean = rep(0,p), sigma = covX)
  
  require(MASS) #simulate from a multivariate normal distribution
  for (i in 1:n)
  {
    
    for (r in 1:k)
    {
      BXir = rep(0,m)
      for (mm in 1:m)
        BXir[mm] = X[i,] %*% beta[,mm,r]
      Y[i,] = Y[i,] + pi[r] * mvrnorm(1,BXir, covY[,,r])
    }
  }
  
  return (list(X=X,Y=Y))
}
