generateIOdefault = function(n, p, m, k){
  rangeX = 100
  meanX = rangeX*(1-matrix(runif(k*p),ncol = p))
  
  covX = array(0, dim=c(p,p,k))
  covY = array(0, dim=c(m,m,k))
  
  for(r in 1:k){
    covX[,,r] = diag(p)
    covY[,,r] = diag(m)
  }
  
  pi = (1/k) * rep(1,k)
  
  beta = array(0, dim=c(p,m,k))
  
  for(j in 1:p){
    nonZeroCount = ceiling(m * runif(1))
    beta[j,1:nonZeroCount,] = matrix(runif(nonZeroCount*k),ncol = k)
  }
  
  generate = generateIO(meanX, covX, covY, pi, beta, n)
  
  return(list(generate[[1]],generate[[2]]))
}