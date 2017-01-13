covariance= function(p,a){
  A = matrix(a, p,p)
  for{i in 1:p}{
    for{k  in 1:p}{
      A[i,]= A[i,]^abs(i-k) 
    }
  }
  return(A=A)
}