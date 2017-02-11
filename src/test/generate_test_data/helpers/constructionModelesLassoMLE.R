constructionModelesLassoMLE = function(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda,X,Y,seuil,tau,A1,A2){
  #get matrix sizes
  n = dim(X)[1];
  p = dim(phiInit)[1]
  m = dim(phiInit)[2]
  k  = dim(phiInit)[3]
  L = length(glambda)
  
  #output parameters
  phi = array(0, dim=c(p,m,k,L))
  rho = array(0, dim=c(m,m,k,L))
  Pi = matrix(0, k, L)
  llh = matrix(0, L, 2) #log-likelihood

  for(lambdaIndex in 1:L){
    a = A1[, 1, lambdaIndex]
    a[a==0] = c()
    if(length(a)==0){
      next
    }
    EMGLLf = EMGLLF(phiInit[a,,],rhoInit,piInit,gamInit,mini,maxi,gamma,0,X[,a],Y,tau)
    
    phiLambda = EMGLLf$phi
    rhoLambda = EMGLLf$rho
    piLambda = EMGLLf$Pi
    
    for(j in 1:length(a)){
      phi[a[j],,,lambdaIndex] = phiLambda[j,,]
    }
    rho[,,,lambdaIndex] = rhoLambda
    Pi[,lambdaIndex] = piLambda
    
    dimension = 0
    for(j in 1:p){
      vec =  c(2, dim(A2)[2])
      b = A2[j,vec,lambdaIndex]
      b[b==0] = c()
      if(length(b) > 0){
        phi[A2[j,1,lambdaIndex],b,,lambdaIndex] = 0
      }
      c = A1[j,vec,lambdaIndex]
      c[c==0] = c()
      dimension = dimension + length(c)
    }
    
    #on veut calculer l'EMV avec toutes nos estimations
		densite = matrix(0, n, L)
		for(i in 1:n){
			for( r in 1:k){
				delta = Y[i,]%*%rho[,,r,lambdaIndex] - (X[i,a]%*%phi[a,,r,lambdaIndex]);
				densite[i,lambdaIndex] = densite[i,lambdaIndex] +	Pi[r,lambdaIndex]*det(rho[,,r,lambdaIndex])/(sqrt(2*pi))^m*exp(-tcrossprod(delta)/2.0)
			}
		}
		llh[lambdaIndex,1] = sum(log(densite[,lambdaIndex]))
		llh[lambdaIndex,2] = (dimension+m+1)*k-1
  }
  return(list(phi=phi, rho=rho, Pi=Pi, llh = llh))
}
