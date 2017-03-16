constructionModelesLassoMLE = function(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda,X,Y,seuil,tau,selected){
  #get matrix sizes
  n = dim(X)[1]
  p = dim(phiInit)[1]
  m = dim(phiInit)[2]
  k = dim(phiInit)[3]
  L = length(glambda)
  
  #output parameters
  phi = array(0, dim = c(p,m,k,L))
  rho = array(dim = c(m,m,k,L))
  pi = array( dim = c(k,L))
  lvraisemblance = array( dim = c(L,2))
  
  for (lambdaIndex in 1:L){
    # Procedure Lasso-MLE  
    a = selected[,1,lambdaIndex]
    a(a==0) = c()
    if (length(a) != 0){
      res_EM = EMGLLF(phiInit[a,,],rhoInit,piInit,gamInit,mini,maxi,gamma,0,X[,a],Y,tau)
      phiLambda = res_EM$phi
      rhoLambda = res_EM$rho
      piLambda = res_EM$pi
      for (j in 1:length(a)){
        phi[a[j],,,lambdaIndex] = phiLambda[j,,]
      }
      rho[,,,lambdaIndex] = rhoLambda
      pi[,lambdaIndex] = piLambda
      
      dimension = 0
      for (j in 1:p){
        b = A2[j,2:end,lambdaIndex]
        b[b==0] = c()
        if (length(b) > 0){
        phi[A2[j,1,lambdaIndex],b,,lambdaIndex] = 0.0
        }
        c = A1[j,2:end,lambdaIndex]
        c[c==0] = c()
        dimension = dimension + length(c)
      }
      
      #on veut calculer l'EMV avec toutes nos estimations
      densite = array(n,L)
      for (i in 1:n){
        for (r in 1:k){
          delta = Y[i,]*rho[,,r,lambdaIndex] - X[i,a]*phi[a,,r,lambdaIndex]
          densite[i,lambdaIndex] = densite[i,lambdaIndex] + pi[r,lambdaIndex]*det(rho[,,r,lambdaIndex])/(sqrt(2*base::pi))^m*exp(-delta %*% delta/2.0)
        }
      }
      
      lvraisemblance(lambdaIndex,1) = sum(log(densite[,lambdaIndex]))
      lvraisemblance(lambdaIndex,2) = (dimension+m+1)*k-1
    }
  }
}