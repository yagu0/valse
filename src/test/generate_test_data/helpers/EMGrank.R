require(MASS)
EMGrank = function(Pi, Rho, mini, maxi, X, Y, tau, rank){
  #matrix dimensions
  n = dim(X)[1]
  p = dim(X)[2]
  m = dim(Rho)[2]
  k = dim(Rho)[3]
  
  #init outputs
  phi = array(0, dim=c(p,m,k))
  Z = rep(1, n)
  Pi = piInit
  LLF = 0
  
  #local variables
  Phi = array(0, dim=c(p,m,k))
  deltaPhi = c(0)
  sumDeltaPhi = 0
  deltaPhiBufferSize = 20
  
  #main loop
  ite = 1
  while(ite<=mini || (ite<=maxi && sumDeltaPhi>tau)){
    #M step: Mise à jour de Beta (et donc phi)
    for(r in 1:k){
      Z_bin = valse:::vec_bin(Z,r)
      Z_vec = Z_bin$vec #vecteur 0 et 1 aux endroits o? Z==r
      Z_indice = Z_bin$indice 
      if(sum(Z_indice) == 0){
        next
      }
      #U,S,V = SVD of (t(Xr)Xr)^{-1} * t(Xr) * Yr
      sv = svd(ginv( crossprod(X[Z_indice,]) )   %*%   crossprod(X[Z_indice,], Y[Z_indice,])     )
      S = diag(sv$d)
      U = sv$u
      V = sv$v
      #Set m-rank(r) singular values to zero, and recompose
      #best rank(r) approximation of the initial product
      if(r==k){
        j_r_1 = length(S)
      }
      else{
        j_r_1 = c(rank[r]+1:length(S))
      }
      S[j_r_1] = 0
      S = diag(S, nrow = ncol(U))
      phi[,,r] = U %*% S %*% t(V) %*% Rho[,,r]
    }
  
  #Etape E et calcul de LLF
  sumLogLLF2 = 0
  for(i in 1:n){
    sumLLF1 = 0
    maxLogGamIR = -Inf
    for(r in 1:k){
      dotProduct = tcrossprod(Y[i,]%*%Rho[,,r]-X[i,]%*%phi[,,r])
      logGamIR = log(Pi[r]) + log(det(Rho[,,r])) - 0.5*dotProduct
      #Z[i] = index of max (gam[i,])
      if(logGamIR > maxLogGamIR){
        Z[i] = r
        maxLogGamIR = logGamIR
      }
    sumLLF1 = sumLLF1 + exp(logGamIR) / (2*pi)^(m/2)
    }
    sumLogLLF2 = sumLogLLF2 + log(sumLLF1)
  }
  
  LLF = -1/n * sumLogLLF2
  
  #update distance parameter to check algorithm convergence (delta(phi, Phi))
  deltaPhi = c(deltaPhi, max(max(max((abs(phi-Phi))/(1+abs(phi))))) )
  if(length(deltaPhi) > deltaPhiBufferSize){
    l_1 = c(2:length(deltaPhi))
    deltaPhi = deltaPhi[l_1]
  }
  sumDeltaPhi = sum(abs(deltaPhi))
  
  #update other local variables
  Phi = phi
  ite = ite+1
  
  }
  return(list(phi=phi, LLF=LLF))
}
