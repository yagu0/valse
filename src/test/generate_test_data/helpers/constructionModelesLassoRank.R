constructionModelesLassoRank = function(Pi,Rho,mini,maxi,X,Y,tau,A1,rangmin,rangmax){
  #get matrix sizes
  n = dim(X)[1]
  p = dim(X)[2]
  m = dim(rho)[2]
  k = dim(rho)[3]
  L = dim(A1)[2]
  
  deltaRank = rangmax - rangmin + 1
  Size = deltaRank^k
  Rank = matrix(0, Size, k)
#  for(r in 1:k) {
#    Rank[,r] = rangmin +  <--- #FIXME:
#  }
  
  phi = array(0, dim=c(p,m,k,L*Size))
  lvraisemblance = matrix(0, L*Size, 2)
  for(lambdaIndex in 1:L){
    #on ne garde que les colonnes actives
    #active sera l'ensemble des variables informatives
    active = A1[, lambdaIndex]
    active[active==0] = c()
    if(length(active)>0){
      for(j in 1:Size){
        EMG_rank = EMGrank(Pi[,lambdaIndex], Rho[,,,lambdaIndex], mini, maxi, X[, active], Y, tau, Rank[j,])
        phiLambda = EMG_rank$phi
        LLF = EMG_rank$LLF
        lvraisemblance[(lambdaIndex-1)*Size+j,] = c(LLF, sum(Rank[j,]^(length(active)- Rank[j,]+m)))
        phi[active,,,(lambdaIndex-1)*Size+j] = phiLambda
      }
    }
  }
  return(list(phi=phi, lvraisemblance = lvraisemblance))
}
