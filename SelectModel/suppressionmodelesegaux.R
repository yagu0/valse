suppressionmodelesegaux = function(B1,B2,glambda,rho,pi){
  ind = c()
  for(j in 1:length(glambda)){
    for(ll in 1:(l-1)){
      if(B1[,,l] == B1[,,ll]){
        ind = c(ind, l)
      }
    }
  }
  ind = unique(ind)
  B1 = B1[,,-ind]
  glambda = glambda[-ind]
  B2 = B2[,,-ind]
  rho = rho[,,,-ind] 
  pi = pi[,-ind]
  
  return(list(B1,B2,glambda,ind,rho,pi))
}