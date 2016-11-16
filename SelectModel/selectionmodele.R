vec_bin = function(X,r){
  Z = c()
  indice = c()
  j=1
  for(i in 1:length(X)){
    if(X[i] == r){
      Z[i] = 1
      indice[j] = i
      j=j+1
    }
    else{
      Z[i] = 0
    }
  }
  return(list(Z,indice))
}

selectionmodele = function(vraisemblance){
  D = vraimsemblance[,2]
  D1 = unique(D)
  
  indice = rep(1, length(D1))
  
  #select argmax MLE
  if(length(D1)>2){
    for(i in 1:length(D1)){
      A = c()
      for(j in 1:length(D)){
        if(D[[j]]==D1[[i]]){
          a = c(a, vraimsemblance[j,1])
        }
      }
      b = max(a)
      indice[i] = which.max(vec_bin(vraimsemblance,b)[[1]]) #retourne le premier indice du vecteur binaire où u_i ==1
    }
  }
  return(list(indice,D1))
}