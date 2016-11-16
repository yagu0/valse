suppressionmodelesegaux2 = function(B1,rho,pi){
  ind = c()
  dim_B1 = dim(B1)
  B2 = array(0,dim=c(dim_B1[1],dim_B1[2],dim_B1[3]))
  nombreLambda=dim_B1[[2]]
  glambda = rep(0,nombreLambda)
  
  #for(j in 1:nombreLambda){
  #  for(ll in 1:(l-1)){
  #    if(B1[,,l] == B1[,,ll]){
  #      ind = c(ind, l)
  #    }
  #  }
  #}
  #ind = unique(ind)
  #B1 = B1[,,-ind]
  #rho = rho[,,,-ind] 
  #pi = pi[,-ind]
  
  suppressmodel = suppressionmodelesegaux(B1,B2,glambda,rho,pi)
  B1 = suppressmodel[[1]]
  ind = suppressmodel[[4]]
  rho = suppressmodel[[5]]
  pi = suppressmodel[[6]]
  return(list(B1,ind,rho,pi))
}