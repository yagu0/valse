basic_Init_Parameters = function(n,p,m,k){
  phiInit = array(0, dim=c(p,m,k))
  
  piInit = (1.0/k)*rep.int(1,k)
  
  rhoInit = array(0, dim=c(m,m,k))
  
  for(i in 1:k){
    rhoInit[,,i] = diag(m)
  }
  
  gamInit = 0.1*array(1, dim=c(n,k))
  
  R = sample(1:k,n, replace= TRUE)
  
  for(i in 1:n){
    gamInit[i,R[i]] = 0.9
  }
  gamInit = gamInit/sum(gamInit[1,])
  
  
  return(list(phiInit, rhoInit, piInit, gamInit))
}

n= 10
p = 10
m = 5
k = 5
list_param = basic_Init_Parameters(n,p,m,k)
