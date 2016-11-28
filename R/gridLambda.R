gridLambda = function(phiInit, rhoInit, piInit, gamInit, X, Y, gamma, mini, maxi, tau){
  n = nrow(X)
  p = dimension(phiInit)[1]
  m = dimension(phiInit)[2]
  k = dimension(phiInit)[3]
  list_EMG = EMGLLF(phiInit,rhoInit,piInit,gamInit,mini,maxi,1,0,X,Y,tau)
  #.C("EMGLLF", phiInit = phiInit, rhoInit = rhoInit, ...)
  phi = list_EMG[[1]]
  rho = list_EMG[[2]]
  pi = list_EMG[[3]]
  S = list_EMG[[5]]
  
  grid = array(0, dim=c(p,m,k))
  for(i in 1:p){
    for(j in 1:m){
      grid[i,j,] = abs(S[i,j,]) / (n*pi^gamma)
    }
  }
  grid = unique(grid)
  grid = grid[grid <=1 ]
  
  return(grid)
}


#test pour voir si formatage à la fin de grid ok
grid= array(mvrnorm(5*5*2,1,1), dim=c(5,5,2))
grid = unique(grid)
grid = grid[grid<= 1 ]
