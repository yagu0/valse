#Return a list of outputs, for each lambda in grid: selected,Rho,Pi
selectiontotale = function(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda,X,Y,thresh,tau, parallel = FALSE){
  if (parallel) {
    require(parallel)
    cl = parallel::makeCluster( parallel::detectCores() / 4) # <-- ça devrait être un argument
    parallel::clusterExport(cl=cl,
                            varlist=c("phiInit","rhoInit","gamInit","mini","maxi","glambda","X","Y","thresh","tau"),
                            envir=environment())
    #Pour chaque lambda de la grille, on calcule les coefficients
    out = parLapply(cl,  1:length(glambda), function(lambdaIndex)
    {
      params = 
        EMGLLF(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda[lambdaIndex],X,Y,tau)
      
      p = dim(phiInit)[1]
      m = dim(phiInit)[2]
      #selectedVariables: list where element j contains vector of selected variables in [1,m]
      selectedVariables = lapply(1:p, function(j) {
        #from boolean matrix mxk of selected variables obtain the corresponding boolean m-vector,
        #and finally return the corresponding indices
        seq_len(m)[ apply( abs(params$phi[j,,]) > thresh, 1, any ) ]
      })
      
      list("selected"=selectedVariables,"Rho"=params$Rho,"Pi"=params$Pi)
    })
    parallel::stopCluster(cl)
  }
  else {
    selectedVariables = list()
    Rho = list()
    Pi = list()
    #Pour chaque lambda de la grille, on calcule les coefficients
    for (lambdaIndex in 1:length(glambda)){
      print(lambdaIndex)
      params = 
        EMGLLF(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda[lambdaIndex],X,Y,tau)
      p = dim(phiInit)[1]
      m = dim(phiInit)[2]
      #selectedVariables: list where element j contains vector of selected variables in [1,m]
      selectedVariables[[lambdaIndex]] = sapply(1:p, function(j) {
        #from boolean matrix mxk of selected variables obtain the corresponding boolean m-vector,
        #and finally return the corresponding indices
        c(seq_len(m)[ apply( abs(params$phi[j,,]) > thresh, 1, any ) ], rep(0, m-length(apply( abs(params$phi[j,,]) > thresh, 1, any ) )))
      })
      Rho[[lambdaIndex]] = params$Rho
      Pi[[lambdaIndex]] = params$Pi
    }
    list("selected"=selectedVariables,"Rho"=Rho,"Pi"=Pi)
  }
}