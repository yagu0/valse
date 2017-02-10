generateRunSaveTest_constructionModelesLassoMLE = function(n=200, p=15, m=10, k=3, mini=5, maxi=10, gamma=1.0, glambda=list(0.0,0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.5,0.7,0.85,0.99)){
  testFolder = "data/"
  dir.create(testFolder, showWarnings=FALSE, mode="0755")
  delimiter = " "
  
  #Generate phiInit,piInit,...
	require(valse)
  params = valse:::basic_Init_Parameters(n, p, m, k)

  #save inputs
  write.table(paste(testFolder,"phiInit",sep=""), params$phiInit, sep=delimiter)
  write.table(paste(testFolder,"rhoInit",sep=""), params$rhoInit, sep=delimiter)
  write.table(paste(testFolder,"piInit",sep=""), params$piInit, sep=delimiter)
  write.table(paste(testFolder,"gamInit",sep=""), params$gamInit, sep=delimiter)
  write.table(paste(testFolder,"mini",sep=""), mini, sep=delimiter)
  write.table(paste(testFolder,"maxi",sep=""), maxi, sep=delimiter)
  write.table(paste(testFolder,"gamma",sep=""), gamma, sep=delimiter)
  write.table(paste(testFolder,"lambda",sep=""), lambda, sep=delimiter)
  write.table(paste(testFolder,"X",sep=""), io$X, sep=delimiter)
  write.table(paste(testFolder,"Y",sep=""), io$Y, sep=delimiter)
  write.table(paste(testFolder,"tau",sep=""), tau, sep=delimiter)
  write.table(paste(testFolder,"dimensions",sep=""), c(n,p,m,k), sep=delimiter)
  
  A2 = array(0, dim=c(p, m+1, L))
  A1 = array(0, dim=c(p, m+1, L))
  for(i in 1:L){
    A2[,1,i] = seq(1,p)
    A2[,2,i] = seq(1,5)
    A1[,1, i] = seq(1,p)
    A1[,2,i] = seq(1,5)
  }
  
  #Generate X and Y
  generateIOdef = generateIOdefault(n, p, m, k)
  construct_LME = constructionModelesLassoMLE(params$phiInit,params$rhoInit,params$piInit,params$gamInit,mini,maxi,gamma,glambda,generateIOdef$X,generateIOdef$Y,seuil,tau,A1,A2)
  phi = construct_LME$phi
  rho = construct_LME$rho
  pi = construct_LME$pi
  lvraisemblance = construct_LME$lvraisemblance
  
  #save output
  write.table(paste(testFolder,"phi",sep=""), construct_LME$phi, sep=delimiter)
  write.table(paste(testFolder,"rho",sep=""), construct_LME$rho, sep=delimiter)
  write.table(paste(testFolder,"pi",sep=""), construct_LME$pi, sep=delimiter)
  write.table(paste(testFolder,"lvraisemblance",sep=""), construct_LME$lvraisemblance, sep=delimiter)
}
