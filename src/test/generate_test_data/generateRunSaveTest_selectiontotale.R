generateRunSaveTest_selectiontotale= function(n=200, p=15, m=10, k=3, mini=5, maxi=10, gamma=1.0, glambda=list(0.0,0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.5,0.7,0.85,0.99)){
  testFolder = "data/"
  dir.create(testFolder, showWarnings=FALSE, mode="0755")
  delimiter = " "
  
	require(valse)

  #Generate phiInit,piInit,...
  params = valse:::basic_Init_Parameters(n, p, m, k)
  
  #Generate X and Y
  generateIOdef = valse:::generateIOdefault(n, p, m, k)
  
  #save inputs
  write.table(paste(testFolder,"phiInit",sep=""), params$phiInit, sep=delimiter)
  write.table(paste(testFolder,"rhoInit",sep=""), params$rhoInit, sep=delimiter)
  write.table(paste(testFolder,"piInit",sep=""), params$piInit, sep=delimiter)
  write.table(paste(testFolder,"gamInit",sep=""), params$gamInit, sep=delimiter)
  write.table(paste(testFolder,"mini",sep=""), mini, sep=delimiter)
  write.table(paste(testFolder,"maxi",sep=""), maxi, sep=delimiter)
  write.table(paste(testFolder,"gamma",sep=""), gamma, sep=delimiter)
  write.table(paste(testFolder,"lambda",sep=""), glambda, sep=delimiter)
  write.table(paste(testFolder,"X",sep=""), io$X, sep=delimiter)
  write.table(paste(testFolder,"Y",sep=""), io$Y, sep=delimiter)
  write.table(paste(testFolder,"tau",sep=""), tau, sep=delimiter)
  write.table(paste(testFolder,"dimensions",sep=""), c(n,p,m,k), sep=delimiter)
  
  
  selec = selectiontotale(params$phiInit,params$rhoInit,params$piInit,params$gamInit,mini,maxi,gamma,glambda,generateIOdef$X,generateIOdef$Y,seuil, tau)
  
  #save output
  write.table(paste(testFolder,"A1",sep=""), selec$A1, sep=delimiter)
  write.table(paste(testFolder,"A2",sep=""), selec$A2, sep=delimiter)
  write.table(paste(testFolder,"rho",sep=""), selec$rho, sep=delimiter)
  write.table(paste(testFolder,"pi",sep=""), selec$pi, sep=delimiter)
}
