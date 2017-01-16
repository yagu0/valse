generateRunSaveTest_constructionModelesLassoRank = function(n=200, p=15, m=10, L=12, mini=5, maxi=10, gamma=1.0, rangmin=3, rangmax=6){
  testFolder = "data/"
  dir.create(testFolder, showWarnings=FALSE, mode="0755")
  delimiter = " "
  
  tau = 1e-6
  
  
  pi = matrix(0, k,L)
  for(i in 1:L){
    pi[,i] = rep(1.0/k, k)
  }
  rho = array(0, dim=c(m,m,k,L))
  for(l in 1:L){
    for(r in 1:k){
      rho[,,r,l] = diag(1,m)
    }
  }
  #Generate X and Y
  generateIOdef = generateIOdefault(n, p, m, k)
  
  A1 = matrix(0,p,L)
  for(i in 1:L){
    A1[,i] = seq(1,p)
  }
  #save inputs
  write.table(paste(testFolder,"rho",sep=""), rho, sep=delimiter)
  write.table(paste(testFolder,"pi",sep=""), pi, sep=delimiter)
  write.table(paste(testFolder,"mini",sep=""), mini, sep=delimiter)
  write.table(paste(testFolder,"maxi",sep=""), maxi, sep=delimiter)
  write.table(paste(testFolder,"X",sep=""), generateIOdef$X sep=delimiter)
  write.table(paste(testFolder,"Y",sep=""), generateIOdef$Y, sep=delimiter)
  write.table(paste(testFolder,"tau",sep=""), tau, sep=delimiter)
  write.table(paste(testFolder,"A1",sep=""), A1, sep=delimiter)
  write.table(paste(testFolder,"rangmin",sep=""), rangmin, sep=delimiter)
  write.table(paste(testFolder,"rangmax",sep=""), rangmax, sep=delimiter)
  write.table(paste(testFolder,"dimensions",sep=""), c(n,p,m,k), sep=delimiter)
  
  construct = constructionModelesLassoRank(pi,rho,mini,maxi,X,Y,tau,A1,rangmin,rangmax))
  
  #save output
  write.table(paste(testFolder,"phi",sep=""), construct$phi, sep=delimiter)
  write.table(paste(testFolder,"lvraisemblance",sep=""), construct$lvraisemblance, sep=delimiter)
  
}