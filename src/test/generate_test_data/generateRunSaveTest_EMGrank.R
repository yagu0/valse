generateRunSaveTest_EMGrank = function(n=200, p=15, m=10, k=3, mini=5, maxi=10, gamma=1.0, rank = c(1,2,4)){
  testFolder = "data/"
  dir.create(testFolder, showWarnings=FALSE, mode="0755")
  delimiter = " "
  
  tau = 1e-6
  
  pi = rep(1.0/k, k)
  rho = array(0, dim=c(m,m,k))
  
  for(i in 1:k){
    rho[,,i] = diag(1,m)
  }

  #Generate X and Y
	require(valse)
  generateIOdef = valse:::generateIOdefault(n, p, m, k)
  
  #save inputs
  write.table(paste(testFolder,"rho",sep=""), rho, sep=delimiter)
  write.table(paste(testFolder,"pi",sep=""), pi, sep=delimiter)
  write.table(paste(testFolder,"mini",sep=""), mini, sep=delimiter)
  write.table(paste(testFolder,"maxi",sep=""), maxi, sep=delimiter)
  write.table(paste(testFolder,"X",sep=""), generateIOdef$X sep=delimiter)
  write.table(paste(testFolder,"Y",sep=""), generateIOdef$Y, sep=delimiter)
  write.table(paste(testFolder,"tau",sep=""), tau, sep=delimiter)
  write.table(paste(testFolder,"rank",sep=""), rank, sep=delimiter)
  write.table(paste(testFolder,"dimensions",sep=""), c(n,p,m,k), sep=delimiter)
  
  EMG_rank = EMG(pi,rho,mini,maxi,X,Y,tau,rank)
  
  #save output
  write.table(paste(testFolder,"phi",sep=""), EMG_rank$phi, sep=delimiter)
  write.table(paste(testFolder,"LLF",sep=""), EMG_rank$LLF, sep=delimiter)
}
