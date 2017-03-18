source("EMGrank.R")

generateRunSaveTest_EMGrank = function(n=200, p=15, m=10, k=3, mini=5, maxi=10, gamma=1.0,
	rank = c(1,2,4))
{
  tau = 1e-6
  pi = rep(1.0/k, k)
  rho = array(dim=c(m,m,k))
  for(i in 1:k)
    rho[,,i] = diag(1,m)
	require(valse)
  xy = valse:::generateXYdefault(n, p, m, k)

  testFolder = "../data/"
  dir.create(testFolder, showWarnings=FALSE, mode="0755")
  #save inputs
  write.table(as.double(rho), paste(testFolder,"rho",sep=""),
		row.names=F, col.names=F)
	write.table(as.double(pi), paste(testFolder,"pi",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(mini), paste(testFolder,"mini",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(maxi), paste(testFolder,"maxi",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(xy$X), paste(testFolder,"X",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(xy$Y), paste(testFolder,"Y",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(tau), paste(testFolder,"tau",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(rank), paste(testFolder,"rank",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(c(n,p,m,k)), paste(testFolder,"dimensions",sep=""),
		row.names=F, col.names=F)

  res = EMGrank_R(pi,rho,mini,maxi,xy$X,xy$Y,tau,rank)

  #save output
  write.table(as.double(res$phi), paste(testFolder,"phi",sep=""), row.names=F,col.names=F)
  write.table(as.double(res$LLF), paste(testFolder,"LLF",sep=""), row.names=F,col.names=F)
}
