source("helper.R")
library(valse)

generateRunSaveTest_EMGrank = function(n=200, p=15, m=10, k=3, mini=5, maxi=10, gamma=1.0, rank = c(1,2,4))
{
  eps = 1e-6
  Pi = rep(1.0/k, k)
  Rho = array(dim=c(m,m,k))
  for(i in 1:k)
    Rho[,,i] = diag(1,m)
  xy = generateXYdefault(n, p, m, k)

  testFolder = "./data/"
  dir.create(testFolder, showWarnings=FALSE, mode="0755")
  #save inputs
	write.table(as.double(Pi), paste(testFolder,"Pi",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(Rho), paste(testFolder,"Rho",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(mini), paste(testFolder,"mini",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(maxi), paste(testFolder,"maxi",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(xy$X), paste(testFolder,"X",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(xy$Y), paste(testFolder,"Y",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(eps), paste(testFolder,"eps",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(rank), paste(testFolder,"rank",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(c(n,p,m,k)), paste(testFolder,"dimensions",sep=""),
		row.names=F, col.names=F)

  res = valse::EMGrank(Pi,Rho,mini,maxi,xy$X,xy$Y,eps,rank,fast=FALSE)

  #save output
  write.table(as.double(res$phi),paste(testFolder,"phi",sep=""),row.names=F,col.names=F)
  write.table(as.double(res$LLF),paste(testFolder,"LLF",sep=""),row.names=F,col.names=F)
}
