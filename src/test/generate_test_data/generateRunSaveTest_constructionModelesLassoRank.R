source("helpers/constructionModelesLassoRank.R")

generateRunSaveTest_constructionModelesLassoRank = function(n=200, p=15, m=10, k=3, L=12, mini=5,
	maxi=10, gamma=1.0, rangmin=3, rangmax=6)
{
  tau = 1e-6
  pi = matrix(1./k, nrow=k, ncol=L)
  rho = array(dim=c(m,m,k,L))
  for (l in 1:L)
	{
    for (r in 1:k)
      rho[,,r,l] = diag(1,m)
  }
  A1 = matrix(seq_len(p), nrow=p, ncol=L)
	require(valse)
  xy = valse:::generateXYdefault(n, p, m, k)

  testFolder = "../data/"
  dir.create(testFolder, showWarnings=FALSE, mode="0755")
  #save inputs
  write.table(as.double(pi), paste(testFolder,"pi",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(rho),paste(testFolder,"rho",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(mini), paste(testFolder,"mini",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(maxi), paste(testFolder,"maxi",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(xy$X), paste(testFolder,"X",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(xy$Y), paste(testFolder,"Y",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(tau),paste(testFolder,"tau",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(A1),paste(testFolder,"A1",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(rangmin),paste(testFolder,"rangmin",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(rangmax),paste(testFolder,"rangmax",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(c(n,p,m,k,L)),paste(testFolder,"dimensions",sep=""),
		row.names=F, col.names=F)

  res = constructionModelesLassoRank(pi,rho,mini,maxi,xy$X,xy$Y,tau,A1,rangmin,rangmax)

  #save output
  write.table(as.double(res$phi), paste(testFolder,"phi",sep=""), row.names=F, col.names=F)
  write.table(as.double(res$llh), paste(testFolder,"llh",sep=""), row.names=F, col.names=F)
}
