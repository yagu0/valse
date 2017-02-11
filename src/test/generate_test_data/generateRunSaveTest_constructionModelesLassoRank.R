generateRunSaveTest_constructionModelesLassoRank = function(n=200, p=15, m=10, L=12, mini=5,
	maxi=10, gamma=1.0, rangmin=3, rangmax=6)
{
  testFolder = "data/"
  dir.create(testFolder, showWarnings=FALSE, mode="0755")

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
	require(valse)
  io = valse:::generateIOdefault(n, p, m, k)
  A1 = matrix(0,p,L)
  for(i in 1:L){
    A1[,i] = seq(1,p)
  }

  #save inputs
  write.table(as.double(rho),paste(testFolder,"rho",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(pi), paste(testFolder,"pi",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(mini), paste(testFolder,"mini",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(maxi), paste(testFolder,"maxi",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(io$X), paste(testFolder,"X",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(io$Y), paste(testFolder,"Y",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(tau),paste(testFolder,"tau",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(A1),paste(testFolder,"A1",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(rangmin),paste(testFolder,"rangmin",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(rangmax),paste(testFolder,"rangmax",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(c(n,p,m,k)),paste(testFolder,"dimensions",sep=""),
		row.names=F, col.names=F)

  res = constructionModelesLassoRank(pi,rho,mini,maxi,X,Y,tau,A1,rangmin,rangmax)

  #save output
  write.table(as.double(res$phi), paste(testFolder,"phi",sep=""), row.names=F, col.names=F)
  write.table(as.double(res$llh), paste(testFolder,"llh",sep=""), row.names=F, col.names=F)
}
