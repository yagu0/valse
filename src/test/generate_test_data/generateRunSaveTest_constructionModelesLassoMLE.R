generateRunSaveTest_constructionModelesLassoMLE = function(n=200, p=15, m=10, k=3, mini=5,
	maxi=10, gamma=1.0, glambda=list(0.0,0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.5,0.7,0.85,0.99))
{
  testFolder = "data/"
  dir.create(testFolder, showWarnings=FALSE, mode="0755")

	require(valse)
  params = valse:::basic_Init_Parameters(n, p, m, k)
  A2 = array(0, dim=c(p, m+1, L))
  A1 = array(0, dim=c(p, m+1, L))
  for (i in 1:L)
	{
    A2[,1,i] = seq_len(p)
    A2[,2,i] = seq_len(5)
    A1[,1,i] = seq_len(p)
    A1[,2,i] = seq_len(5)
  }
  io = generateIOdefault(n, p, m, k)

  #save inputs
  write.table(as.double(params$phiInit), paste(testFolder,"phiInit",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(params$rhoInit), paste(testFolder,"rhoInit",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(params$piInit), paste(testFolder,"piInit",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(params$gamInit), paste(testFolder,"gamInit",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(mini), paste(testFolder,"mini",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(maxi), paste(testFolder,"maxi",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(gamma), paste(testFolder,"gamma",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(lambda), paste(testFolder,"lambda",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(io$X), paste(testFolder,"X",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(io$Y), paste(testFolder,"Y",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(tau), paste(testFolder,"tau",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(c(n,p,m,k)), paste(testFolder,"dimensions",sep=""),
		row.names=F, col.names=F)

  res = constructionModelesLassoMLE(
		params$phiInit,params$rhoInit,params$piInit,params$gamInit,
		mini,maxi,gamma,glambda,io$X,io$Y,seuil,tau,A1,A2)

  #save output
  write.table(as.double(res$phi), paste(testFolder,"phi",sep=""), row.names=F, col.names=F)
  write.table(as.double(res$rho), paste(testFolder,"rho",sep=""), row.names=F, col.names=F)
  write.table(as.double(res$pi), paste(testFolder,"pi",sep=""), row.names=F, col.names=F)
  write.table(as.double(res$llh), paste(testFolder,"llh",sep=""), row.names=F, col.names=F)
}
