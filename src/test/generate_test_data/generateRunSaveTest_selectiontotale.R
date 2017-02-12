generateRunSaveTest_selectiontotale= function(n=200, p=15, m=10, k=3, mini=5, maxi=10, gamma=1.0,
	glambda=list(0.0,0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.5,0.7,0.85,0.99))
{
	tau = 1e-6;
	seuil = 1e-15;
	require(valse)
  params = valse:::basicInitParameters(n, p, m, k)
  xy = valse:::generateXYdefault(n, p, m, k)
	L = length(glambda)

  testFolder = "data/"
  dir.create(testFolder, showWarnings=FALSE, mode="0755")
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
  write.table(as.double(glambda), paste(testFolder,"lambda",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(xy$X), paste(testFolder,"X",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(xy$Y), paste(testFolder,"Y",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(seuil), paste(testFolder,"seuil",sep=""),
		row.names=F, col.names=F)
  write.table(as.double(tau), paste(testFolder,"tau",sep=""),
		row.names=F, col.names=F)
  write.table(as.integer(c(n,p,m,k,L)), paste(testFolder,"dimensions",sep=""),
		row.names=F, col.names=F)

  res = selectiontotale(params$phiInit,params$rhoInit,params$piInit,params$gamInit,
		mini,maxi,gamma,glambda,xy$X,xy$Y,seuil, tau)

  #save output
  write.table(as.double(res$A1), paste(testFolder,"A1",sep=""), row.names=F, col.names=F)
  write.table(as.double(res$A2), paste(testFolder,"A2",sep=""), row.names=F, col.names=F)
  write.table(as.double(res$rho), paste(testFolder,"rho",sep=""), row.names=F, col.names=F)
  write.table(as.double(res$pi), paste(testFolder,"pi",sep=""), row.names=F, col.names=F)
}
