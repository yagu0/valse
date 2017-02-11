generateRunSaveTest_EMGLLF = function(n=200, p=15, m=10, k=3, mini=5, maxi=10,
	gamma=1., lambda=0.5, tau=1e-6)
{
	testFolder = "data/"
	dir.create(testFolder, showWarnings=FALSE, mode="0755")

	require(valse)
	params = valse:::basic_Init_Parameters(n, p, m, k)
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

	res = EMGLLF(params$phiInit,params$rhoInit,params$piInit,params$gamInit,mini,maxi,
		gamma,lambda,io$X,io$Y,tau)

	#save outputs
	write.table(as.double(res$phi), paste(testFolder,"phi",sep=""), row.names=F, col.names=F)
	write.table(as.double(res$rho), paste(testFolder,"rho",sep=""), row.names=F, col.names=F)
	write.table(as.double(res$pi), paste(testFolder,"pi",sep=""), row.names=F, col.names=F)
	write.table(as.double(res$LLF), paste(testFolder,"LLF",sep=""), row.names=F, col.names=F)
	write.table(as.double(res$S), paste(testFolder,"S",sep=""), row.names=F, col.names=F)
}
