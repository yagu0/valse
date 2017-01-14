generateRunSaveTest_EMGLLF = function(n=200, p=15, m=10, k=3, mini=5, maxi=10,
	gamma=1., lambda=0.5, tau=1e-6)
{
	testFolder = "data/"
	dir.create(testFolder, showWarnings=FALSE, mode="0755")
	delimiter = " "

	#save inputs
	params = basicInitParameters(n, p, m, k)
	io = generateIOdefault(n, p, m, k)
	write.table(paste(testFolder,"phiInit",sep=""), params$phiInit, sep=delimiter);
	write.table(paste(testFolder,"rhoInit",sep=""), params$rhoInit, sep=delimiter);
	write.table(paste(testFolder,"piInit",sep=""), params$piInit, sep=delimiter);
	write.table(paste(testFolder,"gamInit",sep=""), params$gamInit, sep=delimiter);
	write.table(paste(testFolder,"mini",sep=""), mini, sep=delimiter);
	write.table(paste(testFolder,"maxi",sep=""), maxi, sep=delimiter);
	write.table(paste(testFolder,"gamma",sep=""), gamma, sep=delimiter);
	write.table(paste(testFolder,"lambda",sep=""), lambda, sep=delimiter);
	write.table(paste(testFolder,"X",sep=""), io$X, sep=delimiter);
	write.table(paste(testFolder,"Y",sep=""), io$Y, sep=delimiter);
	write.table(paste(testFolder,"tau",sep=""), tau, sep=delimiter);
	write.table(paste(testFolder,"dimensions",sep=""), c(n,p,m,k), sep=delimiter);

	res = EMGLLF(params$phiInit,params$rhoInit,params$piInit,params$gamInit,mini,maxi,
		gamma,lambda,io$X,io$Y,tau);

	#save outputs
	write.table(paste(testFolder,"phi",sep=""), res$phi, sep=delimiter);
	write.table(paste(testFolder,"rho",sep=""), res$rho, sep=delimiter);
	write.table(paste(testFolder,"pi",sep=""), res$pi, sep=delimiter);
	write.table(paste(testFolder,"LLF",sep=""), res$LLF, sep=delimiter);
	write.table(paste(testFolder,"S",sep=""), res$S, sep=delimiter);
}
