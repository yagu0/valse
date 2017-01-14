testEMGLLF = function()
{
	testFolder = "data/"
	array_delimiter = " "

	#get dimensions
	dimensions = read.table(paste(testFolder,"dimensions",sep=""), header=FALSE,
		sep=array_delimiter)
	n = dimensions[1]
	p = dimensions[2]
	m = dimensions[3]
	k = dimensions[4]

	#get all input arrays
	phiInit = read.table(paste(testFolder,"phiInit",sep=""), header=FALSE, sep=array_delimiter)
	rhoInit = read.table(paste(testFolder,"rhoInit",sep=""), header=FALSE, sep=array_delimiter)
	piInit = read.table(paste(testFolder,"piInit",sep=""), header=FALSE, sep=array_delimiter)
	gamInit = read.table(paste(testFolder,"gamInit",sep=""), header=FALSE, sep=array_delimiter)
	mini = read.table(paste(testFolder,"mini",sep=""), header=FALSE, sep=array_delimiter)
	maxi = read.table(paste(testFolder,"maxi",sep=""), header=FALSE, sep=array_delimiter)
	gamma = read.table(paste(testFolder,"gamma",sep=""), header=FALSE, sep=array_delimiter)
	lambda = read.table(paste(testFolder,"lambda",sep=""), header=FALSE, sep=array_delimiter)
	X = read.table(paste(testFolder,"X",sep=""), header=FALSE, sep=array_delimiter)
	Y = read.table(paste(testFolder,"Y",sep=""), header=FALSE, sep=array_delimiter)
	tau = read.table(paste(testFolder,"tau",sep=""), header=FALSE, sep=array_delimiter)

	#run EMGLLF.c
	EMG = .Call("EMGLLF_core",phiInit,rhoInit,piInit1,gamInit,mini,maxi,gamma,lambda,X,Y,tau,
		PACKAGE="valse")
	phi = EMG$phi
	rho = EMG$rho
	pi = EMG$pi
	LLF = EMG$LLF
	S = EMG$S

	#get all stored outputs
	ref_phi = read.table(paste(testFolder,"phi",sep=""), header=FALSE, sep=array_delimiter)
	ref_rho = read.table(paste(testFolder,"rho",sep=""), header=FALSE, sep=array_delimiter)
	ref_pi = read.table(paste(testFolder,"pi",sep=""), header=FALSE, sep=array_delimiter)
	ref_LLF = read.table(paste(testFolder,"LLF",sep=""), header=FALSE, sep=array_delimiter)
	ref_S = read.table(paste(testFolder,"S",sep=""), header=FALSE, sep=array_delimiter)

	#check that output correspond to stored output
	tol = 1e-5;
	checkOutput("phi",phi,ref_phi,tol);
	checkOutput("rho",rho,ref_rho,tol);
	checkOutput("pi",pi,ref_pi,tol);
	checkOutput("LLF",LLF,ref_LLF,tol);
	checkOutput("S",S,ref_S,tol);
}
