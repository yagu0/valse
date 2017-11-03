library(valse)
testFolder = "data/"

# NOTE: R typing is really terrible. as.double as.matrix ...and so on; don't remove.

#get inputs
npmk = as.matrix(read.table(paste(testFolder,"dimensions",sep="")))
n = npmk[1]; p=npmk[2]; m=npmk[3]; k=npmk[4]
phiInit = array(as.double(as.matrix(read.table(paste(testFolder,"phiInit",sep="")))), dim=c(p,m,k))
rhoInit = array(as.double(as.matrix(read.table(paste(testFolder,"rhoInit",sep="")))), dim=c(m,m,k))
piInit = as.double(as.matrix(read.table(paste(testFolder,"piInit",sep="")))[,])
gamInit = matrix(as.double(as.matrix(read.table(paste(testFolder,"gamInit",sep="")))), n,k)
mini = as.integer(as.matrix(read.table(paste(testFolder,"mini",sep="")))[1])
maxi = as.integer(as.matrix(read.table(paste(testFolder,"maxi",sep="")))[1])
gamma = as.double(as.matrix(read.table(paste(testFolder,"gamma",sep="")))[1])
lambda = as.double(as.matrix(read.table(paste(testFolder,"lambda",sep="")))[1])
X = matrix(as.double(as.matrix(read.table(paste(testFolder,"X",sep="")))), n,p)
Y = matrix(as.double(as.matrix(read.table(paste(testFolder,"Y",sep="")))), n,m)
eps = as.double(as.matrix(read.table(paste(testFolder,"eps",sep="")))[1])

#get outputs
phi = array(as.double(as.matrix(read.table(paste(testFolder,"phi",sep="")))), dim=c(p,m,k))
rho = array(as.double(as.matrix(read.table(paste(testFolder,"rho",sep="")))), dim=c(m,m,k))
pi = as.double(as.matrix(read.table(paste(testFolder,"pi",sep="")))[,])
llh = as.double(as.matrix(read.table(paste(testFolder,"llh",sep="")))[1])
S = array(as.double(as.matrix(read.table(paste(testFolder,"S",sep="")))), dim=c(p,m,k))
affec = as.double(as.matrix(read.table(paste(testFolder,"affec",sep="")))[,])

res = valse::EMGLLF(
	phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,lambda,X,Y,eps,fast=TRUE)

#compare outputs
nd=7 #number of digits
print( all(round(phi,nd) == round(res$phi,nd)) )
print( all(round(rho,nd) == round(res$rho,nd)) )
print( all(round(pi,nd) == round(res$pi,nd)) )
print( all(round(llh,nd) == round(res$llh,nd)) )
print( all(round(S,nd) == round(res$S,nd)) )
print( all(affec == res$affec) )
