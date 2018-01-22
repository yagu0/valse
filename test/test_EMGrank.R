library(valse)
testFolder = "data/"

# NOTE: R typing is really terrible. as.double as.matrix ...and so on; don't remove.

#get inputs
npmk = as.matrix(read.table(paste(testFolder,"dimensions",sep="")))
n = npmk[1]; p=npmk[2]; m=npmk[3]; k=npmk[4]
Pi = as.double(as.matrix(read.table(paste(testFolder,"Pi",sep="")))[,])
Rho = array(as.double(as.matrix(read.table(paste(testFolder,"Rho",sep="")))), dim=c(m,m,k))
mini = as.integer(as.matrix(read.table(paste(testFolder,"mini",sep="")))[1])
maxi = as.integer(as.matrix(read.table(paste(testFolder,"maxi",sep="")))[1])
X = matrix(as.double(as.matrix(read.table(paste(testFolder,"X",sep="")))), n,p)
Y = matrix(as.double(as.matrix(read.table(paste(testFolder,"Y",sep="")))), n,m)
eps = as.double(as.matrix(read.table(paste(testFolder,"eps",sep="")))[1])
rank = as.double(as.matrix(read.table(paste(testFolder,"rank",sep="")))[,])

#get outputs
phi = array(as.double(as.matrix(read.table(paste(testFolder,"phi",sep="")))), dim=c(p,m,k))
LLF = as.double(as.matrix(read.table(paste(testFolder,"LLF",sep="")))[1])

res = valse::EMGrank(Pi,Rho,mini,maxi,X,Y,eps,rank,fast=TRUE)

#compare outputs
nd=7 #number of digits
print( all(round(phi,nd) == round(res$phi,nd)) )
print( all(round(LLF,nd) == round(res$LLF,nd)) )
