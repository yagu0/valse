### Regression matrices
model = Res
K = dim(model$phi)[3]
valMax = max(abs(model$phi))

require(fields)

if (K<4){
  par(mfrow = c(1,K))
} else op = par(mfrow = c(2, (K+1)/2))

## Phi

for (r in 1:K){
  image.plot(t(abs(model$phi[,,r])),
             col=gray(rev(seq(0,64,length.out=65))/65),breaks=seq(0,valMax,length.out=66))
}
par(mfrow = c(1,K),oma = c(0,0,3,0))
mtext("Regression matrices in each cluster", side=3, line=4, font=2, cex=2, col='red')

par(mfrow = c(1,2), oma=c(0,0,3,0))
for (i in 1:4) 
  plot(runif(20), runif(20), 
       main=paste("random plot (",i,")",sep=''))
par(op)
mtext("Four plots", 
      side=3, line=4, font=2, cex=2, col='red')

### Zoom onto two classes we want to compare 
kSel = c(1,2)
par(mfrow = c(1,3))

for (r in kSel){
  image.plot(t(abs(model$phi[,,r])),xaxt="n",yaxt="n", 
             col=gray(rev(seq(0,64,length.out=65))/65),breaks=seq(0,valMax,length.out=66))
}
image.plot(t(abs(model$phi[,,kSel[1]]-model$phi[,,kSel[2]])), 
           col=gray(rev(seq(0,64,length.out=65))/65),breaks=seq(0,valMax,length.out=66))

### Covariance matrices
par(mfrow = c(K, 1))
for (r in 1:K){
  image.plot(matrix(diag(model$rho[,,r]), ncol= 1), 
             col=gray(rev(seq(0,64,length.out=65))/65),breaks=seq(0,valMax,length.out=66))
}

### proportions
Gam = matrix(0, ncol = K, nrow = n)
gam  = Gam
for (i in 1:n){
  for (r in 1:K){
    sqNorm2 = sum( (Y[i,]%*%model$rho[,,r]-X[i,]%*%model$phi[,,r])^2 )
    Gam[i,r] = model$pi[r] * exp(-0.5*sqNorm2)* det(model$rho[,,r])
  }
  gam[i,] = Gam[i,] / sum(Gam[i,])
}
affec = apply(gam, 1,which.max)
gam2 = matrix(NA, ncol = K, nrow = n)
for (i in 1:n){
  gam2[i, affec[i]] = gam[i, affec[i]]
}
boxplot(gam2)

### Mean in each cluster
XY = cbind(X,Y)
XY_class= list()
meanPerClass= matrix(0, ncol = K, nrow = dim(XY)[2])
for (r in 1:K){
  XY_class[[r]] = XY[affec == r, ]
  meanPerClass[,r] = apply(XY_class[[r]], 2, mean)
}

matplot(meanPerClass, type='l')
