library("mclust")
#library("R.matlab", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
#redrawData = TRUE
#if (redrawData==TRUE){
  ###########
  ## Model
  ###########
  K = 2
  p = 48
  T = seq(0,1.5,length.out = p)
  T2 = seq(0,3, length.out = 2*p)
  n = 100
  x1 = cos(2*base::pi*T) + 0.2*cos(4*2*base::pi*T) + 0.3*c(rep(0,round(length(T)/7)),rep(1,round(length(T)*(1-1/7))))+1
  plot(T,x1)
  lines(T,x1)
  
  sigmaX = 0.12
  sigmaY = 0.12
  beta = list()
  p1= 0.5
  beta[[1]] =diag(c(rep(p1,5),rep(1,5), rep(p1,5), rep(1, p-15)))
  p2 = 1
  beta[[2]] = diag(c(rep(p2,5),rep(1,5), rep(p2,5), rep(1, p-15)))
  ITE = 100
  ARI1 = ARI2 = ARI3 = rep(0,ITE)
  XY = array(0, dim = c(ITE, 2*p,n))
  XYproj = array(0, dim=c(ITE, 96,n))

  affec = list()
  ###########
  ## Iterations
  ###########
  for (ite in c(1:ITE)){
    ###########
    ##Sample
    ###########
    x = x1 + matrix(rnorm(n*p, 0, sigmaX), ncol = n)
    affec[[ite]] = sample(c(1,2), n, replace = TRUE)
    y  = x
    xy = matrix(0,ncol=n, nrow= 2*p)
     for (i in c(1:n)){
       y[,i] = x[,i] %*% beta[[affec[[ite]][i]]] + rnorm(p, 0, sigmaY)
       xy[,i] = c(x[,i],y[,i])
       XY[ite,,i] = xy[,i] - mean(xy[,i])
    #   Dx = dwt(x[,i], filter='haar')@W
    #   Dx = rev(unlist(Dx))
    #   Dx = Dx[2:(1+3+6+12+24)]
    #   Ax = dwt(x[,i], filter='haar')@V
    #   Ax = rev(unlist(Ax))
    #   Ax = Ax[2:(1+3)]
    #   Dy = dwt(y[,i], filter='haar')@W
    #   Dy = rev(unlist(Dy))
    #   Dy = Dy[2:(1+3+6+12+24)]
    #   Ay = dwt(y[,i], filter='haar')@V
    #   Ay = rev(unlist(Ay))
    #   Ay = Ay[2:(1+3)]
    #   XYproj[ite,,i] = c(Ax,Dx,Ay,Dy)
     }
    print(ite)
    #
    #
  }
  xy[c(7,55),] = NA
 # write.table(XY,'data.csv', row.names=FALSE, col.names=FALSE)
matplot(T2,xy[,affec[[ite]]==1],type='l', col='red', lty = 1)
matplot(T2,xy[,affec[[ite]]==2],type='l', col='black', add=TRUE, lty= 1)
abline(v = 1.5)
text(0.75,0,'X', cex = 2 )
text(0.75+1.5,0,'Y', cex = 2 )
#proj =  read.table('dataProj.csv')
#}


#matplot(T,x,type='l', col='black', xlab = '', ylab='', lwd=1.5,lty=1)
#matplot(T,y[,affec[[ite]]==1],type='l', col='red', xlab = '', ylab='', lwd=1.5,lty=1)
#matplot(T,y[,affec[[ite]]==2],type='l', col='black', add=TRUE,lwd=2,lty=1)
# proj2 = array(0,dim=c(ITE,2*p,n))
# for (ite in c(1:ITE)){
#   for (i in c(1:n)){
#     A = proj[ite,(1+(i-1)*96):(i*96)]
#     for (j in 1:96){
#       proj2[ite,j,i] = A[1,j]
#     }
#   }
#   print(ite)
# }
###########
## Iterations
###########
Kmod2 = Kmod1 = rep(0,ITE)
Kmod3 = rep(0,ITE)
for (ite in c(1:ITE)){
  print(ite)
  ###########
  ## k-means 1
  ###########
  mod1 = Mclust(t(XY[ite,,]),G = 1:2, mode='VII')
  ARI1[ite] = adjustedRandIndex(mod1$classification, affec[[ite]])
  Kmod1[ite] = mod1$G
  # ###########
  # ## k-means 2
  # ###########
  # #proj2 =
  # mod2 = Mclust(t(XYproj[ite,,]),G = 1:8, mode='VII')
  # ARI2[ite] = adjustedRandIndex(mod2$classification, affec[[ite]])
  # Kmod2[ite] = mod2$G
  # ###########
  # ## k-means 1
  # ###########
  # #proj3 =
  # mod3 = Mclust(t(XYproj[ite,c(4:12,52:60),]),G = 1:8, mode='VII')
  # ARI3[ite] = adjustedRandIndex(mod3$classification, affec[[ite]])
  # Kmod3[ite] = mod3$G
}
ARI0 = rep(1,ITE)
par(cex.lab=1.5)
par(cex.axis=1.5)
boxplot(ARI0,ARI1, names = c('LassoMLE','K-means'), lwd=1.3)
table(Kmod1)
table(Kmod2)
