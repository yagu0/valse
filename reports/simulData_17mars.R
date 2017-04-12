simulData_17mars = function(ite){
  set.seed = 22021989+ite
  
  ###########
  ## Modele
  ###########
  K = 2
  p = 48
  T = seq(0,1.5,length.out = p)
  T2 = seq(0,3, length.out = 2*p)
  n = 100
  x1 = cos(2*base::pi*T) + 0.2*cos(4*2*base::pi*T) + 0.3*c(rep(0,round(length(T)/7)),rep(1,round(length(T)*(1-1/7))))+1
  sigmaX = 0.12
  sigmaY = 0.12
  beta = list()
  p1= 0.5
  beta[[1]] =diag(c(rep(p1,5),rep(1,5), rep(p1,5), rep(1, p-15)))
  p2 = 1
  beta[[2]] = diag(c(rep(p2,5),rep(1,5), rep(p2,5), rep(1, p-15)))
  ARI1 = ARI2 = ARI3 = 0
  
  ###########
  ## Data + Projection
  ###########
  require(wavelets)
  XY = array(0, dim = c(2*p,n))
  XYproj = array(0, dim=c(96,n))
  x = x1 + matrix(rnorm(n*p, 0, sigmaX), ncol = n)
  affec = sample(c(1,2), n, replace = TRUE)
  y  = x
  xy = matrix(0,ncol=n, nrow= 2*p)
  for (i in c(1:n)){
    y[,i] = x[,i] %*% beta[[affec[i]]] + rnorm(p, 0, sigmaY)
    xy[,i] = c(x[,i],y[,i])
    XY[,i] = xy[,i] - mean(xy[,i])
    Dx = dwt(x[,i], filter='haar')@W
    Dx = rev(unlist(Dx))
    Dx = Dx[2:(1+3+6+12+24)]
    Ax = dwt(x[,i], filter='haar')@V
    Ax = rev(unlist(Ax))
    Ax = Ax[2:(1+3)]
    Dy = dwt(y[,i], filter='haar')@W
    Dy = rev(unlist(Dy))
    Dy = Dy[2:(1+3+6+12+24)]
    Ay = dwt(y[,i], filter='haar')@V
    Ay = rev(unlist(Ay))
    Ay = Ay[2:(1+3)]
    XYproj[,i] = c(Ax,Dx,Ay,Dy)
  }
  
  res_valse = valse(x,y, kmax=2, verbose=TRUE, plot=FALSE, size_coll_mod = 200)
  res_valse_proj = valse(XYproj[1:p,],XYproj[(p+1):(2*p),], kmax=2, verbose=TRUE, plot=FALSE, size_coll_mod = 200)
  
  save(res_valse,file=paste("Res_",ite, ".RData",sep=""))
  save(res_valse_proj,file=paste("ResProj_",ite, ".RData",sep=""))
}
