m=6
p=6

covY = array(0,dim = c(m,m,2))
covY[,,1] = diag(m)
covY[,,2] = diag(m)

Beta = array(0, dim = c(p, m, 2))
Beta[1:4,1:4,1] = 3*diag(4)
Beta[1:4,1:4,2] = -2*diag(4)

Data = generateXY(200, c(0.5,0.5), rep(0,p), Beta, diag(p), covY)
#  
Res = valse(Data$X,Data$Y, fast=TRUE, plot=FALSE, verbose = TRUE, kmax=3, size_coll_mod = 50, selecMod = "DDSE", mini = 50, maxi=100)
plot(Res$tableau[,3], -Res$tableau[,4])
