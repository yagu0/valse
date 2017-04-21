m=11
p=10

covY = array(0,dim = c(m,m,2))
covY[,,1] = diag(m)
covY[,,2] = diag(m)

Beta = array(0, dim = c(p, m, 2))
Beta[1:4,1:4,1] = 3*diag(4)
Beta[1:4,1:4,2] = -2*diag(4)

#Data = generateXY(100, c(0.5,0.5), rep(0,p), Beta, diag(p), covY)
#
#Res = valse(Data$X,Data$Y, fast=FALSE, plot=FALSE, verbose = TRUE, kmax=2, compute_grid_lambda = FALSE, 
#            grid_lambda = seq(0.2,2,length = 50), size_coll_mod = 50)
