m=11
p=10

covY = array(0,dim = c(m,m,2))
covY[,,1] = diag(m)
covY[,,2] = diag(m)

Beta = array(0, dim = c(p, p, 2))
Beta[1:4,1:4,1] = 3*diag(4)
Beta[1:4,1:4,2] = -2*diag(4)

Data = generateXY(100, c(0.5,0.5), rep(0,p), Beta, diag(m), covY)
