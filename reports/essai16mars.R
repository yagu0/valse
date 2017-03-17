p = 10
q = 15
k = 2
D = 20

meanX = rep(0,p)
covX = 0.1*diag(p)

covY = array(dim = c(q,q,k))
covY[,,1] = 0.1*diag(q)
covY[,,2] = 0.2*diag(q)

beta = array(dim = c(p,q,2))
beta[,,2] = matrix(c(rep(2,(D)),rep(0, p*q-D)))
beta[,,1] = matrix(c(rep(1,D),rep(0, p*q-D)))

n = 100

pi = c(0.4,0.6)

data = generateXY(meanX,covX,covY, pi, beta, n)

X = data$X
Y = data$Y

res_valse = valse(X,Y)
