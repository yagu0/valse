library(valse)
n = 50; m = 10; p = 5
beta = array(0, dim=c(p,m,2))
beta[,,1] = 1
beta[,,2] = 2
data = generateXY(n, c(0.4,0.6), rep(0,p), beta, diag(0.5, p), diag(0.5, m))
X = data$X ; Y = data$Y
res = runValse(X, Y, kmax = 5)
