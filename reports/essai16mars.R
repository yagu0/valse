meanX = rep(0,6)
covX = 0.1*diag(6)

covY = array(dim = c(5,5,2))
covY[,,1] = 0.1*diag(5)
covY[,,2] = 0.2*diag(5)

beta = array(dim = c(6,5,2))
beta[,,2] = matrix(c(rep(2,12),rep(0, 18)))
beta[,,1] = matrix(c(rep(1,12),rep(0, 18)))

n = 500

pi = c(0.4,0.6)

source('~/valse/R/generateSampleInputs.R')
data = generateXY(meanX,covX,covY, pi, beta, n)

X = data$X
Y = data$Y

k = 2

n = nrow(Y)
m = ncol(Y)
p = ncol(X)

Zinit1 = array(0, dim=c(n))
betaInit1 = array(0, dim=c(p,m,k))
sigmaInit1 = array(0, dim = c(m,m,k))
phiInit1 = array(0, dim = c(p,m,k))
rhoInit1 = array(0, dim = c(m,m,k))
Gam = matrix(0, n, k)
piInit1 = matrix(0,k)
gamInit1 = array(0, dim=c(n,k))
LLFinit1 = list()

require(MASS) #Moore-Penrose generalized inverse of matrix

  distance_clus = dist(X)
  tree_hier = hclust(distance_clus)
  Zinit1 = cutree(tree_hier, k)
  sum(Zinit1==1)
  
  for(r in 1:k)
  {
    Z = Zinit1
    Z_indice = seq_len(n)[Z == r] #renvoit les indices où Z==r
    if (length(Z_indice) == 1) {
      betaInit1[,,r] = ginv(crossprod(t(X[Z_indice,]))) %*%
        crossprod(t(X[Z_indice,]), Y[Z_indice,])
    } else {
      betaInit1[,,r] = ginv(crossprod(X[Z_indice,])) %*%
        crossprod(X[Z_indice,], Y[Z_indice,])
    }
    sigmaInit1[,,r] = diag(m)
    phiInit1[,,r] = betaInit1[,,r] #/ sigmaInit1[,,r]
    rhoInit1[,,r] = solve(sigmaInit1[,,r])
    piInit1[r] = mean(Z == r)
  }
  
  for(i in 1:n)
  {
    for(r in 1:k)
    {
      dotProduct = tcrossprod(Y[i,]%*%rhoInit1[,,r]-X[i,]%*%phiInit1[,,r])
      Gam[i,r] = piInit1[r]*det(rhoInit1[,,r])*exp(-0.5*dotProduct)
    }
    sumGamI = sum(Gam[i,])
    gamInit1[i,]= Gam[i,] / sumGamI
  }
  
  miniInit = 10
  maxiInit = 101
  
  new_EMG = EMGLLF(phiInit1,rhoInit1,piInit1,gamInit1,miniInit,maxiInit,1,0,X,Y,1e-6)
  
  new_EMG$phi
  new_EMG$pi
  LLFEessai = new_EMG$LLF
  LLFinit1 = LLFEessai[length(LLFEessai)]


b = which.max(LLFinit1)
phiInit = phiInit1[,,,b]
rhoInit = rhoInit1[,,,b]
piInit = piInit1[b,]
gamInit = gamInit1[,,b]
