#' initialization of the EM algorithm
#'
#' @param k number of components
#' @param X matrix of covariates (of size n*p)
#' @param Y matrix of responses (of size n*m)
#' @param tau threshold to stop EM algorithm
#'
#' @return a list with phiInit, rhoInit, piInit, gamInit
#' @export
initSmallEM = function(k,X,Y,tau)
{
	n = nrow(Y)
	m = ncol(Y)
	p = ncol(X)
  
	Zinit1 = array(0, dim=c(n,20)) #doute sur la taille
	betaInit1 = array(0, dim=c(p,m,k,20))
	sigmaInit1 = array(0, dim = c(m,m,k,20))
	phiInit1 = array(0, dim = c(p,m,k,20))
	rhoInit1 = array(0, dim = c(m,m,k,20))
	Gam = matrix(0, n, k)
	piInit1 = matrix(0,20,k)
	gamInit1 = array(0, dim=c(n,k,20))
	LLFinit1 = list()

	require(MASS) #Moore-Penrose generalized inverse of matrix
	require(mclust) # K-means with selection of K
	for(repet in 1:20)
	{
		clusters = Mclust(X,k) #default distance : euclidean  #Mclust(matrix(c(X,Y)),k)
		Zinit1[,repet] = clusters$classification
		
		for(r in 1:k)
		{
			Z = Zinit1[,repet]
			Z_bin = vec_bin(Z,r)
			Z_vec = Z_bin$vec #vecteur 0 et 1 aux endroits o? Z==r
			Z_indice = Z_bin$indice #renvoit les indices o? Z==r
			
			betaInit1[,,r,repet] = ginv( crossprod(X[Z_indice,]) )   %*%   crossprod(X[Z_indice,], Y[Z_indice,]) 
			sigmaInit1[,,r,repet] = diag(m)
			phiInit1[,,r,repet] = betaInit1[,,r,repet]/sigmaInit1[,,r,repet]
			rhoInit1[,,r,repet] = solve(sigmaInit1[,,r,repet])
			piInit1[repet,r] = sum(Z_vec)/n
		}
		
		for(i in 1:n)
		{
			for(r in 1:k)
			{
				dotProduct = (Y[i,]%*%rhoInit1[,,r,repet]-X[i,]%*%phiInit1[,,r,repet]) %*% (Y[i,]%*%rhoInit1[,,r,repet]-X[i,]%*%phiInit1[,,r,repet])
				Gam[i,r] = piInit1[repet,r]*det(rhoInit1[,,r,repet])*exp(-0.5*dotProduct)
			}
			sumGamI = sum(Gam[i,])
			gamInit1[i,,repet]= Gam[i,] / sumGamI
		}
		
		miniInit = 10
		maxiInit = 11
		
		new_EMG = .Call("EMGLLF_core",phiInit1[,,,repet],rhoInit1[,,,repet],piInit1[repet,],
										gamInit1[,,repet],miniInit,maxiInit,1,0,X,Y,tau)
		LLFEessai = new_EMG$LLF
		LLFinit1[repet] = LLFEessai[length(LLFEessai)]
	}

	b = which.max(LLFinit1)
	phiInit = phiInit1[,,,b]
	rhoInit = rhoInit1[,,,b]
	piInit = piInit1[b,]
	gamInit = gamInit1[,,b]

	return (list(phiInit=phiInit, rhoInit=rhoInit, piInit=piInit, gamInit=gamInit))
}
