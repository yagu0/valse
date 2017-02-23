constructionModelesLassoMLE = function(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda,
	X,Y,seuil,tau,A1,A2)
{
  n = dim(X)[1];
  p = dim(phiInit)[1]
  m = dim(phiInit)[2]
  k = dim(phiInit)[3]
  L = length(glambda)

  #output parameters
  phi = array(0, dim=c(p,m,k,L))
  rho = array(0, dim=c(m,m,k,L))
  pi = matrix(0, k, L)
  llh = matrix(0, L, 2) #log-likelihood

  for(lambdaIndex in 1:L)
	{
    a = A1[,1,lambdaIndex]
    a = a[a!=0]
    if(length(a)==0)
      next

    res = EMGLLF(phiInit[a,,],rhoInit,piInit,gamInit,mini,maxi,gamma,0.,X[,a],Y,tau)

		#TODO: supprimer ça et utiliser parLapply(...)
    for (j in 1:length(a))
      phi[a[j],,,lambdaIndex] = res$phi[j,,]
    rho[,,,lambdaIndex] = res$rho
    pi[,lambdaIndex] = res$pi

    dimension = 0
    for (j in 1:p) #TODO: doit pouvoir être fait beaucoup plus simplement
		{
      b = A2[j,2:dim(A2)[2],lambdaIndex]
      b = b[b!=0]
      if (length(b) > 0)
        phi[A2[j,1,lambdaIndex],b,,lambdaIndex] = 0.
      c = A1[j,2:dim(A1)[2],lambdaIndex]
      dimension = dimension + sum(c!=0)
    }

    #on veut calculer l'EMV avec toutes nos estimations
		densite = matrix(0, nrow=n, ncol=L)
		for (i in 1:n) #TODO: pas besoin de cette boucle (vectorize)
		{
			for (r in 1:k)
			{
				delta = Y[i,]%*%rho[,,r,lambdaIndex] - (X[i,a]%*%phi[a,,r,lambdaIndex]);
				densite[i,lambdaIndex] = densite[i,lambdaIndex] + pi[r,lambdaIndex] *
					det(rho[,,r,lambdaIndex])/(sqrt(2*base::pi))^m * exp(-tcrossprod(delta)/2.0)
			}
		}
		llh[lambdaIndex,1] = sum(log(densite[,lambdaIndex]))
		llh[lambdaIndex,2] = (dimension+m+1)*k-1
  }
  return (list("phi"=phi, "rho"=rho, "pi"=pi, "llh" = llh))
}
