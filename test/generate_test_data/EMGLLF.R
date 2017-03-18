EMGLLF_R = function(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,lambda,X,Y,tau)
{
  #matrix dimensions
  n = dim(X)[1]
  p = dim(phiInit)[1]
  m = dim(phiInit)[2]
  k = dim(phiInit)[3]
  
  #init outputs
  phi = phiInit
  rho = rhoInit
  pi = piInit
  LLF = rep(0, maxi)
  S = array(0, dim=c(p,m,k))
  
  gam = gamInit
  Gram2 = array(0, dim=c(p,p,k))
  ps2 = array(0, dim=c(p,m,k))
  b = rep(0, k)
  pen = matrix(0, maxi, k)
  X2 = array(0, dim=c(n,p,k))
  Y2 = array(0, dim=c(n,m,k))
  dist = 0
  dist2 = 0
  ite = 1
  pi2 = rep(0, k)
  ps = matrix(0, m,k)
  nY2 = matrix(0, m,k)
  ps1 = array(0, dim=c(n,m,k))
  Gam = matrix(0, n,k)
  EPS = 1E-15
  
  while(ite <= mini || (ite <= maxi && (dist >= tau || dist2 >= sqrt(tau))))
	{
    Phi = phi
    Rho = rho
    Pi = pi

    #calcul associé à Y et X
    for(r in 1:k)
		{
      for (mm in 1:m)
        Y2[,mm,r] = sqrt(gam[,r]) * Y[,mm]
      for (i in 1:n)
        X2[i,,r] = sqrt(gam[i,r]) * X[i,]
      for (mm in 1:m)
        ps2[,mm,r] = crossprod(X2[,,r],Y2[,mm,r])
      for (j in 1:p)
			{
        for (s in 1:p)
          Gram2[j,s,r] = crossprod(X2[,j,r], X2[,s,r])
      }
    }
    
    ##########
    #Etape M #
    ##########
    
    #pour pi
    for (r in 1:k){
      b[r] = sum(abs(phi[,,r]))}
    gam2 = colSums(gam)
    a = sum(gam %*% log(pi))
    
    #tant que les props sont negatives
    kk = 0
    pi2AllPositive = FALSE
    while (!pi2AllPositive)
		{
      pi2 = pi + 0.1^kk * ((1/n)*gam2 - pi)
      pi2AllPositive = all(pi2 >= 0)
      kk = kk+1
    }

    #t[m] la plus grande valeur dans la grille O.1^k tel que ce soit décroissante ou constante
    while( kk < 1000 && -a/n + lambda * sum(pi^gamma * b) <
			-sum(gam2 * log(pi2))/n + lambda * sum(pi2^gamma * b) )
		{
      pi2 = pi + 0.1^kk * (1/n*gam2 - pi)
      kk = kk + 1
    }
    t = 0.1^kk
    pi = (pi + t*(pi2-pi)) / sum(pi + t*(pi2-pi))
    
    #Pour phi et rho
    for (r in 1:k)
		{
      for (mm in 1:m)
			{
        for (i in 1:n)
				{
          ps1[i,mm,r] = Y2[i,mm,r] * sum(X2[i,,r] * phi[,mm,r])
        }
        ps[mm,r] = sum(ps1[,mm,r])
        nY2[mm,r] = sum(Y2[,mm,r]^2)
        rho[mm,mm,r] = (ps[mm,r]+sqrt(ps[mm,r]^2+4*nY2[mm,r]*gam2[r])) / (2*nY2[mm,r])
      }
    }
    for (r in 1:k)
		{
      for (j in 1:p)
			{
        for (mm in 1:m)
				{
          S[j,mm,r] = -rho[mm,mm,r]*ps2[j,mm,r] + sum(phi[-j,mm,r] * Gram2[j, setdiff(1:p,j),r])
          if (abs(S[j,mm,r]) <= n*lambda*(pi[r]^gamma))
            phi[j,mm,r]=0
          else if(S[j,mm,r] > n*lambda*(pi[r]^gamma))
            phi[j,mm,r] = (n*lambda*(pi[r]^gamma)-S[j,mm,r]) / Gram2[j,j,r]
          else
            phi[j,mm,r] = -(n*lambda*(pi[r]^gamma)+S[j,mm,r]) / Gram2[j,j,r]
        }
      }
    }

    ##########
    #Etape E #
    ##########
    sumLogLLF2 = 0
    for (i in 1:n)
		{
      #precompute sq norms to numerically adjust their values
      sqNorm2 = rep(0,k)
      for (r in 1:k){
        sqNorm2[r] = sum( (Y[i,]%*%rho[,,r]-X[i,]%*%phi[,,r])^2 )}

      #compute Gam(:,:) using shift determined above
      sumLLF1 = 0.0;
      for (r in 1:k)
			{
				Gam[i,r] = pi[r] * exp(-0.5*sqNorm2[r])* det(rho[,,r])
        sumLLF1 = sumLLF1 + Gam[i,r] / (2*base::pi)^(m/2)
      }
      sumLogLLF2 = sumLogLLF2 + log(sumLLF1)
      sumGamI = sum(Gam[i,])
      if(sumGamI > EPS)
        gam[i,] = Gam[i,] / sumGamI
      else
        gam[i,] = rep(0,k)
    }

    sumPen = sum(pi^gamma * b)
    LLF[ite] = -sumLogLLF2/n + lambda*sumPen

    dist = ifelse( ite == 1, LLF[ite], (LLF[ite]-LLF[ite-1]) / (1+abs(LLF[ite])) )

    Dist1 = max( (abs(phi-Phi)) / (1+abs(phi)) )
    Dist2 = max( (abs(rho-Rho)) / (1+abs(rho)) )
    Dist3 = max( (abs(pi-Pi)) / (1+abs(Pi)) )
    dist2 = max(Dist1,Dist2,Dist3)

    ite = ite+1
  }
  
  affec = apply(gam, 1,which.max)
  return(list("phi"=phi, "rho"=rho, "pi"=pi, "LLF"=LLF, "S"=S, "affec" = affec ))
}
