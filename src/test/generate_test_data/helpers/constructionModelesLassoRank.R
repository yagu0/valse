constructionModelesLassoRank = function(pi,rho,mini,maxi,X,Y,tau,A1,rangmin,rangmax)
{
  #get matrix sizes
  n = dim(X)[1]
  p = dim(X)[2]
  m = dim(rho)[2]
  k = dim(rho)[3]
  L = dim(A1)[2]

	# On cherche les rangs possiblement intéressants
  deltaRank = rangmax - rangmin + 1
  Size = deltaRank^k
  Rank = matrix(0, nrow=Size, ncol=k)
  for(r in 1:k)
	{
		# On veut le tableau de toutes les combinaisons de rangs possibles
		# Dans la première colonne : on répète (rangmax-rangmin)^(k-1) chaque chiffre :
		#   ça remplit la colonne
		# Dans la deuxieme : on répète (rangmax-rangmin)^(k-2) chaque chiffre,
		#   et on fait ça (rangmax-rangmin)^2 fois
		# ...
		# Dans la dernière, on répète chaque chiffre une fois,
		#   et on fait ça (rangmin-rangmax)^(k-1) fois.
    Rank[,r] = rangmin + rep(0:(deltaRank-1), deltaRank^(r-1), each=deltaRank^(k-r))
  }

	# output parameters
  phi = array(0, dim=c(p,m,k,L*Size))
  llh = matrix(0, L*Size, 2) #log-likelihood
  for(lambdaIndex in 1:L)
	{
    # on ne garde que les colonnes actives
    # 'active' sera l'ensemble des variables informatives
    active = A1[,lambdaIndex]
    active = active[-(active==0)]
    if (length(active) > 0)
		{
      for (j in 1:Size)
			{
        res = EMGrank(Pi[,lambdaIndex], Rho[,,,lambdaIndex], mini, maxi,
					X[,active], Y, tau, Rank[j,])
        llh[(lambdaIndex-1)*Size+j,] =
					c( res$LLF, sum(Rank[j,] * (length(active)- Rank[j,] + m)) )
        phi[active,,,(lambdaIndex-1)*Size+j] = res$phi
      }
    }
  }
  return (list(phi=phi, llh = llh))
}
