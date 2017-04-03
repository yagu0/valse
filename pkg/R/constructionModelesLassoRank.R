#' constructionModelesLassoRank
#'
#' TODO: description
#'
#' @param ...
#'
#' @return ...
#'
#' export
constructionModelesLassoRank = function(pi, rho, mini, maxi, X, Y, tau, A1, rangmin,
	rangmax, ncores, verbose=FALSE)
{
  n = dim(X)[1]
  p = dim(X)[2]
  m = dim(rho)[2]
  k = dim(rho)[3]
  L = dim(A1)[2]

	# On cherche les rangs possiblement intéressants
  deltaRank = rangmax - rangmin + 1
  Size = deltaRank^k
  Rank = matrix(0, nrow=Size, ncol=k)
  for (r in 1:k)
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

  if (ncores > 1)
	{
    cl = parallel::makeCluster(ncores)
    parallel::clusterExport( cl, envir=environment(),
			varlist=c("A1","Size","Pi","Rho","mini","maxi","X","Y","tau",
			"Rank","m","phi","ncores","verbose") )
	}

	computeAtLambda <- function(lambdaIndex)
	{
		if (ncores > 1)
			require("valse") #workers start with an empty environment

    # on ne garde que les colonnes actives
    # 'active' sera l'ensemble des variables informatives
    active = A1[,lambdaIndex]
    active = active[-(active==0)]
		phi = array(0, dim=c(p,m,k,Size))
		llh = matrix(0, Size, 2) #log-likelihood
    if (length(active) > 0)
		{
      for (j in 1:Size)
			{
        res = EMGrank(Pi[,lambdaIndex], Rho[,,,lambdaIndex], mini, maxi,
					X[,active], Y, tau, Rank[j,])
        llh = rbind(llh,
					c( res$LLF, sum(Rank[j,] * (length(active)- Rank[j,] + m)) ) )
        phi[active,,,] = rbind(phi[active,,,], res$phi)
      }
    }
		list("llh"=llh, "phi"=phi)
	}

	#Pour chaque lambda de la grille, on calcule les coefficients
  out =
		if (ncores > 1)
			parLapply(cl, seq_along(glambda), computeAtLambda)
		else
			lapply(seq_along(glambda), computeAtLambda)

	if (ncores > 1)
    parallel::stopCluster(cl)

	# TODO: this is a bit ugly. Better use bigmemory and fill llh/phi in-place
	# (but this also adds a dependency...)
	llh <- do.call( rbind, lapply(out, function(model) model$llh) )
	phi <- do.call( rbind, lapply(out, function(model) model$phi) )
	list("llh"=llh, "phi"=phi)
}
